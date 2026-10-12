include(joinpath(@__DIR__, "setup.jl"))
using LinearAlgebra, Random, Test, JSON3
import Nemo
const MP=CM.Products
const MN=MP.MN
const F=MP.F

function expand_product_test(cores)
    partial=[ones(F, 1, 1)]
    for core in cores
        d=size(core, 2)
        next=[zeros(F, d*size(partial[1], 1), d*size(partial[1], 2)) for _ in axes(core, 4)]
        for a in axes(core, 4), b in axes(core, 1)
            next[a]+=kron(partial[b], core[b, :, :, a])
        end
        partial=next
    end
    return only(partial)
end

@testset "Interval mixed products and independent exact replay" begin
    rng = MersenneTwister(9087)
    for n in (2, 3, 4), weighted in (false, true), arithmetic in (:rational, :interval)
        da, db = vcat(1, fill(2, n-1), 1), vcat(1, fill(3, n-1), 1)
        A = [F.(rand(rng, -2:2, da[j], 2, 2, da[j+1]))/7 for j in 1:n]
        B = [F.(rand(rng, -3:3, db[j], 2, 2, db[j+1]))/11 for j in 1:n]
        gauges = vcat([ones(F, 1, 1)], [F[1 j//7; -1//5 1] for j in 1:(n-1)], [ones(F, 1, 1)])
        weights = weighted ? F[1, 3] : nothing
        events = Any[]
        proof = MP.product_candidate(
            A,
            B;
            maxdim = 2,
            weights,
            first_gauges = gauges,
            candidate_arithmetic = arithmetic,
            progress = row -> push!(events, row),
        )
        @test length(events) == n && last(events).site == n
        X = expand_product_test(A)*expand_product_test(B)-expand_product_test(proof.cores)
        w = reduce(kron, fill(weighted ? weights : F[1, 1], n))
        actual = sum(X[i, j]^2*w[j]/w[i] for i in axes(X, 1), j in axes(X, 2))
        @test actual <= proof.error^2
        a, b, restored = MP.decode(JSON3.read(JSON3.write(MP.encode(A, B, proof))))
        @test a == A && b == B
        for replay in (:rational, :interval)
            @test MP.verify_product(a, b, restored; residual_arithmetic = replay)
            bad = copy(proof.residual_bounds)
            bad[1] = 0
            @test_throws AssertionError MP.verify_product(
                A,
                B,
                merge(proof, (; residual_bounds = bad));
                residual_arithmetic = replay,
            )
            bad = copy(proof.operator_bounds)
            bad[1] = 0
            @test_throws AssertionError MP.verify_product(
                A,
                B,
                merge(proof, (; operator_bounds = bad));
                residual_arithmetic = replay,
            )
            corrupted = deepcopy(proof.cores)
            corrupted[end][1] += 1
            @test_throws AssertionError MP.verify_product(
                A,
                B,
                merge(proof, (; cores = corrupted));
                residual_arithmetic = replay,
            )
        end
    end
    # A zero product has zero residual bounds. Interval replay must take the
    # exact fallback successfully, while rejecting a changed interior core.
    A = [zeros(F, 1, 2, 2, 1) for _ in 1:3]
    B = fill(reshape(Matrix{F}(I, 2, 2), 1, 2, 2, 1), 3)
    proof = MP.product_candidate(A, B; candidate_arithmetic = :interval)
    @test proof.error == 0
    @test MP.verify_product(A, B, proof; residual_arithmetic = :interval)
    corrupted = deepcopy(proof.cores)
    corrupted[1] .*= 2
    @test_throws AssertionError MP.verify_product(
        A,
        B,
        merge(proof, (; cores = corrupted));
        residual_arithmetic = :interval,
    )
    one = MP.product_candidate(B[1:1], B[1:1]; candidate_arithmetic = :interval)
    @test one.error == 0 && MP.verify_product(B[1:1], B[1:1], one; residual_arithmetic = :interval)
end

@testset "Unequal-rank mixed MPO products" begin
    rng=MersenneTwister(6641)
    for n in (2, 3, 4), weighted in (false, true)
        dimsA=vcat(1, fill(2, n-1), 1)
        dimsB=vcat(1, fill(3, n-1), 1)
        A=[F.(rand(rng, -3:3, dimsA[j], 2, 2, dimsA[j+1]))/8 for j in 1:n]
        B=[F.(rand(rng, -3:3, dimsB[j], 2, 2, dimsB[j+1]))/8 for j in 1:n]
        weights=weighted ? F[1, 3] : nothing
        gauges=MN.rational_gauges(A, randn(rng, sum(dimsA[2:(end-1)])))
        if n == 2 && weighted
            gauges = [ones(F, 1, 1), F[1 1//4; -1//5 1], ones(F, 1, 1)]
        end
        proof=MP.product_candidate(A, B; maxdim = 2, weights, first_gauges = gauges)
        @test MP.verify_product(A, B, proof)
        X=expand_product_test(A)*expand_product_test(B)-expand_product_test(proof.cores)
        w=reduce(kron, fill(weighted ? weights : F[1, 1], n))
        actual=sum(X[i, j]^2*w[j]/w[i] for i in axes(X, 1), j in axes(X, 2))
        @test actual <= proof.error^2
        exact_product=expand_product_test(A)*expand_product_test(B)
        exact_candidate=expand_product_test(proof.cores)
        overlap=sum(
            exact_candidate[i, j]*exact_product[i, j]*w[j]/w[i] for i in axes(X, 1), j in axes(X, 2)
        )
        @test Nemo.contains(MP.product_inner_ball(proof.cores, A, B; weights), overlap)
        @test maximum(size(c, 4) for c in proof.cores)<=2
        a, b, p=MP.decode(JSON3.read(JSON3.write(MP.encode(A, B, proof))))
        @test a==A && b==B && MP.verify_product(a, b, p)
        bad=copy(proof.residual_bounds)
        bad[1]=0
        @test_throws AssertionError MP.verify_product(A, B, merge(proof, (; residual_bounds = bad)))
        bad_cores=deepcopy(proof.cores)
        bad_cores[end][1]+=1
        @test_throws AssertionError MP.verify_product(A, B, merge(proof, (; cores = bad_cores)))
    end
    identity=reshape(Matrix{F}(I, 2, 2), 1, 2, 2, 1)
    long=fill(identity, 256)
    @test Nemo.contains(MP.product_inner_ball(long, long, long; weights = F[1, 3]), big(2)^256)
end
