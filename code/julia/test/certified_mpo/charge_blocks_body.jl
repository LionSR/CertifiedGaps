using LinearAlgebra, Random, Test, JSON3
import Nemo

const CB = MP.MS.CB
const F = MP.F
const MN = MP.MN

function charge_test_mpo(rng, charges, physical)
    return [
        [
            a+s-t == b ? F(rand(rng, (-3, -2, -1, 1, 2, 3)), 7) : F(0) for
            a in charges[j], s in physical, t in physical, b in charges[j+1]
        ] for j in 1:(length(charges)-1)
    ]
end

function charge_test_expand(cores)
    partial = [ones(F, 1, 1)]
    for c in cores
        d = size(c, 2)
        next = [zeros(F, d*size(first(partial), 1), d*size(first(partial), 2)) for _ in axes(c, 4)]
        for a in axes(c, 4), b in axes(c, 1)
            next[a] += kron(partial[b], c[b, :, :, a])
        end
        partial = next
    end
    return only(partial)
end

@testset "Exact charge support and sector products" begin
    rng = MersenneTwister(91563)
    physical = [-1, 0, 1]
    qa, qb = [[0], [-1, 0, 0, 1], [0]], [[0], [-2, -1, 0, 1, 2], [0]]
    A, B = charge_test_mpo(rng, qa, physical), charge_test_mpo(rng, qb, physical)
    @test CB.infer_charges(A, physical) == qa
    @test CB.infer_charges(B, physical) == qb
    @test CB.infer_charges(MN.adjoint_cores(A), physical) == [-q for q in qa]
    @test CB.infer_charges(MN.reflect_cores(A), physical) == [-q for q in reverse(qa)]
    broken = deepcopy(A)
    broken[1][1, 1, 1, 1] = 1
    @test_throws ArgumentError CB.infer_charges(broken, physical)
    # Disconnected auxiliary support and zero MPOs do not justify silently
    # dropping entries or declaring an inconsistent labelling valid.
    zero_mpo = [zeros(F, 1, 3, 3, 2), zeros(F, 2, 3, 3, 1)]
    @test CB.infer_charges(zero_mpo, physical) == [[0], [0, 0], [0]]
    RR = Nemo.ArbField(192)
    # Both bridge contractions, unequal input ranks, nontrivial prefix charge.
    for site in 1:2
        prefix = site == 1 ? [0] : [-2, -1, 0, 0, 1, 2]
        bridge =
            [p == a+b ? F(rand(rng, -3:3), 11) : F(0) for p in prefix, a in qa[site], b in qb[site]]
        exact = MP.MS.product_bridge(bridge, A[site], B[site])
        blocked = CB.product_bridge(
            bridge,
            A[site],
            B[site],
            prefix,
            qa[site],
            qa[site+1],
            qb[site],
            qb[site+1],
            physical,
        )
        @test blocked == exact
        interval = CB.product_bridge(
            bridge,
            A[site],
            B[site],
            prefix,
            qa[site],
            qa[site+1],
            qb[site],
            qb[site+1],
            physical;
            field = RR,
        )
        @test all(Nemo.contains(x, Nemo.QQ(y)) for (x, y) in zip(interval, exact))
    end
    @test_throws ArgumentError CB.matrix_product(
        F[1 1; 0 1],
        Matrix{F}(I, 2, 2),
        [0, 1],
        [0, 1],
        [0, 1],
    )
    @test_throws ArgumentError CB.check_support(
        [RR(1) RR(1)//big(2)^1000; RR(0) RR(1)],
        [0, 1],
        [0, 1],
    )
end

@testset "Global rank allocation and degenerate sectors" begin
    labels = [0, 1, 0, 1]
    M = Matrix(Diagonal([4.0, 4.0, 2.0, 1.0]))
    for rank in 1:4
        block = CB.truncated_svd(M, labels, labels; maxdim = rank)
        @test length(block.S) == rank
        @test block.S ≈ svdvals(M)[1:rank]
        @test CB.check_support(block.U, labels, block.charges)
        @test CB.check_support(block.Vt, block.charges, labels)
        @test norm(M-block.U*Diagonal(block.S)*block.Vt)^2 ≈ sum(abs2, svdvals(M)[(rank+1):end]) atol =
            1e-14
    end
    empty = CB.truncated_svd(zeros(2, 3), [0, 0], [1, 1, 1]; maxdim = 2)
    @test empty.U*Diagonal(empty.S)*empty.Vt == zeros(2, 3)
end

@testset "Charge-preserving candidates, independent general replay" begin
    rng = MersenneTwister(30078)
    for physical in ([0, 1], [-1, 0, 1]), arithmetic in (:rational, :interval)
        d = length(physical)
        qa, qb = [[0], [-1, 0, 1], [-1, 0, 1], [0]], [[0], [0, 1], [0, 1], [0]]
        A, B = charge_test_mpo(rng, qa, physical), charge_test_mpo(rng, qb, physical)
        weights = F.(1:d)
        proof = MP.product_candidate(
            A,
            B;
            physical_charges = physical,
            maxdim = 3,
            weights,
            candidate_arithmetic = arithmetic,
        )
        @test maximum(size(c, 4) for c in proof.cores) <= 3
        @test length(CB.infer_charges(proof.cores, physical)) == 4
        @test MP.verify_product(A, B, proof)
        @test MP.verify_product(A, B, proof; residual_arithmetic = :interval)
        a, b, restored = MP.decode(JSON3.read(JSON3.write(MP.encode(A, B, proof))))
        @test MP.verify_product(a, b, restored)
        residual = charge_test_expand(A)*charge_test_expand(B)-charge_test_expand(proof.cores)
        metric = reduce(kron, fill(weights, 3))
        @test sum(
            residual[i, j]^2*metric[j]/metric[i] for i in axes(residual, 1), j in axes(residual, 2)
        ) <= proof.error^2
        corrupted = deepcopy(proof.cores)
        corrupted[end][1] += 1
        @test_throws AssertionError MP.verify_product(A, B, merge(proof, (; cores = corrupted)))
        for adjoint_first in (false, true), project_bridge in (false, true)
            square = MP.MS.square_candidate(
                A;
                physical_charges = physical,
                maxdim = 3,
                weights,
                candidate_arithmetic = arithmetic,
                adjoint_first,
                project_bridge,
                error_method = :projection,
            )
            @test length(CB.infer_charges(square.cores, physical)) == 4
            @test MP.MS.verify_square_certificate(A, square)
        end
    end
end

@testset "Charge-preserving Schmidt gauges" begin
    rng = MersenneTwister(973)
    physical = [-1, 0, 1]
    q = [[0], [-1, 0, 0, 1], [-1, 0, 0, 1], [0]]
    cores = charge_test_mpo(rng, q, physical)
    proposed = CB.schmidt_gauges([Float64.(c) for c in cores], physical)
    gauges = [MP.MS.relative_quantize(g; bits = 48) for g in proposed]
    @test all(CB.check_support(g, charges, charges) for (g, charges) in zip(gauges, q))
    transformed = MN.gauge_transform(cores, gauges)
    @test CB.infer_charges(transformed, physical) == q
    @test charge_test_expand(transformed) == charge_test_expand(cores)
    proof = MP.product_candidate(
        cores,
        cores;
        physical_charges = physical,
        first_gauges = gauges,
        maxdim = 3,
        candidate_arithmetic = :interval,
    )
    @test MP.verify_product(cores, cores, proof; residual_arithmetic = :interval)
    # Include sectors with no reachable physical prefix and a wide right nullspace.
    null = [zeros(F, 1, 2, 2, 5), zeros(F, 5, 2, 2, 1)]
    ng = CB.schmidt_gauges([Float64.(c) for c in null], [0, 1])
    @test all(isfinite, ng[2]) && rank(ng[2]) == 5
end
