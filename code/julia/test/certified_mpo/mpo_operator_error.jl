include(joinpath(@__DIR__, "setup.jl"))
using Test, LinearAlgebra, Random, JSON3
const OE = CM.OperatorError
const MS = OE.MS
const MN = OE.MN
const F = OE.F
const MP = CM.Products

function dense(cores)
    partial = [ones(F, 1, 1)]
    for core in cores
        next = [zeros(F, 2size(partial[1], 1), 2size(partial[1], 2)) for _ in axes(core, 4)]
        for a in axes(core, 4), b in axes(core, 1)
            next[a] += kron(partial[b], core[b, :, :, a])
        end
        partial = next
    end
    return only(partial)
end

@testset "One-site matrix residual closure" begin
    input = [reshape(Matrix{F}(I, 2, 2), 1, 2, 2, 1)]
    square = MS.square_candidate(input; gauges = [ones(F, 1, 1), ones(F, 1, 1)])
    unit = OE.unit_gauges(input)
    @test isempty(OE.residuals(input, square, unit, unit, unit))
    proof = OE.operator_error_certificate(input, square; arithmetic = :interval)
    @test proof.error == 0 &&
          OE.verify_operator_error(input, square, proof; residual_arithmetic = :interval)
end

function prefixes(cores)
    partial = [ones(F, 1, 1)]
    for core in cores
        next = [zeros(F, 2size(partial[1], 1), 2size(partial[1], 2)) for _ in axes(core, 4)]
        for a in axes(core, 4), b in axes(core, 1)
            next[a] += kron(partial[b], core[b, :, :, a])
        end
        partial = next
    end
    return partial
end

function tails(cores)
    partial = [ones(F, 1, 1)]
    for core in reverse(cores)
        next = [zeros(F, 2size(partial[1], 1), 2size(partial[1], 2)) for _ in axes(core, 1)]
        for a in axes(core, 4), b in axes(core, 1)
            next[b] += kron(core[b, :, :, a], partial[a])
        end
        partial = next
    end
    return partial
end

@testset "Full-matrix operator residual identity and interval replay" begin
    rng = MersenneTwister(74385)
    function matrices(cores)
        dims = vcat(1, [size(c, 4) for c in cores])
        return [
            d == 1 ? ones(F, 1, 1) : Matrix{F}(I, d, d)+F.(rand(rng, -2:2, d, d))/16 for d in dims
        ]
    end
    for n in (2, 3, 4), gram in (false, true)
        dims = n == 4 ? [1, 2, 3, 2, 1] : vcat(1, fill(2, n-1), 1)
        input = [F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1]))/8 for j in 1:n]
        weights = F[1, 3]
        square = MS.square_candidate(
            input;
            maxdim = 2,
            weights,
            adjoint_first = gram,
            gauges = matrices(input),
            candidate_arithmetic = :interval,
        )
        fg, sg, pg = matrices(input), matrices(input), matrices(square.cores)
        first = gram ? MN.adjoint_cores(input; weights) : input
        first_maps = MN.gauge_transform(first, fg)
        second_maps = MN.gauge_transform(input, sg)
        prefix_maps = MN.gauge_transform(square.cores, pg)
        E = OE.residuals(input, square, fg, sg, pg)
        reconstructed = zeros(F, 2^n, 2^n)
        for j in 1:(n-1)
            pre = prefixes(prefix_maps[1:(j-1)])
            a, b = tails(first_maps[(j+1):end]), tails(second_maps[(j+1):end])
            residual = reshape(E[j], length(pre), 2, 2, length(a), length(b))
            for q in eachindex(pre), i in eachindex(a), k in eachindex(b)
                reconstructed += kron(pre[q], kron(residual[q, :, :, i, k], a[i]*b[k]))
            end
        end
        @test reconstructed == dense(first)*dense(input)-dense(square.cores)
        w = reduce(kron, fill(weights, n))
        physical = Float64.(reconstructed) .* transpose(sqrt.(Float64.(w))) ./ sqrt.(Float64.(w))
        for arithmetic in (:rational, :interval)
            proof = OE.operator_error_certificate(
                input,
                square;
                first_gauges = fg,
                second_gauges = sg,
                prefix_gauges = pg,
                arithmetic,
            )
            @test opnorm(physical) <= Float64(proof.error)
            replay = OE.decode(JSON3.read(JSON3.write(OE.encode(proof))))
            @test replay.error == proof.error
            for residual_arithmetic in (:rational, :interval)
                @test OE.verify_operator_error(input, square, replay; residual_arithmetic)
                bad = copy(replay.residual_bounds)
                bad[1] = 0
                @test_throws AssertionError OE.verify_operator_error(
                    input,
                    square,
                    merge(replay, (; residual_bounds = bad));
                    residual_arithmetic,
                )
            end
        end
        bad = deepcopy(pg)
        bad[1][1] = 2
        @test_throws AssertionError OE.residuals(input, square, fg, sg, bad)
    end
end

@testset "Unequal operator-error tails and linear compression" begin
    rng = MersenneTwister(48210)
    for right_dimension in (1, 3)
        n = 3
        a_dims = [1, 2, 2, 1]
        b_dims = [1, right_dimension, right_dimension, 1]
        A = [F.(rand(rng, -3:3, a_dims[j], 2, 2, a_dims[j+1]))/8 for j in 1:n]
        B =
            right_dimension == 1 ? fill(reshape(Matrix{F}(I, 2, 2), 1, 2, 2, 1), n) :
            [F.(rand(rng, -3:3, b_dims[j], 2, 2, b_dims[j+1]))/8 for j in 1:n]
        gauges(c) = [
            Matrix{F}(I, d, d)+(d == 1 ? zeros(F, d, d) : F.(rand(rng, -2:2, d, d))/16) for
            d in vcat(1, [size(x, 4) for x in c])
        ]
        weights = F[1, 3]
        product = MP.product_candidate(A, B; maxdim = 2, weights, first_gauges = gauges(A))
        fg, sg, pg = gauges(A), gauges(B), gauges(product.cores)
        E = OE.residuals(A, product, fg, sg, pg; second_input = B)
        reconstructed = zeros(F, 2^n, 2^n)
        am, bm, cm = MN.gauge_transform(A, fg),
        MN.gauge_transform(B, sg),
        MN.gauge_transform(product.cores, pg)
        for j in 1:(n-1)
            pre, a, b = prefixes(cm[1:(j-1)]), tails(am[(j+1):end]), tails(bm[(j+1):end])
            e = reshape(E[j], length(pre), 2, 2, length(a), length(b))
            for q in eachindex(pre), i in eachindex(a), k in eachindex(b)
                reconstructed += kron(pre[q], kron(e[q, :, :, i, k], a[i]*b[k]))
            end
        end
        @test reconstructed == dense(A)*dense(B)-dense(product.cores)
        for arithmetic in (:rational, :interval)
            proof = OE.operator_error_certificate(
                A,
                product;
                second_input = B,
                first_gauges = fg,
                second_gauges = sg,
                prefix_gauges = pg,
                arithmetic,
            )
            decoded = OE.decode(JSON3.read(JSON3.write(OE.encode(proof))))
            for residual_arithmetic in (:rational, :interval)
                @test OE.verify_operator_error(
                    A,
                    product,
                    decoded;
                    second_input = B,
                    residual_arithmetic,
                )
            end
            w = reduce(kron, fill(weights, n))
            physical =
                Float64.(reconstructed) .* transpose(sqrt.(Float64.(w))) ./ sqrt.(Float64.(w))
            @test opnorm(physical) <= Float64(proof.error)
        end
    end
    A = fill(reshape(Matrix{F}(I, 2, 2), 1, 2, 2, 1), 2)
    B = fill(zeros(F, 1, 2, 2, 1), 2)
    product = MP.product_candidate(A, B)
    proof = OE.operator_error_certificate(A, product; second_input = B, arithmetic = :interval)
    @test proof.error == 0
    @test OE.verify_operator_error(
        A,
        product,
        proof;
        second_input = B,
        residual_arithmetic = :interval,
    )
end

@testset "Independent operator residual gauges" begin
    rng = MersenneTwister(2841)
    gauge(c) = MN.rational_gauges(c, randn(rng, sum(size(x, 4) for x in c[1:(end-1)])))
    for n in (2, 3, 4), gram in (false, true)
        dims = vcat(1, fill(2, n-1), 1)
        input = [F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1]))/8 for j in 1:n]
        weights = F[1, 3]
        square = MS.square_candidate(
            input;
            maxdim = 2,
            weights,
            adjoint_first = gram,
            gauges = gauge(input),
        )
        proof = OE.operator_error_certificate(
            input,
            square;
            first_gauges = gauge(input),
            second_gauges = gauge(input),
            prefix_gauges = gauge(square.cores),
        )
        @test OE.verify_operator_error(input, square, proof)
        @test OE.verify_operator_error(input, square, OE.decode(OE.encode(proof)))
        w = reduce(kron, fill(weights, n))
        physical(A) = Float64.(A) .* transpose(sqrt.(Float64.(w))) ./ sqrt.(Float64.(w))
        B, C = physical(dense(input)), physical(dense(square.cores))
        residual = (gram ? transpose(B)*B : B*B)-C
        @test opnorm(residual) <= Float64(proof.error)
        @test_throws AssertionError OE.verify_operator_error(
            input,
            square,
            merge(proof, (; error = F(0))),
        )
        bad = copy(proof.residual_bounds)
        bad[1] = 0
        @test_throws AssertionError OE.verify_operator_error(
            input,
            square,
            merge(proof, (; residual_bounds = bad)),
        )
        altered = deepcopy(square.cores)
        altered[end][1] += 1
        @test_throws AssertionError OE.verify_operator_error(
            input,
            merge(square, (; cores = altered)),
            proof,
        )
    end
end
