include(joinpath(@__DIR__, "setup.jl"))
using Test, LinearAlgebra, Random, JSON3
const MS = CM.Squaring
const MN = MS.MN
const F = MS.F

function tail_operators(cores)
    partial = [ones(F, 1, 1)]
    for core in reverse(cores)
        l, d, _, r = size(core)
        next = [zeros(F, d*size(first(partial), 1), d*size(first(partial), 2)) for _ in 1:l]
        for b in 1:l, a in 1:r
            next[b] += kron(core[b, :, :, a], partial[a])
        end
        partial = next
    end
    return partial
end

@testset "Separate operator and residual gauges" begin
    input = [zeros(F, 1, 2, 2, 2), zeros(F, 2, 2, 2, 1)]
    input[1][1, :, :, 1] = Matrix{F}(I, 2, 2)/2
    input[1][1, :, :, 2] = F[1 0; 0 -1]/2
    input[2][1, :, :, 1] = Matrix{F}(I, 2, 2)
    input[2][2, :, :, 1] = F[1 0; 0 -1]
    weights = F[1, 3]
    residual_gauges = [ones(F, 1, 1), F[1 0; 0 128], ones(F, 1, 1)]
    operator_gauges = [ones(F, 1, 1), Matrix{F}(I, 2, 2), ones(F, 1, 1)]
    for gram in (false, true)
        proof = MS.square_candidate(
            input;
            gauges = residual_gauges,
            weights,
            maxdim = 2,
            adjoint_first = gram,
            candidate_arithmetic = :interval,
        )
        first_copy = gram ? MN.adjoint_cores(input; weights) : input
        independent = MN.local_product_witness(first_copy, operator_gauges; weights)
        augmented = merge(
            proof,
            (; operator_certificate = (; gauges = operator_gauges, bounds = independent.bounds)),
        )
        previous = MS.verify_square_certificate(input, proof; details = true)
        for arithmetic in (:rational, :interval)
            checked = MS.verify_square_certificate(
                input,
                augmented;
                details = true,
                residual_arithmetic = arithmetic,
            )
            @test 1 <= checked.operator_norm_upper < previous.operator_norm_upper
            @test checked.operator_norm_upper ==
                  min(previous.operator_norm_upper, independent.upper)
            @test checked.triangle_error_upper == previous.triangle_error_upper
            decoded, replay =
                MS.decode_witness(JSON3.read(JSON3.write(MS.witness_encoding(input, augmented))))
            @test decoded == input && replay.cores == proof.cores
            @test MS.verify_square_certificate(
                decoded,
                replay;
                details = true,
                residual_arithmetic = arithmetic,
            ).operator_norm_upper == checked.operator_norm_upper
            bad = merge(
                augmented,
                (;
                    operator_certificate = (;
                        gauges = operator_gauges,
                        bounds = F[0, independent.bounds[2]],
                    )
                ),
            )
            @test_throws AssertionError MS.verify_square_certificate(
                input,
                bad;
                residual_arithmetic = arithmetic,
            )
            bad_gauges = deepcopy(operator_gauges)
            bad_gauges[1][1] = 2
            bad = merge(
                augmented,
                (; operator_certificate = (; gauges = bad_gauges, bounds = independent.bounds)),
            )
            @test_throws AssertionError MS.verify_square_certificate(
                input,
                bad;
                residual_arithmetic = arithmetic,
            )
        end
    end
end

@testset "Full-matrix physical MPO gauges" begin
    rng = MersenneTwister(42022)
    dims = [1, 2, 2, 1]
    cores = [F.(rand(rng, -2:2, dims[j], 2, 2, dims[j+1]))/8 for j in 1:3]
    gauges = [ones(F, 1, 1), F[1 1//3; 2//5 1], F[-2//3 1//7; 1//5 3//2], ones(F, 1, 1)]
    transformed = MN.gauge_transform(cores, gauges)
    enclosed = MN.gauge_transform(cores, gauges; interval_precision = 128)
    @test all(MS.Nemo.contains(x, y) for j in 1:2 for (x, y) in zip(enclosed[j], transformed[j]))
    @test enclosed[3] == transformed[3]
    @test MN.exact_gauge_core(cores[2], gauges[2], gauges[3]) == transformed[2]
    exact = only(tail_operators(cores))
    @test only(tail_operators(transformed)) == exact
    @test MN.decode_gauges(JSON3.read(JSON3.write(MN.encode_gauges(gauges)))) == gauges
    weights = F[1, 3]
    w = kron(weights, weights, weights)
    physical = sqrt.(Float64.(inv.(w))) .* Float64.(exact) .* transpose(sqrt.(Float64.(w)))
    normproof = MN.local_product_witness(cores, gauges; weights)
    @test MN.verify_local_product(cores, gauges, normproof.bounds; weights)
    @test Float64(normproof.upper) >= opnorm(physical)
    for gram in (false, true)
        proof = MS.square_candidate(cores; gauges, weights, maxdim = 2, adjoint_first = gram)
        @test MS.verify_square_certificate(cores, proof)
        @test MS.verify_square_certificate(cores, proof; residual_arithmetic = :interval)
        candidate = only(tail_operators(proof.cores))
        candidate2 =
            sum(candidate[i, j]^2*w[j]/w[i] for i in axes(candidate, 1), j in axes(candidate, 2))
        @test MS.output_hs_upper(proof)^2 >= candidate2
        first = gram ? Diagonal(w)*transpose(exact)*Diagonal(inv.(w)) : exact
        residual = first*exact-only(tail_operators(proof.cores))
        residual2 =
            sum(residual[i, j]^2*w[j]/w[i] for i in axes(residual, 1), j in axes(residual, 2))
        @test proof.error^2 >= residual2
        decoded, replay =
            MS.decode_witness(JSON3.read(JSON3.write(MS.witness_encoding(cores, proof))))
        @test decoded == cores && MS.verify_square_certificate(decoded, replay)
        scaled = MS.rescale_square_gauges(proof, [F[1], F[2, 3//2], F[4//5, 2], F[1]])
        refined = MS.square_candidate(
            cores;
            gauges = scaled.gauges,
            weights,
            adjoint_first = gram,
            fixed_compression = scaled,
            error_method = :projection,
            candidate_arithmetic = :interval,
        )
        @test refined.cores == proof.cores
        @test MS.verify_square_certificate(cores, refined; residual_arithmetic = :rational)
        @test MS.verify_square_certificate(cores, refined; residual_arithmetic = :interval)
    end
end

@testset "Interval residual replay, including tiny and zero operators" begin
    core = reshape(F[1 1//3; -1//5 2], 1, 2, 2, 1)
    for scale in F[0, 1, F(1, big(10)^200)]
        input = [scale*core, scale*core]
        proof = MS.square_candidate(input; maxdim = 1)
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :interval)
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :rational)
    end
    rng = MersenneTwister(3470)
    input = [F.(rand(rng, -3:3, 1, 2, 2, 2))/8, F.(rand(rng, -3:3, 2, 2, 2, 1))/8]
    proof = MS.square_candidate(input; maxdim = 1, weights = F[1, 3], adjoint_first = true)
    @test MS.verify_square_certificate(input, proof; residual_arithmetic = :interval)
    bad_rows = deepcopy(proof.rows)
    bad_rows[1] = merge(bad_rows[1], (; residual_bound = F(0)))
    @test_throws AssertionError MS.verify_square_certificate(
        input,
        merge(proof, (; rows = bad_rows));
        residual_arithmetic = :interval,
    )
end

@testset "Interval candidate generation with exact residual replay" begin
    rng = MersenneTwister(72137)
    dims = [1, 2, 2, 1]
    cores = [F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1]))/8 for j in 1:3]
    weights = F[1, 3]
    gauges = [ones(F, 1, 1), F[1 1//7; 2//5 1], F[-2//3 1//7; 1//5 3//2], ones(F, 1, 1)]
    for scale in F[0, 1, F(1, big(10)^200)], gram in (false, true)
        input = [scale*cores[1], cores[2:end]...]
        proof = MS.square_candidate(
            input;
            maxdim = 2,
            weights,
            gauges,
            adjoint_first = gram,
            candidate_arithmetic = :interval,
        )
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :rational)
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :interval)
        exact = only(tail_operators(input))
        w = kron(weights, weights, weights)
        first = gram ? Diagonal(w)*transpose(exact)*Diagonal(inv.(w)) : exact
        residual = first*exact-only(tail_operators(proof.cores))
        residual2 =
            sum(residual[i, j]^2*w[j]/w[i] for i in axes(residual, 1), j in axes(residual, 2))
        @test proof.error^2 >= residual2
    end
end

@testset "Projection-defect error recurrence" begin
    rows =
        [(; residual_bound = F(1), local_prefix_bound = F(1), projection_bound = F(0)) for _ in 1:4]
    @test 2 <= MS.projection_error_upper(rows, ones(F, 6)) < 2+F(1, big(10)^20)
    parallel = [merge(row, (; projection_bound = F(1))) for row in rows]
    @test MS.projection_error_upper(parallel, ones(F, 6)) == 4
    rng = MersenneTwister(23851)
    dims = [1, 2, 3, 2, 1]
    input = [F.(rand(rng, -4:4, dims[j], 2, 2, dims[j+1]))/8 for j in 1:4]
    weights = F[1, 3]
    w = kron(weights, weights, weights, weights)
    exact = only(tail_operators(input))
    for gram in (false, true), arithmetic in (:rational, :interval)
        proof = MS.square_candidate(
            input;
            weights,
            maxdim = 2,
            bits = 12,
            adjoint_first = gram,
            candidate_arithmetic = arithmetic,
            error_method = :projection,
        )
        checked = MS.verify_square_certificate(input, proof; details = true)
        @test checked.valid
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :interval)
        @test proof.error <= checked.triangle_error_upper
        first = gram ? Diagonal(w)*transpose(exact)*Diagonal(inv.(w)) : exact
        residual = first*exact-only(tail_operators(proof.cores))
        @test proof.error^2 >=
              sum(residual[i, j]^2*w[j]/w[i] for i in axes(residual, 1), j in axes(residual, 2))
        decoded, replay =
            MS.decode_witness(JSON3.read(JSON3.write(MS.witness_encoding(input, proof))))
        @test decoded == input && MS.verify_square_certificate(decoded, replay)
        frozen = MS.square_candidate(
            input;
            weights,
            adjoint_first = gram,
            fixed_compression = proof,
            candidate_arithmetic = arithmetic,
            error_method = :projection,
        )
        @test frozen.cores == proof.cores && frozen.bridges == proof.bridges
        @test MS.verify_square_certificate(input, frozen)
        bad_rows = deepcopy(proof.rows)
        bad_rows[1] = merge(bad_rows[1], (; projection_bound = F(0)))
        @test_throws AssertionError MS.verify_square_certificate(
            input,
            merge(proof, (; rows = bad_rows));
            residual_arithmetic = :interval,
        )
    end
end

@testset "Fitted bridges remain independently certifiable" begin
    rng = MersenneTwister(39254)
    dims = [1, 2, 3, 2, 1]
    weights = F[1, 3]
    w = kron(weights, weights, weights, weights)
    input = [F.(rand(rng, -4:4, dims[j], 2, 2, dims[j+1]))/8 for j in 1:4]
    exact = only(tail_operators(input))
    product = Diagonal(w)*transpose(exact)*Diagonal(inv.(w))*exact
    for bits in (8, 64)
        proof = MS.square_candidate(
            input;
            weights,
            maxdim = 3,
            bits,
            adjoint_first = true,
            candidate_arithmetic = :interval,
            error_method = :projection,
            project_bridge = true,
        )
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :rational)
        @test MS.verify_square_certificate(input, proof; residual_arithmetic = :interval)
        residual = product-only(tail_operators(proof.cores))
        @test proof.error^2 >=
              sum(residual[i, j]^2*w[j]/w[i] for i in axes(residual, 1), j in axes(residual, 2))
    end
end

@testset "Exact compact reflected candidates" begin
    core = reshape(F[1 2; -3 4], 1, 2, 2, 1)
    output = [core, 3core]
    data = Dict(
        "output"=>MN.encode_cores(output),
        "normalization"=>"3//2",
        "reflected_output"=>true,
        "next_candidate_encoding"=>"reflected-output",
    )
    expected = MN.reflect_cores(output)
    expected[1] *= F(2, 3)
    @test MS.flow_next_candidate(data) == expected
    @test MS.flow_next_candidate(data; output_cores = output) == expected
    data["reflected_output"] = false
    @test_throws AssertionError MS.flow_next_candidate(data)
    data["next_candidate"] = MN.encode_cores(expected)
    @test MS.flow_next_candidate(data) == expected
end

@testset "Certified MPO squaring residuals" begin
    rng = MersenneTwister(901)
    for dims in ([1, 2, 1], [1, 2, 3, 1], [1, 2, 2, 2, 1])
        cores = [F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1])) ./ 8 for j in 1:(length(dims)-1)]
        tails = MN.squared_tail_witness(cores)
        for j in 1:length(cores)
            R = tail_operators(cores[j:end])
            multiplication = hcat([vec(a*b) for a in R for b in R]...)
            @test Float64(tails.upper[j]) >= opnorm(Float64.(multiplication))-1e-12
        end
        B = only(tail_operators(cores))
        blocked_cores = vcat([MN.merge_cores(cores[1], cores[2])], cores[3:end])
        @test only(tail_operators(blocked_cores)) == B
        @test MS.hs_squared(cores) == sum(abs2, B)
        for chi in (1, 4)
            certificate = MS.square_candidate(cores; maxdim = chi)
            @test MS.verify_square_certificate(cores, certificate)
            C = only(tail_operators(certificate.cores))
            exact_error2 = sum(abs2, B*B-C)
            @test certificate.error^2 >= exact_error2
            @test maximum(size(c, 4) for c in certificate.cores) <= chi
        end
        gauges = MN.rational_gauges(cores, randn(rng, sum(size(c, 4) for c in cores[1:(end-1)])))
        certificate = MS.square_candidate(cores; maxdim = 2, gauges)
        @test MS.verify_square_certificate(cores, certificate)
        C = only(tail_operators(certificate.cores))
        @test certificate.error^2 >= sum(abs2, B*B-C)
        bad_rows = deepcopy(certificate.rows)
        bad_rows[1] = merge(bad_rows[1], (; residual_bound = F(0)))
        @test_throws AssertionError MS.verify_square_certificate(
            cores,
            merge(certificate, (; rows = bad_rows)),
        )
        bad_cores = deepcopy(certificate.cores)
        bad_cores[end][1] += 1
        @test_throws AssertionError MS.verify_square_certificate(
            cores,
            merge(certificate, (; cores = bad_cores)),
        )
        blocked = MS.square_candidate(cores; maxdim = 2, gauges, blocksize = 2)
        @test MS.verify_square_certificate(cores, blocked)
        @test blocked.error <= certificate.error
        @test blocked.error^2 >= sum(abs2, B*B-only(tail_operators(blocked.cores)))
        weights = F[1, 3]
        global_weights = reduce(kron, fill(weights, length(cores)))
        norm2(X) = sum(
            abs2(X[i, j])*global_weights[j]/global_weights[i] for i in axes(X, 1), j in axes(X, 2)
        )
        weighted = MS.square_candidate(cores; maxdim = 2, gauges, blocksize = 2, weights)
        @test MS.verify_square_certificate(cores, weighted)
        @test MS.hs_squared(cores; weights) == norm2(B)
        @test weighted.error^2 >= norm2(B*B-only(tail_operators(weighted.cores)))
        adjointB = Diagonal(global_weights)*transpose(B)*Diagonal(inv.(global_weights))
        @test only(tail_operators(MN.adjoint_cores(cores; weights))) == adjointB
        @test MN.reflect_cores(MN.reflect_cores(cores)) == cores
        @test MS.hs_squared(MN.reflect_cores(cores); weights) == norm2(B)
        gram = MS.square_candidate(
            cores;
            maxdim = 2,
            gauges,
            blocksize = 2,
            weights,
            adjoint_first = true,
        )
        @test MS.verify_square_certificate(cores, gram)
        @test gram.error^2 >= norm2(adjointB*B-only(tail_operators(gram.cores)))
        decoded, proof = MS.decode_witness(MS.witness_encoding(cores, gram))
        @test decoded == cores && MS.verify_square_certificate(decoded, proof)
    end
    for x in (2.0^-300, 1.0, 2.0^100)
        A = MS.relative_quantize(fill(x, 2, 2))
        @test A == fill(F(x), 2, 2)
    end
    tiny = F(3, big(10)^600)
    bound = MN.relative_upper_sqrt(tiny)
    @test bound^2 >= tiny
    @test bound^2 < tiny*F(1000001, 1000000)
end

@testset "Backward Gram norm bounds" begin
    B = F[1//3 1//10; -1//7 1//5]
    matrices = [B]
    errors = F[]
    kappas = F[1, 2, 3]
    for kappa in kappas
        next = transpose(last(matrices))*last(matrices)/kappa + F[1//100 0; 0 -1//200]
        residual = transpose(last(matrices))*last(matrices)-kappa*next
        push!(errors, MN.relative_upper_sqrt(sum(abs2, residual)))
        push!(matrices, next)
    end
    terminal = MN.relative_upper_sqrt(sum(abs2, last(matrices)))
    bounds = MS.backward_operator_bounds(kappas, errors, terminal)
    @test all(Float64(b) >= opnorm(Float64.(A)) for (b, A) in zip(bounds, matrices))
    propagated = MS.propagate_gram_errors(F(0), kappas, errors, bounds)
    exact = B
    for j in eachindex(kappas)
        exact = transpose(exact)*exact/kappas[j]
        @test propagated[j+1]^2 >= sum(abs2, exact-matrices[j+1])
    end
end

@testset "Rational physical metric" begin
    W = F[0 1; 1 0]
    metric = F[1, 4]
    @test !MN.verify_local_bound(W, F(3, 2); row_metric = metric, col_metric = metric)
    @test MN.verify_local_bound(W, F(201, 100); row_metric = metric, col_metric = metric)
    @test F(2) < MN.local_norm_upper(W; row_metric = metric, col_metric = metric) < F(201, 100)
    for W in (F[1 2 3; 4 5 6], F[1 2; 3 4; 5 6])
        rows, cols = F.(1:size(W, 1)), F.(2:(size(W, 2)+1))
        physical = sqrt.(Float64.(rows)) .* Float64.(W) ./ transpose(sqrt.(Float64.(cols)))
        bound = MN.local_norm_upper(W; row_metric = rows, col_metric = cols)
        @test Float64(bound) >= opnorm(physical)
        @test !MN.verify_local_bound(W, bound*999//1000; row_metric = rows, col_metric = cols)
    end
end

@testset "Bounded-precision rational error propagation" begin
    for x in F[0, 1, 1//3, 3//7, big(3)//big(10)^600, big(10)^600//3]
        upper = MN.relative_rational_upper(x)
        @test upper >= x
        @test upper <= x*(1+F(1, big(2)^127))
    end
    @test_throws AssertionError MN.relative_rational_upper(F(-1))
    steps = 1000
    kappas, deltas, bounds = ones(F, steps), fill(F(1, big(10)^12), steps), fill(F(1, 2), steps+1)
    propagated = MS.propagate_gram_errors(F(1, big(10)^6), kappas, deltas, bounds; round_bits = 128)
    # Every rounded step remains an upper solution of the exact recurrence.
    @test all(propagated[j+1] >= deltas[j]+propagated[j]*(1+propagated[j]) for j in 1:steps)
    @test maximum(ndigits(denominator(x); base = 2) for x in propagated) <= 150
    @test last(propagated) < F(2, big(10)^6)
    short = MS.propagate_gram_errors(F(1, big(10)^6), kappas[1:4], deltas[1:4], bounds[1:5])
    @test all(propagated[1:5] .>= short)
end

@testset "Independent exact and candidate norms in Gram errors" begin
    A = Matrix{F}(I, 4, 4)/4
    B = A+Diagonal(F[1//1000, -1//1000, 1//1000, -1//1000])
    epsilon = F(1, 500)
    a, b, kappa = F(1, 4), F(251, 1000), F(1, 4)
    @test MS.gram_error_coefficient(epsilon, b, a) == a+b
    bound = MS.propagate_gram_errors(epsilon, F[kappa], F[0], F[b, 1]; reference_bounds = F[a, 1])
    @test last(bound)^2 >= sum(abs2, (transpose(A)*A-transpose(B)*B)/kappa)
    # A conditional 256-site normalization example: independent tiny norm
    # bounds prevent the artificial epsilon^2/kappa contribution.
    tiny = F(1, big(5)^128)
    e = F(1, big(10)^12)
    direct = MS.propagate_gram_errors(e, F[tiny], F[0], F[2tiny, 1]; reference_bounds = F[tiny, 1])
    old = MS.propagate_gram_errors(e, F[tiny], F[0], F[tiny, 1])
    @test last(direct) == 3e
    @test last(old) > big(10)^60
    @test_throws AssertionError MS.gram_error_coefficient(e, tiny, F(-1))
end
