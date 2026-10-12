include(joinpath(@__DIR__, "setup.jl"))
using Test, LinearAlgebra, SparseArrays
import Nemo
const MA = CM.MPOAlgebra
const MS = MA.MS
const MN = MA.MN
const PG = CM.PolynomialBounds
const F = MA.F
const WEIGHTS = MA.WEIGHTS

# The spin-2 constructions stay in spin2 (`Spin2PolynomialMPO`). These test-local copies,
# verbatim apart from calling `MA.`, exercise the generic MPO algebra on the same operators.
function local_spins()
    plus, minus = zeros(F, 5, 5), zeros(F, 5, 5)
    for k in 0:3
        plus[k+2, k+1] = 4-k
        minus[k+1, k+2] = k+1
    end
    return (; identity = Matrix{F}(I, 5, 5), z = Matrix(Diagonal(F.(-2:2))), plus, minus)
end

function centered_hamiltonian(n)
    @assert iseven(n) && n >= 4
    s = local_spins()
    W = zeros(F, 5, 5, 5, 5)
    W[1, :, :, 1] = W[5, :, :, 5] = s.identity
    W[1, :, :, 5] = s.identity # Sum one identity per site, giving N I.
    W[1, :, :, 2], W[2, :, :, 5] = s.z, s.z
    W[1, :, :, 3], W[3, :, :, 5] = s.plus, -s.minus/2
    W[1, :, :, 4], W[4, :, :, 5] = s.minus, -s.plus/2
    open_chain = vcat([W[1:1, :, :, :]], fill(W, n-2), [W[:, :, :, 5:5]])
    wraps = Vector{Array{F,4}}[]
    for (a, b) in ((s.z, s.z), (-s.plus/2, s.minus), (-s.minus/2, s.plus))
        push!(
            wraps,
            vcat([reshape(a, 1, 5, 5, 1)], MA.identity_mpo(n-2), [reshape(b, 1, 5, 5, 1)]),
        )
    end
    return MA.linear_combination(vcat([open_chain], wraps), ones(F, 4))
end

function polynomial_mpo(seed)
    K = centered_hamiltonian(seed.n)
    return MA.linear_combination(
        [MA.identity_mpo(seed.n), K, MA.product_mpo(K, K)],
        [seed.identity, seed.linear, seed.quadratic],
    )
end

function normalized_polynomial(seed)
    Q = polynomial_mpo(seed)
    Q[1] = Q[1]/big(5)^div(seed.n, 2)
    return Q
end

function initial_operator_upper(seed)
    radius = F(5seed.n)
    return (seed.identity+abs(seed.linear)*radius+seed.quadratic*radius^2)/big(5)^div(seed.n, 2)
end

function reference_operator_bounds(seed, kappas; precision = 192)
    RR = Nemo.ArbField(precision)
    bounds = F[initial_operator_upper(seed)]
    for kappa in kappas
        @assert kappa > 0
        next = PG.relative_upper(RR(last(bounds))^2/RR(kappa))
        @assert next*kappa >= last(bounds)^2
        push!(bounds, next)
    end
    return bounds
end

@testset "Exact periodic polynomial MPO" begin
    n = 4
    s = local_spins()
    product(operators) = reduce(kron, sparse.(operators))
    reference = n*product(fill(s.identity, n))
    for j in 1:n
        k = mod1(j+1, n)
        for (a, b, c) in
            ((s.z, s.z, F(1)), (s.plus, s.minus, F(-1, 2)), (s.minus, s.plus, F(-1, 2)))
            terms = fill(s.identity, n)
            terms[j], terms[k] = a, b
            reference += c*product(terms)
        end
    end
    K = centered_hamiltonian(n)
    @test maximum(size(c, 4) for c in K) == 8
    @test MA.sparse_operator(K) == reference
    K2 = MA.product_mpo(K, K)
    @test MA.sparse_operator(K2) == reference*reference
    seed = PG.seed(n, 4, 12)
    Q = polynomial_mpo(seed)
    @test MA.sparse_operator(Q) ==
          seed.identity*sparse(Matrix{F}(I, 625, 625))+seed.linear*reference+seed.quadratic*reference^2
    @test maximum(size(c, 4) for c in Q) == 73
    reduced = MA.exact_compress_right(Q)
    @test MA.verify_exact_compression(Q, reduced.cores, reduced.transforms)
    @test MA.sparse_operator(reduced.cores) == MA.sparse_operator(Q)
    reversed = MN.reflect_cores(reduced.cores)
    twice = MA.exact_compress_right(reversed)
    @test maximum(size(c, 4) for c in twice.cores) == 33
    @test MA.sparse_operator(MN.reflect_cores(twice.cores)) == MA.sparse_operator(Q)
    corrupted = deepcopy(reduced.transforms)
    corrupted[2][1, 1] += 1
    @test_throws AssertionError MA.verify_exact_compression(Q, reduced.cores, corrupted)
    @test MS.hs_squared(K; weights = WEIGHTS) == big(5)^n*(n^2+12n)
    @test Nemo.contains(MS.hs_inner_ball(K, K; weights = WEIGHTS), big(5)^n*(n^2+12n))
    q = normalized_polynomial(seed)
    scalar_norm = initial_operator_upper(seed)
    @test scalar_norm > (seed.identity+seed.linear*F(-16)+seed.quadratic*F(256))/25
    norm_bounds = reference_operator_bounds(seed, F[1, 2, 3])
    @test norm_bounds[1] == scalar_norm
    @test all(norm_bounds[j+1]*F(j) >= norm_bounds[j]^2 for j in 1:3)
    altered = deepcopy(q)
    altered[1] *= F(1000001, 1000000)
    checked = MA.initial_error(q, altered)
    exact_norm2 = MS.hs_squared(q; weights = WEIGHTS)
    @test checked.error^2 >= exact_norm2/F(10)^12
    @test checked.error^2 < F(10001, 10000)*exact_norm2/F(10)^12
    tiny_ball = Nemo.ArbField(192)(F(3, big(10)^300))
    @test Nemo.ArbField(192)(PG.relative_upper(tiny_ball)) > tiny_ball
    @test PG.relative_upper(tiny_ball) < F(4, big(10)^300)
    global_weights = reduce(kron, fill(WEIGHTS, n))
    physical =
        Float64.(reference) .* transpose(sqrt.(Float64.(global_weights))) ./
        sqrt.(Float64.(global_weights))
    @test physical ≈ transpose(physical) atol=1e-12
    @test collect(extrema(eigvals(Symmetric(Matrix(physical))))) ≈ [-16, 20] atol=1e-11
    long = centered_hamiltonian(256)
    @test MS.hs_squared(long; weights = WEIGHTS) == big(5)^256*(256^2+12*256)
end

# Not in spin2: the d, weights and radius arguments added in the core.
@testset "MPO algebra in other local dimensions" begin
    @test MA.identity_mpo(3) == fill(reshape(Matrix{F}(I, 5, 5), 1, 5, 5, 1), 3)
    Z = Matrix(Diagonal(F[1, -1]))
    X = F[0 1; 1 0]
    unit = reshape(Matrix{F}(I, 2, 2), 1, 2, 2, 1)
    ising = [reshape(Z, 1, 2, 2, 1), reshape(Z, 1, 2, 2, 1), unit]
    field = [reshape(X, 1, 2, 2, 1), unit, unit]
    H = MA.linear_combination([ising, field, MA.identity_mpo(3; d = 2)], F[1, F(1, 3), 2])
    expected = kron(Z, Z, I(2))+kron(X, I(2), I(2))/3+2I(8)
    @test MA.sparse_operator(H) == expected
    @test MA.linear_combination([ising, field], F[1, 1]; d = 2) ==
          MA.linear_combination([ising, field], F[1, 1])
    weights = F[1, 2]
    reduced = MA.exact_compress_right(H)
    @test maximum(size(c, 4) for c in reduced.cores) == 2
    conditioned = MA.condition_right_exact(reduced.cores; weights, precision = 256, bits = 40)
    @test MA.sparse_operator(conditioned.cores) == expected
    @test MA.decode_transforms(MA.encode_transforms(conditioned.transforms)) ==
          conditioned.transforms
    perturbed = deepcopy(H)
    perturbed[1] = perturbed[1]*F(1001, 1000)
    exact_error2 = MS.hs_squared(H; weights)/F(1000)^2
    checked = MA.initial_error(H, perturbed; weights)
    @test exact_error2 <= checked.error^2 < F(10001, 10000)*exact_error2
    default = PG.taylor_seed(4, 4, 8, 4)
    @test PG.taylor_seed(4, 4, 8, 4; radius = 20) == default
    quadratic, explicit = PG.seed(4, 4, 8), PG.seed(4, 4, 8; radius = 20)
    # `log_inflation` is an Arb ball, and balls never compare equal; every rational must.
    @test all(explicit[k] == quadratic[k] for k in keys(quadratic) if k != :log_inflation)
    smaller = PG.taylor_seed(4, 4, 8, 4; radius = 4)
    @test smaller.spectral_radius == 4 && smaller.remainder < default.remainder
end
