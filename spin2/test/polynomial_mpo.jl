using Spin2Gap, Test, LinearAlgebra, SparseArrays
import Nemo
const PM = Spin2Gap.PolynomialMPO
const F = PM.F

@testset "Exact periodic polynomial MPO" begin
    n = 4
    s = PM.local_spins()
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
    K = PM.centered_hamiltonian(n)
    @test maximum(size(c, 4) for c in K) == 8
    @test PM.sparse_operator(K) == reference
    K2 = PM.product_mpo(K, K)
    @test PM.sparse_operator(K2) == reference*reference
    seed = PM.PG.seed(n, 4, 12)
    Q = PM.polynomial_mpo(seed)
    @test PM.sparse_operator(Q) ==
          seed.identity*sparse(Matrix{F}(I, 625, 625))+seed.linear*reference+seed.quadratic*reference^2
    @test maximum(size(c, 4) for c in Q) == 73
    reduced = PM.exact_compress_right(Q)
    @test PM.verify_exact_compression(Q, reduced.cores, reduced.transforms)
    @test PM.sparse_operator(reduced.cores) == PM.sparse_operator(Q)
    reversed = PM.MN.reflect_cores(reduced.cores)
    twice = PM.exact_compress_right(reversed)
    @test maximum(size(c, 4) for c in twice.cores) == 33
    @test PM.sparse_operator(PM.MN.reflect_cores(twice.cores)) == PM.sparse_operator(Q)
    corrupted = deepcopy(reduced.transforms)
    corrupted[2][1, 1] += 1
    @test_throws AssertionError PM.verify_exact_compression(Q, reduced.cores, corrupted)
    @test PM.MS.hs_squared(K; weights = PM.WEIGHTS) == big(5)^n*(n^2+12n)
    @test Nemo.contains(PM.MS.hs_inner_ball(K, K; weights = PM.WEIGHTS), big(5)^n*(n^2+12n))
    q = PM.normalized_polynomial(seed)
    scalar_norm = PM.initial_operator_upper(seed)
    @test scalar_norm > (seed.identity+seed.linear*F(-16)+seed.quadratic*F(256))/25
    norm_bounds = PM.reference_operator_bounds(seed, F[1, 2, 3])
    @test norm_bounds[1] == scalar_norm
    @test all(norm_bounds[j+1]*F(j) >= norm_bounds[j]^2 for j in 1:3)
    altered = deepcopy(q)
    altered[1] *= F(1000001, 1000000)
    checked = PM.initial_error(q, altered)
    exact_norm2 = PM.MS.hs_squared(q; weights = PM.WEIGHTS)
    @test checked.error^2 >= exact_norm2/F(10)^12
    @test checked.error^2 < F(10001, 10000)*exact_norm2/F(10)^12
    tiny_ball = Nemo.ArbField(192)(F(3, big(10)^300))
    @test Nemo.ArbField(192)(PM.PG.relative_upper(tiny_ball)) > tiny_ball
    @test PM.PG.relative_upper(tiny_ball) < F(4, big(10)^300)
    global_weights = reduce(kron, fill(PM.WEIGHTS, n))
    physical =
        Float64.(reference) .* transpose(sqrt.(Float64.(global_weights))) ./
        sqrt.(Float64.(global_weights))
    @test physical ≈ transpose(physical) atol=1e-12
    @test collect(extrema(eigvals(Symmetric(Matrix(physical))))) ≈ [-16, 20] atol=1e-11
    long = PM.centered_hamiltonian(256)
    @test PM.MS.hs_squared(long; weights = PM.WEIGHTS) == big(5)^256*(256^2+12*256)
end

@testset "Quartic factors in the physical spin-2 basis" begin
    seed = PM.PG.taylor_seed(4, 4, 8, 4)
    saved = Dict(
        "length"=>4,
        "beta"=>"4//1",
        "power"=>8,
        "degree"=>4,
        "construction"=>"quartic-taylor-sos",
    )
    reconstructed = PM.scalar_seed(saved)
    @test reconstructed.coefficients == seed.coefficients
    @test PM.initial_operator_upper(reconstructed) > PM.PG.taylor_polynomial(seed, -16)/25
    @test_throws AssertionError PM.scalar_seed(merge(saved, Dict("degree"=>8)))
    legacy = PM.scalar_seed(Dict("length"=>4, "beta"=>"4//1", "power"=>8))
    @test legacy.remainder == PM.PG.seed(4, 4, 8).remainder
    sos = PM.quartic_sos_factors(seed)
    K = PM.sparse_operator(PM.centered_hamiltonian(4))
    identity = sparse(Matrix{F}(I, 625, 625))
    factors = [c.identity*identity+c.linear*K+c.quadratic*K^2 for c in sos.factors]
    weights = reduce(kron, fill(PM.WEIGHTS, 4))
    for P in factors
        @test P*Diagonal(weights) == Diagonal(weights)*transpose(P)
    end
    from_squares = sum(w*(P*P) for (w, P) in zip(sos.weights, factors))+seed.remainder*identity
    direct = sum(seed.coefficients[k+1]*K^k for k in 0:4)
    @test from_squares == direct

    # Check the normalization separately with scalar identity outputs, where
    # the expected assembly has no dependence on an MPO compression routine.
    sigma = F(1, 25)
    outputs = [PM.identity_mpo(4) for _ in 1:3]
    for (index, C) in enumerate(outputs)
        C[1] = index*sigma^2*C[1]
    end
    assembled = PM.quartic_sos_assembly(seed, outputs)
    expected = sigma*(sum(sos.weights .* (1:3))+seed.remainder)
    @test PM.sparse_operator(assembled) == expected*identity
end
