include(joinpath(@__DIR__, "../scripts/repeated_l2_probe.jl"))
using Test, Random

@testset "Weighted Hilbert--Schmidt and operator-norm witnesses" begin
    rng = MersenneTwister(801)
    states = SECTORS[15]
    S = rand(rng, -3:3, length(states), length(states))
    S += transpose(S)
    M = zeros(F, 625, 625)
    for (i, r) in enumerate(states), (j, c) in enumerate(states)
        M[location(r, c)...] = weight(r)*S[i, j]//100
    end
    @test hermitize(M) == M
    physical(A) =
        [Float64(A[location(r, c)...])*sqrt(weight(c)/weight(r)) for r in states, c in states]
    P = physical(M)
    @test P ≈ transpose(P) atol=1e-14
    @test Float64(norm2(M)) ≈ sum(abs2, P) atol=1e-12
    bound = operator_norm_upper(M)
    witness = operator_norm_upper(M; details = true)
    @test witness.upper == bound
    @test verify_operator_norm_witness(M, witness)
    @test_throws AssertionError verify_operator_norm_witness(M, merge(witness, (; upper = F(0))))
    @test Float64(bound) >= opnorm(Symmetric(P))-1e-12
    @test sqrt_upper(norm2(M))^2 >= norm2(M)
    C, _ = svd_candidate(M, 2)
    full, fullrank = svd_candidate(M, 85)
    @test full == M
    @test fullrank == 625
    _, smallrank = svd_candidate(M, 2)
    @test smallrank == 32
    @test hermitize(C) == C
    @test Float64(norm2(M-C)) ≈ sum(abs2, physical(M-C)) atol=1e-12
    lo, hi = squared_interval(M)
    @test hermitize(lo) == lo && hermitize(hi) == hi
    exact = Nemo.matrix(Nemo.QQ, [M[location(r, c)...] for r in states, c in states])^2
    @test all(
        lo[location(r, c)...] <= exact[i, j] <= hi[location(r, c)...] for
        (i, r) in enumerate(states), (j, c) in enumerate(states)
    )
    # Check the perturbation inequality independently on full physical matrices.
    error = sqrt_upper(norm2(M-C))
    propagation = error*(2operator_norm_upper(C)+error)
    @test opnorm(P*P-physical(C)^2, 2) <= Float64(propagation)+1e-12
    @test norm(P*P-physical(C)^2) <= Float64(propagation)+1e-12
end
