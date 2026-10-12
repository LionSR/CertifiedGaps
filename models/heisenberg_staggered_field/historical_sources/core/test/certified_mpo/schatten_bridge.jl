include(joinpath(@__DIR__, "setup.jl"))
using Test, LinearAlgebra
const SB = CM.SchattenBridge
const F = SB.F

@testset "Operator-error Schatten bridge" begin
    A = F[1//2 1//8; 1//8 1//3]
    B = A+F[1//100 -1//200; 1//300 1//400]
    eta = CM.Norms.relative_upper_sqrt(sum(abs2, A-B))
    for p in (2, 4, 8, 16), t in (F(1, 100), F(1, 10), F(1))
        moment = tr((transpose(B)*B)^div(p, 2))
        bound = SB.trace_upper(moment, eta, 2, p, t)
        @test bound.upper >= tr(A^p)
        @test bound.upper == bound.moment_contribution+bound.floor_contribution
    end
    for p in (16, 32, 64, 128)
        root = big(5)^div(256, p)
        @test root^p == big(5)^256
        @test SB.normalized_upper(root, p, F(1, 10^8)) >= 1
    end
    @test SB.normalized_upper(big(5)^4, 64, F(1, 10^8)) < F(10005, 10000)
    @test_throws AssertionError SB.trace_upper(F(1), F(-1), 2, 4, F(1))
end

@testset "Backward Schatten Gram chain on explicit matrices" begin
    hs(X) = SB.rational_sqrt_upper(sum(abs2, X))
    for n in (2, 3), steps in 0:3
        A = Matrix{F}(I, n, n)/3+fill(F(1, 17), n, n)
        B = A+F.([i-j+1 for i in 1:n, j in 1:n])/1000
        initial = hs(A-B)
        kappas, errors, op_errors = F[], F[], F[]
        for j in 1:steps
            residual = F.([2i-j-l for i in 1:n, l in 1:n])/10000
            C = transpose(B)*B+residual
            kappa = F(j+2, 5)
            push!(kappas, kappa)
            push!(errors, hs(residual))
            push!(op_errors, hs(residual)) # HS dominates operator norm.
            B = reverse(C; dims = (1, 2))/kappa # Unitary spatial reversal.
        end
        actual = tr(A^(big(2)^(steps+1)))
        for op in (fill(nothing, steps), op_errors), initial_hs in (nothing, initial)
            result = SB.gram_chain_upper(
                hs(B),
                kappas,
                errors;
                dimension = n,
                operator_errors = op,
                initial_hs_error = initial_hs,
                initial_operator_error = initial,
            )
            @test result.initial_norm^result.power >= actual
            @test result.dimension_roots[1]^result.power >= n
            @test length(result.candidate_norms) == steps+1
        end
    end
    # Isotropic residuals have a sharp known operator bound. Their large-p
    # Schatten cost must be below the dimension^(1/2) HS cost.
    n, steps = 16, 3
    kappas = fill(F(1), steps)
    errors, op = fill(F(4, 1000), steps), fill(F(1, 1000), steps)
    pure_hs =
        SB.gram_chain_upper(F(4), kappas, errors; dimension = n, initial_hs_error = F(4, 1000))
    mixed = SB.gram_chain_upper(
        F(4),
        kappas,
        errors;
        dimension = n,
        operator_errors = op,
        initial_hs_error = F(4, 1000),
        initial_operator_error = F(1, 1000),
    )
    @test mixed.initial_norm < pure_hs.initial_norm
    @test mixed.local_errors[1] < errors[1]
    @test mixed.local_errors[end] == errors[end]
    # The dimension roots remain cheap even when the ultimate exponent is
    # hundreds of millions; no enormous exact integer power is constructed.
    long = SB.gram_chain_upper(F(1), fill(F(1), 28), fill(F(0), 28); dimension = big(5)^256)
    @test long.power == big(2)^29
    @test long.dimension_roots[1] < F(1000001, 1000000)
    @test long.initial_norm == 1
    @test_throws AssertionError SB.gram_chain_upper(F(1), [F(0)], [F(0)]; dimension = 2)
    @test_throws AssertionError SB.gram_chain_upper(F(1), [F(1)], [F(-1)]; dimension = 2)
    @test_throws AssertionError SB.gram_chain_upper(
        F(1),
        [F(1)],
        [F(0)];
        dimension = 2,
        operator_errors = [F(-1)],
    )
    @test_throws AssertionError SB.gram_chain_upper(
        F(1),
        F[],
        F[];
        dimension = 2,
        initial_hs_error = nothing,
    )
    @test_throws AssertionError SB.gram_chain_upper(
        F(1),
        F[],
        F[];
        dimension = 2,
        root_upper = x -> F(1),
    )
end
