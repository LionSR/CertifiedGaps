include(joinpath(@__DIR__, "setup.jl"))
using LinearAlgebra, Random, Test
const AT = CM.AdjointTrace
const MN = CM.Norms
const F = AT.F

@testset "Adjoint trace telescoping certificate" begin
    rng = MersenneTwister(5112)
    random_symmetric(d) = begin
        x=F.(rand(rng, -3:3, d, d))/100
        (x+transpose(x))/2
    end
    for d in (2, 3, 4), steps in (1, 3)
        kappas=F.(2:(steps+1))
        A=[Matrix{F}(I, d, d)/2+random_symmetric(d)]
        B=[A[1]+random_symmetric(d)]
        for kappa in kappas
            push!(A, last(A)^2/kappa)
            push!(B, last(B)^2/kappa+random_symmetric(d))
        end
        G=Vector{Matrix{F}}(undef, steps+1)
        G[end]=B[end]
        defects=F[]
        for j in steps:-1:1
            G[j]=(B[j]*G[j+1]+G[j+1]*B[j])/kappas[j]+random_symmetric(d)
        end
        errors=[MN.relative_upper_sqrt(sum(abs2, a-b)) for (a, b) in zip(A, B)]
        tau=[
            MN.relative_upper_sqrt(sum(abs2, (B[j]*G[j+1]+G[j+1]*B[j])/kappas[j]-G[j])) for
            j in 1:steps
        ]
        gamma=[MN.relative_upper_sqrt(sum(abs2, g)) for g in G[2:end]]
        center=sum(abs2, B[end]) +
               2tr(G[1]*(A[1]-B[1])) +
               2sum(tr(G[j+1]*(B[j]^2/kappas[j]-B[j+1])) for j in 1:steps)
        bound=AT.trace_interval(
            center,
            center;
            errors,
            dual_defects = tau,
            dual_operator_bounds = gamma,
            normalizations = kappas,
        )
        exact=sum(abs2, A[end])
        @test bound.lower <= exact <= bound.upper
        @test bound.cross >= 0 && bound.quadratic >= 0
        # Independent exact telescoping identity, including both defect terms.
        expanded=center+sum(abs2, A[end]-B[end])+2sum(
            tr(((B[j]*G[j+1]+G[j+1]*B[j])/kappas[j]-G[j])*(A[j]-B[j])) +
            tr(G[j+1]*(A[j]-B[j])^2)/kappas[j] for j in 1:steps
        )
        @test expanded == exact
        @test_throws AssertionError AT.trace_interval(
            center,
            center;
            errors = errors[1:(end-1)],
            dual_defects = tau,
            dual_operator_bounds = gamma,
            normalizations = kappas,
        )
    end
    C=F[1//2 1//10; 1//10 1//3]
    A=F[1//2 1//8; 1//8 1//3]
    eta=MN.relative_upper_sqrt(sum(abs2, A-C))
    center=2tr(C*A)-sum(abs2, C)
    one_step=AT.trace_interval(
        center,
        center;
        errors = F[0, eta],
        dual_defects = F[0],
        dual_operator_bounds = F[1],
        normalizations = F[1],
    )
    @test one_step.upper == center+eta^2
    @test one_step.lower <= sum(abs2, A) <= one_step.upper
end

@testset "Non-Hermitian reflected Gram flow" begin
    rng=MersenneTwister(9711)
    for d in (2, 3, 4)
        steps=3
        random_matrix() = F.(rand(rng, -3:3, d, d))/50
        P=reverse(Matrix{F}(I, d, d); dims = 1)
        reflect(X, j) = isodd(j) ? P*X*transpose(P) : X
        A=[random_matrix()]
        B=[A[1]+random_matrix()]
        kappas=F[2, 3, 4]
        for j in 1:steps
            push!(A, reflect(transpose(A[end])*A[end], j)/kappas[j])
            push!(B, reflect(transpose(B[end])*B[end], j)/kappas[j]+random_matrix())
        end
        G=Vector{Matrix{F}}(undef, steps+1)
        G[end]=B[end]
        T=Vector{Matrix{F}}(undef, steps)
        R=Vector{Matrix{F}}(undef, steps)
        for j in steps:-1:1
            reflected=reflect(G[j+1], j)
            mapped=B[j]*(reflected+transpose(reflected))/kappas[j]
            G[j]=mapped+random_matrix()
            T[j]=mapped-G[j]
            R[j]=reflect(transpose(B[j])*B[j], j)/kappas[j]-B[j+1]
        end
        pair(X, Y) = sum(X .* Y)
        D=A .- B
        center=sum(abs2, B[end])+2pair(G[1], D[1])+2sum(pair(G[j+1], R[j]) for j in 1:steps)
        exact=sum(abs2, A[end])
        expanded=center+sum(abs2, D[end])+2sum(
            pair(T[j], D[j]) + pair(G[j+1], reflect(transpose(D[j])*D[j], j))/kappas[j] for
            j in 1:steps
        )
        @test expanded == exact
        uppernorm(X) = MN.relative_upper_sqrt(sum(abs2, X))
        bound=AT.trace_interval(
            center,
            center;
            errors = uppernorm.(D),
            dual_defects = uppernorm.(T),
            dual_operator_bounds = uppernorm.(G[2:end]),
            normalizations = kappas,
        )
        @test bound.lower <= exact <= bound.upper
    end
end
