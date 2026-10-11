include(joinpath(@__DIR__, "../scripts/mpo_adjoint_flow.jl"))
using LinearAlgebra, Random, Test

function expand_adjoint_test(cores)
    partial = [ones(MN.F, 1, 1)]
    for core in cores
        d = size(core, 2)
        next = [zeros(MN.F, d*size(partial[1], 1), d*size(partial[1], 2)) for _ in axes(core, 4)]
        for a in axes(core, 4), b in axes(core, 1)
            next[a] += kron(partial[b], core[b, :, :, a])
        end
        partial = next
    end
    return only(partial)
end

@testset "Physical adjoint and reflection in backward MPO products" begin
    rng = MersenneTwister(76910)
    w = kron(PM.WEIGHTS, PM.WEIGHTS)
    physical_adjoint(X) = Diagonal(w)*transpose(X)*Diagonal(inv.(w))
    reflect(X) = reshape(permutedims(reshape(X, 5, 5, 5, 5), (2, 1, 4, 3)), 25, 25)
    physical_inner(X, Y) = sum(X[i, j]*Y[i, j]*w[j]/w[i] for i in 1:25, j in 1:25)
    B = [MN.F.(rand(rng, -2:2, 1, 5, 5, 2))/30, MN.F.(rand(rng, -2:2, 2, 5, 5, 1))/30]
    G = [MN.F.(rand(rng, -2:2, 1, 5, 5, 3))/30, MN.F.(rand(rng, -2:2, 3, 5, 5, 1))/30]
    b, g = expand_adjoint_test(B), expand_adjoint_test(G)
    first, second = dual_product_inputs(B, G)
    @test expand_adjoint_test(first) == reflect(b)
    @test expand_adjoint_test(second) == g+physical_adjoint(g)
    kappa = MN.F(3, 5)
    forward_gauges = MN.rational_gauges(B, [0.4, -0.3])
    first_gauges = reflected_dual_gauges(forward_gauges)
    @test expand_adjoint_test(MN.gauge_transform(first, first_gauges)) == reflect(b)
    proof = MP.product_candidate(first, second; maxdim = 4, weights = PM.WEIGHTS, first_gauges)
    @test MP.verify_product(first, second, proof)
    candidate = expand_adjoint_test(dual_product_output(proof, kappa))
    exact = b*(reflect(g)+physical_adjoint(reflect(g)))/kappa
    defect = exact-candidate
    @test physical_inner(defect, defect) <= (proof.error/kappa)^2
    overlap = MP.product_inner_ball(
        MN.reflect_cores(G),
        MN.adjoint_cores(B; weights = PM.WEIGHTS),
        B;
        weights = PM.WEIGHTS,
    )
    @test Nemo.contains(overlap, physical_inner(g, reflect(physical_adjoint(b)*b)))
    matrix_gauges = [ones(MN.F, 1, 1), MN.F[1 1//4; -1//5 1], ones(MN.F, 1, 1)]
    reflected = reflected_dual_gauges(matrix_gauges)
    @test expand_adjoint_test(MN.gauge_transform(first, reflected)) == reflect(b)
    matrix_proof = MP.product_candidate(
        first,
        second;
        maxdim = 4,
        weights = PM.WEIGHTS,
        first_gauges = reflected,
    )
    @test MP.verify_product(first, second, matrix_proof)
    matrix_defect = exact-expand_adjoint_test(dual_product_output(matrix_proof, kappa))
    @test physical_inner(matrix_defect, matrix_defect) <= (matrix_proof.error/kappa)^2
end
