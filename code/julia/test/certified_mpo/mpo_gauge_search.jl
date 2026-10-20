include(joinpath(@__DIR__, "setup.jl"))
using .CM.GaugeSearch, LinearAlgebra, Random, Test
const MG = CM.GaugeSearch
const MN = MG.MN

function expand_gauge_test(cores)
    partial = [ones(MN.F, 1, 1)]
    for core in cores
        next = [zeros(MN.F, 2size(partial[1], 1), 2size(partial[1], 2)) for _ in axes(core, 4)]
        for a in axes(core, 4), b in axes(core, 1)
            next[a] += kron(partial[b], core[b, :, :, a])
        end
        partial = next
    end
    return only(partial)
end

@testset "Schmidt proposals preserve redundant auxiliary directions" begin
    rng = MersenneTwister(9481)
    for dims in ([1, 2, 2, 1], [1, 6, 6, 6, 1])
        exact = [MN.F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1]))/8 for j in 1:(length(dims)-1)]
        cores = [Float64.(c) for c in exact]
        gauges = MG.schmidt_gauge_proposal(cores)
        @test size.(gauges) == [(d, d) for d in dims]
        @test all(g -> minimum(svdvals(g)) > 0, gauges)
        @test isfinite(MG.matrix_gauge_objective(cores, gauges))
        rational = MN.quantize(gauges; bits = 64)
        @test expand_gauge_test(MN.gauge_transform(exact, rational)) == expand_gauge_test(exact)
    end
end
