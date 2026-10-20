# Fast tests of all six studies (about 30 s; no Lean). Each study's tests live in
# test/<study>.jl, in a module of their own, so that study namespaces do not collide.
#
#     julia --project=. test/runtests.jl [study ...]
using Test
include(joinpath(@__DIR__, "..", "src", "LeanSeeds.jl"))
using .LeanSeeds

const SELECTED = isempty(ARGS) ? LeanSeeds.Common.STUDIES : Tuple(ARGS)

@testset "lean_seeds" begin
    for study in SELECTED
        @testset "$study" begin
            include(joinpath(@__DIR__, "$study.jl"))
        end
    end
end
