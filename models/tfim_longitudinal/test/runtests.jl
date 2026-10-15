# Fast tests (~20 s): the fast checks of check.jl on all 30 certificates plus the seed table.
using Test, JSON3
include(joinpath(@__DIR__, "..", "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate
const C=LongitudinalCertificate
@testset "Exact RP, sector count, stored enclosures and recurrence" begin
    for g in C.GS, h in C.HS
        @test C.rp_check(g,h)
    end
    for L in 2:12
        @test sum(C.sector_dimensions(C.orbits(L)))==2^L
    end
    paths=C.result_paths()
    @test C.check_grid(paths)
    for path in paths
        @test C.check_result(path)
    end
end
@testset "Compact seed table matches the certificates" begin
    table=JSON3.read(read(joinpath(C.ROOT,"results","seeds.json"),String),Dict{String,Any})
    @test length(table["points"])==30
    @test [Dict{String,Any}(r) for r in table["points"]]==C.seed_rows()
end
