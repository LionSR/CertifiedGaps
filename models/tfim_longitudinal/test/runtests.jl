using Test
include("../check.jl")
@testset "Exact RP, sector count, stored enclosures and recurrence" begin
    for g in C.GS, h in C.HS
        @test C.rp_check(g,h)
    end
    for L in 2:12
        @test sum(C.sector_dimensions(C.orbits(L)))==2^L
    end
    paths=C.result_paths()
    @test length(paths)==30
    for path in paths
        @test check_result(path)
    end
end
@testset "Compact seed table matches the certificates" begin
    include(joinpath(@__DIR__,"..","seeds.jl"))
    table=JSON3.read(read(joinpath(C.ROOT,"results","seeds.json"),String),Dict{String,Any})
    @test length(table["points"])==30
    @test [Dict{String,Any}(r) for r in table["points"]]==seed_rows()
end
