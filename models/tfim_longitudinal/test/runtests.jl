using Test
include("../check.jl")
@testset "Exact RP, sector count, stored enclosures and recurrence" begin
    for g in C.GS, h in C.HS
        @test C.rp_check(g,h)
    end
    for L in 2:12
        @test sum(C.sector_dimensions(C.orbits(L)))==2^L
    end
    paths=filter(p->endswith(p,".json") && startswith(basename(p),"g"),
        readdir(joinpath(C.ROOT,"results");join=true))
    @test !isempty(paths)
    for path in paths
        @test check_result(path)
    end
end
