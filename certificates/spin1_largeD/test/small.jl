include("../model.jl")
using Test
@testset "Exact RP and dense block completeness" begin
    @test rp_identity(4//1)
    @test rp_identity(13//4)
    for L in 2:6
        D=13//4
        direct=zeros(3^L,3^L)
        for s in 0:3^L-1, (t,h) in transitions(s,L,D)
            direct[t+1,s+1]=Float64(h)
        end
        @test issymmetric(direct)
        data=compute_spectrum(L,D)
        @test data["complete"]
        @test data["energies"] ≈ eigvals(Symmetric(direct)) atol=1e-11
        if L<=4
            cert=compute_spectrum(L,D;rigorous=true)
            @test all(abs.(cert["energies"].-eigvals(Symmetric(direct))) .<= cert["radii"])
        end
    end
end
