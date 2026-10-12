using Test
include("../model.jl")
@testset "Exact transition symmetries and all orbit dimensions" begin
    D=1603//500
    for L in 2:6, s in 0:3^L-1
        h=transitions(s,L,D)
        @test all(charge(t,L)==charge(s,L) for t in keys(h))
        @test transitions(rot(s,L),L,D)==Dict(rot(t,L)=>v for (t,v) in h)
        @test transitions(flip(s,L),L,D)==Dict(flip(t,L)=>v for (t,v) in h)
    end
    for L in 2:12
        count=0
        for q in 0:L
            o=orbit_data(L,q)
            for p in (q==0 ? (1,-1) : (1,)), k in 0:div(L,2)
                mult=(q==0 ? 1 : 2)*(k==0 || 2k==L ? 1 : 2)
                count+=mult*length(active_orbits(o,k,p))
            end
        end
        @test count==3^L
    end
end
