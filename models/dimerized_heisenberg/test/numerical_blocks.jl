include("../numerical_ed.jl")
using CertifiedGaps.Models
@testset "Numerical sector implementation vs core dense ED" begin
    for L in (4,6,8,10)
        d=0.3
        heis=heisenberg(1//2).chain.bonds(L)[1][3]
        chain=Chain(2,L->[(j,mod1(j+1,L),(1+(-1)^j*d)*heis) for j in 1:L];charge=[1,-1])
        spec,blocks=dense_spectrum(d,L)
        reference=spectrum(chain,L;dense_max=10_000)
        @test maximum(abs,spec.E-reference.E)<1e-10
        @test sum(b["dimension"]*b["multiplicity"] for b in blocks)==2^L
        @test length(spec.E)==2^L
        @test abs(sum(spec.E))<1e-10*2^L
        expected=2^L*3/16*(L÷2)*((1+d)^2+(1-d)^2)
        @test abs(sum(abs2,spec.E)-expected)<1e-10*expected
    end
end
