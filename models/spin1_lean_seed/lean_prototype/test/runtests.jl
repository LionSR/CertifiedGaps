# Exact generator regression tests; these do not confer Lean soundness.
using Test
include(joinpath(@__DIR__,"..","scripts","generate.jl"))
@testset "packed generator, independent rows" begin
    for n in 2:6
        basis,lookup,rows=sector(n)
        v=BigInt[mod(s,17) for s in basis]
        full=zeros(BigInt,3^n)
        for (i,s) in enumerate(basis)
            full[s+1]=v[i]
        end
        p=pack(64,full); hp=applypacked(64,n,p); hv=applyrows(rows,v)
        @test all(lane(64,s,hp)==hv[i] for (i,s) in enumerate(basis))
        @test reversepacked(64,3^n,p)==pack(64,reverse(full))
        @test dotpacked(64,3^n,p,p)==dot(v,v)
        @test all(lane(64,i,mask(64,n,j,t))==(i÷3^j%3==t ? big(2)^64-1 : 0)
            for j in 0:n-1 for t in 0:2 for i in 0:3^n-1)
    end
    @test dotpacked(2,2,pack(2,[3,3]),pack(2,[3,3]))==0 # rejects unchecked carry assumption
    @test dot([3,3],[3,3])==18
end
