using Test, CertifiedGaps, LinearAlgebra
using CertifiedGaps.ExactDiagonalization
include("../su2.jl")
include("../dimensions.jl")
BLAS.set_num_threads(1)
@testset "Empty highest-weight sectors" begin
    for (L,S,k) in ((2,1,0),(2,2,1))
        A,B,_=su2_components(L,S,k)
        @test size(A)==(0,0)
        @test isempty(eigvals(Symmetric(A+B)))
    end
end
@testset "Exact highest-weight polynomial on small charge spaces" begin
    for L in 2:4,M in 0:L
        states=[s for s in 0:3^L-1 if charge(s,L)==M]
        next=[s for s in 0:3^L-1 if charge(s,L)==M+1]
        index=Dict(s=>i for (i,s) in enumerate(next))
        C=zeros(Rational{BigInt},length(next),length(states))
        for (j,s) in enumerate(states),site in 0:L-1
            s÷3^site%3==2 && continue
            C[index[s+3^site],j]+=1
        end
        casimir=2C'*C
        P=Matrix{Rational{BigInt}}(I,length(states),length(states))
        for J in M+1:L
            P=(I-casimir/(J*(J+1)-M*(M+1)))*P
        end
        @test P^2==P
        @test iszero(C*P)
        @test tr(P)==length(states)-length(next)
    end
end
@testset "SU(2) highest weights versus full dense spectra" begin
    X=sz_bond()
    for L in 2:7
        levels=[Float64[] for _ in 1:3];dim=0
        for S in 0:L,k in 0:L÷2
            o=orbits(L,S)
            count(p->mod(k*p,L)==0,o.periods)==0 && continue
            A,B,diag=su2_components(L,S,k)
            @test size(A,1)==dimension(L,S,k)-dimension(L,S+1,k)
            mult=(2S+1)*((k==0 || 2k==L) ? 1 : 2)
            dim+=mult*size(A,1)
            for (t,theta) in enumerate((0.0,-0.3,-0.7))
                ev=eigvals(Symmetric(cos(theta)*A+sin(theta)*B))
                for _ in 1:mult;append!(levels[t],ev);end
            end
        end
        @test dim==3^L
        for (t,theta) in enumerate((0.0,-0.3,-0.7))
            h=cos(theta)*X+sin(theta)*X^2
            expected=spectrum(Chain(3,ring_bonds(h);charge=[-1,0,1]),L;
                dense_max=10000).E
            @test sort(levels[t])≈expected atol=1e-8 rtol=1e-9
        end
        println("validated SU2 L=$L dimension=$dim");flush(stdout)
    end
end
