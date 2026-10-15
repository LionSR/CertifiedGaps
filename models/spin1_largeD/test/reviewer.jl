# Independent checks from the stage-1 review (supplementing, not replacing, the enclosures):
# an exact Cartesian-to-normalised-Sz basis bridge for the bond, and 256-bit diagnostics of
# the residual and orthogonality bounds. The BigFloat comparisons are not new enclosures.
using Test
include("../model.jl")

@testset "Independent Cartesian to normalized Sz bridge" begin
    G = [[0 0 0; 0 0 1; 0 -1 0], [0 0 -1; 0 0 0; 1 0 0],
         [0 1 0; -1 0 0; 0 0 0]]
    S = [-im*g for g in G]
    # Columns have m=-1,0,+1. The phase i on the middle column makes
    # the ladder matrix elements positive. Norm squares are exactly 2,1,2.
    V = Complex{Int}[im 0 -im; 1 0 1; 0 im 0]
    ns = [2,1,2]
    @test V'*V == Diagonal(ns)
    @test S[3]*V == V*Diagonal([-1,0,1])
    h = sum(kron(s,s) for s in S)
    B = kron(V,V)'*h*kron(V,V)
    normalized = zeros(Rational{Int},9,9)
    for a in 1:3, b in 1:3, c in 1:3, d in 1:3
        i=3(a-1)+b; j=3(c-1)+d
        n=ns[a]*ns[b]*ns[c]*ns[d]; r=isqrt(n)
        @test isreal(B[i,j])
        if r*r==n
            normalized[i,j]=real(B[i,j])//r
        else
            @test B[i,j]==0
        end
    end
    expected=zeros(Rational{Int},9,9)
    for a in 0:2,b in 0:2
        j=3a+b+1
        expected[j,j]=(a-1)*(b-1)
        a<2 && b>0 && (expected[3(a+1)+(b-1)+1,j]=1)
        a>0 && b<2 && (expected[3(a-1)+(b+1)+1,j]=1)
    end
    @test normalized==expected
    @test V'*(S[3]^2)*V == Diagonal(ns.*[1,0,1])
    # At L=2 the periodic sum contains the same two-site bond twice.
    for D in (4//1,1603//500)
        direct=zeros(Rational{Int},9,9)
        for s in 0:8,(t,v) in transitions(s,2,D)
            direct[t+1,s+1]=v
        end
        onsite=Diagonal([D*((a-1)^2+(b-1)^2) for a in 0:2 for b in 0:2])
        @test direct==2expected+onsite
        @test rp_identity(D)
    end
end

@testset "Residual and orthogonality bounds against 256-bit accumulation" begin
    for D in (4//1,1603//500),L in (3,4),q in (0,1)
        o=orbit_data(L,q)
        for p in (q==0 ? (1,-1) : (1,)),k in 0:div(L,2)
            H,c=block(o,k,p,D;verified=true)
            isempty(H) && continue
            A=Matrix(H); F=eigen(Hermitian(A))
            _,eta,info=enclose_eigen(H,c,L,D)
            setprecision(BigFloat,256) do
                Ab=Complex{BigFloat}.(A); Qb=Complex{BigFloat}.(F.vectors)
                lb=BigFloat.(F.values); m=length(lb)
                residual=sqrt(sum(abs2,Ab*Qb-Qb*Diagonal(lb)))
                orthogonal=sqrt(sum(abs2,Qb'*Qb-Matrix{Complex{BigFloat}}(I,m,m)))
                @test residual < BigFloat(info["residual_upper"])
                @test orthogonal < BigFloat(info["orthogonality_upper"])
                @test eta>=c && info["orthogonality_upper"]<1
            end
        end
    end
end
