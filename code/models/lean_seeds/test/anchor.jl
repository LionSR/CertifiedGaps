# Study `anchor`: Hamiltonian identities, exact-checker rejection, the comparison decomposition
# and the missing-level counterexample (about 10 s).
module TestAnchor
using Test, LinearAlgebra, SparseArrays
using Main.LeanSeeds.AnchorStudy
const A=Main.LeanSeeds.AnchorStudy

@testset "Integer spin-one Hamiltonian" begin
    h=Matrix(A.bond(2,collect(0:8),1,2))
    @test (h+2I)*(h+I)*(h-I)==zeros(Int,9,9)
    @test tr(h)==0
    @test tr(h*h)==12
    @test length(A.basis(12,0))==73789
    @test sum(length(A.basis(6,m)) for m in -6:6)==729
    @test eigmin(Symmetric(Matrix{Float64}(A.ring(4,0)))) ≈ -6
    for m in 0:3
        @test eigvals(Symmetric(Matrix{Float64}(A.ring(4,m;twist=-1)))) ≈
            eigvals(Symmetric(Matrix{Float64}(A.ring(4,-m;twist=-1))))
    end
end

@testset "Exact checker detects failure" begin
    @test exact_psd([1 1;1 1]).ok
    @test exact_psd([0 0;0 2]).ok
    @test !exact_psd([1 2;2 1]).ok
    @test !exact_psd([0 1;1 0]).ok
    @test !exact_psd([-1 0;0 1]).ok
    @test_throws AssertionError exact_psd([1 1;0 1])
    mat=QQ[1 1;1 2]
    cert=exact_spectrum_enclosure(mat)
    vals=eigvals(Symmetric(Float64.(mat)))
    @test all(Float64.(cert.lower).<=vals.<=Float64.(cert.upper))
    @test cert.radius<1//10^6
end

@testset "Four-site comparison decomposition" begin
    n=8; l=4; m=6
    b=A.basis(n,m)
    h=[A.bond(n,b,i,mod1(i+1,n)) for i in 1:n]
    u,v,c=-1//3,1//10,-5//1
    B(i)=u*h[i]+v*h[i]^2
    H0=sum(h[i] for i in [1,2,3,5,6,7])+B(1)+B(3)+B(5)+B(7)+2c*I
    R1=h[4]-B(3)-B(5)-c*I
    R2=h[8]-B(7)-B(1)-c*I
    @test sum(h)-H0==R1+R2
    for (j,X) in enumerate(A.LOCAL)
        R=QQ.(X.H)-u*X.L-v*X.Q-c*I
        @test exact_psd(R).ok
        m=j-1; b=A.basis(4,m)
        signs=[(-1)^(((x÷9)%3)+((x÷27)%3)) for x in b]
        Rp=Matrix(A.bond(4,b,2,3;twist=-1))-u*X.L-v*X.Q-c*I
        @test Diagonal(signs)*R*Diagonal(signs)==Rp
    end
end

@testset "Tensor sum multiplicities" begin
    spec=A.comparison_spectra([4,4],0.,0.)
    for m in 0:8
        @test length(spec[m])==length(A.basis(8,m))
    end
end

@testset "Tail scalar certificates and missing-level obstruction" begin
    @test exp_negative_upper(0//1)==1
    @test exp_negative_upper(1//1)<368//1000
    @test exp_negative_upper(1//1)>367//1000
    @test exp_negative_upper(5//1;degree=32)>=exp_negative_upper(5//1;degree=64)
    # Residual-zero Ritz vectors and an upper count of two below rho=0 do not
    # imply that the first Ritz value encloses the first true eigenvalue.
    H=Diagonal([-2,-1,10]);V=[0 0;1 0;0 1];θ=Diagonal([-1,10])
    @test H*V==V*θ
    @test H[3,3]>=0
    @test H[1,1]<θ[1,1]
end
end
