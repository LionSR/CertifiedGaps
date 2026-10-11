using Test, CertifiedGaps, JSON3, LinearAlgebra
using CertifiedGaps.Results
BLAS.set_num_threads(1)
include("../rp.jl")
root=normpath(joinpath(@__DIR__,"../../.."))
parseq(s)=parse(Rational{BigInt},s)
@testset "Exact reflection positivity" begin
    for J in (4//1,7//2,3//1,5//2,2//1,3//2)
        @test exact_rp(J)
    end
end
@testset "Stored rational certificate" begin
    c=JSON3.read(read(joinpath(@__DIR__,"../results/certificate.json"),String))
    @test iseven(c.n0)
    b=certified_gap_bound(parseq(c.eps_s),parseq(c.eps_t),parseq(c.beta))
    @test b.bound !== nothing
    @test b.bound == c.bound
    @test b.bound > 0
    for L in (Int(c.n0),2Int(c.n0))
        s=JSON3.read(read(joinpath(@__DIR__,"../results/verified_L$(L).json"),String))
        @test sum(r.dimension*r.multiplicity for r in s.sectors)==4^L==s.dimension
        @test length(s.energies)==length(s.radii)==4^L
        @test all(parseq(r.orthogonality_upper)<1 for r in s.sectors)
    end
    for (path,hash) in pairs(c.metadata.source_sha256)
        @test sha256_file(joinpath(root,string(path)))==hash
    end
end
include("../recheck_seeds.jl")
@testset "Fresh full spectral verification" begin
    cert=JSON3.read(read(joinpath(@__DIR__,"../results/certificate.json"),String))
    n=Int(cert.n0); beta=parseq(cert.beta)
    a=full_spectrum(n,4//1;verified=true)
    b=full_spectrum(2n,4//1;verified=true)
    la=enclosed_logZ(a,beta); lb=enclosed_logZ(b,beta); lc=enclosed_logZ(b,2beta)
    @test Rational{BigInt}(sup(-expm1(lb-2la))) <= parseq(cert.eps_s)
    @test Rational{BigInt}(sup(-expm1(lc-2lb))) <= parseq(cert.eps_t)
end
