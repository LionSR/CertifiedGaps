using SparseArrays, Test, LinearAlgebra, CertifiedGaps, CertifiedGaps.Results
include(joinpath(@__DIR__, "../src/Ladder.jl"))
using .Ladder
include(joinpath(@__DIR__, "../src/Provenance.jl"))
using .Provenance
BLAS.set_num_threads(1)
@testset "Exact ladder RP and full sector resolution" begin
    for g in (2//1,5//2,3//1,4//1)
        @test rp_identity(g)
    end
    for L in 2:4
        bs=blocks(L)
        @test sum(length(b.reps) for b in bs)==4^L
        # Independent full Hamiltonian in the X basis, including the L=2 doubled bonds.
        N=4^L; H=zeros(N,N)
        for s in 0:N-1
            H[s+1,s+1]=-3*(2L-2count_ones(s))
            for j in 0:L-1
                H[(s ⊻ (1<<(2j)) ⊻ (1<<(2j+1)))+1,s+1]-=1
                for a in 0:1
                    H[(s ⊻ (1<<(2j+a)) ⊻ (1<<(2mod(j+1,L)+a)))+1,s+1]-=1
                end
            end
        end
        ev=Float64[]
        for b in bs
            M,A=block_matrix(b,3//1)
            append!(ev,eigvals(Symmetric(M)))
        end
        @test sort(ev) ≈ eigvals(Symmetric(H)) atol=1e-11
    end
end
@testset "Stored rational seed and outward recurrence" begin
    paths=filter(p->startswith(basename(p),"certificate_g") && endswith(p,".json"),
        readdir(joinpath(@__DIR__,"../results");join=true))
    @test !isempty(paths)
    for p in paths
        c=read_exact_json(p)
        @test rp_identity(parse(Rational{BigInt},c["g"]))
        es=parse(Rational{BigInt},c["epsilon_s"])
        et=parse(Rational{BigInt},c["epsilon_t"])
        beta=parse(Rational{BigInt},c["beta"])
        n=parse(Int,c["n0"]); @test iseven(n)
        @test certified_gap_bound(es,et,beta).bound == c["bound"]
        loaded=[(;sectors=[(;eigenvalues=Float64.(sec["eigenvalues"]),
            verification=(;radius=sec["verification"]["radius"])) for sec in sp["sectors"]])
            for sp in c["spectra"]]
        ln=logz_interval(loaded[1],beta); l2=logz_interval(loaded[2],beta)
        l22=logz_interval(loaded[2],2beta)
        ls=Ladder.down(()->l2[1]-2ln[2]); lt=Ladder.down(()->l22[1]-2l2[2])
        @test Rational{BigInt}(-Ladder.down(()->expm1(ls))) == es
        @test Rational{BigInt}(-Ladder.down(()->expm1(lt))) == et
        for (rel,hash) in c["metadata"]["source_sha256"]
            snapshot=basename(p)=="certificate_g3_1.json" ? HISTORICAL_SNAPSHOT : nothing
            @test source_sha256(joinpath(@__DIR__,".."),rel;snapshot)==hash
        end

        for spec in c["spectra"]
            L=parse(Int,spec["L"])
            @test sum(parse(Int,s["dimension"]) for s in spec["sectors"])==4^L
            @test all(length(s["eigenvalues"])==parse(Int,s["dimension"]) for s in spec["sectors"])
        end
    end
end
@testset "Enclosure arithmetic on an exactly soluble block" begin
    # Known spectrum {−sqrt(5),sqrt(5)} for an exact rational symmetric matrix.
    H=[1.0 2.0;2.0 -1.0]
    A=2Int.(H); b=(sizes=[1,1],)
    F=eigen(Symmetric(H)); v=verify_eigen(H,sparse(A),b,F)
    eta=BigFloat(parse(Rational{BigInt},v.radius))
    @test abs(BigFloat(F.values[1])+sqrt(big(5.0)))<=eta
    @test abs(BigFloat(F.values[2])-sqrt(big(5.0)))<=eta
end
