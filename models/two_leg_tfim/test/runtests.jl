using SparseArrays, Test, LinearAlgebra, SHA, CertifiedGaps, CertifiedGaps.Results
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
# Certificate files passed as arguments (check.jl --full passes a fresh one) are checked too.
@testset "Stored rational seed and outward recurrence" begin
    paths=filter(p->startswith(basename(p),"certificate_g") &&
        (endswith(p,".json") || endswith(p,".json.gz")),
        readdir(joinpath(@__DIR__,"../results");join=true))
    @test !isempty(paths)
    append!(paths,abspath.(ARGS))
    for p in paths
        c=read_certificate(p)
        @test rp_identity(parse(Rational{BigInt},c["g"]))
        es=parse(Rational{BigInt},c["epsilon_s"])
        et=parse(Rational{BigInt},c["epsilon_t"])
        beta=parse(Rational{BigInt},c["beta"])
        n=parse(Int,c["n0"]); @test iseven(n)
        @test certified_gap_bound(es,et,beta).bound == c["bound"]
        loaded=stored_spectra(c)
        seeds=seed_defects(loaded[1],loaded[2],beta)
        @test seeds.es == es
        @test seeds.et == et
        # The original certificate's source hashes describe the audited historical snapshot
        # (when this directory was certificates/two_leg_tfim and the core was ThermalMPOGap).
        # They are provenance, not part of the mathematical check: verify them when that
        # commit is in the local object store, and say loudly when it is not.
        historical=basename(p)=="certificate_g3_1.json.gz"
        if historical
            @test certificate_bytes(p) |> sha256 |> bytes2hex == ORIGINAL_CERTIFICATE_SHA256
        end
        if historical && !snapshot_available(joinpath(@__DIR__,".."))
            @warn "Historical snapshot $(HISTORICAL_SNAPSHOT) not in the Git object store; " *
                "skipping $(length(c["metadata"]["source_sha256"])) provenance hashes"
            @test_skip false
        else
            for (rel,hash) in c["metadata"]["source_sha256"]
                # A current certificate enforces only model-local sources; a path outside
                # this directory (the shared core) is provenance and is not compared.
                !historical && startswith(normpath(rel),"..") && continue
                snapshot=historical ? HISTORICAL_SNAPSHOT : nothing
                @test source_sha256(joinpath(@__DIR__,".."),rel;snapshot)==hash
            end
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
