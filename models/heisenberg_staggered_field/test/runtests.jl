# Fast checker: exact RP, small-ring spectra against the core ED, sector counts, log Z and the
# rational seeds recomputed from the stored eigenvalues and radii, the recurrence, the published
# bounds, source provenance, and the exact parts of the shared-core cross-check certificate.
using Test, SHA
include(joinpath(@__DIR__, "..", "src", "common.jl"))
@testset "exact RP and rotation" begin
    for h in FIELDS; @test exact_rp(h);end
end
@testset "independent small-ring spectra and sector counts" begin
    for L in (2,4,6)
        a=all_spectra(L,3//4;verified=true)
        b=spectrum(CertifiedGaps.Models.staggered_heisenberg(3//4).chain,L;dense_max=1000)
        @test a.spec.N==2^L
        @test maximum(abs.(a.spec.E-b.E))<1e-12
        @test sum(r["dimension"]*r["multiplicity"] for r in a.records)==2^L
        for beta in (1//2,3//1)
            z=enclosed_logz(a.records,beta)
            @test inf(z)<=log_partition(b.E,beta)<=sup(z)
        end
    end
end
@testset "stored rational seeds and bound" begin
    files=filter(f->startswith(f,"certificate_") && endswith(f,".json"),readdir(joinpath(ROOT,"results")))
    @test !isempty(files)
    for file in files
        c=JSON3.read(read(joinpath(ROOT,"results",file),String),Dict{String,Any})
        @test exact_rp(parse_rational(c["h"]))
        for (key,L) in (("small_sectors",c["n0"]),("large_sectors",2c["n0"]))
            records=c[key]
            @test all(r["dense"] && length(r["eigenvalues"])==r["dimension"] for r in records)
            @test sum(r["dimension"]*r["multiplicity"] for r in records)==2^L
            @test all(0<=r["orthogonality_upper"]<1 && r["radius"]>0 for r in records)
        end
        beta=parse_rational(c["beta"])
        s=parse_rational(c["epsilon_s"]);t=parse_rational(c["epsilon_t"])
        seeds=enclosed_seeds(c["small_sectors"],c["large_sectors"],beta)
        # Recomputed defects are valid upper bounds; the stored seeds may only be larger.
        @test seeds.s<=s && seeds.t<=t
        bound=certified_gap_bound(s,t,beta)
        @test bound.bound==c["bound"]
        @test bound.bound>0
    end
end
@testset "published bounds follow from both stored and replayed seeds" begin
    summary=JSON3.read(read(joinpath(ROOT,"results","spectral_replay_summary.json"),String),Dict{String,Any})
    @test summary["status"]=="complete" && length(summary["records"])==length(FIELDS)
    for r in summary["records"]
        h=parse_rational(r["h"]);beta=parse_rational(r["beta"])
        c=JSON3.read(read(joinpath(ROOT,"results","certificate_$(tag(h)).json"),String),Dict{String,Any})
        @test r["n0"]==c["n0"] && iseven(r["n0"]) && beta==parse_rational(c["beta"])
        @test r["original_epsilon_s"]==c["epsilon_s"] && r["original_epsilon_t"]==c["epsilon_t"]
        published=parse_rational(r["published_bound_rational"])
        # Each rigorous seed pair (original certificate and fresh replay) proves the claim.
        for (s,t) in ((c["epsilon_s"],c["epsilon_t"]),(r["epsilon_s"],r["epsilon_t"]))
            gap=certified_gap_bound(parse_rational(s),parse_rational(t),beta).bound
            @test gap!==nothing && Rational{BigInt}(gap)>=published>0
        end
        @test published==floor(BigInt,min(Rational{BigInt}(r["bound"]),
            Rational{BigInt}(c["bound"]))*10^6)//big(10)^6
    end
end
@testset "provenance of the stored certificates" begin
    # The certificates were made when this directory was certificates/heisenberg_staggered_field
    # and the core was named ThermalMPOGap; their source hashes refer to the snapshot of those
    # sources kept in historical_sources/ at commit HISTORICAL. They attest provenance, not
    # correctness: checked when that commit is in the local object store, skipped otherwise.
    HISTORICAL = "cb635a1cdf2ab06630b5b0838fac4356d3c33509"
    available = success(pipeline(`git -C $REPO cat-file -e $HISTORICAL^\{commit\}`;
        stderr=devnull))
    for file in filter(f->startswith(f,"certificate_"),readdir(joinpath(ROOT,"results")))
        c=JSON3.read(read(joinpath(ROOT,"results",file),String),Dict{String,Any})
        hashes=c["metadata"]["source_sha256"]
        @test haskey(c["metadata"],"ThermalMPOGap") && !isempty(hashes)
        if !available
            @warn "Historical snapshot $HISTORICAL unavailable; skipping provenance of $file"
            @test_skip false
            continue
        end
        for (path,hash) in hashes
            spec="$HISTORICAL:certificates/heisenberg_staggered_field/historical_sources/$path"
            @test bytes2hex(sha256(read(`git -C $REPO show $spec`)))==hash
        end
    end
end
module CoreCrossCheck
    # The model definition used by scripts/core_certificate.jl (its main block does not run).
    include(joinpath(@__DIR__, "..", "scripts", "core_certificate.jl"))
end
@testset "shared-core cross-check at h = 1 (exact parts)" begin
    path=joinpath(ROOT,"results","core_certificate_1_1.json")
    cert=CertifiedGaps.Results.read_exact_json(path)
    # The stored model must be the intended Hamiltonian, not merely self-consistent.
    model=CertifiedGaps.RingModels.decode_model(cert["model"])
    reference=CoreCrossCheck.staggered(1)
    @test model.d==reference.d && model.bond==reference.bond && model.onsite==reference.onsite
    seed=cert["seed"]
    c=JSON3.read(read(joinpath(ROOT,"results","certificate_1_1.json"),String),Dict{String,Any})
    @test parse(Int,string(seed["n0"]))==c["n0"]
    @test parse_rational(seed["beta0"])==parse_rational(c["beta"])
    # Exact RP identity and the exact rational chain; the spectra are re-verified by
    # `check.jl --full` (about 6 minutes).
    recomputed=CertifiedGaps.GapCertificates.check_certificate(cert; spectra=false)
    @test recomputed["gap_lower_bound"]==cert["gap_lower_bound"]
    @test parse_rational(cert["gap_lower_bound"])>1147699//1000000
end
