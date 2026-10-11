# Full bounded spectral replay from exact model parameters in saved certificates.
# No scan is read or run; original certificate files are not overwritten.
include("common.jl")
include("historical_provenance.jl")
start=time(); records=Any[]
prior=joinpath(ROOT,"results","spectral_replay_attempt2.json")
if "--resume-attempt2" in ARGS
    saved=JSON3.read(read(prior,String),Dict{String,Any})
    @assert saved["source_sha256"]["common.jl"]==sha256_file(joinpath(ROOT,"historical_sources","common.jl"))
    @assert saved["source_sha256"]["spectral_replay.jl"]==sha256_file(joinpath(ROOT,"spectral_replay_attempt2.jl"))
    for r in saved["records"]
        original=joinpath(ROOT,"results","certificate_$(tag(rat(r["h"]))).json")
        @assert sha256_file(original)==r["original_certificate_sha256"]
        s=enclosed_seeds(r["small_sectors"],r["large_sectors"],rat(r["beta"]))
        @assert s.s==rat(r["epsilon_s"]) && s.t==rat(r["epsilon_t"])
        @assert certified_gap_bound(s.s,s.t,rat(r["beta"])).bound==r["bound"]>0
        r["supports_original_nine_decimal_bound"]=r["supports_published_bound"]
        r["original_nine_decimal_bound_rational"]=r["published_bound_rational"]
        safe=min(Rational{BigInt}(r["bound"]),Rational{BigInt}(r["original_bound"]))
        r["published_bound_rational"]=string(floor(BigInt,safe*10^6)//big(10)^6)
        r["supports_published_bound"]=true
        r["spectral_execution_provenance"]=Dict("attempt"=>"spectral_replay_attempt2.json",
            "sha256"=>sha256_file(prior),"note"=>"Full fresh dense spectral stage executed in preceding bounded attempt; only enclosure propagation replayed during resume")
        push!(records,r)
    end
end
function save_replay(records,start,status)
    write_result(joinpath(ROOT,"results","migration_spectral_replay.json"),Dict(
        "status"=>status,
        "scope"=>"Reconstructed exact model blocks; recomputed complete dense eigenvectors, residuals, orthogonality defects, rigorous radii, logZ, rational seeds and recurrence. No parameter search; no saved eigenvalues or radii reused as verification inputs.",
        "implementation"=>"Fresh execution of existing verifier, not a second independently implemented spectral verifier",
        "acceptance"=>"Every fresh outward bound must be positive. Published bounds are rounded DOWN to six decimals from min(original, fresh). Support of the original nine-decimal claim and binary equality are recorded separately, not assumed. Eigenvalue interval overlap is only a consistency check.",
        "source_sha256"=>file_hashes(ROOT,["common.jl","spectral_replay.jl"]),
        "julia"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
        "seconds"=>time()-start,"records"=>records))
end
for h in FIELDS
    any(r->rat(r["h"])==h,records) && continue
    path=joinpath(ROOT,"results","certificate_$(tag(h)).json")
    c=JSON3.read(read(path,String),Dict{String,Any})
    for (name,hash) in c["metadata"]["source_sha256"]
        source=joinpath(ROOT,"historical_sources",name)
        @assert sha256_file(source)==hash
    end
    hh=rat(c["h"]);n=c["n0"];beta=rat(c["beta"])
    @assert hh==h && exact_rp(hh)
    t=time()
    small=all_spectra(n,hh;verified=true)
    large=all_spectra(2n,hh;verified=true)
    seeds=enclosed_seeds(small.records,large.records,beta)
    bound=certified_gap_bound(seeds.s,seeds.t,beta)
    displayed=floor(BigInt,Rational{BigInt}(c["bound"])*10^9)//big(10)^9
    supports=bound.bound!==nothing && Rational{BigInt}(bound.bound)>=displayed
    positive=bound.bound!==nothing && bound.bound>0
    safe=positive ? min(Rational{BigInt}(bound.bound),Rational{BigInt}(c["bound"])) : 0//big(1)
    published=floor(BigInt,safe*10^6)//big(10)^6
    exact_fields=true;overlap=true;comparisons=Any[]
    for (old,new,L) in ((c["small_sectors"],small.records,n),
                       (c["large_sectors"],large.records,2n))
        @assert length(old)==length(new)
        for (a,b) in zip(old,new)
            @assert all(a[key]==b[key] for key in ("p","k","dimension","multiplicity"))
            same=all(a[key]==b[key] for key in ("eigenvalues","radius","entry_error",
                "residual_upper","orthogonality_upper"))
            joint=all(sup(abs(iv(x)-iv(y)))<=inf(iv(a["radius"])+iv(b["radius"]))
                for (x,y) in zip(a["eigenvalues"],b["eigenvalues"]))
            exact_fields &= same;overlap &= joint
            push!(comparisons,Dict("L"=>L,"p"=>a["p"],"k"=>a["k"],
                "all_spectral_fields_identical"=>same,"all_eigenvalue_intervals_overlap"=>joint,
                "max_eigenvalue_difference_numerical"=>maximum(abs.(a["eigenvalues"].-b["eigenvalues"])),
                "radius_difference_numerical"=>b["radius"]-a["radius"]))
        end
    end
    identical=exact_fields && seeds.s==rat(c["epsilon_s"]) &&
        seeds.t==rat(c["epsilon_t"]) && bound.bound==c["bound"]
    push!(records,Dict("h"=>string(hh),"n0"=>n,"beta"=>string(beta),
        "original_certificate_sha256"=>sha256_file(path),"seconds"=>time()-t,
        "small_sectors"=>small.records,"large_sectors"=>large.records,
        "epsilon_s"=>string(seeds.s),"epsilon_t"=>string(seeds.t),
        "original_epsilon_s"=>c["epsilon_s"],"original_epsilon_t"=>c["epsilon_t"],
        "epsilon_s_difference"=>string(seeds.s-rat(c["epsilon_s"])),
        "epsilon_t_difference"=>string(seeds.t-rat(c["epsilon_t"])),
        "logZ_intervals"=>seeds.logs,"bound"=>bound.bound,"original_bound"=>c["bound"],
        "original_nine_decimal_bound_rational"=>string(displayed),
        "supports_original_nine_decimal_bound"=>supports,
        "published_bound_rational"=>string(published),"supports_published_bound"=>positive,
        "spectral_execution_provenance"=>Dict("attempt"=>"current full dense execution",
            "script_sha256"=>sha256_file(joinpath(ROOT,"spectral_replay.jl"))),
        "all_spectral_fields_and_seeds_identical"=>identical,
        "all_original_and_fresh_intervals_overlap"=>overlap,"comparisons"=>comparisons))
    # Preserve fresh witness data and differences even if acceptance fails.
    save_replay(records,start,positive && overlap ?
        (length(records)==length(FIELDS) ? "complete" : "in_progress") : "failed")
    @assert positive "fresh certificate did not yield a positive bound"
    @assert overlap "original and fresh spectral intervals are inconsistent"
    println("SPECTRAL REPLAY PASSED h=$hh bound=$(bound.bound) binary_identical=$identical");flush(stdout)
end
