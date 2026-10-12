# Full spectral replay from the exact model parameters stored in the certificates.
# No scan is read or run; the certificate files are not overwritten. Every 8/16-site block
# is rebuilt, densely re-diagonalised and re-enclosed, then propagated to logZ, rational
# seeds and the recurrence. This is a fresh execution of the same verifier in common.jl,
# not an independent second implementation. Writes results/spectral_replay_summary.json;
# with --full also writes the per-block fresh witnesses (about 5.7 MB, not committed) to
# results/spectral_replay_full.json.
include("common.jl")
start=time(); records=Any[]; full=Any[]
for h in FIELDS
    path=joinpath(ROOT,"results","certificate_$(tag(h)).json")
    c=JSON3.read(read(path,String),Dict{String,Any})
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
    # The published value is floor_6(min(original, fresh)), so either seed pair proves it.
    safe=positive ? min(Rational{BigInt}(bound.bound),Rational{BigInt}(c["bound"])) : 0//big(1)
    published=floor(BigInt,safe*10^6)//big(10)^6
    exact_fields=true;overlap=true;maxdiff=0.0;maxrad=-Inf;blocks=0
    for (old,new) in ((c["small_sectors"],small.records),(c["large_sectors"],large.records))
        @assert length(old)==length(new)
        for (a,b) in zip(old,new)
            @assert all(a[key]==b[key] for key in ("p","k","dimension","multiplicity"))
            exact_fields &= all(a[key]==b[key] for key in ("eigenvalues","radius",
                "entry_error","residual_upper","orthogonality_upper"))
            # Overlap is a consistency check only; the fresh enclosure is its own proof.
            overlap &= all(sup(abs(iv(x)-iv(y)))<=inf(iv(a["radius"])+iv(b["radius"]))
                for (x,y) in zip(a["eigenvalues"],b["eigenvalues"]))
            maxdiff=max(maxdiff,maximum(abs.(a["eigenvalues"].-b["eigenvalues"])))
            maxrad=max(maxrad,b["radius"]-a["radius"]);blocks+=1
        end
    end
    identical=exact_fields && seeds.s==rat(c["epsilon_s"]) &&
        seeds.t==rat(c["epsilon_t"]) && bound.bound==c["bound"]
    push!(records,Dict("h"=>string(hh),"n0"=>n,"beta"=>string(beta),"seconds"=>time()-t,
        "epsilon_s"=>string(seeds.s),"epsilon_t"=>string(seeds.t),
        "original_epsilon_s"=>c["epsilon_s"],"original_epsilon_t"=>c["epsilon_t"],
        "epsilon_s_difference"=>string(seeds.s-rat(c["epsilon_s"])),
        "epsilon_t_difference"=>string(seeds.t-rat(c["epsilon_t"])),
        "logZ_intervals"=>seeds.logs,"bound"=>bound.bound,"original_bound"=>c["bound"],
        "supports_original_nine_decimal_bound"=>supports,
        "published_bound_rational"=>string(published),
        "all_spectral_fields_and_seeds_identical"=>identical,
        "all_original_and_fresh_intervals_overlap"=>overlap,"blocks_compared"=>blocks,
        "max_eigenvalue_difference_numerical"=>maxdiff,
        "max_radius_increase_numerical"=>maxrad))
    "--full" in ARGS && push!(full,Dict("h"=>string(hh),"small_sectors"=>small.records,
        "large_sectors"=>large.records))
    @assert positive "fresh certificate did not yield a positive bound"
    @assert overlap "original and fresh spectral intervals are inconsistent"
    println("SPECTRAL REPLAY PASSED h=$hh bound=$(bound.bound) binary_identical=$identical")
    flush(stdout)
end
write_result(joinpath(ROOT,"results","spectral_replay_summary.json"),Dict(
    "status"=>"complete",
    "scope"=>"Reconstructed exact model blocks; recomputed complete dense eigenvectors, residuals, orthogonality defects, rigorous radii, logZ, rational seeds and recurrence. No parameter search; no saved eigenvalues or radii reused as verification inputs.",
    "implementation"=>"Fresh execution of existing verifier, not a second independently implemented spectral verifier",
    "acceptance"=>"Every fresh outward bound must be positive. Published bounds are rounded DOWN to six decimals from min(original, fresh).",
    "source_sha256"=>file_hashes(ROOT,["common.jl","spectral_replay.jl"]),
    "julia"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
    "seconds"=>time()-start,"records"=>records))
"--full" in ARGS && write_result(joinpath(ROOT,"results","spectral_replay_full.json"),
    Dict("records"=>full))
