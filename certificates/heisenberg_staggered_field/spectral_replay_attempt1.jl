# Full bounded spectral replay from exact model parameters in saved certificates.
# No scan is read or run; original certificate files are not overwritten.
include("common.jl")
start=time(); records=Any[]
for h in FIELDS
    path=joinpath(ROOT,"results","certificate_$(tag(h)).json")
    c=JSON3.read(read(path,String),Dict{String,Any})
    for (name,hash) in c["metadata"]["source_sha256"]
        source=startswith(name,"core/") ? joinpath(REPO,"julia",name[6:end]) : joinpath(ROOT,name)
        @assert sha256_file(source)==hash
    end
    hh=rat(c["h"]);n=c["n0"];beta=rat(c["beta"])
    @assert exact_rp(hh)
    t=time()
    small=all_spectra(n,hh;verified=true)
    large=all_spectra(2n,hh;verified=true)
    seeds=enclosed_seeds(small.records,large.records,beta)
    bound=certified_gap_bound(seeds.s,seeds.t,beta)
    @assert seeds.s==rat(c["epsilon_s"])
    @assert seeds.t==rat(c["epsilon_t"])
    @assert bound.bound==c["bound"]
    for (old,new) in ((c["small_sectors"],small.records),(c["large_sectors"],large.records))
        @assert length(old)==length(new)
        for (a,b) in zip(old,new)
            for key in ("p","k","dimension","multiplicity","eigenvalues","radius",
                        "entry_error","residual_upper","orthogonality_upper")
                @assert a[key]==b[key] "spectral replay mismatch: $key"
            end
        end
    end
    push!(records,Dict("h"=>string(hh),"n0"=>n,"beta"=>string(beta),
        "original_certificate_sha256"=>sha256_file(path),"seconds"=>time()-t,
        "small_sectors"=>small.records,"large_sectors"=>large.records,
        "epsilon_s"=>string(seeds.s),"epsilon_t"=>string(seeds.t),
        "logZ_intervals"=>seeds.logs,"bound"=>bound.bound,
        "all_spectral_fields_and_seeds_identical"=>true))
    println("SPECTRAL REPLAY PASSED h=$hh bound=$(bound.bound)");flush(stdout)
    write_result(joinpath(ROOT,"results","spectral_replay.json"),Dict(
        "status"=>length(records)==length(FIELDS) ? "complete" : "in_progress",
        "scope"=>"Reconstructed exact model blocks; recomputed complete dense eigenvectors, residuals, orthogonality defects, rigorous radii, logZ, rational seeds and recurrence. No parameter search; no saved eigenvalues or radii reused as verification inputs.",
        "implementation"=>"Fresh execution of existing verifier, not a second independently implemented spectral verifier",
        "source_sha256"=>file_hashes(ROOT,["common.jl","spectral_replay.jl"]),
        "julia"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
        "seconds"=>time()-start,"records"=>records))
end
