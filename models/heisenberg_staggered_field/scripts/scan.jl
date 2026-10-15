# Numerical seed scan (not certified): spectra of rings up to maxL sites (partial Krylov spectra
# above dimension 2500), cached in results/spectrum_<h>_<L>.json (gitignored), and a β grid
# written to results/scan_<h>_<maxL>.json.
#     julia --project=models/heisenberg_staggered_field \
#         models/heisenberg_staggered_field/scripts/scan.jl 16 [h]
include(joinpath(dirname(@__DIR__), "src", "common.jl"))
maxL=length(ARGS)>0 ? parse(Int,ARGS[1]) : 16
fields=length(ARGS)>1 ? [rat(ARGS[2])] : FIELDS
for h in fields
    spectra=Dict{Int,PartialSpectrum}(); timings=Dict(); sectors=Dict()
    lengths=sort!(unique(vcat(collect(2:2:maxL÷2),collect(4:4:maxL))))
    for L in lengths
        path=joinpath(ROOT,"results","spectrum_$(tag(h))_$(L).json")
        if isfile(path)
            d=JSON3.read(read(path,String),Dict{String,Any})
            spectra[L]=PartialSpectrum(Float64.(d["eigenvalues"]),2^L,d["complete"])
            timings[L]=d["seconds"];sectors[L]=d["sectors"]
        else
            a=all_spectra(L,h;dense_max=2500)
            spectra[L]=a.spec;timings[L]=a.seconds;sectors[L]=a.records
            write_result(path,Dict("metadata"=>provenance(Dict("h"=>string(h),"L"=>L)),
                "eigenvalues"=>a.spec.E,"complete"=>a.spec.exact,"seconds"=>a.seconds,
                "sectors"=>a.records))
        end
    end
    rows=scan_seeds(spectral_defects(spectra),2:2:maxL÷2,0.25:0.25:32)
    best=best_seed(rows)
    table=[Dict("n0"=>n,"max_min_purity"=>maximum(min(1-r.εₛ,1-r.εₜ) for r in rows if r.n₀==n),
        "delta"=>spectra[2n].E[1]-2spectra[n].E[1]) for n in 2:2:maxL÷2]
    bestdict=best===nothing ? nothing : Dict("n0"=>best.n₀,"beta"=>best.β₀,
        "epsilon_s"=>best.εₛ,"epsilon_t"=>best.εₜ,"bound"=>best.bound)
    data=Dict("metadata"=>provenance(Dict("h"=>string(h),"maxL"=>maxL)),
        "h"=>string(h),"best"=>bestdict,"per_n"=>table,"timings"=>timings,
        "reference_gap"=>spectra[maximum(lengths)].E[2]-spectra[maximum(lengths)].E[1],
        "grid"=>[Dict("n0"=>r.n₀,"beta"=>r.β₀,"es"=>r.εₛ,"et"=>r.εₜ,"bound"=>r.bound) for r in rows])
    write_result(joinpath(ROOT,"results","scan_$(tag(h))_$(maxL).json"),data)
    println("BEST h=$h $bestdict");flush(stdout)
end
