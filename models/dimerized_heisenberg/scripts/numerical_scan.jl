# Numerical seed scan (NOT certified): complete Float64 spectra of the rings of 4, 6, 8, 12 and
# 16 spins at six values of delta, and the core seed scan over n0 = 2, 3, 4 cells and
# beta0 = 0.25:0.25:32. Writes results/numerical_delta_<p>_<q>.json (gitignored, about 2.5 min);
# scripts/summarize_numerics.jl condenses them into the tracked results/numerical_summary.json.
include(joinpath(@__DIR__, "..", "src", "sectors.jl"))
using CertifiedGaps.Results, CertifiedGaps.Seeds

const MODEL_DIR = normpath(joinpath(@__DIR__, ".."))
const REPO_ROOT = normpath(joinpath(MODEL_DIR, "..", ".."))

"A seed-scan row as a JSON dictionary (`nothing` stays `nothing`)."
function seedrow(r)
    r===nothing && return nothing
    Dict("n0_cells"=>r.n₀,"beta0"=>r.β₀,"epsilon_s"=>r.εₛ,
        "epsilon_t"=>r.εₜ,"numerical_bound"=>r.bound)
end

"""
    run_numerics(; maxseed=4)

Scan seeds n0 = 2, …, `maxseed` cells for the six tested deltas and write one JSON per delta.
"""
function run_numerics(; maxseed=4)
    # n0=3 is retained as a diagnostic but excluded from theorem-eligible best seeds.
    ns=collect(2:maxseed); lengths=sort!(unique(vcat(ns,2 .* ns)))
    deltas=[big(1)//2,big(3)//10,big(1)//4,big(1)//5,big(3)//20,big(1)//10]
    for d in deltas
        delta=Float64(d); specs=Dict{Int,PartialSpectrum}(); sizes=Any[]
        total=@elapsed begin
            for cells in lengths
                elapsed=@elapsed spec,blocks=dense_spectrum(delta,2*cells)
                specs[cells]=spec
                push!(sizes,Dict("cells"=>cells,"spins"=>2*cells,"dimension"=>spec.N,
                    "complete_dense_spectrum"=>true,"ground_energy"=>spec.E[1],
                    "reference_gap"=>spec.E[2]-spec.E[1],"seconds"=>elapsed,"sectors"=>blocks))
                println("delta=$d cells=$cells dimension=$(spec.N) ",
                    "gap=$(spec.E[2]-spec.E[1]) seconds=$elapsed"); flush(stdout)
            end
        end
        grid=scan_seeds(spectral_defects(specs),ns,0.25:0.25:32.0)
        diagnostics=Any[]
        for n in ns
            rows=filter(r->r.n₀==n,grid)
            peak=argmax(r->min(1-r.εₛ,1-r.εₜ),rows)
            push!(diagnostics,Dict("n0_cells"=>n,"even_seed_eligible"=>iseven(n),
                "max_min_purity"=>min(1-peak.εₛ,1-peak.εₜ),"peak_beta"=>peak.β₀,
                "energy_difference"=>specs[2n].E[1]-2*specs[n].E[1],
                "best_seed"=>seedrow(best_seed(rows))))
        end
        best=best_seed(filter(r->iseven(r.n₀),grid))
        smallest=minimum((r.n₀ for r in grid if iseven(r.n₀)&&r.bound!==nothing);
            init=typemax(Int))
        refgap=specs[maximum(lengths)].E[2]-specs[maximum(lengths)].E[1]
        hashes=source_hashes(REPO_ROOT;dirs=("julia/src",),files=(
            "models/dimerized_heisenberg/src/sectors.jl",
            "models/dimerized_heisenberg/scripts/numerical_scan.jl"))
        data=Dict("metadata"=>run_metadata(Dict("delta"=>string(d),"maxseed"=>maxseed,
                "beta_grid"=>"0.25:0.25:32.0");
                source_sha256=hashes,packages=(CertifiedGaps,),
                status="numerical_not_certified",certified=false),
            "delta"=>string(d),"seconds_spectra"=>total,"sizes"=>sizes,
            "diagnostics"=>diagnostics,"best_even_seed"=>seedrow(best),
            "smallest_successful_even_n0"=>(smallest==typemax(Int) ? nothing : smallest),
            "reference_gap"=>refgap,"reference_spins"=>2*maximum(lengths),
            "bound_to_finite_ring_gap_ratio"=>(best===nothing ? nothing : best.bound/refgap),
            "scan"=>[seedrow(r) for r in grid],
            "rigor_limit"=>"Full dense spectra with dimension checks, but no eigenvalue " *
                "enclosures or certified purity seeds; all gap values numerical.")
        path=joinpath(MODEL_DIR,"results",
            "numerical_delta_$(numerator(d))_$(denominator(d)).json")
        write_result(path,data)
        println("saved $path best=$(seedrow(best))"); flush(stdout)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_numerics()
end
