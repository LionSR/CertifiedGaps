# Numerical seed scan (NOT certified): complete Float64 spectra (core `float_ring_spectrum`) of
# the rings of 4, 6, 8, 12 and 16 spins at six values of delta, and the core seed scan over
# n0 = 2, 3, 4 cells and beta0 = 0.25:0.25:32. Writes results/numerical_delta_<p>_<q>.json
# (gitignored); scripts/summarize_numerics.jl condenses them into the tracked
# results/numerical_summary.json.
using CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds, LinearAlgebra
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
include(joinpath(@__DIR__, "..", "src", "DimerizedHeisenberg.jl"))
using .DimerizedHeisenberg
BLAS.set_num_threads(1)

const REPO_ROOT = normpath(joinpath(MODEL_DIR, "..", ".."))

"`seed_record` of a scan row, or `nothing`."
seedrow(r) = r === nothing ? nothing : seed_record(r)

"""
    run_numerics(; maxseed = 4)

Scan seeds n0 = 2, …, `maxseed` cells for the six grid deltas and write one JSON per delta.
"""
function run_numerics(; maxseed = 4)
    # n0 = 3 is kept as a diagnostic but excluded from the theorem-eligible best seeds.
    ns = collect(2:maxseed)
    lengths = sort!(unique(vcat(ns, 2 .* ns)))
    for d in DELTAS
        specs = Dict{Int,PartialSpectrum}()
        sizes = Any[]
        total = @elapsed for cells in lengths
            elapsed = @elapsed spec = float_ring_spectrum(rotated_model(d), cells)
            specs[cells] = spec
            push!(sizes, Dict("cells" => cells, "spins" => 2cells, "dimension" => spec.N,
                "complete_dense_spectrum" => spec.exact, "ground_energy" => spec.E[1],
                "reference_gap" => spec.E[2] - spec.E[1], "seconds" => elapsed))
            println("delta=$d cells=$cells gap=$(spec.E[2] - spec.E[1]) seconds=$elapsed")
        end
        grid = scan_seeds(spectral_defects(specs), ns, 0.25:0.25:32.0)
        diagnostics = map(purity_peaks(grid)) do peak
            n = peak["n0"]
            merge(peak, Dict("even_seed_eligible" => iseven(n),
                "energy_difference" => specs[2n].E[1] - 2specs[n].E[1]))
        end
        best = best_seed(filter(r -> iseven(r.n₀), grid))
        smallest = minimum((r.n₀ for r in grid if iseven(r.n₀) && r.bound !== nothing);
            init = typemax(Int))
        top = maximum(lengths)
        refgap = specs[top].E[2] - specs[top].E[1]
        hashes = source_hashes(REPO_ROOT; dirs = ("julia/src",), files = (
            "models/dimerized_heisenberg/src/DimerizedHeisenberg.jl",
            "models/dimerized_heisenberg/scripts/numerical_scan.jl"))
        data = Dict("metadata" => run_metadata(Dict("delta" => string(d),
                "maxseed" => maxseed, "beta_grid" => "0.25:0.25:32.0");
                source_sha256 = hashes, packages = (CertifiedGaps,),
                status = "numerical_not_certified", certified = false),
            "delta" => string(d), "seconds_spectra" => total, "sizes" => sizes,
            "diagnostics" => diagnostics, "best_even_seed" => seedrow(best),
            "smallest_successful_even_n0" => (smallest == typemax(Int) ? nothing : smallest),
            "reference_gap" => refgap, "reference_spins" => 2top,
            "bound_to_finite_ring_gap_ratio" => (best === nothing ? nothing : best.bound / refgap),
            "scan" => [seed_record(r) for r in grid],
            "rigor_limit" => "Full dense spectra with dimension checks, but no eigenvalue " *
                "enclosures or certified purity seeds; all gap values numerical.")
        path = joinpath(MODEL_DIR, "results",
            "numerical_delta_$(numerator(d))_$(denominator(d)).json")
        write_result(path, data)
        println("saved $path best=$(seedrow(best))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_numerics()
end
