# Numerical seed scan (not certified) at δ = d/10: floating-point spectra of rings of 2, 3, 4, 6
# cells from the core (`VerifiedSpectra.float_ring_spectrum` of the certificate model, whose
# spectrum is that of the physical chain), cached in results/spectrum_d<d>_c<cells>.json
# (gitignored), and the β grid in results/scan_d<d>.json.
#     julia --project=julia models/spin1_dimerized/scripts/scan.jl 7 2,3,4,6
include(joinpath(dirname(@__DIR__), "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps.Results, CertifiedGaps.Seeds, JSON3
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
const D = DimerizedCertificate
# Separate invocations per parameter keep the runtime independently bounded.
delta = length(ARGS) > 0 ? parse(Int, ARGS[1]) // 10 : 9 // 10
sizes = length(ARGS) > 1 ? parse.(Int, split(ARGS[2], ",")) : [2, 3, 4, 6]
tag = Int(10delta)
spectrum_path(cells) = joinpath(D.ROOT, "results", "spectrum_d$(tag)_c$cells.json")
for cells in sizes
    isfile(spectrum_path(cells)) && continue
    seconds = @elapsed spec = float_ring_spectrum(D.dimerized(delta), cells)
    config = Dict("delta" => D.str(delta), "cells" => cells)
    write_result(spectrum_path(cells), Dict("metadata" => model_metadata(D.ROOT, config),
        "energies" => spec.E, "dimension" => spec.N, "seconds" => seconds,
        "ground_energy" => spec.E[1], "gap" => spec.E[2] - spec.E[1]))
end
spectra = Dict{Int,PartialSpectrum}()
for cells in (2, 3, 4, 6)
    isfile(spectrum_path(cells)) || continue
    data = JSON3.read(read(spectrum_path(cells), String))
    spectra[cells] = PartialSpectrum(collect(Float64, data.energies), Int(data.dimension))
end
ns = [n for n in (2, 3) if haskey(spectra, n) && haskey(spectra, 2n)]
results = scan_seeds(spectral_defects(spectra), ns, 0.25:0.25:32)
best = best_seed(results)
# Per n₀ the purity peak, plus the ground-energy doubling mismatch E₀(2n) − 2E₀(n).
summary = [merge(peak, Dict("energy_defect" => spectra[2n].E[1] - 2spectra[n].E[1]))
           for (n, peak) in zip(ns, purity_peaks(results))]
write_result(joinpath(D.ROOT, "results", "scan_d$tag.json"), Dict(
    "metadata" => model_metadata(D.ROOT,
        Dict("delta" => D.str(delta), "beta_grid" => "1//4:1//4:32//1")),
    "rows" => seed_record.(results), "summary" => summary,
    "best" => best === nothing ? nothing : seed_record(best)))
println("SUMMARY delta=$delta $summary BEST $best")
