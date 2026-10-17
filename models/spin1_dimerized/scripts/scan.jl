# Numerical seed scan (not certified) at δ = d/10: dense symmetry-block spectra of rings of
# 2, 3, 4, 6 cells, cached in results/spectrum_d<d>_c<cells>.json (gitignored; six cells take
# about 9 min), and the β grid in results/scan_d<d>.json.
#     julia --project=julia models/spin1_dimerized/scripts/scan.jl 7 2,3,4,6
include(joinpath(dirname(@__DIR__), "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps.Results, CertifiedGaps.Seeds, JSON3
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
const D = DimerizedCertificate
# Separate invocations per parameter keep the runtime independently bounded.
delta = length(ARGS) > 0 ? parse(Int, ARGS[1]) // 10 : 9 // 10
sizes = length(ARGS) > 1 ? parse.(Int, split(ARGS[2], ",")) : [2, 3, 4, 6]
tag = Int(10delta)
for cells in sizes
    path = joinpath(D.ROOT, "results", "spectrum_d$(tag)_c$cells.json")
    isfile(path) && continue
    spec, sectors, seconds = D.dense_spectrum(cells, delta)
    config = Dict("delta" => D.str(delta), "cells" => cells)
    write_result(path, Dict("metadata" => model_metadata(D.ROOT, config; D.SOURCES...),
        "energies" => spec.E, "dimension" => spec.N, "sectors" => sectors, "seconds" => seconds,
        "ground_energy" => spec.E[1], "gap" => spec.E[2] - spec.E[1]))
end
spectra = Dict{Int,PartialSpectrum}()
for cells in (2, 3, 4, 6)
    path = joinpath(D.ROOT, "results", "spectrum_d$(tag)_c$cells.json")
    if isfile(path)
        data = JSON3.read(read(path, String))
        spectra[cells] = PartialSpectrum(collect(Float64, data.energies), Int(data.dimension), true)
    end
end
ns = [n for n in (2, 3) if haskey(spectra, n) && haskey(spectra, 2n)]
results = scan_seeds(spectral_defects(spectra), ns, 0.25:0.25:32)
rows = [Dict("n0" => r.n₀, "beta0" => r.β₀, "epsilon_s" => r.εₛ, "epsilon_t" => r.εₜ,
    "bound" => r.bound, "eligible_positive_transfer_seed" => true) for r in results]
best = best_seed(results)
summary = [Dict("n0" => n,
    "max_min_purity" => maximum(min(1 - r.εₛ, 1 - r.εₜ) for r in results if r.n₀ == n),
    "energy_defect" => spectra[2n].E[1] - 2spectra[n].E[1]) for n in ns]
write_result(joinpath(D.ROOT, "results", "scan_d$tag.json"), Dict(
    "metadata" => model_metadata(D.ROOT,
        Dict("delta" => D.str(delta), "beta_grid" => "1//4:1//4:32//1"); D.SOURCES...),
    "rows" => rows, "summary" => summary,
    "best" => best === nothing ? nothing : Dict("n0" => best.n₀, "beta0" => best.β₀,
        "epsilon_s" => best.εₛ, "epsilon_t" => best.εₜ, "bound" => best.bound)))
println("SUMMARY delta=$delta $summary BEST $best")
