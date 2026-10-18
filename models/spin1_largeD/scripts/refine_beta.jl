# Floating-point beta refinement near the threshold (not a certificate): n0 = 6 and
# beta = 2.75:0.002:3.75, from the cached spectra of scripts/scan.jl.
#
#   julia --project=julia models/spin1_largeD/scripts/refine_beta.jl D ...
#
# Writes the gitignored results/beta_refinement_D*.json, summarized by scripts/report.jl.
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
using CertifiedGaps.Seeds: scan_seeds, spectral_defects, best_seed, seed_record
using CertifiedGaps.Results, JSON3
include(joinpath(dirname(@__DIR__), "src", "model.jl"))
using .LargeDModel
const ROOT = dirname(@__DIR__)
for arg in ARGS
    D = parse(Rational{Int}, arg)
    spectra = Dict{Int,PartialSpectrum}()
    for L in (6, 12)
        data = JSON3.read(read(joinpath(ROOT, "results", "scan_$(stem(D))_L$L.json"), String),
            Dict{String,Any})
        spectra[L] = PartialSpectrum(Float64.(data["energies"]), Int(data["dimension"]))
    end
    results = scan_seeds(spectral_defects(spectra), [6], 2.75:0.002:3.75)
    best = best_seed(results)
    data = Dict("D" => string(D), "certified" => false,
        "beta_grid" => "2.75:0.002:3.75", "n0" => 6,
        "best" => best === nothing ? nothing : merge(seed_record(best),
            Dict("beta_exact" => string(round(Int, 1000best.β₀) // 1000))),
        "rows" => seed_record.(results))
    write_result(joinpath(ROOT, "results", "beta_refinement_$(stem(D)).json"), data)
    println("FINE D=$D best=$best"); flush(stdout)
end
