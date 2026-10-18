# Numerical seed scan (not certified): floating-point spectra of the rings of 2, 3, 4, 6, 8
# rungs at field g from the core's exact block reduction (`VerifiedSpectra.float_ring_spectrum`)
# and a β grid; writes results/spectra_g<g>.json (gitignored) and results/scan_g<g>.json.
#     julia --project=julia models/two_leg_tfim/scripts/scan.jl 3//1
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "src", "Ladder.jl"))
using .Ladder, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
using LinearAlgebra
BLAS.set_num_threads(1)

g = parse(Rational{Int}, get(ARGS, 1, "2//1"))
tag = "$(numerator(g))_$(denominator(g))"
metadata(config) = model_metadata(ROOT, config; files = ("scripts/scan.jl",),
    packages = (CertifiedGaps,))
specs, records = Dict{Int,Any}(), []
for L in [2, 3, 4, 6, 8]
    seconds = @elapsed specs[L] = float_ring_spectrum(ladder(g), L)
    push!(records, (; L, g = string(g), energies = specs[L].E, seconds))
    println("L=$L g=$g seconds=$(round(seconds, digits = 2))"); flush(stdout)
    write_result(joinpath(ROOT, "results", "spectra_g$tag.json"),
        (; status = "numerical, not certified", records, metadata = metadata((; g = string(g))));
        compact = true)
end
rs = scan_seeds(spectral_defects(specs), 2:4, 1/8:1/8:32)
best = best_seed(filter(r -> iseven(r.n₀), rs))  # odd n₀ would need a positivity argument
# Per n: the purity peak of the scan and the ground-energy doubling mismatch δ(n).
summary = [merge(peak, Dict("delta" => specs[2peak["n0"]].E[1] - 2specs[peak["n0"]].E[1]))
           for peak in purity_peaks(rs)]
write_result(joinpath(ROOT, "results", "scan_g$tag.json"),
    (; g = string(g), best = best === nothing ? nothing : seed_record(best), summary,
        grid = seed_record.(rs),
        metadata = metadata((; g = string(g), ns = "2:4", beta = "1//8:1//8:32//1",
            max_rungs = 8))); compact = true)  # a one-line file, as the grid is long
println("BEST ", best); println(summary)
