# Floating-point seed scan (not a certificate): for every grid point, the complete spectra of
# the rings of 2 to 12 sites (core `float_ring_spectrum`), the seed defects over
# n₀ ∈ {2, 4, 6} and β₀ ∈ {1/4, …, 64}, and the seed with the best `certified_gap_bound`.
# Writes results/scan.json and reports whether each best β₀ agrees with the table `BETA0` that
# certify.jl uses.
#     julia --project=julia models/tfim_longitudinal/scripts/scan.jl
include(joinpath(@__DIR__, "..", "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, Printf, CertifiedGaps
using CertifiedGaps.Results, CertifiedGaps.Seeds
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
const C = LongitudinalCertificate

start = time(); ns = [2, 4, 6]; betas = collect(0.25:0.25:64); rows = Any[]; agree = 0
for (g, h) in C.points()
    model = C.longitudinal(g, h)
    spectra = Dict(L => float_ring_spectrum(model, L) for L in sort(unique(vcat(ns, 2ns))))
    best = best_seed(scan_seeds(spectral_defects(spectra), ns, betas))
    gap12 = spectra[12].E[2] - spectra[12].E[1]
    same = best !== nothing && best.n₀ == C.N0 && Rational{Int}(best.β₀) == C.BETA0[(g, h)]
    global agree += same
    push!(rows, merge(seed_record(best), Dict("g" => string(g), "h" => string(h),
        "gap_L12" => gap12, "matches_BETA0" => same)))
    @printf("g=%s h=%s best (n0, beta0) = (%d, %g) bound %.4f gap(12) %.4f %s\n", g, h,
        best.n₀, best.β₀, best.bound, gap12, same ? "" : "differs from BETA0")
    flush(stdout)
end
write_result(joinpath(C.ROOT, "results", "scan.json"), Dict(
    "status" => "floating-point diagnostics, not a certificate",
    "beta_grid" => ["1//4", "1//4", "64//1"], "n0_grid" => ns, "points" => rows,
    "matches_BETA0" => agree, "seconds" => time() - start,
    "metadata" => model_metadata(C.ROOT, Dict("kind" => "floating-point seed scan");
        packages = (CertifiedGaps,))))
println("$agree of $(length(rows)) best seeds match BETA0")
