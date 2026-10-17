# Floating-point seed scan (not a certificate): for every grid point, the complete spectra of
# the rings of 2 to 12 sites from the core momentum blocks (`VerifiedSpectra.ring_blocks`,
# block centers diagonalized in Float64), the seed defects over n₀ ∈ {2, 4, 6} and
# β₀ ∈ {1/4, …, 64}, and the seed with the best `certified_gap_bound`. Writes
# results/scan.json and reports whether each best β₀ agrees with the table `BETA0` that
# certify.jl uses.
#     julia --project=julia models/tfim_longitudinal/scripts/scan.jl
include(joinpath(@__DIR__, "..", "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, LinearAlgebra, Printf
using CertifiedGaps.Results, CertifiedGaps.Seeds
using CertifiedGaps.VerifiedSpectra: ring_blocks
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
const C = LongitudinalCertificate

"Floating complete spectrum of the ring of L sites from the core momentum blocks."
function float_spectrum(model, L)
    E = Float64[]
    ring_blocks(model, L) do M, _, mult, _
        for _ in 1:mult
            append!(E, eigvals(Hermitian(M)))
        end
    end
    PartialSpectrum(E, 2^L)
end

start = time(); ns = [2, 4, 6]; betas = collect(0.25:0.25:64); rows = Any[]; agree = 0
for (g, h) in C.points()
    model = C.longitudinal(g, h)
    spectra = Dict(L => float_spectrum(model, L) for L in sort(unique(vcat(ns, 2ns))))
    best = best_seed(scan_seeds(spectral_defects(spectra), ns, betas))
    gap12 = spectra[12].E[2] - spectra[12].E[1]
    same = best !== nothing && best.n₀ == C.N0 && Rational{Int}(best.β₀) == C.BETA0[(g, h)]
    global agree += same
    push!(rows, Dict("g" => C.str(g), "h" => C.str(h), "gap_L12" => gap12,
        "best_n0" => best.n₀, "best_beta0" => best.β₀, "epsilon_s" => best.εₛ,
        "epsilon_t" => best.εₜ, "float_bound" => best.bound, "matches_BETA0" => same))
    @printf("g=%s h=%s best (n0, beta0) = (%d, %g) bound %.4f gap(12) %.4f %s\n", C.str(g),
        C.str(h), best.n₀, best.β₀, best.bound, gap12, same ? "" : "differs from BETA0")
    flush(stdout)
end
write_result(joinpath(C.ROOT, "results", "scan.json"), Dict(
    "status" => "floating-point diagnostics, not a certificate",
    "beta_grid" => ["1//4", "1//4", "64//1"], "n0_grid" => ns, "points" => rows,
    "matches_BETA0" => agree, "seconds" => time() - start,
    "metadata" => C.metadata(Dict("kind" => "floating-point seed scan"))))
println("$agree of $(length(rows)) best seeds match BETA0")
