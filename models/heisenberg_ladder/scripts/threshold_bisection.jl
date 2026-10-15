# Floating-point bisection in J⊥ for the finite-grid seed threshold (NOT a certificate;
# about 20 s). Each trial computes unverified complete spectra of 2, 4 and 8 rungs and scans
# n₀ ∈ {2, 4}, β₀ = 1/4:1/4:16. Sixteen bisections of [3/2, 4] give the rational endpoints
# 326133/131072 (fails) and 163069/65536 (succeeds); the certificate suite then re-decides
# both endpoints rigorously. Monotonicity in J⊥ is not assumed or proved.
# Writes results/scans/threshold_bisection.json.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder \
#         models/heisenberg_ladder/scripts/threshold_bisection.jl
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
using .HeisenbergLadder: BASE, full_spectrum, provenance
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
using CertifiedGaps.Seeds: scan_seeds, spectral_defects, best_seed
using CertifiedGaps.Results: write_result

"Bisect between a failing and a succeeding coupling, recording every trial."
function threshold(lo = 3 // 2, hi = 4 // 1; steps = 16)
    start = time(); trials = Any[]
    function trial(J)
        specs = Dict(L => PartialSpectrum(full_spectrum(L, J).energies, 4^L) for L in (2, 4, 8))
        best = best_seed(scan_seeds(spectral_defects(specs), (2, 4), 0.25:0.25:16))
        push!(trials, Dict("J" => string(J), "succeeds" => best !== nothing,
            "best" => best === nothing ? nothing : Dict("n0" => best.n₀, "beta" => best.β₀,
                "eps_s" => best.εₛ, "eps_t" => best.εₜ, "bound" => best.bound)))
        best !== nothing
    end
    @assert !trial(lo) && trial(hi)
    for _ in 1:steps
        mid = (lo + hi) / 2
        trial(mid) ? (hi = mid) : (lo = mid)
    end
    write_result(joinpath(BASE, "results", "scans", "threshold_bisection.json"), Dict(
        "metadata" => provenance(Dict("ns" => [2, 4], "beta_grid" => "1//4:1//4:16//1")),
        "last_failure" => string(lo), "first_success" => string(hi),
        "interpretation" => "Numerical local transition bracket; monotonicity is not proved.",
        "trials" => trials, "seconds" => time() - start))
    println("NUMERICAL BRACKET ", lo, " ", hi)
end

threshold()
