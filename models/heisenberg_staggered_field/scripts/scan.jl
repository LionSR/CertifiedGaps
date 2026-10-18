# Numerical seed scan (floating point, not certified) on rings up to maxL sites, with a β grid
# written to results/scan_<h>_<maxL>.json (gitignored). Up to 16 sites the spectra are complete:
# unverified LAPACK eigenvalues of the certificate's blocks (`float_ring_spectrum` of the blocked
# cell model; the one-site model for the 2-site ring, a single cell). Larger rings, whose dense
# blocks are too big, use the core exact diagonalization of the physical chain in Sᶻ sectors
# (the lowest 80 Lanczos levels above dimension 2500). The certified points are in src/common.jl.
#     julia --project=julia models/heisenberg_staggered_field/scripts/scan.jl 16 [h]
include(joinpath(dirname(@__DIR__), "src", "common.jl"))
using CertifiedGaps.Seeds
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
using CertifiedGaps.ExactDiagonalization: spectrum
maxL = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 16
fields = length(ARGS) > 1 ? [parse_rational(ARGS[2])] : FIELDS
for h in fields
    cell = cell_model(site_model(h), BLOCKING.cell; charge = BLOCKING.charge)
    chain = CertifiedGaps.Models.staggered_heisenberg(h).chain
    ring(L) = L == 2 ? float_ring_spectrum(site_model(h), 2) :
              L <= 16 ? float_ring_spectrum(cell, L ÷ 2) :
              spectrum(chain, L; dense_max = 2500, k = 80)
    lengths = sort!(unique(vcat(collect(2:2:maxL÷2), collect(4:4:maxL))))
    spectra, timings = Dict{Int,Any}(), Dict{Int,Float64}()
    for L in lengths
        timings[L] = @elapsed spectra[L] = ring(L)
    end
    rows = scan_seeds(spectral_defects(spectra), 2:2:maxL÷2, 0.25:0.25:32)
    best = best_seed(rows)
    # Per n₀, the purity peak of the grid and the ground-energy doubling mismatch δ(n₀).
    table = map(purity_peaks(rows)) do r
        n = r["n0"]
        merge(r, Dict("delta" => spectra[2n].E[1] - 2spectra[n].E[1]))
    end
    top = spectra[maximum(lengths)].E
    write_result(joinpath(ROOT, "results", "scan_$(tag(h))_$(maxL).json"), Dict(
        "status" => "floating-point scan, not certified", "h" => string(h),
        "best" => best === nothing ? nothing : seed_record(best), "per_n" => table,
        "timings" => timings, "reference_gap" => top[2] - top[1],
        "grid" => seed_record.(rows)))
    println("BEST h=$h ", best === nothing ? nothing : seed_record(best))
    flush(stdout)
end
