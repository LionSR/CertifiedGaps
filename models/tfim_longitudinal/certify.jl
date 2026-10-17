# Produce the per-point core certificates results/certificate_g<g>_h<h>.json.gz and
# results/seeds.json.
#
# For each grid point, `GapCertificates.certify_gap` verifies the exact RP identity
# −h = Z ⊗ Z, encloses the complete spectra of the rings of 6 and 12 sites (core
# `verified_ring_spectrum`), bounds the seed defects at (n₀, β₀) = (6, BETA0[(g, h)]) and
# runs the exact rational doubling chain. The seeds β₀ come from the floating-point scan
# `scripts/scan.jl`. Usage (about 1 minute for all 30 points):
#     julia --project=julia models/tfim_longitudinal/certify.jl \
#         [point indices 1-30, g-major order]
include(joinpath(@__DIR__, "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, CertifiedGaps.Results, Printf
using CertifiedGaps.ExactArithmetic: parse_rational
const C = LongitudinalCertificate

points = C.points()
!isempty(ARGS) && (points = points[parse.(Int, ARGS)])
for (g, h) in points
    seconds = @elapsed cert = C.certify_point(g, h)
    write_result(C.certificate_path(g, h), cert; level = 9)
    s = cert["seed"]
    @printf("g=%s h=%s beta0=%s eps_s=%.6g eps_t=%.6g gap > %s (interval %.6g) %.1fs\n",
        C.str(g), C.str(h), s["beta0"], Float64(parse_rational(s["eps_s"])),
        Float64(parse_rational(s["eps_t"])), cert["gap_lower_bound_decimal"],
        cert["interval_bound_decimal"], seconds)
    flush(stdout)
end
# The seed table is rebuilt from every stored certificate, not only the points just run.
length(C.result_paths()) == length(C.points()) && C.write_seed_table()
