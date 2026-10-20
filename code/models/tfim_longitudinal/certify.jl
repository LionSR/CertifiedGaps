# Produce the per-point certificates results/certificate_g<g>_h<h>.json.gz and
# results/seeds.json.
#
# For each grid point, the core `certify_point` verifies the exact RP identity −h = Z ⊗ Z,
# encloses the complete spectra of the rings of 6 and 12 sites (momentum blocks, real by the
# lattice reflection), bounds the seed defects at (n₀, β₀) = (6, BETA0[(g, h)]), runs the exact
# rational doubling chain and writes the certificate with its params (g, h) and run metadata.
# The seeds β₀ come from the floating-point scan `scripts/scan.jl`. Usage (about 1 minute for
# all 30 points):
#     julia --project=julia models/tfim_longitudinal/certify.jl \
#         [point indices 1-30, g-major order]
include(joinpath(@__DIR__, "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, Printf
using CertifiedGaps.ExactArithmetic: parse_rational
const C = LongitudinalCertificate

points = C.points()
!isempty(ARGS) && (points = points[parse.(Int, ARGS)])
for (g, h) in points
    seconds = @elapsed cert = C.certify(g, h)
    s = cert["seed"]
    @printf("g=%s h=%s beta0=%s eps_s=%.6g eps_t=%.6g gap > %s (interval %.6g) %.1fs\n",
        g, h, s["beta0"], Float64(parse_rational(s["eps_s"])),
        Float64(parse_rational(s["eps_t"])), cert["gap_lower_bound_decimal"],
        cert["interval_bound_decimal"], seconds)
    flush(stdout)
end
# The seed table is rebuilt from every stored certificate, not only the points just run.
length(C.result_paths()) == length(C.points()) && C.write_seed_table()
