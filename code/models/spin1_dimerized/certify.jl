# Produce the core-format certificate results/certificate_d<d>.json of δ = d/10.
#
# `GapCertificates.certify_point` verifies the exact reflected RP identity, encloses the complete
# spectra of the rings of n₀ and 2n₀ cells (charge sectors ±q paired by the spin flip, real
# momentum blocks by the lattice reflection), bounds the seed defects at the point's (n₀, β₀),
# runs the exact rational doubling chain and writes the certificate. The odd seed n₀ = 3 of
# δ = 7/10 is passed with the `positive_transfer` reference (X_β = K K* ≥ 0, METHOD.md). Usage:
#     julia --project=julia models/spin1_dimerized/certify.jl 9   # about 25 s
#     julia --project=julia models/spin1_dimerized/certify.jl 7   # 6-cell ring, see README.md
# Without arguments every point of `POINTS` is certified.
include(joinpath(@__DIR__, "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, Printf
using CertifiedGaps.ExactArithmetic: parse_rational
const D = DimerizedCertificate

if abspath(PROGRAM_FILE) == @__FILE__
    deltas = isempty(ARGS) ? sort(collect(keys(D.POINTS)); rev = true) :
             [parse(Int, a) // 10 for a in ARGS]
    for δ in deltas
        cert = D.certify(δ; progress = true)
        s = cert["seed"]
        @printf("delta=%s n0=%s beta0=%s eps_s=%.12g eps_t=%.12g gap > %s (interval %.10g) %.1fs\n",
            D.str(δ), s["n0"], s["beta0"], Float64(parse_rational(s["eps_s"])),
            Float64(parse_rational(s["eps_t"])), cert["gap_lower_bound_decimal"],
            cert["interval_bound_decimal"], cert["metadata"]["config"]["seconds"])
    end
end
