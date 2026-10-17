# Produce the core-format certificate results/certificate_d<d>.json of δ = d/10.
#
# `GapCertificates.certify_gap` verifies the exact reflected RP identity, encloses the complete
# spectra of the rings of n₀ and 2n₀ cells (core `verified_ring_spectrum`, charge sectors ±q
# paired by the spin flip), bounds the seed defects at the point's (n₀, β₀) and runs the exact
# rational doubling chain. The odd seed n₀ = 3 of δ = 7/10 is passed with the core
# `positive_transfer` reference (X_β = K K* ≥ 0, METHOD.md). Usage:
#     julia --project=julia models/spin1_dimerized/certify.jl 9   # about 25 s
#     julia --project=julia models/spin1_dimerized/certify.jl 7   # 6-cell ring
# Without arguments every point of `POINTS` is certified.
include(joinpath(@__DIR__, "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps.Results, Printf
using CertifiedGaps.ExactArithmetic: parse_rational
const D = DimerizedCertificate

if abspath(PROGRAM_FILE) == @__FILE__
    deltas = isempty(ARGS) ? sort(collect(keys(D.POINTS)); rev = true) :
             [parse(Int, a) // 10 for a in ARGS]
    for δ in deltas
        seconds = @elapsed cert = D.certify_point(δ; progress = true)
        cert["seconds"] = seconds
        write_result(D.certificate_path(δ), cert)
        s = cert["seed"]
        @printf("delta=%s n0=%s beta0=%s eps_s=%.10g eps_t=%.10g gap > %s (interval %.10g) %.1fs\n",
            D.str(δ), s["n0"], s["beta0"], Float64(parse_rational(s["eps_s"])),
            Float64(parse_rational(s["eps_t"])), cert["gap_lower_bound_decimal"],
            cert["interval_bound_decimal"], seconds)
    end
end
