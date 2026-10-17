# Checker for the dimerized spin-1 certificates results/certificate_d9.json, certificate_d7.json.
#
#     julia --project=julia models/spin1_dimerized/check.jl [--full]
#
# Fast mode (default): per certificate, the recorded model-local source hashes, the stored
# model, RP data, seed and positive-transfer reference against the point definition, the Lean
# seed values, the exact RP identity and the exact rational chain (core
# `check_certificate(cert; spectra = false)`). `--full` also re-encloses the complete spectra of
# the rings of n₀ and 2n₀ cells with the core (for δ = 7/10 the 6-cell ring, about an hour and
# 15 GB) and requires the fresh seed defects not to exceed the stored εₛ, εₜ. Any failure throws.
include(joinpath(@__DIR__, "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate
const D = DimerizedCertificate

if abspath(PROGRAM_FILE) == @__FILE__
    full = "--full" in ARGS
    paths = D.result_paths()
    basename.(paths) == [basename(D.certificate_path(δ)) for δ in sort(collect(keys(D.POINTS)))] ||
        error("the stored certificates are not exactly the certified points")
    for path in paths
        r = D.check_point(path; full)
        println(basename(path), ": gap > ", r["gap_lower_bound"],
            full ? "  fresh defects ≤ $(Float64(r["defect_s_upper"])), " *
                   "$(Float64(r["defect_t_upper"]))" : "")
    end
    println("All $(length(paths)) certificates verified", full ? " (full replay)." : ".")
end
