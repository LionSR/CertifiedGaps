# Checker for the dimerized spin-1 certificates results/certificate_d9.json, certificate_d7.json.
#
#     julia --project=julia models/spin1_dimerized/check.jl [--full]
#
# Fast mode (default): per certificate, `GapCertificates.check_point` against the point
# definition (recorded model-local source hashes, stored model, RP data, seed, params and
# positive-transfer reference, the Lean seed defects of δ = 9/10), the exact RP identity and the
# exact rational chain. `--full` also re-encloses the complete spectra of the rings of n₀ and
# 2n₀ cells with the core (for δ = 7/10 the 6-cell ring, as long as `certify.jl 7`) and
# requires the fresh seed defects not to exceed the stored εₛ, εₜ. Any failure throws.
include(joinpath(@__DIR__, "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate
const D = DimerizedCertificate

if abspath(PROGRAM_FILE) == @__FILE__
    full = "--full" in ARGS
    paths = D.result_paths()
    basename.(paths) == [basename(D.certificate_path(δ)) for δ in sort(collect(keys(D.POINTS)))] ||
        error("the stored certificates are not exactly the certified points")
    for path in paths
        r = D.check(path; full)
        println(basename(path), ": gap > ", r["gap_lower_bound"], " (interval ",
            r["interval_bound_decimal"], ")", full ? "; fresh defects ≤ " *
            "$(Float64(r["defect_s_upper"])), $(Float64(r["defect_t_upper"]))" : "")
    end
    println("All $(length(paths)) certificates verified", full ? " (full replay)." : ".")
end
