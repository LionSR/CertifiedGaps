# Checker for the 30 longitudinal-field TFIM certificates in results/certificate_g*_h*.json.gz.
#
#     julia --project=julia models/tfim_longitudinal/check.jl [--full]
#
# Fast mode (default): per certificate, the recorded model-local source hashes, the stored
# model, RP factor and seed against the grid definition, the exact RP identity and the exact
# rational chain (core `check_certificate(cert; spectra = false)`); then the seed table.
# `--full` also re-encloses the complete spectra of the rings of 6 and 12 sites with the core
# and requires the fresh seed defects not to exceed the stored εₛ, εₜ
# (`check_certificate(cert; spectra = true)`). Any failure throws.
include(joinpath(@__DIR__, "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, CertifiedGaps.Results
const C = LongitudinalCertificate

full = "--full" in ARGS
paths = C.result_paths()
C.check_grid(paths)
for path in paths
    r = C.check_point(path; full)
    println(basename(path), ": gap > ", r["gap_lower_bound"],
        full ? "  fresh defects ≤ $(Float64(r["defect_s_upper"])), " *
               "$(Float64(r["defect_t_upper"]))" : "")
end
table = read_result(joinpath(C.ROOT, "results", "seeds.json"))
table["points"] == parse_exact_json(Results.JSON3.write(C.seed_rows())) ||
    error("results/seeds.json does not match the certificates")
println("All $(length(paths)) certificates verified", full ? " (full replay)." : ".")
