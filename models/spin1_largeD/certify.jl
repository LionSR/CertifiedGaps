# Stage 1: core gap certificates (CertifiedGaps.GapCertificates.certify_gap) at the two points
# of `POINTS` with seed length n0 = 6, and the verified spectra of the rings L = 6, 12 that
# stage 2 (scripts/transfer_certify.jl) reuses. Writes results/certificate_D*.json and
# results/spectrum_D*_L{6,12}.json.gz (their SHA-256 recorded in the certificate).
#
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/certify.jl
#   ... certify.jl 1603//500          # one point only
using CertifiedGaps, CertifiedGaps.GapCertificates, CertifiedGaps.Results, IntervalArithmetic
using CertifiedGaps.VerifiedSpectra: write_spectrum
include(joinpath(@__DIR__, "src", "model.jl"))
using .LargeDModel

const ROOT = @__DIR__
const SELECTED = isempty(ARGS) ? [p.D for p in POINTS] : [parse(Rational{Int}, a) for a in ARGS]

for (D, beta0) in POINTS
    D in SELECTED || continue
    started = time()
    model = ring_model(D)
    seed = seed_bounds(model, N0, beta0; progress = true)
    cert = certify_gap(model; n0 = N0, β0 = beta0, rp = rp_data(model), seed)
    files = Dict(N0 => spectrum_file(D, N0), 2N0 => spectrum_file(D, 2N0))
    write_spectrum(joinpath(ROOT, files[N0]), seed.small)
    write_spectrum(joinpath(ROOT, files[2N0]), seed.large)
    cert["spectrum_sha256"] = Dict(f => sha256_file(joinpath(ROOT, f)) for f in values(files))
    # Model-local producers only (STAGE1_SOURCES), so stage 2 and the checkers can change
    # without invalidating these hour-long runs.
    cert["metadata"] = model_metadata(ROOT, Dict("D" => string(D), "n0" => N0,
            "beta0" => string(beta0), "wall_seconds" => time() - started);
        STAGE1_SOURCES..., packages = (CertifiedGaps, IntervalArithmetic), certified = true)
    write_result(joinpath(ROOT, certificate_file(D)), cert)
    println("D=$D: eps_s=$(cert["seed"]["eps_s"]) eps_t=$(cert["seed"]["eps_t"]) ",
        "gap > $(cert["interval_bound_decimal"]) (chain $(cert["gap_lower_bound_decimal"])) ",
        "in $(round(time() - started; digits = 1)) s")
end
