# Stage 1: core gap certificates (GapCertificates.certify_point) at the two points of `POINTS`
# with seed length n0 = 6, and the verified spectra of the rings L = 6, 12 that stage 2
# (scripts/transfer_certify.jl) reuses. Writes results/certificate_D*.json and
# results/spectrum_D*_L{6,12}.json.gz (their SHA-256 recorded in the certificate).
#
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/certify.jl
#   ... certify.jl 1603//500          # one point only
using CertifiedGaps, CertifiedGaps.GapCertificates, IntervalArithmetic
using CertifiedGaps.Results: sha256_file
using CertifiedGaps.VerifiedSpectra: write_spectrum
include(joinpath(@__DIR__, "src", "model.jl"))
using .LargeDModel

const ROOT = @__DIR__
const SELECTED = isempty(ARGS) ? [p.D for p in POINTS] : [parse(Rational{Int}, a) for a in ARGS]

for (D, beta0) in POINTS
    D in SELECTED || continue
    model = ring_model(D)
    # The spectra are computed here rather than inside certify_point so that they can be stored
    # for stage 2 and the checker, and their hashes recorded in the certificate.
    seconds = @elapsed seed = seed_bounds(model, N0, beta0; progress = true)
    files = Dict(N0 => spectrum_file(D, N0), 2N0 => spectrum_file(D, 2N0))
    write_spectrum(joinpath(ROOT, files[N0]), seed.small)
    write_spectrum(joinpath(ROOT, files[2N0]), seed.large)
    record = Dict{String,Any}("spectrum_seconds" => seconds, "spectrum_sha256" =>
        Dict(f => sha256_file(joinpath(ROOT, f)) for f in values(files)))
    cert = certify_point(joinpath(ROOT, certificate_file(D)), model; n0 = N0, β0 = beta0,
        params = Dict("D" => D), model_dir = ROOT, sources = STAGE1_SOURCES,
        packages = (IntervalArithmetic,), record, rp = rp_data(model), seed)
    println(certificate_row(cert), " (spectra in $(round(seconds; digits = 1)) s)")
end
