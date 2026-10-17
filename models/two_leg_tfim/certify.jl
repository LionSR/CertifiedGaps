# Core gap certificate of the two-leg TFIM ladder at every point of `Ladder.POINTS` (or at the
# field given as the first argument, with that point's seed): CertifiedGaps.certify_gap on the
# rotated-frame RingModel (exact RP identity, verified spectra of the rings of n₀ and 2n₀ rungs,
# log Z enclosures, seed defects, exact rational chain). Writes results/certificate_g<g>.json
# with `Results.model_metadata` (model-local and core source hashes) under "metadata".
#     julia --project=julia models/two_leg_tfim/certify.jl [g]     (about 9 min)
include(joinpath(@__DIR__, "src", "Ladder.jl"))
using .Ladder, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.GapCertificates
using LinearAlgebra
BLAS.set_num_threads(1)
const ROOT = @__DIR__

points = isempty(ARGS) ? POINTS : filter(p -> p.g == parse(Rational{Int}, ARGS[1]), POINTS)
isempty(points) && error("no certified point at g = $(ARGS[1]); see Ladder.POINTS")
for p in points
    model = ladder(p.g)
    seconds = @elapsed cert = certify_gap(model; n0 = p.n0, β0 = p.β0, rp = ladder_rp(model),
        progress = true)
    cert["point"] = Dict("g" => string(p.g), "lengths" => "even L >= $(p.n0) rungs")
    cert["metadata"] = model_metadata(ROOT, (; g = string(p.g), n0 = p.n0,
            beta0 = string(p.β0), blas = string(BLAS.get_config()), seconds);
        PROVENANCE..., packages = (CertifiedGaps,), certified = true)
    tag = replace(string(p.g), "//" => "_")
    path = write_result(joinpath(ROOT, "results", "certificate_g$tag.json"), cert)
    println(path, ": g = $(p.g), εₛ ≤ $(Float64(parse(Rational{BigInt}, cert["seed"]["eps_s"]))),",
        " εₜ ≤ $(Float64(parse(Rational{BigInt}, cert["seed"]["eps_t"]))),",
        " gap > $(cert["gap_lower_bound_decimal"]), $(round(seconds; digits = 1)) s")
end
