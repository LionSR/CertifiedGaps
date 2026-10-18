# Core gap certificate of the two-leg TFIM ladder at every point of `Ladder.POINTS` (or at the
# field given as the first argument, with that point's seed): GapCertificates.certify_point on
# the rotated-frame RingModel (exact RP identity, verified spectra of the rings of n₀ and 2n₀
# rungs, log Z enclosures, seed defects, exact rational chain, params and source hashes).
# Writes results/certificate_g<g>.json.
#     julia --project=julia models/two_leg_tfim/certify.jl [g]     (run time in README.md)
include(joinpath(@__DIR__, "src", "Ladder.jl"))
using .Ladder, CertifiedGaps.GapCertificates
using LinearAlgebra
BLAS.set_num_threads(1)
const ROOT = @__DIR__

points = isempty(ARGS) ? POINTS : filter(p -> p.g == parse(Rational{Int}, ARGS[1]), POINTS)
isempty(points) && error("no certified point at g = $(ARGS[1]); see Ladder.POINTS")
decimal(x) = Float64(parse(Rational{BigInt}, x))
for p in points
    model = ladder(p.g)
    path = point_path(ROOT, p.g)
    cert = certify_point(path, model; n0 = p.n0, β0 = p.β0, params = (; g = p.g),
        model_dir = ROOT, rp = ladder_rp(model), progress = true)
    row = certificate_row(cert)
    println(path, ": g = $(p.g), εₛ ≤ $(decimal(row["eps_s"])), εₜ ≤ $(decimal(row["eps_t"])),",
        " gap > $(row["gap_lower_bound_decimal"]),",
        " interval bound $(row["interval_bound_decimal"]),",
        " $(round(cert["metadata"]["config"]["seconds"]; digits = 1)) s")
end
