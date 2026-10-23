# Core gap certificates of the two-leg TFIM ladder (rotated-frame RingModel `Ladder.ladder(g)`).
#
#   * `Ladder.POINTS` (complete-spectrum seeds, n₀ = 4): GapCertificates.certify_point (exact RP
#     identity, verified spectra of the rings of n₀ and 2n₀ rungs, log Z enclosures, seed
#     defects, exact rational chain, params and source hashes); results/certificate_g<g>.json.
#   * `RitzSeeds.RITZ_POINTS` (coupled seeds, n₀ = 6): the complete verified spectrum of the
#     6-rung ring and rigorous Rayleigh–Ritz upper bounds on the low levels of every block of
#     the 12-rung ring (core `RitzBounds`), combined by `coupled_seed_bounds` and passed to
#     certify_point as its seed; results/certificate_g<g>_n6.json and the stored bounds
#     results/partial_g<g>_L12.json. They also store `finite_bridge`: the verified gap of the
#     4-rung ring, so that together with the criterion (every even L ≥ 6) they cover every
#     even L ≥ 4.
#
# Every certificate pins the model-local sources `RitzSeeds.SOURCES` (src/Ladder.jl,
# src/RitzSeeds.jl, certify.jl, check.jl); the list of coupled points, src/ritz_points.jl, is not
# pinned, so that adding a point leaves the other certificates valid.
#
#     julia --project=code/julia code/models/two_leg_tfim/certify.jl [g] [--ritz | --complete]
#
# With a field g only that field's points are certified; `--ritz` or `--complete` restricts to
# one kind (run times in README.md).
include(joinpath(@__DIR__, "src", "Ladder.jl"))
include(joinpath(@__DIR__, "src", "RitzSeeds.jl"))
using .Ladder, .RitzSeeds, CertifiedGaps.GapCertificates
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
using CertifiedGaps.RitzBounds: write_partial_spectrum
using CertifiedGaps.Results: sha256_file
using LinearAlgebra
BLAS.set_num_threads(1)
const ROOT = @__DIR__

fields = filter(a -> !startswith(a, "--"), ARGS)
selected(p) = isempty(fields) || p.g == parse(Rational{Int}, only(fields))
complete = filter(selected, "--ritz" in ARGS ? [] : POINTS)
coupled = filter(selected, "--complete" in ARGS ? [] : RITZ_POINTS)
isempty(complete) && isempty(coupled) && error("no point selected; see Ladder.POINTS and " *
                                              "RitzSeeds.RITZ_POINTS")
decimal(x) = Float64(parse(Rational{BigInt}, x))
function report(path, p, cert)
    row = certificate_row(cert)
    println(path, ": g = $(p.g), (n0, β0) = ($(p.n0), $(p.β0)), εₛ ≤ $(decimal(row["eps_s"])),",
        " εₜ ≤ $(decimal(row["eps_t"])), gap > $(row["gap_lower_bound_decimal"]),",
        " interval bound $(row["interval_bound_decimal"]),",
        " $(round(cert["metadata"]["config"]["seconds"]; digits = 1)) s")
end
for p in complete
    model = ladder(p.g)
    path = point_path(ROOT, p.g)
    cert = certify_point(path, model; n0 = p.n0, β0 = p.β0, params = (; g = p.g),
        model_dir = ROOT, sources = SOURCES, rp = ladder_rp(model), progress = true)
    report(path, p, cert)
end
for p in coupled
    model = ladder(p.g)
    path = ritz_point_path(ROOT, p.g)
    seconds = @elapsed begin
        small = verified_ring_spectrum(model, p.n0)
        large = doubled_ring_bounds(p.g; progress = true)
    end
    partial = partial_spectrum_path(ROOT, p.g)
    write_partial_spectrum(partial, large)
    seed = coupled_seed_bounds(model, p.n0, p.β0; small, large)
    cert = certify_point(path, model; n0 = p.n0, β0 = p.β0, params = (; g = p.g),
        model_dir = ROOT, sources = SOURCES, rp = ladder_rp(model), seed,
        record = Dict{String,Any}("spectra_seconds" => seconds, "partial_spectrum" => Dict(
            "path" => basename(partial), "sha256" => sha256_file(partial)),
            "finite_bridge" => finite_bridge(p.g)))
    bound = max(parse(Rational{BigInt}, cert["gap_lower_bound"]),
        Rational{BigInt}(cert["interval_bound_decimal"]))
    all(b -> parse(Rational{BigInt}, b["gap_lower"]) >= bound, cert["finite_bridge"]) ||
        error("g = $(p.g): a bridged ring has a gap below the certified bound")
    report(path, p, cert)
end
