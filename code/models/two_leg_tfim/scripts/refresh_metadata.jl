# Refreshes the provenance metadata of the stored complete-spectrum certificates
# (`Ladder.POINTS`, results/certificate_g<g>.json) without changing any certified number.
#
#     julia --project=code/julia code/models/two_leg_tfim/scripts/refresh_metadata.jl [g]
#
# Why: the certificates pin the model-local sources `RitzSeeds.SOURCES` by SHA-256, so editing
# certify.jl, check.jl or src/ invalidates them, yet a fresh certify.jl run does not reproduce
# the stored numbers bit for bit (the verified radii of LAPACK eigenpairs vary from run to run,
# and with them the log Z enclosures and the seed defects). The repository rule is that a
# certificate's numbers stay identical, so this script reuses the stored seed (εₛ, εₜ, the
# log Z enclosures and spectrum summaries of the original certify.jl run, which
# `check.jl --full` recomputes and compares) and re-derives everything else through the core:
# `certify_gap(...; seed)` re-verifies the exact RP identity and recomputes the model and RP
# records, the exact rational chain, the gap bound and the interval bound. It requires that
# every value outside "metadata" equals the stored one exactly (parsed with
# `parse_exact_json`, so integers and rationals compare as strings), and only then rewrites the
# metadata: `Results.model_metadata` with the same `SOURCES` (`dirs`, `files`) as certify.jl,
# the original run configuration (`config`, including its run time) and, under
# `metadata_refresh`, the creation time and versions of the original run.
include(joinpath(@__DIR__, "..", "src", "Ladder.jl"))
include(joinpath(@__DIR__, "..", "src", "RitzSeeds.jl"))
using .Ladder, .RitzSeeds, CertifiedGaps, CertifiedGaps.GapCertificates
using CertifiedGaps.Results: read_result, write_result, model_metadata, parse_exact_json
using CertifiedGaps.ExactArithmetic: parse_rational
using Dates, JSON3
const ROOT = dirname(@__DIR__)

# The JSON text of `data` as write_result stores it, parsed back exactly.
exact(data) = parse_exact_json(sprint(JSON3.pretty, data))
without_metadata(d) = Dict(k => v for (k, v) in d if k != "metadata")

points = isempty(ARGS) ? POINTS : filter(p -> p.g == parse(Rational{Int}, only(ARGS)), POINTS)
isempty(points) && error("no complete-spectrum point selected; see Ladder.POINTS")
for p in points
    path = point_path(ROOT, p.g)
    stored = read_result(path, Dict{String,Any})   # integers stay integers on rewrite
    s = stored["seed"]
    # A complete-spectrum seed without blocking or positive transfer, as certify.jl writes it.
    (haskey(s, "route") || haskey(s, "positive_transfer") || haskey(stored, "blocking")) &&
        error("$(basename(path)): unexpected seed record")
    seed = (; n0 = p.n0, β0 = parse_rational(s["beta0"]), blocking = nothing,
        positive_transfer = "", eps_s = parse_rational(s["eps_s"]),
        eps_t = parse_rational(s["eps_t"]), log_partition = stored["log_partition"],
        spectra = stored["spectra"])
    model = ladder(p.g)
    cert = certify_gap(model; n0 = p.n0, β0 = p.β0, params = (; g = p.g), rp = ladder_rp(model),
        seed)
    exact(cert) == without_metadata(exact(stored)) ||
        error("$(basename(path)): the re-derived certificate differs from the stored one")
    old = stored["metadata"]
    metadata = model_metadata(ROOT, old["config"]; SOURCES...,
        packages = (CertifiedGaps,), certified = true)
    metadata["metadata_refresh"] = Dict{String,Any}(
        "script" => "scripts/refresh_metadata.jl",
        "note" => "numbers of the original certify.jl run, unchanged; only the source and " *
                  "core hashes, versions and creation time were refreshed",
        "original_created_utc" => old["created_utc"],
        "original_julia" => old["julia"],
        "original_CertifiedGaps" => old["CertifiedGaps"])
    cert["metadata"] = metadata
    write_result(path, cert)
    without_metadata(exact(read_result(path, Dict{String,Any}))) ==
        without_metadata(exact(stored)) || error("$(basename(path)): rewrite changed a number")
    println(basename(path), ": metadata refreshed (", length(metadata["source_sha256"]),
        " model-local sources); every value outside \"metadata\" is unchanged; gap > ",
        cert["gap_lower_bound_decimal"], ", interval bound ", cert["interval_bound_decimal"])
end
