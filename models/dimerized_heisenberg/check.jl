# Checker for the dimerized Heisenberg certificates.
#
# Fast (about 10 s): from the stored files it
#   1. re-derives every exact RP identity of results/rp_audit.json at all six deltas;
#   2. for every δ of results/summary.json: checks the frame identity and the model-local source
#      hashes; for each core certificate checks its SHA-256, that its model and RP data are
#      `rotated_model(δ)` and `rotated_rp(δ)` with the seed of `SEEDS`, and runs
#      `GapCertificates.check_certificate(cert; spectra = false)` (exact RP identity and exact
#      rational chain); checks every all-β exclusion from its stored lower bounds (tiling of
#      [1, B], the tail, every box rejected by the recurrence), the dimer perturbation bound and
#      the bound-to-gap ratio.
# `check.jl --full` recomputes the verified spectra of the rings of 2, 4 and 8 cells from the
# model and repeats everything from them alone: `check_certificate(cert; spectra = true)`
# re-derives the seed defects, every exclusion bound is recomputed, and the stored 16-spin gap
# interval must contain the fresh enclosure.
# Every check throws on failure, so the tests can feed it corrupted records.
using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, LinearAlgebra, Test
using CertifiedGaps.RingModels, CertifiedGaps.VerifiedSpectra, CertifiedGaps.GapCertificates
using CertifiedGaps.ExactArithmetic: parse_rational
isdefined(@__MODULE__, :DimerizedHeisenberg) ||
    include(joinpath(@__DIR__, "src", "DimerizedHeisenberg.jl"))
isdefined(@__MODULE__, :DimerizedRP) || include(joinpath(@__DIR__, "src", "rp.jl"))
using .DimerizedHeisenberg, .DimerizedRP
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
asint(x) = parse(Int, string(x))
require(ok, what) = ok === true || error("check failed: $what")

# Exact model-local source set and SHA-256 (the core hashes are provenance only).
check_source_hashes(metadata) = verify_model_metadata(metadata, MODEL_DIR; SOURCES...) > 0

"Re-derive the exact RP audit (reflected form holds, identical-factor form fails)."
function check_audit()
    data = read_result(joinpath(RESULTS_DIR, "rp_audit.json"))
    require(data["exact_algebraic_audit"] === true && data["gap_certificate"] === false,
        "audit scope")
    require([parse_rational(row["delta"]) for row in data["parameters"]] == DELTAS, "audit grid")
    for row in data["parameters"]
        a = audit_rp(parse_rational(row["delta"]))
        require(a.defect_norm2 == parse_rational(row["swap_defect_frobenius_squared"]) > 0,
            "swap defect")
        require(a.real_rhs_error_norm2 ==
                parse_rational(row["proposed_real_rhs_error_frobenius_squared"]), "real form")
        require(row["conjugate_identity_exact"] === a.conjugate_identity === true &&
                row["reflected_real_identity_exact"] === a.reflected_real_identity === true &&
                row["theta_invariant_onsite_exact"] === true &&
                row["real_identical_factor_identity"] === false, "audit flags")
    end
    return check_source_hashes(data["metadata"])
end

"""
    check_delta(row; full = false) -> true

All claims of one summary row (see the header). With `full = true` the spectra are recomputed
from the model and no stored defect, bound or energy is used.
"""
function check_delta(row; full = false)
    δ = parse_rational(row["delta"])
    model, rp = rotated_model(δ), rotated_rp(δ)
    require(frame_identity(δ), "frame identity")
    spectra = full ? Dict(L => verified_ring_spectrum(model, L) for L in (2, 4, 8)) : nothing
    seeds = [(asint(c["n0"]), parse_rational(c["beta0"])) for c in row["certificates"]]
    require(seeds == get(SEEDS, δ, Tuple{Int,Rational{BigInt}}[]), "seed list")
    for c in row["certificates"]
        path = joinpath(@__DIR__, c["file"])
        require(sha256_file(path) == c["sha256"], "certificate hash")
        cert = read_result(path)
        check_source_hashes(cert["metadata"])
        require(encode_model(decode_model(cert["model"])) == encode_model(model), "model")
        require(encode_rp(decode_rp(cert["reflection_positivity"])) == encode_rp(rp), "RP data")
        require((asint(cert["seed"]["n0"]), parse_rational(cert["seed"]["beta0"])) ==
                (asint(c["n0"]), parse_rational(c["beta0"])), "seed (n0, beta0)")
        require(cert["seed"]["eps_s"] == c["eps_s"] && cert["seed"]["eps_t"] == c["eps_t"] &&
                cert["gap_lower_bound"] == c["gap_lower_bound"], "summary row")
        # --full hands the core the spectra just computed above from the same model.
        n0 = asint(c["n0"])
        fresh = full ? Dict(n0 => spectra[n0], 2n0 => spectra[2n0]) : false
        require(check_certificate(cert; spectra = fresh)["gap_lower_bound"] ==
                c["gap_lower_bound"], "core certificate")
        screen = certified_gap_bound(parse_rational(c["eps_s"]), parse_rational(c["eps_t"]),
            parse_rational(c["beta0"])).bound
        require(screen == cert["interval_bound_decimal"] == c["interval_bound_decimal"],
            "interval bound")
    end
    excluded = [asint(e["n0_cells"]) for e in row["negative_results"]]
    require(sort(vcat(unique(first.(seeds)), excluded)) == [2, 4], "each n0 classified once")
    for e in row["negative_results"]
        n = asint(e["n0_cells"])
        full ? check_cover(e; small = spectra[n], large = spectra[2n]) : check_cover(e)
    end
    gl, gh = parse_rational.(row["reference_gap_16_interval"])
    require(0 < gl <= gh, "16-spin gap interval")
    if full
        lo, hi = spectral_gap_enclosure(spectra[8])
        require(gl <= lo <= hi <= gh, "fresh 16-spin gap")
    end
    if isempty(seeds)
        require(row["bound_to_gap16_interval"] === nothing, "ratio without a bound")
    else
        best = maximum(Rational{BigInt}(c["interval_bound_decimal"]) for c in row["certificates"])
        rl, rh = parse_rational.(row["bound_to_gap16_interval"])
        require(rl <= best / gh <= best / gl <= rh, "bound-to-gap ratio")
    end
    dimer, d = row["dimer_perturbation"], dimer_bounds(δ, 4, DIMER_BETA(δ))
    endpoints(x) = string.([Rational{BigInt}(inf(x)), Rational{BigInt}(sup(x))])
    require(dimer["norm_small"] == string(d.norm_small) &&
            dimer["eps_s_enclosure"] == endpoints(d.s) &&
            dimer["eps_t_enclosure"] == endpoints(d.t), "dimer perturbation")
    return true
end

"""
    check_summary(; full = false, summary = stored) -> true

Check results/summary.json and every certificate it indexes, from the stored files or
(`full = true`) from freshly verified spectra.
"""
function check_summary(; full = false,
    summary = read_result(joinpath(RESULTS_DIR, "summary.json")))
    check_source_hashes(summary["metadata"])
    require(summary["format"] == "dimerized-heisenberg-summary-v3", "summary format")
    require([parse_rational(r["delta"]) for r in summary["parameters"]] == DELTAS, "grid")
    for row in summary["parameters"]
        seconds = @elapsed check_delta(row; full)
        bounds = [c["interval_bound_decimal"] for c in row["certificates"]]
        println("checked delta=$(row["delta"]) full=$full bounds=$bounds ",
            "exclusions=$(length(row["negative_results"])) seconds=$(round(seconds; digits=1))")
    end
    # Every smaller grid delta has no seed and its n0 = 2, 4 exclusions were checked above.
    positive = [parse_rational(r["delta"]) for r in summary["parameters"]
                if !isempty(r["certificates"])]
    require(minimum(positive) == parse_rational(summary["smallest_certified_grid_delta"]) &&
            all(δ -> δ in positive, filter(>=(minimum(positive)), DELTAS)), "grid minimum")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Dimerized Heisenberg certificates" begin
        @test check_audit()
        @test check_summary()
        "--full" in ARGS && @test check_summary(; full = true)
    end
end
