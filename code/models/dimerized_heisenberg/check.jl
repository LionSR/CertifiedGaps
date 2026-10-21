# Checker for the dimerized Heisenberg certificates.
#
# Fast (default): from the stored files it
#   1. re-derives every exact RP identity of results/rp_audit.json at all six deltas;
#   2. for every δ of results/summary.json: checks the frame identity and the model-local source
#      hashes; for each certificate checks its SHA-256 and summary row and runs the core
#      `check_point` (exactly the model `rotated_model(δ)`, RP data `rotated_rp(δ)`, seed of
#      `SEEDS` and params, the defects of the Lean point when the Lean file has it, the exact RP
#      identity, the exact rational chain and the interval bound); checks every all-β exclusion
#      with the core `check_exclusion` (tiling of [0, B], the tail, every window rejected), the
#      dimer perturbation bound and the bound-to-gap ratio.
# `check.jl --full` recomputes the verified spectra of the rings of 2, 4 and 8 cells from the
# model and repeats everything from them alone: the seed defects, every exclusion bound, and the
# 16-spin gap interval, which must contain the fresh enclosure.
# Every check throws on failure, so the tests can feed it corrupted records.
using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, LinearAlgebra, Test
using CertifiedGaps.RingModels, CertifiedGaps.VerifiedSpectra, CertifiedGaps.GapCertificates
using CertifiedGaps.Exclusions: check_exclusion
using CertifiedGaps.ExactArithmetic: parse_rational
isdefined(@__MODULE__, :DimerizedHeisenberg) ||
    include(joinpath(@__DIR__, "src", "DimerizedHeisenberg.jl"))
isdefined(@__MODULE__, :DimerizedRP) || include(joinpath(@__DIR__, "src", "rp.jl"))
using .DimerizedHeisenberg, .DimerizedRP
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const LEAN = joinpath(REPO, "lean", "CertifiedGaps", "Models", "DimerizedHeisenberg.lean")
asint(x) = parse(Int, string(x))
require(ok, what) = ok === true || error("check failed: $what")

"The Lean seed defects `(; s, t)` of the point (δ, β₀), or `nothing` if the Lean file lacks it."
function lean_seed(δ, β0)
    hits = [p for p in lean_points(LEAN) if p.params == [δ, β0]]
    return isempty(hits) ? nothing : only(hits)
end

"Every Lean point is a certified seed (δ, 4, β₀) of `SEEDS`."
check_lean_points() = all(p -> length(p.params) == 2 &&
    (4, p.params[2]) in get(SEEDS, p.params[1], ()), lean_points(LEAN))

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
    return verify_model_metadata(data["metadata"], MODEL_DIR) > 0
end

# A summary row in a form comparable with `certificate_row` of the stored certificate.
_row(r) = Dict(k => string(v) for (k, v) in r if !(k in ("file", "sha256")))

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
    for (c, (n0, β0)) in zip(row["certificates"], seeds)
        path = joinpath(@__DIR__, c["file"])
        require(sha256_file(path) == c["sha256"], "certificate hash")
        require(_row(certificate_row(read_result(path))) == _row(c), "summary row")
        # --full hands the core the spectra just computed above from the same model.
        fresh = full ? Dict(n0 => spectra[n0], 2n0 => spectra[2n0]) : false
        check_point(path, model; n0, β0, params = point_params(δ), model_dir = MODEL_DIR, rp,
            lean = lean_seed(δ, β0), spectra = fresh)
    end
    excluded = [asint(e["n0"]) for e in row["negative_results"]]
    require(sort(vcat(unique(first.(seeds)), excluded)) == [2, 4], "each n0 classified once")
    for (e, n) in zip(row["negative_results"], excluded)
        full ? check_exclusion(e; small = spectra[n], large = spectra[2n]) : check_exclusion(e)
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
    verify_model_metadata(summary["metadata"], MODEL_DIR)
    require(summary["format"] == "dimerized-heisenberg-summary-v4", "summary format")
    require([parse_rational(r["delta"]) for r in summary["parameters"]] == DELTAS, "grid")
    require(check_lean_points(), "every Lean point has a certificate")
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
