# Checker for the two-leg TFIM certificates (results/certificate_g*.json, core format).
#
#     julia --project=code/julia code/models/two_leg_tfim/check.jl          # fast
#     julia --project=code/julia code/models/two_leg_tfim/check.jl --full   # + spectra
#
# For every point of `Ladder.POINTS` (seed n₀ = 4 from complete spectra) and of
# `RitzSeeds.RITZ_POINTS` (coupled seed n₀ = 6), `GapCertificates.check_point` requires the
# model-local source hashes, exactly the point's rotated-frame model `Ladder.ladder(g)` (charge,
# leg-exchange symmetry, reflection), its RP factors, seed and params, and then re-verifies the
# exact RP identity, the exact rational chain and the interval bound. The n₀ = 4 points must
# also satisfy the Lean seed hypotheses of lean/CertifiedGaps/Models/TwoLegIsing.lean (whose
# theorem uses those seeds). Every point checks the exact Hadamard frame change from the
# physical (Lean) ladder.
#
# Coupled seeds: the fast check recomputes the seed defects from a freshly verified spectrum
# of the 6-rung ring and the stored Rayleigh–Ritz upper bounds of the 12-rung ring
# (results/partial_g*_L12.json, SHA-256 pinned in the certificate); `--full` recomputes those
# bounds from scratch as well. Without `--full` the stored bounds are trusted as produced by
# certify.jl. Complete seeds: `--full` recomputes every verified spectrum. Fresh defects must
# not exceed the stored εₛ, εₜ.
# The coupled certificates also bridge the even lengths below n₀ = 6 (L = 4): the stored and the
# freshly verified gap of the 4-rung ring must be at least the certified bound.
include(joinpath(@__DIR__, "src", "Ladder.jl"))
include(joinpath(@__DIR__, "src", "RitzSeeds.jl"))
using .Ladder, .RitzSeeds, CertifiedGaps.GapCertificates, CertifiedGaps.RingModels
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
using CertifiedGaps.RitzBounds: read_partial_spectrum
using CertifiedGaps.Results: read_result, sha256_file
using LinearAlgebra
BLAS.set_num_threads(1)
const ROOT = @__DIR__
const LEAN = joinpath(ROOT, "..", "..", "lean", "CertifiedGaps", "Models", "TwoLegIsing.lean")

full = "--full" in ARGS
lean = lean_points(LEAN)
# Every stored certificate must belong to a listed point, so none goes unchecked.
stored = filter(f -> startswith(f, "certificate_"), readdir(joinpath(ROOT, "results")))
listed = [[basename(point_path(ROOT, p.g)) for p in POINTS];
          [basename(ritz_point_path(ROOT, p.g)) for p in RITZ_POINTS]]
Set(stored) == Set(listed) || error("results/ holds certificates $stored, listed $listed")
for p in [POINTS; RITZ_POINTS]
    is_frame_change(physical_ladder(p.g), ladder(p.g), rung_hadamard()) ||
        error("g = $(p.g): the Hadamard frame change fails")
end
for p in POINTS
    model = ladder(p.g)
    leanpoint = only(filter(q -> q.params == [p.g, p.β0], lean))
    path = point_path(ROOT, p.g)
    fresh = check_point(path, model; n0 = p.n0, β0 = p.β0, params = (; g = p.g),
        model_dir = ROOT, sources = SOURCES, rp = ladder_rp(model),
        lean = (s = leanpoint.s, t = leanpoint.t), spectra = full)
    println(basename(path), ": OK, gap > ", fresh["gap_lower_bound"],
        ", interval bound ", fresh["interval_bound_decimal"],
        full ? ", fresh defect upper ends $(fresh["defect_s_upper"]), $(fresh["defect_t_upper"])" :
        "")
end
for p in RITZ_POINTS
    model = ladder(p.g)
    path = ritz_point_path(ROOT, p.g)
    cert = read_result(path)
    record = cert["partial_spectrum"]
    partial = partial_spectrum_path(ROOT, p.g)
    record["path"] == basename(partial) && record["sha256"] == sha256_file(partial) ||
        error("$(basename(partial)) is not the partial spectrum pinned by $(basename(path))")
    large = full ? doubled_ring_bounds(p.g; progress = true) : read_partial_spectrum(partial)
    spectra = Dict(p.n0 => verified_ring_spectrum(model, p.n0), 2p.n0 => large)
    fresh = check_point(path, model; n0 = p.n0, β0 = p.β0, params = (; g = p.g),
        model_dir = ROOT, sources = SOURCES, rp = ladder_rp(model), spectra)
    # Finite bridge: the stored and the freshly verified gaps of the rings below n₀ are at least
    # the certified bound, so the certificate covers every even L ≥ 4.
    bound = max(parse(Rational{BigInt}, fresh["gap_lower_bound"]),
        Rational{BigInt}(fresh["interval_bound_decimal"]))
    stored_bridge = cert["finite_bridge"]
    [parse(Int, b["L"]) for b in stored_bridge] == collect(BRIDGE_LENGTHS) ||
        error("$(basename(path)): the bridged lengths are not $(BRIDGE_LENGTHS)")
    for (b, f) in zip(stored_bridge, finite_bridge(p.g))
        parse(Rational{BigInt}, b["gap_lower"]) >= bound &&
            parse(Rational{BigInt}, f["gap_lower"]) >= bound ||
            error("$(basename(path)): the L = $(f["L"]) gap is not above the certified bound")
    end
    println(basename(path), ": OK, gap > ", fresh["gap_lower_bound"],
        ", interval bound ", fresh["interval_bound_decimal"],
        ", defect upper ends from ", full ? "fresh" : "stored", " Ritz bounds ",
        "$(fresh["defect_s_upper"]), $(fresh["defect_t_upper"])")
end
