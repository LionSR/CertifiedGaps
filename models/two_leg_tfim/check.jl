# Checker for the two-leg TFIM certificates (results/certificate_g*.json, core format).
#
#     julia --project=julia models/two_leg_tfim/check.jl          # fast
#     julia --project=julia models/two_leg_tfim/check.jl --full   # + spectra
#
# For every point of `Ladder.POINTS`, `GapCertificates.check_point` requires the model-local
# source hashes, exactly the point's rotated-frame model `Ladder.ladder(g)` (charge, leg-exchange
# symmetry, reflection), its RP factors, seed and params, the Lean seed hypotheses of
# lean/CertifiedGaps/Models/TwoLegIsing.lean, and then re-verifies the exact RP identity, the
# exact rational chain and the interval bound. It also checks the exact Hadamard frame change
# from the physical (Lean) ladder. `--full` also recomputes every verified spectrum and the seed
# defects, which must not exceed the stored εₛ, εₜ.
include(joinpath(@__DIR__, "src", "Ladder.jl"))
using .Ladder, CertifiedGaps.GapCertificates, CertifiedGaps.RingModels
using LinearAlgebra
BLAS.set_num_threads(1)
const ROOT = @__DIR__
const LEAN = joinpath(ROOT, "..", "..", "lean", "CertifiedGaps", "Models", "TwoLegIsing.lean")

full = "--full" in ARGS
lean = lean_points(LEAN)
# Every stored certificate must belong to a listed point, so none goes unchecked.
stored = filter(f -> startswith(f, "certificate_"), readdir(joinpath(ROOT, "results")))
Set(stored) == Set(basename(point_path(ROOT, p.g)) for p in POINTS) ||
    error("results/ holds certificates of unlisted points: $stored")
for p in POINTS
    model = ladder(p.g)
    is_frame_change(physical_ladder(p.g), model, rung_hadamard()) ||
        error("g = $(p.g): the Hadamard frame change fails")
    leanpoint = only(filter(q -> q.params == [p.g, p.β0], lean))
    path = point_path(ROOT, p.g)
    fresh = check_point(path, model; n0 = p.n0, β0 = p.β0, params = (; g = p.g),
        model_dir = ROOT, rp = ladder_rp(model), lean = (s = leanpoint.s, t = leanpoint.t),
        spectra = full)
    println(basename(path), ": OK, gap > ", fresh["gap_lower_bound"],
        ", interval bound ", fresh["interval_bound_decimal"],
        full ? ", fresh defect upper ends $(fresh["defect_s_upper"]), $(fresh["defect_t_upper"])" :
        "")
end
