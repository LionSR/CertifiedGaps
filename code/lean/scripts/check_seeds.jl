# Cross-check the seed defects quoted in the Lean model files against the model certificates.
#
# Every model file `lean/CertifiedGaps/Models/<Model>.lean` lists its certified points as
# `⟨param₁, …, β₀, { s := [εₛ, …], t := [εₜ, …], … }⟩`. The Lean seed statement for a point is
# `1 - εₛ ≤ P_s(n₀, β₀)` and `1 - εₜ ≤ P_τ(2n₀, β₀)`, so it is implied by a core certificate
# (`GapCertificates.certify_gap`) in `models/<slug>/results` exactly when the Lean decimals are no
# smaller than the certified defects at the same parameters (the certificate's `params`) and the
# same seed `(n₀, β₀)`. Both sides are parsed as exact rationals, the Lean side by the core
# `GapCertificates.lean_points`. Each model's `check.jl` makes the same comparison for its own
# points (`check_point(...; lean)`); this script checks every Lean file at once. It does not edit
# anything: a failure is reported, not fixed.
#
# Usage (from the repository root):
#     julia --project=julia lean/scripts/check_seeds.jl
# Exit status 1 on any mismatch or missing record.

using CertifiedGaps.Results: read_result
using CertifiedGaps.ExactArithmetic: parse_rational
using CertifiedGaps.GapCertificates: lean_points

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REPO = normpath(joinpath(ROOT, ".."))

# Lean file, model directory, the certificate `params` keys in the order of the Lean point
# fields before β₀, and the seed length n₀ of the Lean `seeds` statement.
const MODELS = [
    ("SpinOneLargeD.lean", "spin1_largeD", ["D"], 6),
    ("SpinOneDimerized.lean", "spin1_dimerized", ["delta"], 2),
    ("HeisenbergLadder.lean", "heisenberg_ladder", ["J"], 4),
    ("TwoLegIsing.lean", "two_leg_tfim", ["g"], 4),
    ("StaggeredHeisenberg.lean", "heisenberg_staggered_field", ["h"], 8),
    ("DimerizedHeisenberg.lean", "dimerized_heisenberg", ["delta"], 4),
    ("IsingLongitudinal.lean", "tfim_longitudinal", ["g", "h"], 6),
]

"Core certificates of a model as (params, n₀, β₀, εₛ, εₜ, source), all exact."
function records(slug, keys)
    dir = joinpath(REPO, "models", slug, "results")
    out = NamedTuple[]
    for f in sort(readdir(dir))
        occursin(r"^certificate.*\.json(\.gz)?$", f) || continue
        x = read_result(joinpath(dir, f), Dict{String,Any})
        startswith(string(get(x, "format", "")), "CertifiedGaps ring gap certificate") || continue
        haskey(x, "params") || error("$slug/$f: core certificate without params")
        s = x["seed"]
        push!(out, (params = [parse_rational(string(x["params"][k])) for k in keys],
            n₀ = parse(Int, string(s["n0"])), β₀ = parse_rational(string(s["beta0"])),
            εs = parse_rational(s["eps_s"]), εt = parse_rational(s["eps_t"]),
            source = "models/$slug/results/$f"))
    end
    return out
end

"Check every Lean point against the certified records; return the number of failures."
function check()
    failures = 0
    for (file, slug, keys, n₀) in MODELS
        recs = records(slug, keys)
        for p in lean_points(joinpath(ROOT, "CertifiedGaps", "Models", file))
            params, β₀ = p.params[1:end-1], p.params[end]
            label = "$file $(join(string.(params), ", ")) at (n₀, β₀) = ($n₀, $β₀)"
            same = [r for r in recs if r.params == params && r.n₀ == n₀ && r.β₀ == β₀]
            if isempty(same)
                println("MISSING  $label: no certified record with these parameters and seed")
                failures += 1
                continue
            end
            ok = [r for r in same if r.εs <= p.s && r.εt <= p.t]
            if isempty(ok)
                r = same[1]
                println("FAIL     $label: Lean εₛ = $(Float64(p.s)), εₜ = $(Float64(p.t)) ",
                        "< certified $(Float64(r.εs)), $(Float64(r.εt)) ($(r.source))")
                failures += 1
            else
                r = ok[1]
                println("ok       $label: εₛ $(Float64(r.εs)) ≤ $(Float64(p.s)), ",
                        "εₜ $(Float64(r.εt)) ≤ $(Float64(p.t)) ($(r.source))")
            end
        end
    end
    return failures
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    n = check()
    println(n == 0 ? "all Lean seeds are implied by the certificates" : "$n failure(s)")
    exit(n == 0 ? 0 : 1)
end
