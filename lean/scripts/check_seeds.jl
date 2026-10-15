# Cross-check the seed defects quoted in the Lean model files against the model certificates.
#
# Every model file `lean/CertifiedGaps/Models/<Model>.lean` lists its certified points as
# `⟨param₁, …, β₀, { s := [εₛ, …], t := [εₜ, …], … }⟩`. The Lean seed statement for a point is
# `1 - εₛ ≤ P_s(n₀, β₀)` and `1 - εₜ ≤ P_τ(2n₀, β₀)`, so it is implied by the certificate in
# `models/<slug>/results` exactly when the Lean decimals are no smaller than the certified
# defects at the same parameters and the same seed `(n₀, β₀)`. This script parses both sides as
# exact rationals and checks that. It does not edit anything: a failure is reported, not fixed.
#
# Usage (from the repository root; JSON3 comes from the core project):
#     julia --project=julia lean/scripts/check_seeds.jl
# Exit status 1 on any mismatch or missing record.

using JSON3

include(joinpath(@__DIR__, "chain_certificate.jl"))   # `parse_rational` (literal parser)

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const LEAN_MODELS = joinpath(ROOT, "lean", "CertifiedGaps", "Models")
const RESULTS(slug) = joinpath(ROOT, "models", slug, "results")

"Exact rational from a JSON value: an integer, or a string `a//b` / decimal."
exact(x::Integer) = Rational{BigInt}(x)
exact(x::AbstractString) = parse_rational(x)
exact(x) = error("not an exact JSON number: $(repr(x))")

"Read a JSON file, decompressing `.gz` with the system `gzip`."
readjson(path) = JSON3.read(endswith(path, ".gz") ? read(`gzip -dc $path`, String) :
                            read(path, String))

"Files in `dir` whose names match `rx`."
matching(dir, rx) = [joinpath(dir, f) for f in sort(readdir(dir)) if occursin(rx, f)]

"""
    lean_points(file) -> Vector{NamedTuple}

The points of a Lean model file: `params` (all anonymous-constructor fields before `β₀`), `β₀`
and the first entries `εₛ = s[0]`, `εₜ = t[0]` of the certificate.
"""
function lean_points(file)
    text = read(joinpath(LEAN_MODELS, file), String)
    body = split(text, r"def points : List Point := \["; limit = 2)
    length(body) == 2 || error("$file: no `def points : List Point := [`")
    pts = NamedTuple[]
    rx = r"⟨([^⟨{]*?),\s*\{\s*s := \[\s*([0-9e.+-]+)\s*,.*?t := \[\s*([0-9e.+-]+)\s*,"s
    for m in eachmatch(rx, body[2])
        fields = parse_rational.(split(replace(m[1], " " => ""), ","))
        push!(pts, (params = fields[1:end-1], β₀ = fields[end],
                    εs = parse_rational(m[2]), εt = parse_rational(m[3])))
    end
    isempty(pts) && error("$file: no points parsed")
    return pts
end

# Certified records per model: (params, n₀, β₀, εₛ, εₜ, source), all exact.
rec(params, n0, β, εs, εt, src) =
    (params = exact.(params), n₀ = Int(exact(n0)), β₀ = exact(β), εs = exact(εs),
     εt = exact(εt), source = relpath(src, ROOT))

function records_spin1_largeD()
    map(matching(RESULTS("spin1_largeD"), r"^certificate_D.*\.json$")) do f
        x = readjson(f)
        rec([x.D], x.n0, x.beta, x.eps_s, x.eps_t, f)
    end
end

function records_spin1_dimerized()
    map(matching(RESULTS("spin1_dimerized"), r"^certificate_d.*\.json$")) do f
        x = readjson(f)
        rec([x.delta], x.n0, x.beta0, x.epsilon_s, x.epsilon_t, f)
    end
end

function records_heisenberg_ladder()
    f = joinpath(RESULTS("heisenberg_ladder"), "certificate.json")
    [rec([p.J], p.best.n0, p.best.beta0, p.best.eps_s_upper, p.best.eps_t_upper, f)
     for p in readjson(f).parameters if p.best !== nothing]
end

function records_two_leg_tfim()
    map(matching(RESULTS("two_leg_tfim"), r"^certificate_g.*\.json(\.gz)?$")) do f
        x = readjson(f)
        rec([x.g], x.n0, x.beta, x.epsilon_s, x.epsilon_t, f)
    end
end

function records_heisenberg_staggered_field()
    map(matching(RESULTS("heisenberg_staggered_field"), r"^certificate_\d+_\d+\.json$")) do f
        x = readjson(f)
        rec([x.h], x.n0, x.beta, x.epsilon_s, x.epsilon_t, f)
    end
end

function records_dimerized_heisenberg()
    f = joinpath(RESULTS("dimerized_heisenberg"), "certificate.json")
    out = []
    for p in readjson(f).parameters
        seeds = collect(p.seeds)
        haskey(p, :survey_seed) && p.survey_seed !== nothing && push!(seeds, p.survey_seed)
        for s in seeds
            push!(out, rec([p.delta], s.n0_cells, s.beta0, s.eps_s, s.eps_t, f))
        end
    end
    return out
end

function records_tfim_longitudinal()
    f = joinpath(RESULTS("tfim_longitudinal"), "seeds.json")
    [rec([p.g, p.h], p.n0, p.beta0, p.epsilon_s, p.epsilon_t, f) for p in readjson(f).points]
end

# Lean file, seed length n₀ (as in its `seeds` statement), certified records.
const MODELS = [
    ("SpinOneLargeD.lean", 6, records_spin1_largeD),
    ("SpinOneDimerized.lean", 2, records_spin1_dimerized),
    ("HeisenbergLadder.lean", 4, records_heisenberg_ladder),
    ("TwoLegIsing.lean", 4, records_two_leg_tfim),
    ("StaggeredHeisenberg.lean", 8, records_heisenberg_staggered_field),
    ("DimerizedHeisenberg.lean", 4, records_dimerized_heisenberg),
    ("IsingLongitudinal.lean", 6, records_tfim_longitudinal),
]

"Check every Lean point against the certified records; return the number of failures."
function check()
    failures = 0
    for (file, n₀, records) in MODELS
        recs = records()
        for p in lean_points(file)
            label = "$file $(join(string.(p.params), ", ")) at (n₀, β₀) = ($n₀, $(p.β₀))"
            same = [r for r in recs if r.params == p.params && r.n₀ == n₀ && r.β₀ == p.β₀]
            if isempty(same)
                println("MISSING  $label: no certified record with these parameters and seed")
                failures += 1
                continue
            end
            ok = [r for r in same if r.εs <= p.εs && r.εt <= p.εt]
            if isempty(ok)
                r = same[1]
                println("FAIL     $label: Lean εₛ = $(Float64(p.εs)), εₜ = $(Float64(p.εt)) ",
                        "< certified $(Float64(r.εs)), $(Float64(r.εt)) ($(r.source))")
                failures += 1
            else
                r = ok[1]
                println("ok       $label: εₛ $(Float64(r.εs)) ≤ $(Float64(p.εs)), ",
                        "εₜ $(Float64(r.εt)) ≤ $(Float64(p.εt)) ($(r.source))")
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
