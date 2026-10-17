# Generate the Lean chain certificates of the seeded leaderboard records:
#
#     julia --project=tools/leaderboard tools/leaderboard/lean_chains.jl
#
# For every row of home_page/_data/records.yml whose Lean theorem `Leaderboard.Records.<M>.<x>`
# is listed in CHAINED below, this reads the exact seed defects (εₛ, εₜ) from the row's
# certificate (or `lean.certificate`) with the validator's own reader, computes a chain with
# lean/scripts/chain_certificate.jl, and writes it as `Leaderboard.Chains.<M>.<x>` to
# leaderboard/Leaderboard/Chains/<M>.lean, where `decide +kernel` re-checks every inequality.
# It also writes the theorem's `r` and the decimal of log r / β₀ (rounded down, at the
# precision of the row's headline decimal) into the row's `lean:` line. This script only
# searches; Lean and the validator check everything it writes.
#
# The chain has J = 14 doublings and 8 significant digits. The Lean tail condition costs about
# 2⁻ᴶ log 5 / β₀ against the computer bound, which follows the recurrence to J = 48 in
# logarithms; at J = 14 the two agree to about 10⁻⁵ relative.

using LeaderboardCheck
using LeaderboardCheck: certificate_candidates, load_certificate, match_certificate, cert_exact,
    log_enclosure, params_agree
using YAML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "lean", "scripts", "chain_certificate.jl"))

const J = 14
const DIGITS = 8

# Records whose seed is proved in Lean rather than read from a certificate: the defects of the
# kernel-checked seed theorems.
const KERNEL_SEEDS = Dict(
    "Leaderboard.Records.IsingLongitudinal.g_1_4_h_1_20" => (big(9064) // 10^9, big(3141) // 10^8),
)

const HEADER = """
/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import CertifiedGaps.ChainCheck
"""

"""
    lean_seed(rec, name)

The exact seed `(n₀, β₀, εₛ, εₜ)` of a row's Lean theorem: the row's seed in its certificate, or,
when the row names `lean.certificate`, the seed of that certificate at `β₀ = lean.beta` (which
must have a single seed length).
"""
function lean_seed(rec, name)
    seed = rec["seed"]
    n0, β0 = parse_exact(string(seed["n0"])), parse_exact(string(seed["beta0"]))
    haskey(KERNEL_SEEDS, name) && return (n0, β0, KERNEL_SEEDS[name]...)
    params = something(get(rec, "params", nothing), Dict())
    path = get(rec["lean"], "certificate", rec["certificate"])
    if haskey(rec["lean"], "certificate")
        β0 = parse_exact(rec["lean"]["beta"])
        data = load_certificate(joinpath(ROOT, path))
        lengths = unique(cert_exact(c.n0) for c in certificate_candidates(data)
                         if cert_exact(c.beta) == β0 && params_agree(c.scope, params))
        length(lengths) == 1 || error("$(rec["id"]): $path has seed lengths $lengths at β₀ = $β0")
        n0 = lengths[1]
    end
    c = match_certificate(joinpath(ROOT, path), n0, β0, params)
    (c === nothing || c.s === nothing || c.t === nothing) &&
        error("$(rec["id"]): $path states no exact defects for the seed ($n0, $β0)")
    return (n0, β0, cert_exact(c.s), cert_exact(c.t))
end

"The decimal of log r / β rounded down at `k` places, or `nothing` if the enclosure straddles."
function floor_decimal(r, β, k)
    lo, hi = log_enclosure(r)
    lo, hi = lo / β, hi / β
    q = floor(BigInt, lo * big(10)^k) // big(10)^k
    q + big(1) // big(10)^k > hi || return nothing
    s = string(numerator(q * big(10)^k))
    s = lpad(s, k + 1, '0')
    return s[1:end-k] * "." * s[end-k+1:end]
end

"""
    fast_r(s, t, ε)

The largest `DIGITS`-digit decimal `r` with `valid(s, t, r, ε)`. The rate conditions bound `r`
by `((1 - zⱼ)/zⱼ)^(2⁻ʲ)` and `(5ε)^(-2⁻ᴶ)`; their minimum is found in `BigFloat`, and only the
rounded candidate is checked exactly (the bisection of `best_r` checks every step exactly, which
is slow at J = 14).
"""
function fast_r(s, t, ε)
    setprecision(BigFloat, 512) do
        logs = [log((1 - BigFloat(z(s[j], t[j]))) / BigFloat(z(s[j], t[j]))) / big(2)^(j - 1)
                for j in 1:length(s)-1]
        push!(logs, -log(5 * BigFloat(ε)) / big(2)^(length(s) - 1))
        r = round_down(Rational{BigInt}(exp(minimum(logs)) * (1 - big(10)^-40)), DIGITS)
        e = floor(Int, log10(BigFloat(r))) - DIGITS + 1
        step = e >= 0 ? big(10)^e // big(1) : big(1) // big(10)^(-e)
        while r > 1 && !valid(s, t, r, ε)
            r -= step
        end
        return r > 1 ? r : nothing
    end
end

function lean_chain(εs, εt)
    c = chain(εs, εt, J, DIGITS)
    c === nothing && error("the recurrence leaves z < 1 within $J steps")
    s, t = c
    ε = max(s[end], t[end])
    ε <= 1 // 100 || error("defects at scale $J are $(Float64(ε)) > 1/100")
    r = fast_r(s, t, ε)
    r === nothing && error("no r > 1 satisfies the rate conditions")
    return s, t, r, ε
end

"""
    lean_literal(q)

The Lean literal `me-k` of a chain entry `q = m / 10^k` with at most `DIGITS` significant digits.
`lean` of chain_certificate.jl finds `k` by trying `k = 0, 1, …`, which is quadratic for the
entries near 10⁻⁹⁰⁰⁰⁰ of a deep chain; here `k` follows from the size of `q`.
"""
function lean_literal(q::Rational{BigInt})
    k = max(DIGITS - 1 - floor(Int, log10(BigFloat(q))), 0)
    m = q * big(10)^k
    denominator(m) == 1 || error("$q is not a $DIGITS-digit decimal")
    m = numerator(m)
    while k > 0 && m % 10 == 0
        m, k = div(m, 10), k - 1
    end
    return k == 0 ? string(m) : string(m, "e-", k)
end

"`  name := [a, b, …]` of chain entries, wrapped at 100 columns (as `lean_list`)."
function literal_list(name, xs)
    lines, cur = String[], "  $name := ["
    for (i, x) in enumerate(xs)
        tok = lean_literal(x) * (i < length(xs) ? "," : "]")
        cand = endswith(cur, "[") ? cur * tok : cur * " " * tok
        if length(cand) <= 100
            cur = cand
        else
            push!(lines, cur)
            cur = "      " * tok
        end
    end
    push!(lines, cur)
    return join(lines, "\n")
end

function rational_string(q)
    return denominator(q) == 1 ? string(numerator(q)) : string(numerator(q), "/", denominator(q))
end

function main()
    yaml_path = joinpath(ROOT, "home_page", "_data", "records.yml")
    data = YAML.load_file(yaml_path)
    text = read(yaml_path, String)
    modules = Dict{String,Vector{String}}()
    for rec in data["records"]
        leanrow = get(rec, "lean", nothing)
        leanrow === nothing && continue
        name = leanrow["name"]
        parts = split(name, '.')
        mod, thm = parts[3], parts[4]
        n0, β0, εs, εt = lean_seed(rec, name)
        s, t, r, ε = lean_chain(εs, εt)
        k = max(length(split(rec["bound"]["decimal"], '.')[2]), 6)
        dec = something(floor_decimal(r, β0, k), floor_decimal(r, β0, k - 1))
        src = haskey(KERNEL_SEEDS, name) ? "the kernel-checked seed of `SeedKernel.Ising.Seed`" :
              "`" * get(leanrow, "certificate", rec["certificate"]) * "`"
        block = """
            /-- `$(rec["id"])`: seed `(n₀, β₀) = ($(numerator(n0)), $(rational_string(β0)))` from
            $src; `log r / β₀ ≥ $dec`. -/
            def $thm : ChainCertificate where
            $(literal_list("s", s))
            $(literal_list("t", t))
              r := $(lean_literal(r))
              ε := $(lean_literal(ε))

            theorem $(thm)_valid : $thm.Valid := by decide +kernel
            """
        push!(get!(modules, mod, String[]), block)
        # The row's `lean:` line: replace r and decimal, keep everything else.
        line = match(Regex("lean: \\{name: " * replace(name, "." => "\\.") * ",[^\\n]*"), text)
        line === nothing && error("$(rec["id"]): no single-line `lean:` entry for $name")
        new = replace(line.match, r"r: \"[^\"]*\"" => "r: \"$(rational_string(r))\"",
            r"decimal: \"[^\"]*\"" => "decimal: \"$dec\"")
        text = replace(text, line.match => new)
        println(rpad(rec["id"], 32), " r = ", rational_string(r), "  log r / β₀ ≥ ", dec,
            "  (row: ", rec["bound"]["decimal"], ")")
    end
    write(yaml_path, text)
    dir = joinpath(ROOT, "leaderboard", "Leaderboard", "Chains")
    mkpath(dir)
    for (mod, blocks) in sort(collect(modules))
        open(joinpath(dir, mod * ".lean"), "w") do io
            print(io, HEADER, "\n")
            print(io, """
                /-!
                # Chain certificates: `Leaderboard.Records.$mod`

                Generated by `tools/leaderboard/lean_chains.jl` from the seed defects of the
                certificates named below; do not edit. Each `x` is the chain of the record
                `Leaderboard.Records.$mod.x`, and `x_valid` checks it in the kernel.
                -/

                namespace Leaderboard.Chains.$mod

                open CertifiedGaps

                """)
            print(io, join(blocks, "\n"))
            print(io, "\nend Leaderboard.Chains.$mod\n")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
