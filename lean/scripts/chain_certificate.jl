# Rational chain data for `CertifiedGaps.ChainCertificate` (lean/CertifiedGaps/ChainCheck.lean).
#
# Given seed defect bounds εₛ, εₜ and the seed inverse temperature β₀ (exact rationals), this
# computes a chain s_j, t_j of short decimal upper bounds for the exact doubling recurrence
#     z(a, b) = 1 - (1 - b)(1 - a)²,   f(z) = z² / (2(1 - z)²),
#     s_{j+1} ≥ f(z(s_j, t_j)),        t_{j+1} ≥ f(z(t_j, s_{j+1})),
# a hand-over defect ε ≥ s_J, t_J with ε ≤ 1/100, and the largest short decimal r with
#     r^{2^j} z_j < 1 - z_j   (j < J),     r^{2^J} · 5ε ≤ 1.
# Lean re-checks every inequality by `decide`; this script only searches. The gap bound is
# log r / β₀. Everything here is exact `Rational{BigInt}` arithmetic.
#
# Usage: julia lean/scripts/chain_certificate.jl εₛ εₜ β₀ [J] [digits]
#   e.g. julia lean/scripts/chain_certificate.jl 577405872923865//9007199254740992 \
#          45632040552745//1125899906842624 9//4 8
# (this reproduces the `D = 4` point of `Models/SpinOneLargeD.lean`). The arguments are parsed
# as literals, never evaluated: an integer, `a//b`, `a/b`, or a decimal such as `0.064` or
# `6411e-5` (converted exactly in base ten).

"""
    parse_rational(x) -> Rational{BigInt}

Parse an exact rational from a string without evaluating it as code: an integer (`-12`), a
fraction `a//b` or `a/b`, or a decimal with optional exponent (`0.064`, `6411e-5`, `1.5E+3`).
Decimals are converted exactly, never through floating point.
"""
function parse_rational(x::AbstractString)
    str = strip(x)
    m = match(r"^([+-]?\d+)\s*//?\s*(\d+)$", str)
    if m !== nothing
        den = parse(BigInt, m[2])
        iszero(den) && throw(ArgumentError("zero denominator in $(repr(x))"))
        return parse(BigInt, m[1]) // den
    end
    m = match(r"^([+-]?)(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$", str)
    frac = m === nothing || m[3] === nothing ? "" : m[3]
    (m === nothing || isempty(m[2] * frac)) &&
        throw(ArgumentError("not an exact rational literal: $(repr(x))"))
    mant = parse(BigInt, m[2] * frac)
    k = (m[4] === nothing ? 0 : parse(Int, m[4])) - length(frac)
    q = k >= 0 ? (mant * big(10)^k) // big(1) : mant // big(10)^(-k)
    return m[1] == "-" ? -q : q
end

z(a, b) = 1 - (1 - b) * (1 - a)^2
f(x) = x^2 / (2 * (1 - x)^2)

"Smallest decimal `m / 10^k` ≥ `x` with `digits` significant digits."
function round_up(x::Rational{BigInt}, digits)
    x <= 0 && return zero(x)
    k = digits - 1 - floor(Int, log10(BigFloat(x)))
    scale = big(10)^max(k, 0) // big(10)^max(-k, 0)
    return ceil(BigInt, x * scale) // scale
end

"Largest decimal `m / 10^k` ≤ `x` with `digits` significant digits."
function round_down(x::Rational{BigInt}, digits)
    k = digits - 1 - floor(Int, log10(BigFloat(x)))
    scale = big(10)^max(k, 0) // big(10)^max(-k, 0)
    return floor(BigInt, x * scale) // scale
end

function chain(εs, εt, J, digits)
    s, t = [round_up(εs, digits)], [round_up(εt, digits)]
    for _ in 1:J
        z(s[end], t[end]) < 1 || return nothing
        push!(s, round_up(f(z(s[end], t[end])), digits))
        z(t[end], s[end]) < 1 || return nothing
        push!(t, round_up(f(z(t[end], s[end])), digits))
    end
    return s, t
end

function valid(s, t, r, ε)
    J = length(s) - 1
    all(j -> r^(big(2)^(j - 1)) * z(s[j], t[j]) < 1 - z(s[j], t[j]), 1:J) || return false
    return ε <= 1 // 100 && s[end] <= ε && t[end] <= ε && r^(big(2)^J) * 5ε <= 1
end

"Largest `r` (bisection in floating point, then a `digits`-digit decimal) with `valid`."
function best_r(s, t, ε, digits)
    lo, hi = big"1.0", big"1e6"
    ok(x) = valid(s, t, Rational{BigInt}(x), ε)
    ok(lo) || return nothing
    for _ in 1:200
        mid = (lo + hi) / 2
        ok(mid) ? (lo = mid) : (hi = mid)
    end
    r = round_down(Rational{BigInt}(lo), digits)
    return valid(s, t, r, ε) ? r : nothing
end

"Decimal literal `me-k` for a rational with denominator dividing a power of ten (exact in Lean)."
function lean(q::Rational)
    k = 0
    while denominator(q * big(10)^k) != 1
        k += 1
    end
    m = numerator(q * big(10)^k)
    return k == 0 ? string(m) : string(m, "e-", k)
end

"`  name := [a, b, …]`, wrapped at 100 columns with a 6-space continuation indent."
function lean_list(name, xs)
    lines, cur = String[], "  $name := ["
    for (i, x) in enumerate(xs)
        tok = lean(x) * (i < length(xs) ? "," : "]")
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

function main(args)
    εs, εt, β₀ = parse_rational.(args[1:3])
    J = length(args) ≥ 4 ? parse(Int, args[4]) : 6
    digits = length(args) ≥ 5 ? parse(Int, args[5]) : 4
    c = chain(εs, εt, J, digits)
    c === nothing && error("the recurrence leaves z < 1 within $J steps")
    s, t = c
    ε = max(s[end], t[end])
    ε <= 1 // 100 || error("defects at scale J = $J are $(Float64(ε)) > 1/100; increase J")
    r = best_r(s, t, ε, digits)
    r === nothing && error("no r ≥ 1 satisfies the rate conditions")
    println("-- gap bound log r / β₀ = ", Float64(log(BigFloat(r)) / BigFloat(β₀)))
    println(lean_list("s", s))
    println(lean_list("t", t))
    println("  r := ", lean(r))
    println("  ε := ", lean(ε))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS)
end
