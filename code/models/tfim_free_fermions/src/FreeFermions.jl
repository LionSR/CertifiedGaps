# Closed-form partition functions of the periodic transverse-field Ising ring
#
#     H_L = -Σ_j Z_j Z_{j+1} - g Σ_j X_j     (L sites, Z_{L+1} = Z_1)
#
# in BigFloat interval arithmetic, for g > 1. The Jordan–Wigner transformation (with X_j as
# the occupation parity) maps H_L to free fermions whose boundary condition depends on the
# fermion parity P = Π X_j: antiperiodic (Neveu–Schwarz, k = 2π(m + 1/2)/L) on P = +1 and
# periodic (Ramond, k = 2πm/L) on P = -1. After the Bogoliubov rotation every mode has energy
# ε_k = 2√(1 + g² - 2g cos k) ≥ 2(g - 1) > 0, and projecting each sector onto its parity gives
#
#     Z_L(β) = ½ [Π_NS 2cosh(βε_k/2) + Π_NS 2sinh(βε_k/2) + Π_R 2cosh(βε_k/2) - Π_R 2sinh(βε_k/2)].
#
# For g < 1 the Ramond zero mode carries the signed energy 2(g - 1) < 0; this file only
# handles g > 1, where every ε_k is positive (checked below). The formula is the classical
# result of Lieb–Schultz–Mattis (1961) and Pfeuty (1970); METHOD.md derives it and
# test/runtests.jl compares it with exact diagonalization of rings of 4 to 10 sites.

using IntervalArithmetic
using CertifiedGaps.ExactArithmetic: RationalBig

"""
    exact(x::BigFloat) -> Rational{BigInt}

The exact value of `x`. `Rational{BigInt}(x)` alone rounds to the *global* BigFloat precision,
not to that of `x`, so it is converted at `precision(x)` bits.
"""
exact(x::BigFloat) = setprecision(() -> Rational{BigInt}(x), BigFloat, precision(x))

"BigFloat interval enclosing the exact rational `x`."
ball(x::Rational) = interval(BigFloat, numerator(x)) / interval(BigFloat, denominator(x))
ball(x::Integer) = interval(BigFloat, x)

"""
    mode_energies(L, g, sector) -> Vector{Interval{BigFloat}}

Enclosures of the L single-fermion energies ε_k = 2√(1 + g² − 2g cos k) of the ring of L sites,
with k = 2π(m + 1/2)/L for `sector = :ns` and k = 2πm/L for `sector = :r` (m = 0, …, L − 1).
Requires a rational g > 1, so that every ε_k ≥ 2(g − 1) > 0.
"""
function mode_energies(L::Integer, g::Rational, sector::Symbol)
    g > 1 || throw(ArgumentError("the closed form here needs g > 1, got g = $g"))
    shift = sector === :ns ? 1 // 2 : sector === :r ? 0 // 1 : throw(ArgumentError("sector"))
    G = ball(g)
    twopi = 2 * interval(BigFloat, π)
    energies = map(0:(L-1)) do m
        k = twopi * ball(big(m) + shift) / ball(L)
        2 * sqrt(1 + G^2 - 2G * cos(k))
    end
    all(e -> inf(e) > 0, energies) || error("an enclosure of ε_k is not positive")
    return energies
end

# log(2 cosh(x/2)) and log(2 sinh(x/2)) for x > 0, written so that large x does not overflow
# and small x loses no relative accuracy: x/2 + log1p(±e^{−x}).
_logcosh2(x) = x / 2 + log1p(exp(-x))
_logsinh2(x) = x / 2 + log(-expm1(-x))

"""
    logZ(L, g, β; precision = 512) -> Interval{BigFloat}

Rigorous enclosure of log Z_L(β) = log Tr e^{−βH_L} for the periodic transverse-field Ising ring
of L sites at rational g > 1 and rational β > 0, from the free-fermion closed form (see the file
header). Every factor is positive, and the four products are combined relative to the largest,
the Neveu–Schwarz cosh product, so no exponent overflows:

    log Z = C + log((1 + e^{S − C} + e^{C_R − C}(1 − e^{S_R − C_R})) / 2)

with C, S, C_R, S_R the logarithms of the four products.
"""
function logZ(L::Integer, g::Rational, β::Rational; precision::Integer = 512)
    L >= 2 && β > 0 || throw(ArgumentError("need L ≥ 2 and β > 0"))
    setprecision(BigFloat, precision) do
        B = ball(β)
        ns = [B * e for e in mode_energies(L, g, :ns)]
        r = [B * e for e in mode_energies(L, g, :r)]
        C, S = sum(_logcosh2, ns), sum(_logsinh2, ns)
        CR, SR = sum(_logcosh2, r), sum(_logsinh2, r)
        # 1 − e^{S_R − C_R} = 1 − Π tanh(βε_k/2) > 0, enclosed by −expm1 for accuracy.
        bracket = 1 + exp(S - C) + exp(CR - C) * (-expm1(SR - CR))
        inf(bracket) > 0 || error("the enclosure of Z_L is not positive")
        C + log(bracket / 2)
    end
end
