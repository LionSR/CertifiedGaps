# Closed-form partition functions of the periodic cluster chain in a transverse field
#
#     H = -Σ_j Z_{j-1} X_j Z_{j+1} - g Σ_j X_j     (N = 2L sites, L two-site cells)
#
# in BigFloat interval arithmetic, for rational g > 0 with g ≠ 1. Jordan–Wigner with X_j as the
# occupation parity turns Z_{j-1} X_j Z_{j+1} into a next-nearest-neighbor Majorana coupling;
# the boundary terms pick up the fermion parity, so even parity sees antiperiodic
# (Neveu–Schwarz, k = 2π(m + 1/2)/N) and odd parity periodic (Ramond, k = 2πm/N) fermions, with
#
#     ε_k = 2√(1 + g² − 2g cos 2k),
#     Z_N(β) = ½ [Π_NS 2cosh(βε_k/2) + Π_NS 2sinh(βε_k/2) + Π_R 2cosh(βε_k/2) − Π_R 2sinh(βε_k/2)].
#
# On even rings the Ramond modes k = 0, π share the signed energy 2(g − 1), so their signs
# cancel in the sinh product and the unsigned roots serve for every g ≠ 1. test/runtests.jl
# fixes these conventions against exact diagonalization of rings of 4 to 12 sites. The interval
# code follows code/models/tfim_free_fermions (see its METHOD.md for the BigFloat precision traps
# handled by `exact`).

using IntervalArithmetic

"BigFloat interval enclosing the exact rational `x`."
ball(x::Rational) = interval(BigFloat, numerator(x)) / interval(BigFloat, denominator(x))
ball(x::Integer) = interval(BigFloat, x)

"""
    exact(x::BigFloat) -> Rational{BigInt}

The exact value of `x`. `Rational{BigInt}(x)` alone rounds to the *global* BigFloat precision,
not to that of `x`, so it is converted at `precision(x)` bits.
"""
exact(x::BigFloat) = setprecision(() -> Rational{BigInt}(x), BigFloat, precision(x))

"""
    mode_energies(N, g, sector) -> Vector{Interval{BigFloat}}

Enclosures of the N single-fermion energies ε_k = 2√(1 + g² − 2g cos 2k) of the ring of N sites
(N even), with k = 2π(m + 1/2)/N for `sector = :ns` and k = 2πm/N for `sector = :r`. Requires
rational g > 0, g ≠ 1, so that every ε_k ≥ 2|1 − g| > 0.
"""
function mode_energies(N::Integer, g::Rational, sector::Symbol)
    g > 0 && g != 1 || throw(ArgumentError("need g > 0 and g ≠ 1, got g = $g"))
    iseven(N) || throw(ArgumentError("need an even number of sites"))
    shift = sector === :ns ? 1 // 2 : sector === :r ? 0 // 1 : throw(ArgumentError("sector"))
    G = ball(g)
    twopi = 2 * interval(BigFloat, π)
    energies = map(0:(N-1)) do m
        k = twopi * ball(big(m) + shift) / ball(N)
        2 * sqrt(1 + G^2 - 2G * cos(2k))
    end
    all(e -> inf(e) > 0, energies) || error("an enclosure of ε_k is not positive")
    return energies
end

# log(2 cosh(x/2)) and log(2 sinh(x/2)) for x > 0 without overflow: x/2 + log1p(±e^{−x}).
_logcosh2(x) = x / 2 + log1p(exp(-x))
_logsinh2(x) = x / 2 + log(-expm1(-x))

"""
    logZ(L, g, β; precision = 512) -> Interval{BigFloat}

Rigorous enclosure of log Z_L(β) for the periodic cluster chain of L two-site cells (N = 2L
sites), rational g > 0 with g ≠ 1 and rational β > 0, from the closed form in the file header.
The four products are combined relative to the Neveu–Schwarz cosh product:

    log Z = C + log((1 + e^{S − C} + e^{C_R − C}(1 − e^{S_R − C_R})) / 2).
"""
function logZ(L::Integer, g::Rational, β::Rational; precision::Integer = 512)
    L >= 2 && β > 0 || throw(ArgumentError("need L ≥ 2 cells and β > 0"))
    setprecision(BigFloat, precision) do
        b = ball(β)
        ns = [b * e for e in mode_energies(2L, g, :ns)]
        r = [b * e for e in mode_energies(2L, g, :r)]
        C, S = sum(_logcosh2, ns), sum(_logsinh2, ns)
        CR, SR = sum(_logcosh2, r), sum(_logsinh2, r)
        bracket = 1 + exp(S - C) + exp(CR - C) * (-expm1(SR - CR))
        inf(bracket) > 0 || error("the enclosure of Z_L is not positive")
        C + log(bracket / 2)
    end
end
