# Closed-form partition functions of the periodic XY ring in a field
#
#     H_L = Σ_j (J_x X_j X_{j+1} + J_y Y_j Y_{j+1}) - h Σ_j Z_j     (L sites)
#
# in BigFloat interval arithmetic, in the paramagnetic phase h > |J_x + J_y| with J_x J_y < 0.
# Jordan–Wigner (Z_j as the occupation parity) gives free fermions with hopping ∝ J_x + J_y and
# pairing ∝ J_x − J_y; the boundary condition is antiperiodic (Neveu–Schwarz,
# k = 2π(m + 1/2)/L) on even parity and periodic (Ramond, k = 2πm/L) on odd parity. Every mode
# has energy
#
#     ε_k = 2√((h − (J_x + J_y) cos k)² + (J_x − J_y)² sin² k),
#
# where the unpaired Ramond modes k = 0, π carry the signed energies 2(h ∓ (J_x + J_y)), which
# are positive here, so the square root gives them too. Projecting each sector onto its parity,
#
#     Z_L(β) = ½ [Π_NS 2cosh(βε_k/2) + Π_NS 2sinh(βε_k/2) + Π_R 2cosh(βε_k/2) − Π_R 2sinh(βε_k/2)].
#
# test/runtests.jl compares this with exact diagonalization of rings of 4 to 10 sites and the
# interval evaluation with Arb. The interval code follows code/models/tfim_free_fermions
# (see its METHOD.md for the BigFloat precision traps handled by `exact`).

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

"Whether (J_x, J_y, h) lies in the certified regime J_x J_y < 0, h > |J_x + J_y|."
in_regime(Jx, Jy, h) = Jx * Jy < 0 && h > abs(Jx + Jy)

"""
    mode_energies(L, Jx, Jy, h, sector) -> Vector{Interval{BigFloat}}

Enclosures of the L single-fermion energies ε_k of the ring of L sites, with k = 2π(m + 1/2)/L
for `sector = :ns` and k = 2πm/L for `sector = :r` (m = 0, …, L − 1), for rational parameters in
the regime of `in_regime`, where every ε_k ≥ 2(h − |J_x + J_y|) > 0.
"""
function mode_energies(L::Integer, Jx::Rational, Jy::Rational, h::Rational, sector::Symbol)
    in_regime(Jx, Jy, h) || throw(ArgumentError("need J_x J_y < 0 and h > |J_x + J_y|"))
    shift = sector === :ns ? 1 // 2 : sector === :r ? 0 // 1 : throw(ArgumentError("sector"))
    A, B, H = ball(Jx + Jy), ball(Jx - Jy), ball(h)
    twopi = 2 * interval(BigFloat, π)
    energies = map(0:(L-1)) do m
        k = twopi * ball(big(m) + shift) / ball(L)
        2 * sqrt((H - A * cos(k))^2 + B^2 * sin(k)^2)
    end
    all(e -> inf(e) > 0, energies) || error("an enclosure of ε_k is not positive")
    return energies
end

# log(2 cosh(x/2)) and log(2 sinh(x/2)) for x > 0 without overflow: x/2 + log1p(±e^{−x}).
_logcosh2(x) = x / 2 + log1p(exp(-x))
_logsinh2(x) = x / 2 + log(-expm1(-x))

"""
    logZ(L, Jx, Jy, h, β; precision = 512) -> Interval{BigFloat}

Rigorous enclosure of log Z_L(β) for the periodic XY ring of L sites (see the file header), for
rational parameters in the regime of `in_regime` and rational β > 0. The four products are
combined relative to the Neveu–Schwarz cosh product, every term being positive:

    log Z = C + log((1 + e^{S − C} + e^{C_R − C}(1 − e^{S_R − C_R})) / 2).
"""
function logZ(L::Integer, Jx::Rational, Jy::Rational, h::Rational, β::Rational;
    precision::Integer = 512)
    L >= 2 && β > 0 || throw(ArgumentError("need L ≥ 2 and β > 0"))
    setprecision(BigFloat, precision) do
        b = ball(β)
        ns = [b * e for e in mode_energies(L, Jx, Jy, h, :ns)]
        r = [b * e for e in mode_energies(L, Jx, Jy, h, :r)]
        C, S = sum(_logcosh2, ns), sum(_logsinh2, ns)
        CR, SR = sum(_logcosh2, r), sum(_logsinh2, r)
        bracket = 1 + exp(S - C) + exp(CR - C) * (-expm1(SR - CR))
        inf(bracket) > 0 || error("the enclosure of Z_L is not positive")
        C + log(bracket / 2)
    end
end
