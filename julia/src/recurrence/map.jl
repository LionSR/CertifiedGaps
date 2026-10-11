# The exact-constant doubling recurrence (Theorem 2.4) and its scalar constants (§2.6).
#
# Seed pair (s, t) = (εₛ, εₜ). One doubling step is the map R(s, t) = (s', t') with
#     z₁ = z(s, t),  s' = f(z₁),  z₂ = z(t, s'),  t' = f(z₂),
#     z(a, b) = 1 - (1 - b)(1 - a)²,   f(z) = z² / (2(1 - z)²).
# At scale j the thermal purity of every even L ∈ [n_j, 2n_j] is at least 1 - z_j, where
# z_j = z(s_j, t_j), so Δ_L ≥ b_j = log((1 - z_j)/z_j) / β_j whenever z_j < 1/2.

export squaring_defect, zdefect, doubling_map, doubling_orbit, gap_rates, gap_bound,
    diagonal_defect, diagonal_fixed_point, fixed_point_bracket_holds, majorant_loss

"Squaring bound f(z) = z²/(2(1-z)²): the defect after squaring a spectrum of defect z."
squaring_defect(z) = z^2 / (2 * (1 - z)^2)

"""
    diagonal_defect(ε) = ((1−ε)^{−3} − 1)² / 2

First component of one doubling step from the diagonal seed (ε, ε). Its fixed point is ε⋆
(`diagonal_fixed_point`); the scope note writes it δ(ε).
"""
diagonal_defect(ε) = squaring_defect(1 - (1 - ε)^3)

"z(a, b) = 1 - (1-b)(1-a)², written without cancellation."
zdefect(a, b) = -expm1(log1p(-b) + 2 * log1p(-a))

"One step of the exact recurrence; `nothing` when an argument of f leaves [0, 2-√2)."
function doubling_map(s, t)
    zstar = 2 - sqrt(oftype(float(s), 2))
    z1 = zdefect(s, t)
    z1 < zstar || return nothing
    s′ = squaring_defect(z1)
    z2 = zdefect(t, s′)
    z2 < zstar || return nothing
    return (s′, squaring_defect(z2))
end

"""
    doubling_orbit(s, t; steps = 60)

Orbit `[(s₀,t₀), (s₁,t₁), …]` of the exact map, in the number type of `s` and `t`.
Use `BigFloat` for more than a handful of steps: the defects fall doubly exponentially.
Returns `nothing` if the orbit leaves the domain of the map.
"""
function doubling_orbit(s, t; steps = 60)
    orbit = [(s, t)]
    for _ in 1:steps
        next = doubling_map(orbit[end]...)
        next === nothing && return nothing
        push!(orbit, next)
    end
    return orbit
end

"Per-scale rates b_j = log((1-z_j)/z_j)/β_j along an orbit (−Inf when z_j ≥ 1/2)."
function gap_rates(orbit, β₀)
    map(enumerate(orbit)) do (j, (s, t))
        z = zdefect(s, t)
        z < 1 // 2 ? log((1 - z) / z) / (β₀ * 2^(j - 1)) : oftype(z, -Inf)
    end
end

"""
    gap_bound(εₛ, εₜ, β₀; steps = 60, precision = 4096)

Uniform gap bound min_j b_j from seed defects, in `BigFloat` arithmetic (not outward rounded;
see `certified_gap_bound` for the rigorous version). For seeds in the admissible region the
rates are nonincreasing (Theorem 2.20), so this returns b_steps, which exceeds the limit b_∞ by at
most 2^{-steps} log 5 / β₀ once r_steps ≤ 1/10 (Lemma 2.21(ii)).
"""
function gap_bound(εₛ, εₜ, β₀; steps = 60, precision = 4096)
    setprecision(BigFloat, precision) do
        orbit = doubling_orbit(BigFloat(εₛ), BigFloat(εₜ); steps)
        orbit === nothing && return nothing
        Float64(minimum(gap_rates(orbit, BigFloat(β₀))))
    end
end

"""
    diagonal_fixed_point(; precision = 256)

The diagonal threshold ε⋆ (Theorem 2.17): the positive root of (1-ε)^{-3} = 1 + √(2ε),
found by bisection. Diagonal seeds ε < ε⋆ contract to zero; seeds ≥ (ε⋆, ε⋆) do not.
"""
function diagonal_fixed_point(; precision = 256)
    setprecision(BigFloat, precision) do
        G(ε) = ((1 - ε)^-3 - 1)^2 - 2ε
        lo, hi = BigFloat("0.1"), BigFloat("0.2")
        for _ in 1:precision
            mid = (lo + hi) / 2
            G(mid) < 0 ? (lo = mid) : (hi = mid)
        end
        (lo + hi) / 2
    end
end

"Exact rational check of the bracket 0.12731 < ε⋆ < 0.12732: G changes sign on it."
function fixed_point_bracket_holds(lo = 12731 // 100000, hi = 12732 // 100000)
    G(ε) = ((1 - ε)^-3 - 1)^2 - 2ε
    G(big(lo)) < 0 < G(big(hi))
end

"""
    majorant_loss(R₀, P₀; J = 40)

The majorant recursion Φ of Lemma 2.21(i), started at (R₀, P₀) ≥ (r₀, ρ₀). Returns
`(Λ_J, tail_ok)`, where Λ_J bounds the loss Λ in β₀ b_∞ = log(1/r₀) − Λ, and `tail_ok` is the
hypothesis Φ(R_{J-1}, P_{J-1}) ≤ (R_{J-1}, P_{J-1}) under which the bound is valid.
Works for any number type that supports the arithmetic, including intervals.
"""
function majorant_loss(R₀, P₀; J = 40)
    function Φ(R, P)
        U = R * (2P + R / 2)
        Q = (2P + R / 2)^2 / (2 * (1 - U)^2)
        T = Q * R^2
        η = (1 + Q) / (1 - R^2 - T) - 1
        P′ = min(one(P), Q / ((1 - R^2 / 4) + Q * (1 - R^2 / 2)^2))
        return R^2 * (1 + η), P′, η
    end
    R, P = R₀, P₀
    Λ = zero(R₀)
    η = zero(R₀)
    Rprev, Pprev = R, P
    for j in 0:(J - 1)
        Rprev, Pprev = R, P
        R, P, η = Φ(R, P)
        Λ += log(1 + η) / 2^(j + 1)
    end
    Λ += log(1 + η) / 2^J
    # Tail hypothesis: Φ(R_{J-1}, P_{J-1}) = (R_J, P_J) ≤ (R_{J-1}, P_{J-1}).
    tail_ok = _le(R, Rprev) && _le(P, Pprev)
    return Λ, tail_ok
end

_le(a::Real, b::Real) = a <= b
_le(a::Interval, b::Interval) = sup(a) <= inf(b)
