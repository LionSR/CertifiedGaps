# Outward-rounded evaluation of the doubling recurrence.
#
# The survey's `tmlib.uniform_gap_bound` returns b_60 in floating point, an upper bound on b_∞,
# and notes that the rigorous tail correction "is not yet implemented". Here the orbit is enclosed
# in interval arithmetic on the logarithms log s_j, log t_j, which stay representable even when the
# defects are far below the double-precision range. The certificate combines:
#   * admissibility: some s_j, t_j < 0.12731 < ε⋆ puts the orbit in Q⋆ ⊂ A (Theorem 2.18(1));
#   * monotonicity: on A the rates b_j are nonincreasing and b_∞ = inf_j b_j (Theorem 2.20);
#   * the tail: r_J ≤ 1/10 gives b_J − b_∞ ≤ 2^{-J} log 5 / β₀ (Lemma 2.21(ii)).
# R is order-preserving (Lemma 2.16), so upper bounds on the seed defects give a valid bound.

export certified_gap_bound, to_interval

const LOG_HALF = log(interval(1) / 2)
const LOG_ZSTAR = log(interval(2) - sqrt(interval(2)))       # f(z) < 1 iff z < 2 − √2
const LOG_EPS_LOW = log(interval(12731) / 100000)            # 0.12731 < ε⋆ (Theorem 2.17)
const LOG_TENTH = log(interval(1) / 10)

# log z(a, b) from log a, log b. For tiny defects z = 2a + b − (a² + 2ab − a²b) is enclosed by
# log(2a + b) + [log(1 − a − b), 0]; otherwise use the exact cancellation-free formula.
# Deep in the orbit the log-intervals are wide in absolute terms (about 10³ at 10¹⁶), so
# log(2a + b), which is increasing in both arguments, is evaluated separately at the lower and
# the upper endpoints rather than on the intervals themselves.
function _logz(la::Interval, lb::Interval)
    if sup(max(la, lb)) < -50
        lower = inf(_log2a_plus_b(inf(la), inf(lb)))
        upper = sup(_log2a_plus_b(sup(la), sup(lb)))
        # log(1 − a − b) > −1e-20 since a + b < 4e-22, which is below one ulp of |log z| > 50
        return interval(prevfloat(lower), upper)
    end
    return log(zdefect(exp(la), exp(lb)))
end

# Enclosure of log(2eˣ + eʸ) for floats x, y. The shift c is the upper end of the larger exponent,
# so every argument of exp is ≤ 0 even when an ulp of x is large.
function _log2a_plus_b(x::Float64, y::Float64)
    u, v = log(interval(2)) + x, interval(y)
    c = sup(max(u, v))
    return c + log(exp(u - c) + exp(v - c))
end

# log of a seed defect. A zero defect is raised to floatmin: R is order-preserving, so a larger
# seed still gives a valid (weaker) bound, and log stays finite.
_log_defect(ε) = log(max(to_interval(ε), interval(floatmin(Float64))))

_logf(lz::Interval) = 2lz - log(interval(2)) - 2log1p(-exp(lz))    # log f(z)
_logr(lz::Interval) = lz - log1p(-exp(lz))                         # log(z/(1 − z))

"""
    certified_gap_bound(εₛ, εₜ, β₀; J = 48)

Rigorous lower bound on Δ_L for every even L ≥ n₀, given Pₛ(n₀,β₀) ≥ 1 − εₛ and
P_τ(2n₀,β₀) ≥ 1 − εₜ under the moment form (Theorem 2.25). Pass the defects as exact rationals,
intervals or floats (a float is taken to be exactly the binary number it denotes).

Returns `(bound, upper, j_admissible, log_rJ)`: `bound` is a certified lower bound (`nothing` if
certification fails), `upper` the enclosure's upper end of b_J ≥ b_∞, `j_admissible` the first
scale with the orbit inside Q⋆, and `log_rJ` an upper bound on log r_J.

The default J = 48 makes the tail allowance 2^{-J} log 5/β₀ < 6e-15/β₀ while |log r_J| stays
well inside the range where a Float64 ulp is small.
"""
function certified_gap_bound(εₛ, εₜ, β₀; J = 48)
    failed = (bound = nothing, upper = nothing, j_admissible = nothing, log_rJ = nothing)
    ls, lt = _log_defect(εₛ), _log_defect(εₜ)
    j_adm = nothing
    for j in 0:J
        lz = _logz(ls, lt)
        sup(lz) < inf(LOG_HALF) || return failed
        if j_adm === nothing && sup(max(ls, lt)) < inf(LOG_EPS_LOW)
            j_adm = j
        end
        if j == J
            return _certify_tail(lz, to_interval(β₀), J, j_adm)
        end
        # One step of R; each argument of f must stay below 2 − √2.
        sup(lz) < inf(LOG_ZSTAR) || return failed
        ls′ = _logf(lz)
        lz₂ = _logz(lt, ls′)
        sup(lz₂) < inf(LOG_ZSTAR) || return failed
        ls, lt = ls′, _logf(lz₂)
    end
end

function _certify_tail(lz, β, J, j_adm)
    lr = _logr(lz)
    βJ = β * interval(2)^J
    bJ = -lr / βJ
    tail = log(interval(5)) / βJ
    ok = j_adm !== nothing && sup(lr) <= inf(LOG_TENTH)
    return (bound = ok ? inf(bJ - tail) : nothing, upper = sup(bJ), j_admissible = j_adm,
        log_rJ = sup(lr))
end

"Float64 interval enclosing `x` (exact for rationals, integers and floats)."
to_interval(x::Interval{Float64}) = x
to_interval(x::Interval{BigFloat}) = interval(Float64(inf(x), RoundDown), Float64(sup(x), RoundUp))
to_interval(x::Rational) = to_interval(interval(BigFloat, numerator(x)) / interval(BigFloat, denominator(x)))
to_interval(x::BigFloat) = to_interval(interval(x))
to_interval(x::Real) = interval(Float64(x))
