"""
    CertifiedGaps.Seeds

Seed search: seed defects from exact spectra, and the certified bound over a grid of seed
lengths and inverse temperatures; and the rigorous step from log Z enclosures to exact rational
seed defects.
"""
module Seeds

using IntervalArithmetic
using ..ExactArithmetic: RationalBig
using ..ExactDiagonalization: logZ_bounds
using ..Recurrence: certified_gap_bound

export SeedResult, seed_defects, spectral_defects, scan_seeds, best_seed

"""
    seed_defects(spectra, n, β) -> (εₛ, εₜ)

Upper bounds on 1 − Pₛ(n, β) and 1 − P_τ(2n, β) from the spectra of the rings of lengths n and
2n (`spectra[L]` is a `PartialSpectrum`), using the partition-function enclosures.
"""
function seed_defects(spectra, n::Integer, β::Real)
    _, lnhi = logZ_bounds(spectra[n], β)
    l2lo, l2hi = logZ_bounds(spectra[2n], β)
    l22lo, _ = logZ_bounds(spectra[2n], 2β)
    return -expm1(l2lo - 2lnhi), -expm1(l22lo - 2l2hi)
end

"""
    seed_defects(logZn, logZ2n, logZ2n_2β) -> (; eps_s, eps_t, defect_s, defect_t)

Rigorous seed defects from interval enclosures (`Interval{Float64}` or `Interval{BigFloat}`)
of log Z_n(β), log Z_{2n}(β) and log Z_{2n}(2β):

    1 − Pₛ(n, β) = 1 − exp(log Z_{2n}(β) − 2 log Z_n(β)),
    1 − P_τ(2n, β) = 1 − exp(log Z_{2n}(2β) − 2 log Z_{2n}(β)).

`defect_s`, `defect_t` are interval enclosures, and `eps_s`, `eps_t` exact rationals (the
binary upper endpoints, converted without rounding) with εₛ ≥ 1 − Pₛ and εₜ ≥ 1 − P_τ. Each
is the smaller of the upper ends of `1 − exp(x)` and `−expm1(x)`, both valid enclosures, so it
never exceeds the plain `1 − exp` evaluation on the same inputs and, with correctly rounded
`expm1` (the default), equals a directed-rounding `−expm1`. No sign is imposed: callers that
need 0 ≤ ε < 1 check it.
"""
function seed_defects(logZn::Interval, logZ2n::Interval, logZ2n_2β::Interval)
    defect_s, eps_s = _defect(logZ2n - 2logZn)
    defect_t, eps_t = _defect(logZ2n_2β - 2logZ2n)
    return (; eps_s, eps_t, defect_s, defect_t)
end

# Enclosure of 1 − eˣ, the intersection of the two enclosures 1 − exp(x) and −expm1(x), and
# its exact rational upper bound.
function _defect(x::Interval)
    d = intersect_interval(1 - exp(x), -expm1(x))
    return d, RationalBig(sup(d))
end

"Seed defects as a function of (n, β), from ring spectra keyed by length."
spectral_defects(spectra) = (n, β) -> seed_defects(spectra, n, β)

"Seed defects at (n₀, β₀) and the resulting gap bound (`nothing` if the criterion fails)."
struct SeedResult
    n₀::Int
    β₀::Float64
    εₛ::Float64
    εₜ::Float64
    bound::Union{Nothing,Float64}
end

"""
    scan_seeds(defects, ns, βs)

`defects(n, β)` returns upper bounds (εₛ, εₜ). Each bound is `certified_gap_bound`, so the
recurrence is outward rounded; the defects are only as reliable as their computation.
"""
function scan_seeds(defects, ns, βs)
    [let (εₛ, εₜ) = defects(n, β)
         SeedResult(n, β, εₛ, εₜ, _bound(εₛ, εₜ, β))
     end for n in ns for β in βs]
end

_bound(εₛ, εₜ, β) = (0 <= εₛ < 1 && 0 <= εₜ < 1) ? certified_gap_bound(εₛ, εₜ, β).bound : nothing

"The grid point with the largest bound, or `nothing` if no seed satisfies the criterion."
function best_seed(results)
    ok = filter(r -> r.bound !== nothing, results)
    isempty(ok) ? nothing : argmax(r -> r.bound, ok)
end

end
