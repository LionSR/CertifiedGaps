"""
    ThermalMPOGap.Seeds

Seed search: seed defects from exact spectra, and the certified bound over a grid of seed
lengths and inverse temperatures.
"""
module Seeds

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
