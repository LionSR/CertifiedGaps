# Seeds at n₀ = 6 cells from the coupling identity: the complete verified spectrum of the ring of
# 6 cells (12 spins) and rigorous Rayleigh–Ritz upper bounds on the low levels of the ring of 12
# cells (24 spins), from the core `RitzBounds`. See IDEAS.md (route 7) and METHOD.md.
module RitzSeeds

using CertifiedGaps.ExactArithmetic: RationalBig, round_up
using CertifiedGaps.RitzBounds: partial_ring_spectrum
using ..DimerizedHeisenberg: rotated_model

export RITZ_N0, RITZ_CHARGES, RITZ_SEEDS, ritz_levels, doubled_ring_bounds, padded_seed

const Q = RationalBig

"Seed length (cells) of the coupled seeds; the doubled ring has 2n₀ = 12 cells (24 spins)."
const RITZ_N0 = 6

"Charges (total Sᶻ) of the doubled ring whose blocks are bounded; −q is covered by the flip."
const RITZ_CHARGES = 0:2

"Levels bounded per block at |Sᶻ| = 0, 1, 2 (more levels change the defects below 10⁻⁵)."
ritz_levels(label) = (16, 16, 10)[abs(label["charge"])+1]

"""
The coupled seeds β₀ at n₀ = 6 cells per dimerization: the best β₀ on the grid 1/8 of a scan
of the interval bound over the stored Ritz bounds (`scripts/ritz_scan.jl`).
"""
const RITZ_SEEDS = Dict{Q,Q}(1 // 2 => 8 // 1, 3 // 10 => 65 // 8, 1 // 4 => 17 // 2,
    1 // 5 => 19 // 2)

"""
    doubled_ring_bounds(δ; progress = false) -> Vector{PartialBlock}

Rigorous upper bounds on the low levels of every block of the ring of 2n₀ = 12 cells with
charge in `RITZ_CHARGES` (and their flips), `ritz_levels` levels per block, with the core's
deterministic Lanczos start vectors. Blocks of dimension ≤ 1000 are enclosed completely.
"""
doubled_ring_bounds(δ; progress = false) = partial_ring_spectrum(rotated_model(δ), 2RITZ_N0;
    howmany = ritz_levels, charges = RITZ_CHARGES, dense_max = 1000, progress)

"Grid of the stored seed defects: a replay may differ from the stored bounds in the last bits."
const PAD = big(10)^9

"""
    padded_seed(seed) -> NamedTuple

The seed with εₛ, εₜ raised by 10⁻⁹ and rounded up to the grid 10⁻⁹ (still upper bounds), so
that the defects of a `check.jl --full` replay stay below the stored ones: the LAPACK
eigenvectors behind `verified_eigen` radii and the Lanczos vectors behind the Ritz bounds vary
between runs (observed: 4·10⁻¹¹ in εₛ at n₀ = 4), while the core grid is 10⁻¹⁵. Used for every
seed of this model.
"""
padded_seed(seed) = merge(seed, (; eps_s = round_up(seed.eps_s + 1 // PAD, PAD),
    eps_t = round_up(seed.eps_t + 1 // PAD, PAD)))

end
