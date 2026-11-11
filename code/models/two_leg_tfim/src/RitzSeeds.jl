"""
    RitzSeeds

Coupled seeds of the ladder at n₀ = 6 rungs: the complete verified spectrum of the 6-rung ring
(core `verified_ring_spectrum`) and rigorous Rayleigh–Ritz upper bounds on the low levels of
every block of the 12-rung ring (24 spins; core `RitzBounds.partial_ring_spectrum`), combined
by the core `GapCertificates.coupled_seed_bounds`:

    εₛ ≥ 1 − Z₁₂(β₀)⁻ / Z₆(β₀)²,    εₜ ≥ 1 − Z₁₂(2β₀)⁻ / Z₆(β₀)⁴.

Only lower bounds on Z₁₂ are needed, so states and blocks may be omitted. In the rotated frame
the Z₂ charge is the global spin-flip parity of the physical ladder: charge 0 holds the ground
state and the two-quasiparticle states, charge 1 the one- and three-quasiparticle states. The
floating-point scan (`scripts/ritz_scan.jl`) shows that at the seed temperature the bound is
converged to 10⁻⁸ with 4 levels per even and 2 per odd block; a few more are bounded here.

The certified points are listed in `src/ritz_points.jl` (`RITZ_POINTS`), outside the hash-pinned
sources of the certificates (`SOURCES`), so that adding a point does not touch the
pinned code of the others; each certificate stores its own seed, which `check.jl` compares.
"""
module RitzSeeds

using CertifiedGaps.RitzBounds: partial_ring_spectrum
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, spectral_gap_enclosure
using ..Ladder: ladder

export RITZ_N0, RITZ_POINTS, SOURCES, ritz_levels, doubled_ring_bounds, ritz_point_path,
    partial_spectrum_path, BRIDGE_LENGTHS, finite_bridge

"Seed length (rungs) of the coupled seeds; the doubled ring has 2n₀ = 12 rungs (24 spins)."
const RITZ_N0 = 6

include("ritz_points.jl")

"""
Model-local sources pinned by every certificate of this directory: the model, this module
and the drivers (`src/ritz_points.jl`, the list of points, is deliberately not pinned).
"""
const SOURCES = (dirs = (), files = ("src/Ladder.jl", "src/RitzSeeds.jl", "certify.jl",
    "check.jl"))

"Levels bounded per block of the 12-rung ring: more in the even charge (see the module docs)."
ritz_levels(label) = label["charge"] == 0 ? 6 : 3

"""
    doubled_ring_bounds(g; progress = false) -> Vector{PartialBlock}

Rigorous upper bounds on the low levels of every block of the 12-rung ring at field g
(charge, momentum, leg-exchange parity; real blocks from the lattice reflection), from
Lanczos vectors with deterministic start vectors turned into bounds by the core.
"""
doubled_ring_bounds(g; progress = false) = partial_ring_spectrum(ladder(g), 2RITZ_N0;
    howmany = ritz_levels, dense_max = 1000, tol = 1e-8, krylovdim = 32, progress)

"""
Even ring lengths below n₀ covered by `finite_bridge`, so that the coupled certificates cover
every even L ≥ 4 like the n₀ = 4 certificates (and the Lean theorem, whose L₀ is 4).
"""
const BRIDGE_LENGTHS = 4:2:(RITZ_N0-2)

"""
    finite_bridge(g) -> Vector{Dict}

For every L in `BRIDGE_LENGTHS`, exact rational bounds `gap_lower` ≤ E₁ − E₀ ≤ `gap_upper` on
the gap of the L-rung ring at field g from its complete verified spectrum
(`spectral_gap_enclosure`); `gap_lower > 0` proves a unique ground state. The certificate
stores them, and `check.jl` recomputes them and requires `gap_lower` ≥ the certified bound.
"""
function finite_bridge(g)
    map(BRIDGE_LENGTHS) do L
        lo, hi = spectral_gap_enclosure(verified_ring_spectrum(ladder(g), L))
        grid = big(10)^15                        # outward to 10⁻¹⁵, as the seed defects
        lo, hi = floor(BigInt, lo * grid) // grid, ceil(BigInt, hi * grid) // grid
        Dict{String,Any}("L" => L, "gap_lower" => string(lo), "gap_upper" => string(hi),
            "gap_lower_decimal" => Float64(lo, RoundDown))
    end
end

_tag(g) = "g$(numerator(g))_$(denominator(g))"

"Certificate file of the coupled seed at field g (results/certificate_g3_1_n6.json)."
ritz_point_path(root, g) = joinpath(root, "results", "certificate_$(_tag(g))_n6.json")

"Stored Ritz upper bounds of the 12-rung ring at field g (results/partial_g3_1_L12.json)."
partial_spectrum_path(root, g) = joinpath(root, "results", "partial_$(_tag(g))_L12.json")

end
