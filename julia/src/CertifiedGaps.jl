"""
    CertifiedGaps

Shared numerics for the thermal-MPO gap criterion: a ring with spatial and thermal purities

    Pₛ(n, β) = Z_{2n}(β) / Z_n(β)²,        P_τ(n, β) = Z_n(2β) / Z_n(β)²

has a uniform gap Δ_L ≥ b_∞ for every even L ≥ n₀ when its transfer operator is self-adjoint
(reflection positivity) and the seed defects at (n₀, β₀) lie in the admissible region.

Submodules, each used on its own (`using CertifiedGaps.ExactDiagonalization`):

| submodule | content |
|---|---|
| `Recurrence` | the doubling map and the uniform bound, in float, interval and rational arithmetic (re-exported here) |
| `ExactArithmetic` | outward-rounded rational helpers: logarithm and power enclosures, grids, exact tokens |
| `Results` | atomic JSON output, source fingerprints, run metadata, exact JSON parsing |
| `Spins` | spin operators in the standard, integer-ladder and Cartesian bases |
| `ReflectionPositivity` | the realignment test and decomposition −h = Σ C ⊗ C̄ |
| `ExactDiagonalization` | chains on any bond graph, charge sectors, spectra, log Z enclosures |
| `RingModels` | exact ring models and the exact RP certificate −h = Σ w C ⊗ C̄ |
| `VerifiedSpectra` | rigorous complete spectra of rings by symmetry blocks, and log Z enclosures |
| `GapCertificates` | `certify_gap` and `check_certificate`: end-to-end certificates from small seeds |
| `Models` | TFIM (with its exact log Z), spin-s Heisenberg, staggered-field Heisenberg |
| `Seeds` | seed defects from spectra and scans over (n₀, β₀) |
| `Purification` | thermal purification with ITensors (load ITensors and ITensorMPS to enable) |
| `CertifiedMPO` | certified MPO norms, squares, products and operator errors (rationals, Arb) |

Theorem numbers in docstrings refer to `survey/theory/STRENGTHENING.md`.
"""
module CertifiedGaps

include("ExactArithmetic.jl")
include("Results.jl")
include("Spins.jl")
include("ReflectionPositivity.jl")
include("ExactDiagonalization.jl")
include("Recurrence.jl")
include("Seeds.jl")
include("RingModels.jl")
include("VerifiedSpectra.jl")
include("GapCertificates.jl")
include("Models.jl")
include("Purification.jl")
include("CertifiedMPO/CertifiedMPO.jl")

using .Recurrence
for name in names(Recurrence)
    name === :Recurrence || @eval export $name
end

end
