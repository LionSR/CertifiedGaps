"""
    HeisenbergLadder

Gap certificates for the periodic spin-1/2 two-leg Heisenberg ladder (leg coupling 1, rung
coupling `J`, one rung per cell), built on the shared core pipeline of `CertifiedGaps`:

- `model.jl`: the ladder as an exact core `RingModel` in the leg-rotated frame (reflection
  positive with the cell reflection σʸ ⊗ σʸ, conserved Sᶻ, charge flip and leg exchange), its
  RP certificate and the exact identities relating it to the physical ladder.
- `suite.jl`: per coupling, core verified spectra and seed defects on the (n₀, β₀) grid, a core
  `certify_gap` certificate at the best seed, and the model-specific rejection witnesses,
  all-β exclusions, coupling continuation, `produce` and `check`.
"""
module HeisenbergLadder

using CertifiedGaps, LinearAlgebra, JSON3, IntervalArithmetic
using CertifiedGaps.Results, CertifiedGaps.RingModels, CertifiedGaps.VerifiedSpectra
using CertifiedGaps.GapCertificates
using CertifiedGaps.ExactArithmetic: round_up, round_down
using CertifiedGaps.Spins: exact_spin_half

"This model directory."
const BASE = normpath(joinpath(@__DIR__, ".."))

"Exact rationals used for every stored proof number."
const Q = Rational{BigInt}

# The error allowances hold for any BLAS summation order; one thread keeps run times and
# memory comparable with the recorded ones (AGENTS.md: OPENBLAS_NUM_THREADS=1).
__init__() = BLAS.set_num_threads(1)

include("model.jl")
include("suite.jl")

end
