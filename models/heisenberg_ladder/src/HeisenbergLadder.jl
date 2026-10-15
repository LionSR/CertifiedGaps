"""
    HeisenbergLadder

Gap certificates for the periodic spin-1/2 two-leg Heisenberg ladder (leg coupling 1, rung
coupling `J`, one rung per cell). Three layers:

- `rp.jl`: exact checkerboard reflection-positivity identities and the core `RingModel`.
- `spectra.jl`: symmetry-resolved complete spectra, either verified (interval matrix, Weyl
  residual enclosure) or numerical (for the floating-point scans only).
- `suite.jl`: the certificate: grid classification, rejection witnesses, all-beta seed
  exclusions, coupling continuation, `produce` and `check`.
"""
module HeisenbergLadder

using CertifiedGaps, LinearAlgebra, SparseArrays, Printf, JSON3, IntervalArithmetic
using CertifiedGaps.Results, CertifiedGaps.ExactDiagonalization, CertifiedGaps.Seeds
using CertifiedGaps.RingModels
using CertifiedGaps.ExactArithmetic: round_up, round_down

"Repository root, against which source hashes are keyed."
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
"This model directory."
const BASE = normpath(joinpath(@__DIR__, ".."))

"Exact rationals used for every stored proof number."
const Q = Rational{BigInt}

# The error allowances hold for any BLAS summation order; one thread keeps run times and
# memory comparable with the recorded ones (AGENTS.md: OPENBLAS_NUM_THREADS=1).
__init__() = BLAS.set_num_threads(1)

include("rp.jl")
include("spectra.jl")
include("suite.jl")

end
