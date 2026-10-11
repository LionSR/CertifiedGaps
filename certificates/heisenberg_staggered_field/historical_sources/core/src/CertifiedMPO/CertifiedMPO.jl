"""
    ThermalMPOGap.CertifiedMPO

Certified matrix-product-operator toolkit in exact rationals and Arb balls (Nemo). Cores are
`Array{RationalBig,4}` indexed (left, output, input, right); numerical routines only propose
candidates, and every accepted bound is replayed in exact or outward-rounded arithmetic.

| submodule | content |
|---|---|
| `ChargeBlocks` | opt-in exact U(1) support, sector products, SVDs and Schmidt gauges |
| `Norms` | local and global operator/HS norm witnesses, gauges, exact core encodings |
| `Squaring` | certified compression of B², Gram flows, squaring witnesses |
| `Products` | certified zip-up compression of A·B and three-layer overlaps |
| `OperatorError` | telescoping operator-norm error of a compressed square or product |
| `GaugeSearch` | numerical gauge proposals (Schmidt balancing, log-gauge descent) |
| `AdjointTrace` | conditional trace enclosure for approximate Gram/square sequences |
| `SchattenBridge` | operator error and trace moment to a Schatten-p trace bound |
| `PositiveBounds` | exact entrywise nonnegative majorants M ≤ U V |
| `PolynomialBounds` | Taylor thermal majorants with an explicit spectral radius |
| `MPOAlgebra` | exact sums, products, exact compression and sparse expansion of MPOs |

Moved verbatim from spin2 (`Spin2MPONorms`, `Spin2MPOSquaring`, `Spin2MPOProducts`,
`Spin2MPOOperatorError`, `Spin2MPOGaugeSearch`, `Spin2AdjointTrace`, `Spin2SchattenBridge`,
`Spin2PositiveBounds`, `Spin2PolynomialBounds`, and the generic part of `Spin2PolynomialMPO`).
Each submodule keeps the spin2 aliases (`F`, `MN`, `MS`) so that spin2 call sites only change
the module they bind. Where spin2 hard-coded the spin-2 local dimension, weights or spectral
radius, the core takes an argument whose default is the spin-2 value.
"""
module CertifiedMPO

import Nemo
using ..ExactArithmetic: RationalBig, encode_fraction, parse_fraction

"Exact rationals, as in every spin2 certificate (`Rational{BigInt}`)."
const F = RationalBig

"""
    nemo_array(A, field = Nemo.QQ)

Entries of a Nemo matrix as a Julia matrix: `Rational{BigInt}` over `Nemo.QQ`, and the ring's
own elements (Arb balls) over any other field. Shared by every submodule; it replaces eight
identical local converters of the spin2 sources.
"""
function nemo_array(A, field = Nemo.QQ)
    field === Nemo.QQ && return [F(A[i, j]) for i in 1:Nemo.nrows(A), j in 1:Nemo.ncols(A)]
    return [A[i, j] for i in 1:Nemo.nrows(A), j in 1:Nemo.ncols(A)]
end

include("ChargeBlocks.jl")
include("Norms.jl")
include("Squaring.jl")
include("Products.jl")
include("OperatorError.jl")
include("GaugeSearch.jl")
include("AdjointTrace.jl")
include("SchattenBridge.jl")
include("PositiveBounds.jl")
include("PolynomialBounds.jl")
include("MPOAlgebra.jl")

end
