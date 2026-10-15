"""
    LeanSeeds

Toward a Lean-checked seed: five studies sharing one environment, each a submodule.

| submodule | study | content |
|---|---|---|
| `Common` | all | paths, source hashes, Sᶻ basis, Lake driver and axiom audit |
| `SeedDesign` | `seed_design` | smaller seed inputs; exact designs to `ChainCertificate.Valid` |
| `Spectrum` | `spectrum` | Float64 seed spectra, filtering precision, kernel cost model |
| `AnchorStudy` | `anchor` | exact local PSD comparison anchors |
| `Spin1Kernel` | `spin1_kernel` | packed-`Nat` witnesses checked in `lean/SeedKernel/Spin1` |
| `IsingSeed` | `ising` | longitudinal-Ising seed arithmetic checked in `lean/SeedKernel/Ising` |

Scripts load this file and then `using .LeanSeeds.<Submodule>` (each submodule exports all its
names); `certify.jl`, `check.jl` and the tests use qualified names.
"""
module LeanSeeds

include("Common.jl")
include("seed_design.jl")
include("spectrum.jl")
include("anchor.jl")
include("spin1_kernel.jl")
include("ising.jl")

end
