"""
    Spin2Gap

The spin-2 Heisenberg chain: the S=2 ITensors site type, the model and its exact, DMRG and
thermal-purification references (`model.jl`, `exact.jl`, `dmrg.jl`, `thermal.jl`), and the
spin-2 certificate submodules

| submodule | content |
|---|---|
| `Trials` | SU(2) trial MPS templates, numerical searches, exact integer energy certificates |
| `Variational` | channel-reduced variational energies and optimisers |
| `Certificate` | reduced (charge-block) exact and ball certificates of trial energies |
| `ThermalBounds` | rational Trotter-layer thermal upper bounds |
| `PolynomialMPO` | the spin-2 polynomial (Taylor, quartic SOS) thermal seeds |

Shared numerics come from `CertifiedGaps` (`Spins`, `ExactDiagonalization`, `Results`,
`Purification`, `CertifiedMPO`).
"""
module Spin2Gap

using ITensors, ITensorMPS
using LinearAlgebra, SparseArrays, Random, Dates, SHA, Serialization
using JSON3
import CertifiedGaps
import CertifiedGaps.Spins
import CertifiedGaps.ExactDiagonalization as ED
import CertifiedGaps.Purification
import CertifiedGaps.Results
using CertifiedGaps.ExactDiagonalization: log_partition
using CertifiedGaps.Results: write_result

include("model.jl")
include("exact.jl")
include("dmrg.jl")
include("thermal.jl")

export spin_matrices, real_generators, bonds, hamiltonian, sector_hamiltonian
export exact_spectrum, log_partition, thermal_curve, seed_estimates, dmrg_levels
export write_result, provenance

function provenance(config)
    root = dirname(@__DIR__)
    hashes = Results.source_hashes(root; files = ("Project.toml", "Manifest.toml"))
    return Results.run_metadata(
        config;
        source_sha256 = hashes,
        packages = (ITensors, ITensorMPS),
        status = "numerical exploration; no rigorous error enclosure",
    )
end

include("Trials.jl")
include("Variational.jl")
include("Certificate.jl")
include("ThermalBounds.jl")
include("PolynomialMPO.jl")

end
