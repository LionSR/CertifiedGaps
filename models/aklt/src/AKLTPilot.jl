module AKLTPilot

using ITensors, ITensorMPS, LinearAlgebra
using CertifiedGaps.Spins: spin_matrices
using CertifiedGaps.ExactDiagonalization: Chain, graph_bonds, spectrum, log_partition
using CertifiedGaps.Purification: purified_log_partition
using CertifiedGaps.Results: write_result, run_metadata, source_hashes, file_hashes

include("model.jl")
include("exact.jl")
include("thermal.jl")

export honeycomb, bond_projector, full_spectrum
export thermal_curve, log_partition, diagnostics, write_result, provenance

"""
    provenance(config)

Run metadata for an output: `config`, package versions, and the SHA-256 of this package's
sources and of the active environment (which need not be `models/aklt`), with `certified=false`.
"""
function provenance(config)
    root = dirname(@__DIR__)
    environment = dirname(Base.active_project())
    hashes = merge(
        source_hashes(root),
        file_hashes(environment, ("Project.toml", "Manifest.toml"); prefix = "environment/"),
    )
    return run_metadata(
        config;
        source_sha256 = hashes,
        packages = (ITensors, ITensorMPS),
        status = "Floating-point exploration; no rigorous error enclosure",
    )
end

end
