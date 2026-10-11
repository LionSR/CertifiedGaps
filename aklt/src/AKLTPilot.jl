module AKLTPilot

using ITensors, ITensorMPS, LinearAlgebra
using ThermalMPOGap.Spins: spin_matrices
using ThermalMPOGap.ExactDiagonalization: Chain, graph_bonds, spectrum, log_partition
using ThermalMPOGap.Purification: purified_log_partition
using ThermalMPOGap.Results: write_result, run_metadata, source_hashes, file_hashes

include("model.jl")
include("exact.jl")
include("thermal.jl")

export honeycomb, bond_projector, full_spectrum
export thermal_curve, log_partition, diagnostics, write_result, provenance

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
