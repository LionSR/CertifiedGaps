# Source-hash checks for the stored results (reading them is `CertifiedGaps.Results.read_result`).
using Test
using CertifiedGaps.Results: verify_source_hashes

"Repository root (three levels above this src/ directory)."
repository_root() = normpath(joinpath(@__DIR__, "..", "..", ".."))

"""
    check_source_hashes(metadata)

Every SHA-256 in `metadata["source_sha256"]` must match the file at that repository path, and
every such path must lie in this model directory. All results here were regenerated from the
current sources, so a mismatch means that a model-local source changed after its result was
written (rerun the matching `certify.jl` step). The shared core is recorded separately under
`metadata["core_sha256"]` as provenance only and is not compared with the checkout.
"""
function check_source_hashes(metadata)
    key(name) = only(k for k in keys(metadata) if String(k) == name)
    hashes = metadata[key("source_sha256")]
    @test !isempty(hashes)
    @test all(startswith(String(file), "models/dimerized_heisenberg/") for file in keys(hashes))
    @test verify_source_hashes(hashes, repository_root();
        scope = "models/dimerized_heisenberg/") == length(hashes)
end
