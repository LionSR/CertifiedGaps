# Shared loaders and source-hash checks for the stored results.
using SHA, Test

"Read a JSON document as text; `.json.gz` files are decompressed with the system `gzip`."
read_json_text(path) = endswith(path, ".gz") ? read(pipeline(`gzip -dc $path`), String) :
    read(path, String)

"Repository root (three levels above this src/ directory)."
repository_root() = normpath(joinpath(@__DIR__, "..", "..", ".."))

"""
    check_source_hashes(metadata)

Every SHA-256 in `metadata["source_sha256"]` must match the file at that repository path. All
results in this directory were regenerated from the current sources, so a mismatch means that a
source changed after its result was written (rerun the matching `certify.jl` step).
"""
function check_source_hashes(metadata)
    key(name) = only(k for k in keys(metadata) if String(k) == name)
    hashes = metadata[key("source_sha256")]
    @test !isempty(hashes)
    root = repository_root()
    for (file, digest) in hashes
        @test bytes2hex(sha256(read(joinpath(root, String(file))))) == String(digest)
    end
end
