# Shared loaders and historical-provenance checks for the stored results.
#
# The results were produced before the core was renamed (ThermalMPOGap -> CertifiedGaps) and
# before this directory moved from certificates/ to models/. Their `source_sha256` records
# therefore describe the producing sources, which are no longer the files in the checkout.
# We verify them against the git objects of the producing commits instead of rewriting them.
using SHA, Test

"Core commit whose `julia/` and `outputs/` trees produced the stored results."
const PRODUCER_CORE_COMMIT = "ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b"
"Commit holding the certificate sources (under `certificates/`) that produced the results."
const PRODUCER_CERTIFICATE_COMMIT = "fd3c3d5592a91c0d37f9f4dfb85048ca5d9cfba2"

"Read a JSON document as text; `.json.gz` files are decompressed with the system `gzip`."
read_json_text(path) = endswith(path, ".gz") ? read(pipeline(`gzip -dc $path`), String) :
    read(path, String)

"Repository root (two levels above this directory)."
repository_root() = normpath(joinpath(@__DIR__, "..", ".."))

function _git_object_available(root, rev)
    return success(pipeline(`git -C $root cat-file -e $(rev * "^{commit}")`;
        stdout=devnull, stderr=devnull))
end

"""
    check_source_hashes(metadata)

Verify `metadata.source_sha256`. Results produced under the former core name (their metadata
lists `ThermalMPOGap`) are checked with `check_historical_hashes`; results regenerated with
`CertifiedGaps` from this directory are checked against the files in the checkout.
"""
function check_source_hashes(metadata)
    key(name) = only(k for k in keys(metadata) if String(k) == name)
    hashes = metadata[key("source_sha256")]
    if any(k -> String(k) == "ThermalMPOGap", keys(metadata))
        check_historical_hashes(hashes)
    else
        root = repository_root()
        for (file, digest) in hashes
            @test bytes2hex(sha256(read(joinpath(root, String(file))))) == String(digest)
        end
    end
end

"""
    check_historical_hashes(source_sha256)

Compare each recorded SHA-256 with the bytes of that path at the producing commit. Paths under
`certificates/` are read from `PRODUCER_CERTIFICATE_COMMIT`, all others from
`PRODUCER_CORE_COMMIT`. A mismatch fails the enclosing testset; if a commit is not present in
the local clone (shallow clone, or the branch was squashed away) the check is skipped with a
warning, since the mathematical content is re-verified independently of these hashes.
"""
function check_historical_hashes(source_sha256)
    root = repository_root()
    for (file, digest) in source_sha256
        path = String(file)
        rev = startswith(path, "certificates/") ? PRODUCER_CERTIFICATE_COMMIT :
            PRODUCER_CORE_COMMIT
        if !_git_object_available(root, rev)
            @warn "Provenance commit $rev not available; skipping hash of $path"
            continue
        end
        bytes = read(`git -C $root show $(rev * ":" * path)`)
        @test bytes2hex(sha256(bytes)) == String(digest)
    end
end
