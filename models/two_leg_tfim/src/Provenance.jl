module Provenance
using SHA
using CertifiedGaps.Results: parse_exact_json
export HISTORICAL_SNAPSHOT, HISTORICAL_DIR, ORIGINAL_CERTIFICATE, ORIGINAL_CERTIFICATE_SHA256,
    snapshot_available, source_sha256, certificate_bytes, read_certificate

# Published, tree-identical equivalent of the complete pre-rename audited snapshot.
const HISTORICAL_SNAPSHOT = "8a1aad63e1a285392c97bf5cf17e8da03f0b8314"
# Where this directory lived in that snapshot (it now lives in models/two_leg_tfim).
const HISTORICAL_DIR = "certificates/two_leg_tfim"
# The original certificate is stored gzip-compressed; its uncompressed bytes are unchanged.
const ORIGINAL_CERTIFICATE = joinpath("results", "certificate_g3_1.json.gz")
const ORIGINAL_CERTIFICATE_SHA256 =
    "aba98d9865f1b496a96e7a7e17dec0239ec07b34e012b73714d373c1e8a74928"

"Uncompressed bytes of a JSON result, decompressing `.gz` files with the system `gzip`."
certificate_bytes(path) = endswith(path, ".gz") ? read(`gzip -dc $path`) : read(path)

"`parse_exact_json` of a (possibly gzip-compressed) result file."
read_certificate(path) = parse_exact_json(String(certificate_bytes(path)))

"Whether the pinned historical commit is present in the local Git object store."
function snapshot_available(root)
    repo = normpath(joinpath(root, "../.."))
    object = HISTORICAL_SNAPSHOT * "^{commit}"
    return success(pipeline(`git -C $repo cat-file -e $object`; stderr=devnull))
end

"""
Hash current bytes, or bytes from an explicitly pinned historical Git snapshot.

`path` is relative to this certificate directory. Historical paths are resolved relative to
`HISTORICAL_DIR`, the directory's location when the hashes were recorded.
"""
function source_sha256(root, path; snapshot=nothing)
    if snapshot === nothing
        return bytes2hex(sha256(read(normpath(joinpath(root, path)))))
    end
    repo = normpath(joinpath(root, "../.."))
    relative = normpath(joinpath(HISTORICAL_DIR, path))
    startswith(relative, "..") && error("Historical source must be inside the repository")
    object = string(snapshot, ":", relative)
    return bytes2hex(sha256(read(`git -C $repo show $object`)))
end
end
