module Provenance
using SHA
export HISTORICAL_SNAPSHOT, source_sha256

# Published, tree-identical equivalent of the complete pre-rename audited snapshot.
const HISTORICAL_SNAPSHOT = "8a1aad63e1a285392c97bf5cf17e8da03f0b8314"

"Hash current bytes, or bytes from an explicitly pinned historical Git snapshot."
function source_sha256(root, path; snapshot=nothing)
    full=normpath(joinpath(root,path))
    if snapshot===nothing
        return bytes2hex(sha256(read(full)))
    end
    repo=normpath(joinpath(root,"../.."))
    relative=relpath(full,repo)
    startswith(relative,"..") && error("Historical source must be inside the repository")
    object=string(snapshot,":",relative)
    return bytes2hex(sha256(read(`git -C $repo show $object`)))
end
end
