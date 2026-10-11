# Fetch the source paper's finite-verification evidence, pinned by commit and checked by hash.
#
# spin1/ reproduces the certificates of "The periodic spin-one Haldane gap" (OpenAI, September
# 24, 2026) and checks its results against the paper's own evidence files. Those files are not
# stored in this repository: this script downloads exactly the files spin1 verifies from
# github.com/openai/math at a fixed commit, into references/periodic-gap-source/ (git-ignored).
#
# Integrity chain: the downloaded source-manifest.json must have the SHA-256 that spin1 pins
# (MANIFEST_SHA in models/spin1/src/Spin1Certificates.jl); every other file must match the hash the
# manifest lists for it. Any mismatch aborts.
#
# Run from the repository root (no project environment needed):
#     julia tools/fetch_upstream.jl

using Downloads, SHA

const REPOSITORY = "openai/math"
const COMMIT = "adc7f1241b42e322a6451854ab7e4b4c146bf78a"
const PAPER = "preprints/The-periodic-spin-one-Haldane-gap-September-24-2026"
const MANIFEST_SHA = "2d2461e8f87ef9c3af62ce226412b106fcfb34f01d5139683a90708859e2dc92"
const TARGET = normpath(joinpath(@__DIR__, "..", "references", "periodic-gap-source"))

raw_url(path) = "https://raw.githubusercontent.com/$REPOSITORY/$COMMIT/$PAPER/$path"
sha256_hex(path) = bytes2hex(sha256(read(path)))

"Download `path` (relative to the paper directory) unless a copy with the expected hash exists."
function fetch(path, expected)
    destination = joinpath(TARGET, path)
    isfile(destination) && sha256_hex(destination) == expected && return :present
    mkpath(dirname(destination))
    Downloads.download(raw_url(path), destination)
    actual = sha256_hex(destination)
    actual == expected || error("hash mismatch for $path: expected $expected, got $actual")
    return :downloaded
end

# The manifest is pinned by spin1; it lists the hash of every other file.
manifest_path = "verification/computations/source-manifest.json"
fetch(manifest_path, MANIFEST_SHA)
manifest_text = read(joinpath(TARGET, manifest_path), String)

# The manifest is a flat JSON object of {"name": "sha256"} groups; read the two groups spin1
# verifies, and the hash of the top-level verifier, without a JSON dependency.
function group(name)
    block = match(Regex("\"$name\"\\s*:\\s*\\{([^}]*)\\}"), manifest_text)
    block === nothing && error("manifest has no group $name")
    return [(m[1], m[2]) for m in eachmatch(r"\"([^\"]+)\"\s*:\s*\"([0-9a-f]{64})\"", block[1])]
end
verify_py = match(r"\"verify_py_sha256\"\s*:\s*\"([0-9a-f]{64})\"", manifest_text)[1]

entries = vcat(
    [("verification/computations/" * name, hash) for (name, hash) in group("files")],
    [("verification/computations/" * name, hash)
     for (name, hash) in group("supplementary_certificate_files")],
    [("verification/verify.py", verify_py)],
)
counts = Dict(:present => 0, :downloaded => 0)
for (path, hash) in entries
    counts[fetch(path, hash)] += 1
end
println("upstream evidence at $REPOSITORY@$(COMMIT[1:12]): ", length(entries) + 1,
    " files verified ($(counts[:downloaded]) downloaded, $(counts[:present]) already present)")
