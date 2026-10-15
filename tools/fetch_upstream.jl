# Fetch files of "The periodic spin-one Haldane gap" (OpenAI, September 24, 2026) from
# github.com/openai/math at a pinned commit, checked by hash. Run from the repository root (no
# project environment needed); network access required.
#
#     julia tools/fetch_upstream.jl                  # evidence for models/spin1
#     julia tools/fetch_upstream.jl --paper [NAME…]  # paper sources for models/lean_seeds
#
# Default: the finite-verification evidence that models/spin1 reproduces, into
# references/periodic-gap-source/ (git-ignored). Integrity chain: the downloaded
# source-manifest.json must have the SHA-256 that spin1 pins (MANIFEST_SHA in
# models/spin1/src/Spin1Certificates.jl, read from there so that the pin lives in one place);
# every other file must match the hash the manifest lists for it. Any mismatch aborts.
#
# `--paper`: the LaTeX build sources `build/<NAME>.tex` cited by the lean_seeds studies, into
# models/lean_seeds/upstream/ (git-ignored). NAME defaults to initialization, transfer,
# thermal72, trials and bootstrap (all five are read by
# models/lean_seeds/scripts/seed_design/summarize.jl). Every file listed in
# models/lean_seeds/results/spin1_kernel/paper-source-sha256.txt must match that hash; the
# SHA-256 of every fetched file is printed.

using Downloads, SHA

const REPOSITORY = "openai/math"
const COMMIT = "adc7f1241b42e322a6451854ab7e4b4c146bf78a"
const PAPER = "preprints/The-periodic-spin-one-Haldane-gap-September-24-2026"
const REPO = normpath(joinpath(@__DIR__, ".."))
const TARGET = joinpath(REPO, "references", "periodic-gap-source")
const LEAN_SEEDS = joinpath(REPO, "models", "lean_seeds")
const PAPER_NAMES = ("initialization", "transfer", "thermal72", "trials", "bootstrap")

"The manifest hash pinned by models/spin1 (hash-pinned code, so it is read, not copied)."
function pinned_manifest_sha()
    source = read(joinpath(REPO, "models", "spin1", "src", "Spin1Certificates.jl"), String)
    found = match(r"const MANIFEST_SHA = \"([0-9a-f]{64})\"", source)
    found === nothing && error("MANIFEST_SHA not found in Spin1Certificates.jl")
    return String(found[1])
end

raw_url(path) = "https://raw.githubusercontent.com/$REPOSITORY/$COMMIT/$PAPER/$path"
sha256_hex(path) = bytes2hex(sha256(read(path)))

"""
    fetch(path, expected; target = TARGET, destination = joinpath(target, path))

Download `path` (relative to the paper directory) unless a copy with the `expected` hash exists;
`expected = nothing` accepts any content. Returns `:present` or `:downloaded`.
"""
function fetch(path, expected; target = TARGET, destination = joinpath(target, path))
    expected !== nothing && isfile(destination) && sha256_hex(destination) == expected &&
        return :present
    mkpath(dirname(destination))
    Downloads.download(raw_url(path), destination)
    actual = sha256_hex(destination)
    expected === nothing || actual == expected ||
        error("hash mismatch for $path: expected $expected, got $actual")
    return :downloaded
end

"Fetch the evidence files spin1 verifies, all checked against the pinned manifest."
function fetch_evidence()
    manifest_path = "verification/computations/source-manifest.json"
    fetch(manifest_path, pinned_manifest_sha())
    manifest_text = read(joinpath(TARGET, manifest_path), String)
    # The manifest is a flat JSON object of {"name": "sha256"} groups; read the two groups spin1
    # verifies, and the hash of the top-level verifier, without a JSON dependency.
    function group(name)
        block = match(Regex("\"$name\"\\s*:\\s*\\{([^}]*)\\}"), manifest_text)
        block === nothing && error("manifest has no group $name")
        return [(m[1], m[2])
                for m in eachmatch(r"\"([^\"]+)\"\s*:\s*\"([0-9a-f]{64})\"", block[1])]
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
end

"Fetch `build/<name>.tex` for each name into models/lean_seeds/upstream/, checking pinned hashes."
function fetch_paper(names)
    pinned = Dict{String,String}()
    for line in eachline(joinpath(LEAN_SEEDS, "results", "spin1_kernel", "paper-source-sha256.txt"))
        hash, path = split(line)
        pinned[path] = hash
    end
    for name in names
        relative = "upstream/$name.tex"
        destination = joinpath(LEAN_SEEDS, relative)
        fetch("build/$name.tex", get(pinned, relative, nothing); destination)
        println(name, ".tex ", sha256_hex(destination),
            haskey(pinned, relative) ? " (matches the pinned hash)" : " (no pinned hash)")
    end
end

if isempty(ARGS)
    fetch_evidence()
elseif ARGS[1] == "--paper"
    fetch_paper(length(ARGS) > 1 ? ARGS[2:end] : PAPER_NAMES)
else
    error("usage: julia tools/fetch_upstream.jl [--paper [NAME ...]]")
end
