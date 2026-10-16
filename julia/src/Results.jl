"""
    CertifiedGaps.Results

Writing and reading results and recording where they came from: atomic JSON output (pretty or
compact, optionally gzip-compressed), SHA-256 fingerprints of source trees and environments,
model-local provenance and its verification, run metadata, and JSON parsing that never routes a
large integer through `Float64`.

Compression uses the system `gzip` (no extra package dependency, so no model manifest changes):
`gzip -n` stores neither file name nor timestamp, so equal JSON gives equal compressed bytes.
"""
module Results

using Dates, JSON3, LinearAlgebra, SHA

export write_result, read_result, read_result_bytes, sha256_file, file_hashes, source_hashes,
    model_provenance, verify_source_hashes, run_metadata, parse_exact_json, read_exact_json

"""
    write_result(path, data; compact = false, gzip = endswith(path, ".gz"), level = 6) -> path

Write `data` as JSON with a trailing newline, atomically: a crash leaves either the old file or
the new one, never a partial file. The JSON is pretty-printed (`JSON3.pretty`) by default and
single-line (`JSON3.write`) with `compact = true`. With `gzip = true` (the default exactly when
`path` ends in `.gz`, which is then required) the JSON is compressed by `gzip -n -<level>`;
`read_result` restores it.
"""
function write_result(path, data; compact::Bool = false, gzip::Bool = endswith(path, ".gz"),
    level::Integer = 6)
    gzip == endswith(path, ".gz") ||
        throw(ArgumentError("gzip = $gzip needs a path that $(gzip ? "ends" : "does not end") in .gz"))
    1 <= level <= 9 || throw(ArgumentError("gzip level must be in 1:9"))
    mkpath(dirname(abspath(path)))
    temporary = path * ".tmp"
    plain = gzip ? temporary * ".json" : temporary
    open(plain, "w") do io
        compact ? JSON3.write(io, data) : JSON3.pretty(io, data)
        println(io)
    end
    if gzip
        run(pipeline(`gzip -n -$level -c -- $plain`; stdout = temporary))
        rm(plain)
    end
    mv(temporary, path; force = true)
    return path
end

"""
    read_result_bytes(path) -> Vector{UInt8}

The uncompressed bytes of a result file: `.gz` files are decompressed with the system `gzip`,
others are read as they are. Hash these bytes to fingerprint a compressed result by content.
"""
read_result_bytes(path) = endswith(path, ".gz") ? read(`gzip -dc -- $path`) : read(path)

"""
    read_result(path; exact = true)
    read_result(path, T)

Read a `.json` or `.json.gz` result. With `exact = true` (the default) it is parsed by
`parse_exact_json`, so every bare integer comes back as a string; with `exact = false` by
`JSON3.read` (a `JSON3.Object`). `read_result(path, T)` is `JSON3.read(text, T)`, e.g.
`T = Dict{String,Any}`.
"""
function read_result(path; exact::Bool = true)
    text = String(read_result_bytes(path))
    return exact ? parse_exact_json(text) : JSON3.read(text)
end
read_result(path, ::Type{T}) where {T} = JSON3.read(String(read_result_bytes(path)), T)

"Hex SHA-256 of a file's bytes."
sha256_file(path) = bytes2hex(sha256(read(path)))

"""
    file_hashes(root, names; prefix = "") -> Dict

SHA-256 of each `root/name`, keyed by `prefix * name`.
"""
file_hashes(root, names; prefix = "") =
    Dict(prefix * name => sha256_file(joinpath(root, name)) for name in names)

"""
    source_hashes(root; dirs = ("src", "scripts", "test"), files = (), skip_missing = true)

SHA-256 of every `.jl` file under `root/dir` for each `dir`, plus the listed `files`, keyed by
path relative to `root`. Missing directories are skipped unless `skip_missing = false`.
"""
function source_hashes(root; dirs = ("src", "scripts", "test"), files = (),
    skip_missing = true)
    paths = String[joinpath(root, name) for name in files]
    for dir in dirs
        top = joinpath(root, dir)
        isdir(top) || (skip_missing ? continue : error("missing source directory $top"))
        for (folder, _, names) in walkdir(top), name in names
            endswith(name, ".jl") && push!(paths, joinpath(folder, name))
        end
    end
    return Dict(relpath(path, root) => sha256_file(path) for path in sort(paths))
end

"""
    model_provenance(model_dir; dirs = ("src",), files = ("Project.toml", "Manifest.toml"),
        repo = normpath(joinpath(model_dir, "..", ".."))) -> (; source_sha256, core_sha256)

The repository's source-hash policy in one call. `source_sha256` holds the model-local sources
(every `.jl` file under `model_dir/dir` for `dir` in `dirs`, plus `files`, all relative to
`model_dir`); a checker enforces these, e.g. with `verify_source_hashes`. `core_sha256` holds the
shared core (`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`) as provenance only. Both
are keyed by path relative to the repository root `repo` (`models/<slug>/src/...`,
`julia/src/...`).
"""
function model_provenance(model_dir; dirs = ("src",), files = ("Project.toml", "Manifest.toml"),
    repo = normpath(joinpath(model_dir, "..", "..")))
    rel = relpath(normpath(model_dir), repo)
    source_sha256 = source_hashes(repo; dirs = Tuple(joinpath(rel, d) for d in dirs),
        files = Tuple(joinpath(rel, f) for f in files), skip_missing = false)
    core_sha256 = source_hashes(repo; dirs = ("julia/src",),
        files = ("julia/Project.toml", "julia/Manifest.toml"), skip_missing = false)
    return (; source_sha256, core_sha256)
end

"""
    verify_source_hashes(recorded, root; scope = "") -> Int

Compare recorded SHA-256 digests (`path => digest`, paths relative to `root`) with the files in
the checkout and return the number compared. Only paths inside `root` that start with `scope`
(after `normpath`) are compared; the others, e.g. the shared core recorded as provenance, are
skipped. Throws if any compared file is missing or differs, or if no path is in scope.
"""
function verify_source_hashes(recorded, root; scope::AbstractString = "")
    checked, bad = 0, String[]
    for (path, digest) in recorded
        p = normpath(String(path))
        (startswith(p, "..") || isabspath(p) || !startswith(p, scope)) && continue
        file = joinpath(root, p)
        checked += 1
        (isfile(file) && sha256_file(file) == String(digest)) || push!(bad, p)
    end
    isempty(bad) || error("source hash mismatch (or missing file): $(join(bad, ", "))")
    checked > 0 || error("no recorded source hash in scope \"$scope\" under $root")
    return checked
end

"""
    run_metadata(config; source_sha256, packages = (), status, certified = false)

The provenance record stored with a numerical result: creation time, Julia version, the
version of each package module in `packages`, thread counts, the run configuration and the
source fingerprint.
"""
function run_metadata(config; source_sha256, packages = (), status, certified = false)
    record = Dict{String,Any}(
        "created_utc" => string(now(UTC)),
        "julia" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "config" => config,
        "source_sha256" => source_sha256,
        "certified" => certified,
        "status" => status,
    )
    for package in packages
        record[string(nameof(package))] = string(pkgversion(package))
    end
    return record
end

"Parse JSON with every bare integer token read as a string, so big integers stay exact."
function parse_exact_json(raw::AbstractString)
    tokens = r"\"(?:\\.|[^\"\\])*\"|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"
    exact =
        replace(raw, tokens => s -> startswith(s, '"') || occursin(r"[.eE]", s) ? s : '"' * s * '"')
    return JSON3.read(exact, Any)
end

"`parse_exact_json` of a file (plain JSON; `read_result` also handles `.json.gz`)."
read_exact_json(path) = parse_exact_json(read(path, String))

end
