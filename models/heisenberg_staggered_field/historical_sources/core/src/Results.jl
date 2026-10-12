"""
    ThermalMPOGap.Results

Writing results and recording where they came from: atomic JSON output, SHA-256 fingerprints
of source trees and environments, run metadata, and JSON parsing that never routes a large
integer through `Float64`.
"""
module Results

using Dates, JSON3, LinearAlgebra, SHA

export write_result, sha256_file, file_hashes, source_hashes, run_metadata,
    parse_exact_json, read_exact_json

"""
    write_result(path, data) -> path

Write `data` as pretty JSON with a trailing newline, atomically: a crash leaves either the old
file or the new one, never a partial file.
"""
function write_result(path, data)
    mkpath(dirname(abspath(path)))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        JSON3.pretty(io, data)
        println(io)
    end
    mv(temporary, path; force = true)
    return path
end

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

"`parse_exact_json` of a file."
read_exact_json(path) = parse_exact_json(read(path, String))

end
