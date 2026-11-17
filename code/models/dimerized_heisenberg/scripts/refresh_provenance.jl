# Provenance-only refresh of the n₀ = 4 results (no recomputation).
#
#   julia --project=code/julia code/models/dimerized_heisenberg/scripts/refresh_provenance.jl REV
#   julia --project=code/julia code/models/dimerized_heisenberg/scripts/refresh_provenance.jl \
#       --verify REV
#
# Adding model-local sources (src/RitzSeeds.jl, the coupled seeds in certify.jl and check.jl)
# changes the source hashes recorded by results/rp_audit.json, results/summary.json and the
# n₀ = 4 certificates results/certificate_delta_*_beta_*.json. Regenerating them with certify.jl
# would not reproduce their numbers: the `verified_eigen` radii depend on the LAPACK eigenvectors,
# which vary between runs and platforms (about 4·10⁻¹¹ in εₛ), so every enclosure, defect and
# bound would move in its last digits. This script instead reads each file as it is at the git
# revision REV (e.g. the upstream `main` the branch starts from) and rewrites only provenance:
#   - `metadata.source_sha256` and `metadata.core_sha256`: the current `Results.model_provenance`
#     with the `dirs`/`files` of `certify_point` and of `metadata` in certify.jl (the defaults of
#     `Results.model_metadata`), so `Results.verify_model_metadata` in check.jl accepts them;
#   - `metadata.provenance_refresh`: this script, REV and the time of the refresh;
#   - in summary.json, the `sha256` of every certificate row (the certificates' bytes changed).
# The original run record (created_utc, versions, threads, config and its seconds) is kept, since
# it describes the run that produced the numbers. Afterwards (and with --verify, alone) every
# JSON value of each file is compared with REV; any difference outside the keys above throws.
using CertifiedGaps.Results: read_result_bytes, write_result, model_provenance, sha256_file
using JSON3, Dates

const MODEL_DIR = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(MODEL_DIR, "results")
const REPO = normpath(joinpath(MODEL_DIR, "..", ".."))
const DATA = "data/models/dimerized_heisenberg"         # path of results/ in git
const SCRIPT = "models/dimerized_heisenberg/scripts/refresh_provenance.jl"

"The n₀ = 4 certificates (every certificate of the model except the coupled `_n0_6_` ones)."
certificates() = sort!(filter(f -> occursin(r"^certificate_delta_.*\.json$", f) &&
                                   !occursin("_n0_", f), readdir(RESULTS_DIR)))

"""
A mutable JSON object that keeps its key order, so a rewritten file keeps the layout of REV and
its `git diff` shows only the provenance lines (a `Dict` would reorder the keys).
"""
struct Ordered <: AbstractDict{String,Any}
    order::Vector{String}
    values::Dict{String,Any}
end
Base.length(d::Ordered) = length(d.order)
Base.haskey(d::Ordered, k) = haskey(d.values, k)
Base.getindex(d::Ordered, k) = d.values[k]
Base.get(d::Ordered, k, default) = get(d.values, k, default)
Base.iterate(d::Ordered, i = 1) =
    i > length(d.order) ? nothing : (d.order[i] => d.values[d.order[i]], i + 1)
function Base.setindex!(d::Ordered, v, k)
    haskey(d.values, k) || push!(d.order, k)
    d.values[k] = v
    return d
end

ordered(x::JSON3.Object) = Ordered([String(k) for k in keys(x)],
    Dict{String,Any}(String(k) => ordered(v) for (k, v) in x))
ordered(x::JSON3.Array) = Any[ordered(v) for v in x]
ordered(x) = x

"The text of the file `name` of results/ at the git revision `rev`."
revision_text(rev, name) = read(Cmd(`git show $rev:$DATA/$name`; dir = REPO), String)

"The file `name` of results/ at `rev`, parsed with nested `Dict{String,Any}` (for comparison)."
at_revision(rev, name) = JSON3.read(revision_text(rev, name), Dict{String,Any})

"The file `name` of results/ at `rev`, key order kept (for rewriting)."
ordered_at_revision(rev, name) = ordered(JSON3.read(revision_text(rev, name)))

current(name) = JSON3.read(String(read_result_bytes(joinpath(RESULTS_DIR, name))),
    Dict{String,Any})

# The keys whose values may differ from REV: provenance only.
function provenance_key(path)
    length(path) >= 2 && path[1] == "metadata" &&
        path[2] in ("source_sha256", "core_sha256", "provenance_refresh") && return true
    # summary.json: parameters/<i>/certificates/<j>/sha256
    return length(path) == 5 && path[1] == "parameters" && path[3] == "certificates" &&
           path[5] == "sha256"
end

"Every leaf path where `a` and `b` differ (missing keys included)."
function differences!(out, a, b, path = String[])
    if a isa AbstractDict && b isa AbstractDict
        for k in union(keys(a), keys(b))
            p = [path; String(k)]
            haskey(a, k) && haskey(b, k) ? differences!(out, a[k], b[k], p) : push!(out, p)
        end
    elseif a isa AbstractVector && b isa AbstractVector
        length(a) == length(b) || return push!(out, path)
        for (i, (x, y)) in enumerate(zip(a, b))
            differences!(out, x, y, [path; string(i)])
        end
    else
        # Same type and value: an integer 4 and a float 4.0 or a string "4" count as different.
        (typeof(a) == typeof(b) && isequal(a, b)) || push!(out, path)
    end
    return out
end

leaves(x) = x isa AbstractDict ? sum(leaves, values(x); init = 0) :
    x isa AbstractVector ? sum(leaves, x; init = 0) : 1

"Compare every file with `rev`; print the differing keys; throw unless all are provenance."
function verify(rev)
    for name in [certificates(); "summary.json"; "rp_audit.json"]
        old, new = at_revision(rev, name), current(name)
        diff = differences!(Vector{String}[], old, new)
        bad = filter(!provenance_key, diff)
        changed = sort!(unique(join(p[1:min(end, p[1] == "metadata" ? 2 : 5)], ".")
                               for p in diff))
        println(rpad(name, 40), "values = ", leaves(old), ", equal to $rev except ",
            length(diff), " provenance leaves (", join(changed, ", "), ")")
        isempty(bad) || error("$name: non-provenance values differ from $rev: " *
                              join(join.(bad, "."), ", "))
    end
    println("every non-provenance value equals $rev")
end

function refresh!(data, rev, sha)
    prov = model_provenance(MODEL_DIR; dirs = ("src",), files = ("certify.jl", "check.jl"))
    meta = data["metadata"]
    meta["source_sha256"], meta["core_sha256"] = prov.source_sha256, prov.core_sha256
    meta["provenance_refresh"] = Dict{String,Any}("script" => SCRIPT, "from_revision" => sha,
        "utc" => string(now(UTC)), "scope" => "only source_sha256, core_sha256 (and in " *
        "summary.json the certificate sha256) were rewritten; every other value, including " *
        "this run record, is that of the original run at from_revision")
    return data
end

function refresh(rev)
    sha = strip(read(Cmd(`git rev-parse $rev`; dir = REPO), String))
    for name in certificates()
        data = refresh!(ordered_at_revision(rev, name), rev, sha)
        write_result(joinpath(RESULTS_DIR, name), data)
    end
    summary = refresh!(ordered_at_revision(rev, "summary.json"), rev, sha)
    for row in summary["parameters"], c in row["certificates"]
        c["sha256"] = sha256_file(joinpath(MODEL_DIR, c["file"]))
    end
    write_result(joinpath(RESULTS_DIR, "summary.json"), summary)
    write_result(joinpath(RESULTS_DIR, "rp_audit.json"),
        refresh!(ordered_at_revision(rev, "rp_audit.json"), rev, sha))
    println("refreshed the provenance of ", length(certificates()) + 2, " files from $rev ($sha)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = filter(!=("--verify"), ARGS)
    length(args) == 1 || error("usage: refresh_provenance.jl [--verify] REV")
    "--verify" in ARGS || refresh(only(args))
    verify(only(args))
end
