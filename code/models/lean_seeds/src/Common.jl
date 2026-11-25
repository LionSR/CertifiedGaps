"""
    Common

Paths and helpers shared by the six studies: directory layout, source fingerprints, the Sᶻ
product basis of a spin-1 chain, and the Lake/axiom-audit driver for `lean/SeedKernel`.
"""
module Common

using CertifiedGaps.Results: sha256_file

"This model directory (`models/lean_seeds`)."
const MODEL = dirname(@__DIR__)
"Repository root."
const REPO = normpath(joinpath(MODEL, "..", ".."))
"The repository's Lake package, which hosts the `SeedKernel` library."
const LEAN = joinpath(REPO, "lean")
"Stored certificate data of the spin-1 Heisenberg model (read only)."
const SPIN1_RESULTS = joinpath(REPO, "models", "spin1", "results")
"Paper sources fetched by `scripts/fetch_paper.jl` (gitignored)."
const UPSTREAM = joinpath(MODEL, "upstream")
"The six studies, in the order the checker runs them."
const STUDIES = ("seed_design", "spectrum", "anchor", "spin1_kernel", "ising", "ising_open")

"`results/<study>` of this model."
results_dir(study) = joinpath(MODEL, "results", study)

"""
    study_hashes(study; extra = ())

SHA-256 of the sources a study's results depend on, keyed by path relative to MODEL: the shared
drivers (`certify.jl`, `check.jl`, `src/Common.jl`, `src/LeanSeeds.jl`), `src/<study>.jl`,
every `.jl` under `scripts/<study>`, `test/<study>.jl`, and `extra` (relative paths, which may
leave the directory, e.g. `../../julia/src/Results.jl`).
"""
function study_hashes(study; extra=())
    files = String["certify.jl", "check.jl", "src/Common.jl", "src/LeanSeeds.jl",
        "src/$study.jl", "test/$study.jl"]
    dir = joinpath(MODEL, "scripts", study)
    if isdir(dir)
        for (folder, _, names) in walkdir(dir), name in names
            endswith(name, ".jl") && push!(files, relpath(joinpath(folder, name), MODEL))
        end
    end
    append!(files, extra)
    return Dict(f => sha256_file(joinpath(MODEL, f)) for f in sort(files))
end

"""
    sz_basis(n, m)

Sorted Sᶻ = m product basis of an n-site spin-1 chain, encoded as Σⱼ (sⱼ + 1) 3ʲ.
"""
function sz_basis(n, m)
    states = Int[]
    function visit(site, remaining, word)
        if site == n
            remaining == 0 && push!(states, word)
            return
        end
        for s in -1:1
            abs(remaining - s) <= n - site - 1 || continue
            visit(site + 1, remaining - s, word + (s + 1) * 3^site)
        end
    end
    visit(0, m, 0)
    return sort!(states)
end

"""
    export_all(mod)

Export every name defined in `mod`, so that `using .LeanSeeds.<Study>` gives the scripts the
flat namespace they were written for. Only the top-level drivers mix studies, and they qualify
names.
"""
function export_all(mod::Module)
    skip = (:eval, :include, nameof(mod))
    for name in names(mod; all=true)
        s = string(name)
        (startswith(s, '#') || startswith(s, '@') || name in skip) && continue
        Core.eval(mod, Expr(:export, name))
    end
end

# ---------------------------------------------------------------------------------------------
# Lake driver
# ---------------------------------------------------------------------------------------------

"The `lake` executable (override with LAKE_BIN)."
const LAKE = get(ENV, "LAKE_BIN", "lake")
"The only axioms a SeedKernel declaration may use."
const STANDARD_AXIOMS = Set(["propext", "Classical.choice", "Quot.sound"])

"""
    lean_output(cmd; dir = LEAN)

Run `cmd` in `dir` and return its combined output; throw (with the output) if it fails or
reports an error.
"""
function lean_output(cmd; dir=LEAN)
    buffer = IOBuffer()
    process = run(pipeline(ignorestatus(Cmd(cmd; dir)); stdout=buffer, stderr=buffer))
    text = String(take!(buffer))
    process.exitcode == 0 || error("$(cmd) failed:\n$text")
    occursin(r"(?m)^error|\berror:", text) && error("$(cmd) reported an error:\n$text")
    return text
end

"""
    check_axioms(text; expected)

Every `#print axioms` report in `text` must name only the standard axioms (no `sorryAx`, no
`Lean.ofReduceBool`, i.e. no `sorry` and no `native_decide`); every name in `expected` must have
a report. Returns the reported declaration names.
"""
function check_axioms(text; expected)
    occursin("sorryAx", text) && error("an audited declaration depends on sorryAx")
    occursin("ofReduceBool", text) && error("an audited declaration uses native_decide")
    seen = String[]
    for m in eachmatch(r"'([^']+)' depends on axioms: \[([^\]]*)\]", text)
        axioms = Set(strip.(split(m[2], ',')))
        issubset(axioms, STANDARD_AXIOMS) || error("$(m[1]) uses nonstandard axioms $(axioms)")
        push!(seen, m[1])
    end
    for m in eachmatch(r"'([^']+)' does not depend on any axioms", text)
        push!(seen, m[1])
    end
    missing_names = setdiff(expected, seen)
    isempty(missing_names) || error("no axiom report for $(missing_names)")
    return seen
end

"Declarations whose axioms `SeedKernel/Audit.lean` must report (all of its lines)."
audited_declarations() = [m[1] for m in eachmatch(r"(?m)^#print axioms (\S+)",
    read(joinpath(LEAN, "SeedKernel", "Audit.lean"), String))]

"""
    build_seed_kernel(; fresh = false) -> (seconds, audited names)

`lake build SeedKernel`, then the axiom audit (`lake env lean SeedKernel/Audit.lean`). With
`fresh`, the SeedKernel build products are deleted first so that the kernel re-checks every
`decide +kernel` proof (CertifiedGaps and Mathlib are reused).
"""
function build_seed_kernel(; fresh=false)
    if fresh
        build = joinpath(LEAN, ".lake", "build")
        for dir in ("lib/lean", "ir"), suffix in ("", ".olean", ".ilean", ".c", ".trace",
                ".hash", ".olean.hash", ".ilean.hash", ".c.hash", ".setup.json")
            path = joinpath(build, dir, "SeedKernel" * suffix)
            ispath(path) && rm(path; recursive=true)
        end
    end
    seconds = @elapsed lean_output(`$LAKE build SeedKernel`)
    audit = lean_output(`$LAKE env lean SeedKernel/Audit.lean`)
    expected = audited_declarations()
    names = check_axioms(audit; expected)
    length(names) == length(expected) || error("unexpected axiom reports: $(names)")
    return seconds, names
end

end
