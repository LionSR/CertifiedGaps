# End-to-end validation: Lake build, axiom audit of every audited declaration, Julia tests,
# generator hash and toolchain/Mathlib pins. Needs `lake` on PATH (or LAKE_EXE). Logs go to
# the ignored `.cache/`; nothing is written to results/.
using JSON, SHA
const HERE = normpath(joinpath(@__DIR__, ".."))
const LEAN = joinpath(HERE, "lean")
const REPO_LEAN = normpath(joinpath(HERE, "..", "..", "lean"))
const ALLOWED_AXIOMS = ["propext", "Classical.choice", "Quot.sound"]
const AUDITED = 17

"Run `cmd` in `dir`, logging to `.cache/<name>.log`; return the log text."
function timed(name, cmd, dir)
    mkpath(joinpath(HERE, ".cache"))
    log = joinpath(HERE, ".cache", "$name.log")
    seconds = @elapsed open(log, "w") do io
        run(pipeline(Cmd(cmd; dir=dir); stdout=io, stderr=io))
    end
    println(rpad(name, 12), round(seconds; digits=1), " s")
    return read(log, String)
end

function main()
    lake = get(ENV, "LAKE_EXE", "lake")
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    timed("build", `$lake build`, LEAN)
    audit = timed("axioms", `$lake env lean Audit.lean`, LEAN)
    timed("tests", `$julia --project=$HERE $(joinpath(HERE, "test", "runtests.jl"))`, HERE)
    # Reject sorry, native reduction (Lean.ofReduceBool) and anything beyond the standard three.
    @assert !occursin("sorryAx", audit) && !occursin("ofReduceBool", audit)
    lines = filter(l -> occursin("LeanSeedSmall.", l), split(audit, '\n'))
    @assert length(lines) == AUDITED
    for l in lines
        m = match(r"\[(.*)\]", l)
        m === nothing && continue # "does not depend on any axioms"
        @assert all(x -> strip(x) in ALLOWED_AXIOMS, split(m[1], ','))
    end
    cert = JSON.parsefile(joinpath(HERE, "results", "certificate.json"))
    @assert cert["source_hashes"]["certify.jl"] ==
            bytes2hex(sha256(read(joinpath(HERE, "certify.jl"))))
    # Same toolchain and Mathlib commit as the repository's Lean package.
    @assert read(joinpath(LEAN, "lean-toolchain"), String) ==
            read(joinpath(REPO_LEAN, "lean-toolchain"), String)
    rev(p) = only(x["rev"] for x in JSON.parsefile(p)["packages"] if x["name"] == "mathlib")
    @assert rev(joinpath(LEAN, "lake-manifest.json")) ==
            rev(joinpath(REPO_LEAN, "lake-manifest.json"))
    println("Validation passed: $AUDITED declarations use only $(join(ALLOWED_AXIOMS, ", ")).")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
