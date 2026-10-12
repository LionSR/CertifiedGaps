# Historical provenance is immutable. Replay against its exact audited core snapshot.
const PIN = "ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b"
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "../.."))
const REPLAY = joinpath(HERE, ".replay", PIN)
mkpath(REPLAY)
run(pipeline(`git -C $ROOT archive $PIN julia outputs/thermal-mpo-gap-note.tex outputs/thermal-mpo-gap-scope.tex`,
             `tar -x -C $REPLAY`))
mkpath(joinpath(REPLAY, "certificates"))
link = joinpath(REPLAY, "certificates", "dimerized_heisenberg")
ispath(link) || symlink(HERE, link)
if !isempty(ARGS) && first(ARGS) == "--setup"
    run(`$(Base.julia_cmd()) --project=$(joinpath(REPLAY,"julia")) -e 'using Pkg; Pkg.instantiate()'`)
else
    script = isempty(ARGS) ? "test/runtests.jl" : first(ARGS)
    target = normpath(joinpath(link, script))
    startswith(target, link * "/") || error("script must remain inside this certificate")
    rest = isempty(ARGS) ? String[] : ARGS[2:end]
    run(Cmd(`$(Base.julia_cmd()) --project=$(joinpath(REPLAY,"julia")) $target $rest`; dir=REPLAY))
end
