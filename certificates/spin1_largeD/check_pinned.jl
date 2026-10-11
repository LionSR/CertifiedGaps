# Replay with the exact core used by the computation, including after its upstream rename.
# This launcher uses only Julia Base and never changes the checkout or certificate files.
const certificate_root = @__DIR__
const repository_root = normpath(joinpath(certificate_root,"../.."))
const core_revision = "ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b"
@assert VERSION==v"1.12.4"
check = isempty(ARGS) ? "runtests.jl" : only(ARGS)
check in ("runtests.jl","small.jl","symmetry.jl") || error("unknown checker: $check")
mktempdir() do workspace
    # git archive reads immutable committed objects; no checkout, reset or network access.
    run(pipeline(`git -C $repository_root archive $core_revision julia`,
        `tar -x -C $workspace`))
    target=joinpath(workspace,"certificates","spin1_largeD")
    mkpath(dirname(target)); cp(certificate_root,target)
    legacy=joinpath(certificate_root,"provenance","legacy")
    if isdir(legacy)
        for (folder,_,files) in walkdir(legacy), file in files
            destination=joinpath(target,relpath(folder,legacy),file)
            mkpath(dirname(destination)); cp(joinpath(folder,file),destination;force=true)
        end
    end
    project=joinpath(workspace,"julia")
    script=joinpath(target,"test",check)
    println("Replaying $check against core $core_revision")
    withenv("OPENBLAS_NUM_THREADS"=>"1") do
        run(`$(Base.julia_cmd()) --project=$project $script`)
    end
end
