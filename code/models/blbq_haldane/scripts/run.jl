# Launcher for the regeneration drivers: runs `julia --project=julia ARGS...` from the
# repository root with one BLAS thread. If a private Julia 1.12.4 runtime exists under
# models/blbq_haldane/runtime/ (as in the original 4-CPU run), it and its depot are used.
#
#     julia models/blbq_haldane/scripts/run.jl [julia options] SCRIPT [script args]
#     e.g. julia models/blbq_haldane/scripts/run.jl -t 4 models/blbq_haldane/scripts/scan.jl 2 3

"""
    blbq_command(args) -> Cmd

The `julia --project=julia args` command used by every regeneration driver, run from the
repository root with `OPENBLAS_NUM_THREADS=1` (and the private runtime when it is present).
"""
function blbq_command(args)
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    runtime = joinpath(repo, "models", "blbq_haldane", "runtime")
    private = joinpath(runtime, "julia-1.12.4", "bin", "julia")
    env = Dict("OPENBLAS_NUM_THREADS" => "1")
    julia = Base.julia_cmd()
    if Sys.isexecutable(private)
        env["JULIA_DEPOT_PATH"] = joinpath(runtime, "depot")
        julia = `$private`
    end
    return addenv(Cmd(`$julia --project=julia $args`; dir = repo), env)
end

abspath(PROGRAM_FILE) == (@__FILE__) && run(blbq_command(ARGS))
