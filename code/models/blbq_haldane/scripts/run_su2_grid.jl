# All 120 (S, k) highest-weight sectors of the 14-site ring, resumable (finished sectors are
# skipped), each in its own 4-thread process with a 45-minute timeout. Stops at the first
# failure or timeout.
#
#     julia models/blbq_haldane/scripts/run_su2_grid.jl
include(joinpath(@__DIR__, "run.jl"))

const TIMEOUT = 2700.0 # seconds per sector process
const MODEL = joinpath(@__DIR__, "..")

for spin in 0:14, momentum in 0:7
    output = joinpath(MODEL, "results", "su2_L14_S$(spin)_k$(momentum).json")
    isfile(output) && continue
    process = run(blbq_command(["-t", "4", "models/blbq_haldane/scripts/su2_sector.jl",
        "14", string(spin), string(momentum)]); wait = false)
    if timedwait(() -> process_exited(process), TIMEOUT) == :timed_out
        kill(process)
        error("sector S=$spin k=$momentum exceeded $(TIMEOUT) s")
    end
    success(process) || error("sector S=$spin k=$momentum failed")
end
