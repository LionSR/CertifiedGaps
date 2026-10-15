# After run_su2_grid.jl: collect the L14 sectors, gzip the spectra, rebuild scaling.json and the
# RP/environment records, then run the tests. Stops at the first failing step.
#
#     julia models/blbq_haldane/scripts/finalize_su2.jl
include(joinpath(@__DIR__, "run.jl"))

for step in (["models/blbq_haldane/scripts/collect_su2.jl", "14"],
    ["models/blbq_haldane/scripts/archive_spectra.jl"],
    ["models/blbq_haldane/scripts/analyze.jl"],
    ["models/blbq_haldane/scripts/audit.jl"],
    ["models/blbq_haldane/test/runtests.jl"])
    run(blbq_command(step))
end
