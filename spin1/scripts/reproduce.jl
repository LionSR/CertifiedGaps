using Spin1Certificates
const SC = Spin1Certificates
mode = get(ARGS, 1, "all")
directory = abspath(get(ARGS, 2, joinpath(@__DIR__, "../results")))
if mode == "all"
    # An incomplete run cannot retain an earlier top-level success record.
    SC.write_result(
        joinpath(directory, "verification.json"),
        Dict("status" => "running", "execution_mode" => "fresh-recomputation"),
    )
    SC.run_arithmetic(joinpath(directory, "arithmetic.json"))
    SC.run_trials(joinpath(directory, "trials.json"))
    SC.run_thermal(1, joinpath(directory, "source1.json"))
    SC.run_thermal(2, joinpath(directory, "source2.json"))
    SC.run_scalars(joinpath(directory, "scalars.json"))
    SC.ITensorsChecks.run_itensors_checks(joinpath(directory, "itensors.json"))
    SC.audit_results(directory; execution_mode = "fresh-recomputation")
elseif mode == "audit"
    SC.audit_results(directory)
else
    error(
        "Usage: julia --project=spin1 --threads=4 spin1/scripts/reproduce.jl [all|audit] [output-directory]",
    )
end
