using Spin1Certificates
Spin1Certificates.ITensorsChecks.run_itensors_checks(
    get(ARGS, 1, joinpath(@__DIR__, "../results/itensors.json")),
)
