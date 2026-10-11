using Spin1Certificates
Spin1Certificates.run_trials(get(ARGS, 1, joinpath(@__DIR__, "../results/trials.json")))
