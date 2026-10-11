using Spin1Certificates
source = parse(Int, get(ARGS, 1, "1"))
source in (1, 2) || error("source must be 1 or 2")
output = get(ARGS, 2, joinpath(@__DIR__, "../results/source$source.json"))
Spin1Certificates.run_thermal(source, output)
