# Twenty-spin numerical diagnostics (NOT certified; see src/Diagnostics.jl).
#
#   julia scripts/scan20.jl p q      one delta = p/q, about 1.5-2 min
#                                    -> results/diagnostic20_p_q.json (gitignored)
#   julia scripts/scan20.jl summary  condense the six runs -> results/diagnostic20_summary.json
#
# Run after certify.jl: each run reads the 8/16-spin ground energies of its seed certificate.
include(joinpath(@__DIR__, "..", "src", "Diagnostics.jl"))
using LinearAlgebra
BLAS.set_num_threads(1)

const GRID = [big(1)//2, big(3)//10, big(1)//4, big(1)//5, big(3)//20, big(1)//10]

if length(ARGS) == 2
    Diagnostics.run_diagnostic(parse(BigInt, ARGS[1]) // parse(BigInt, ARGS[2]))
elseif ARGS == ["summary"]
    Diagnostics.summarize_diagnostics(GRID)
else
    error("usage: scan20.jl p q | scan20.jl summary")
end
