# Twenty-spin numerical diagnostics (NOT certified; see src/Diagnostics.jl).
#
#   julia scripts/scan20.jl p q      one delta = p/q, about 1.5-2 min
#                                    -> results/diagnostic20_p_q.json (gitignored)
#   julia scripts/scan20.jl summary  condense the six runs -> results/diagnostic20_summary.json
#
# The complete spectra of the rings up to 16 spins are the core `float_ring_spectrum`.
include(joinpath(@__DIR__, "..", "src", "DimerizedHeisenberg.jl"))
include(joinpath(@__DIR__, "..", "src", "Diagnostics.jl"))
using LinearAlgebra, .DimerizedHeisenberg
BLAS.set_num_threads(1)

if length(ARGS) == 2
    Diagnostics.run_diagnostic(parse(BigInt, ARGS[1]) // parse(BigInt, ARGS[2]))
elseif ARGS == ["summary"]
    Diagnostics.summarize_diagnostics(DELTAS)
else
    error("usage: scan20.jl p q | scan20.jl summary")
end
