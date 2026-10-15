# Checker for the two-leg TFIM certificate.
#
#     julia --project=models/two_leg_tfim models/two_leg_tfim/check.jl          # fast, ~35 s
#     julia --project=models/two_leg_tfim models/two_leg_tfim/check.jl --full   # + ~11 min
#
# The fast check (test/runtests.jl) recomputes log Z, the rational seeds and the bound from the
# stored eigenvalues and radii. `--full` first reruns certify.jl at g = 3, which rebuilds every
# block, eigensystem and radius from the model and writes a fresh certificate that the fast
# check then verifies as well. scripts/replay_audit.jl does the same and compares bit for bit.
if "--full" in ARGS
    filter!(!=("--full"), ARGS)
    isempty(ARGS) && push!(ARGS, "3//1")
    include("certify.jl")
end
include("test/runtests.jl")
