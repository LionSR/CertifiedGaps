# Checker for the stored survey: there is no gap certificate to verify, so this runs the test
# suite, which re-proves the exact RP identities, confirms that every stored candidate seed is
# rejected by the recurrence, checks the hashes of the stored references and recomputes the
# small-ring rows of results/scaling.json.
include(joinpath(@__DIR__, "test", "runtests.jl"))
