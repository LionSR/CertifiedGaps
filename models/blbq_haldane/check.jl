# Checker for the stored survey: there is no gap certificate to verify, so this runs the test
# suite, which re-proves the exact RP identities, confirms that every stored candidate seed is
# rejected by the recurrence, checks the hashes of the stored references and recomputes the
# small-ring rows (n₀ = 2, 3) of results/scaling.json.
#
#     julia --project=models/blbq_haldane models/blbq_haldane/check.jl          # fast
#     julia --project=models/blbq_haldane models/blbq_haldane/check.jl --full   # + n₀ = 4
#
# `--full` also recomputes the n₀ = 4 rows from complete dense spectra of the 4- and 8-site
# rings. Larger rows need the regeneration pipeline in scripts/ (hours; see README).
include(joinpath(@__DIR__, "test", "runtests.jl"))
if "--full" in ARGS
    @testset "n₀ = 4 rows of the summary tables reproduced from scratch" begin
        check_small_rows(4:4, (4, 8))
    end
end
