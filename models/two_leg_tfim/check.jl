# Checker for the two-leg TFIM certificate.
#
#     julia --project=models/two_leg_tfim models/two_leg_tfim/check.jl          # fast, ~35 s
#     julia --project=models/two_leg_tfim models/two_leg_tfim/check.jl --full   # + ~11 min
#
# The fast check (test/runtests.jl) recomputes log Z, the rational seeds and the bound from the
# stored eigenvalues and radii. `--full` first reruns certify.jl at g = 3 into a temporary file
# (never into results/), which rebuilds every block, eigensystem and radius from the model; the
# fast checks then verify that fresh certificate as well, with the same seed (g, n0, beta) and a
# bound at least the stored one. scripts/replay_audit.jl does the same and also reports
# sector-by-sector agreement.
using Test

full = "--full" in ARGS
empty!(ARGS)
if full
    fresh = tempname() * ".json"
    certify = joinpath(@__DIR__, "certify.jl")
    project = dirname(Base.active_project())
    run(`$(Base.julia_cmd()) --startup-file=no --project=$project $certify 3//1 $fresh`)
    push!(ARGS, fresh)
end
# Defines the Ladder and Provenance modules used below.
include(joinpath(@__DIR__, "test", "runtests.jl"))
if full
    # Same seed, and a fresh bound at least the stored one (as replay_audit.jl's
    # `original_bound_supported`). The radii depend on the LAPACK/BLAS run, so the enclosures
    # need not be bit-identical: a run on a loaded machine gave 0.636473... vs 0.636386....
    @testset "Fresh certificate supports the stored one" begin
        new = read_result(fresh)
        old = read_result(joinpath(@__DIR__, ORIGINAL_CERTIFICATE))
        for key in ("g", "n0", "beta")
            @test new[key] == old[key]
        end
        @test new["bound"] >= old["bound"]
        identical = all(new[k] == old[k] for k in ("epsilon_s", "epsilon_t", "bound", "logZ"))
        println("fresh bound $(new["bound"]) vs stored $(old["bound"]); ",
            identical ? "seeds and log Z identical" : "enclosures differ (radii)")
    end
    rm(fresh)
end
