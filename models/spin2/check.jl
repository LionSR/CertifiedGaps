# Checker for the three tracked spin-2 finite-energy trial certificates in results/ (finite
# variational energies only; no gap is certified in this directory).
#
#     julia --project=models/spin2 models/spin2/check.jl          # fast
#     julia --project=models/spin2 models/spin2/check.jl --full   # + exact replays, ~10 min
#
# Fast: the test files that read the tracked certificates. `test/trials.jl` rebuilds the D = 24
# trial from its integer parameters and checks the stored integers, block sums and source
# hashes; `test/certificate.jl` re-contracts it with the reduced SU(2) verifier;
# `test/one_trace.jl` checks the D = 137 and D = 295 intervals against the one-trace targets.
# `--full` also replays both integer witnesses from scratch into temporary files:
# `scripts/certify_trial.jl` asserts that every U, V and block record of
# `results/trial-certificate.json` is reproduced, and `scripts/trials/certify_reduced_trial.jl`
# re-derives the D = 137 intervals of `results/trial-137-certificate.json`.
using Test, JSON3

const MODEL = @__DIR__

"Run `julia --project=MODEL args` in a fresh process; did it succeed?"
julia_ok(args...) = success(pipeline(
    `$(Base.julia_cmd()) --startup-file=no --project=$MODEL $args`; stdout, stderr))

@testset "tracked trial certificates" begin
    for name in ("trials.jl", "certificate.jl", "one_trace.jl")
        @test julia_ok(joinpath(MODEL, "test", name))
    end
end

if "--full" in ARGS
    @testset "exact replays of the integer witnesses" begin
        mktempdir() do directory
            @test julia_ok(joinpath(MODEL, "scripts", "certify_trial.jl"),
                joinpath(directory, "trial-replay.json"),
                joinpath(MODEL, "results", "trial-certificate.json"))
            stored = joinpath(MODEL, "results", "trial-137-certificate.json")
            fresh = joinpath(directory, "trial-137-replay.json")
            @test julia_ok(joinpath(MODEL, "scripts", "trials", "certify_reduced_trial.jl"),
                fresh, stored, "137", "replay")
            # The replay itself does not compare with the stored record; do it here.
            intervals(path) = [(c["length"], c["energy_per_bond_interval"])
                for c in JSON3.read(read(path, String))["certificates"]]
            @test JSON3.read(read(fresh, String))["complete"] == true
            @test intervals(fresh) == intervals(stored)
        end
    end
end
