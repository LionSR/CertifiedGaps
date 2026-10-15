# Verify the stored spin-1 certificate.
#
#     julia --project=models/spin1 models/spin1/check.jl          # fast: stored-output audit
#     julia --threads=4 --project=models/spin1 models/spin1/check.jl --full   # ~45 min
#
# Fast mode audits a temporary copy of `results/` (the audit rewrites `verification.json` in the
# directory it checks, and the experiments pin the SHA-256 of the committed file), checks that
# the audit reproduces every component hash and the implementation fingerprint recorded in the
# committed `verification.json`, runs the corrupted-artifact tests on that copy, and replays the
# two `experiments/` verifiers. `--full` instead recomputes every thermal and trial contraction
# and the ITensors controls into a temporary directory, audits that fresh run, and compares its
# component outputs with the committed ones.
#
# This file sits outside the implementation fingerprint (`src/`, `scripts/`, `test/` and the
# three environment files), so editing it does not invalidate `results/verification.json`.
# Both modes need the upstream evidence: run `julia tools/fetch_upstream.jl` first.

using Test
import Spin1Certificates as SC

const MODEL = @__DIR__
const RESULTS = joinpath(MODEL, "results")
const COMPONENTS = ("source1", "source2", "trials", "scalars", "arithmetic", "itensors")

isdir(SC.SOURCE) ||
    error("upstream evidence missing at $(SC.SOURCE); run `julia tools/fetch_upstream.jl`")

"Run a Julia file of this project in a fresh process (so its own `ARGS` apply); did it pass?"
run_file(path, args...) = success(pipeline(
    `$(Base.julia_cmd()) --startup-file=no --project=$MODEL $path $args`; stdout, stderr))

"The committed audit record, the reference for every replay."
committed() = SC.read_exact_json(joinpath(RESULTS, "verification.json"))

"Copy the six component outputs of `results/` into a fresh directory."
function copy_components(directory)
    for name in COMPONENTS
        cp(joinpath(RESULTS, "$name.json"), joinpath(directory, "$name.json"))
    end
    return directory
end

function check_fast()
    reference = committed()
    @testset "stored-output audit of a copy of results/" begin
        mktempdir() do directory
            copy_components(directory)
            SC.audit_results(directory)
            audit = SC.read_exact_json(joinpath(directory, "verification.json"))
            @test audit["status"] == "passed"
            @test audit["execution_mode"] == "stored-output-audit"
            @test audit["thermal_certificates"] == reference["thermal_certificates"]
            @test SC.integer(audit["thermal_certificates"]) == 27
            @test audit["trial_certificates"] == reference["trial_certificates"]
            @test SC.integer(audit["trial_certificates"]) == 2
            # Same implementation as the committed run, and byte-identical component outputs.
            @test audit["implementation_sha256"] == reference["implementation_sha256"]
            for name in COMPONENTS
                @test audit["components"][name]["sha256"] ==
                      reference["components"][name]["sha256"]
            end
            @test run_file(joinpath(MODEL, "test", "audit_failures.jl"), directory)
        end
    end
    @testset "experiments/ verifiers" begin
        @test run_file(joinpath(MODEL, "experiments", "test_bounds.jl"))
        @test run_file(joinpath(MODEL, "experiments", "test_temperature_lift.jl"))
    end
end

function check_full()
    reference = committed()
    @testset "fresh recomputation into a temporary directory" begin
        mktempdir() do directory
            reproduce = joinpath(MODEL, "scripts", "reproduce.jl")
            threads = Threads.nthreads()
            run(`$(Base.julia_cmd()) --startup-file=no --threads=$threads --project=$MODEL
                $reproduce all $directory`)
            fresh = SC.read_exact_json(joinpath(directory, "verification.json"))
            @test fresh["status"] == "passed"
            @test fresh["execution_mode"] == "fresh-recomputation"
            @test fresh["implementation_sha256"] == reference["implementation_sha256"]
            # Timings differ between runs; the certified integers and intervals must not.
            for (name, fields) in (
                "source1" => ("certificates",),
                "source2" => ("certificates",),
                "trials" => ("certificates", "multiplicity_dimensions"),
                "scalars" => ("source1", "source2", "bootstrap"),
                "arithmetic" => ("matrices", "triad", "bounds"),
            )
                old = SC.read_exact_json(joinpath(RESULTS, "$name.json"))
                new = SC.read_exact_json(joinpath(directory, "$name.json"))
                for field in fields
                    @test old[field] == new[field]
                end
            end
        end
    end
end

"--full" in ARGS ? check_full() : check_fast()
