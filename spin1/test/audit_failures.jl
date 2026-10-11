using Test, JSON3
import Spin1Certificates as SC
const RESULTS = abspath(get(ARGS, 1, joinpath(@__DIR__, "../results")))

function mutated_audit(mutate!)
    mktempdir() do directory
        for name in ("source1", "source2", "trials", "scalars", "arithmetic", "itensors")
            cp(joinpath(RESULTS, "$name.json"), joinpath(directory, "$name.json"))
        end
        SC.write_result(
            joinpath(directory, "verification.json"),
            Dict("status" => "passed", "stale" => true),
        )
        mutate!(directory)
        @test_throws AssertionError SC.audit_results(directory)
        @test SC.read_exact_json(joinpath(directory, "verification.json"))["status"] == "failed"
    end
end

@testset "certificate audit rejects corrupted or incomplete artifacts" begin
    mutated_audit() do directory
        path = joinpath(directory, "source1.json")
        data = SC.read_exact_json(path)
        pop!(data["certificates"])
        SC.write_result(path, data)
    end
    mutated_audit() do directory
        path = joinpath(directory, "source2.json")
        data = SC.read_exact_json(path)
        c = first(data["certificates"])
        c["total"] = string(SC.integer(c["total"]) + 1)
        SC.write_result(path, data)
    end
    mutated_audit() do directory
        path = joinpath(directory, "trials.json")
        data = SC.read_exact_json(path)
        first(data["certificates"])["positive_margin"] = "0"
        SC.write_result(path, data)
    end
    mutated_audit() do directory
        path = joinpath(directory, "scalars.json")
        data = SC.read_exact_json(path)
        data["source1"]["rounded_defects"][end][1] = "0"
        SC.write_result(path, data)
    end
end
