using CertifiedGaps.Results
elapsed=@elapsed include("test/runtests.jl")
@assert elapsed < 120 "Fast checker exceeded the requested two-minute budget"
write_result(joinpath(@__DIR__,"results","tests.json"),(;status="passed",seconds=elapsed,
    julia=string(VERSION),source_hashes=source_hashes(@__DIR__;files=["run_tests.jl"])))
println("All certificate checks passed in $elapsed seconds.")
