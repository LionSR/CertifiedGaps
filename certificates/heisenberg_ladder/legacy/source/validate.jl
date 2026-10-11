include("ladder.jl")
checks=Any[]
for script in ("smoke.jl","test/runtests.jl")
    logname=script=="smoke.jl" ? "smoke.log" : "tests.log"
    seconds=@elapsed open(joinpath(@__DIR__,"results",logname),"w") do io
        run(pipeline(`$(Base.julia_cmd()) --project=$(joinpath(ROOT,"julia")) $(joinpath(@__DIR__,script))`,stdout=io,stderr=io))
    end
    script=="test/runtests.jl" && @assert seconds<120
    push!(checks,Dict("script"=>script,"seconds"=>seconds,"passed"=>true))
end
write_result(joinpath(@__DIR__,"results","validation.json"),Dict(
    "metadata"=>provenance(Dict("checks"=>["independent core Hamiltonian", "full dense small spectra",
        "exact RP", "stored rational seeds", "source hashes", "complete dimensions"]);certified=true),
    "checks"=>checks))
println(checks)
