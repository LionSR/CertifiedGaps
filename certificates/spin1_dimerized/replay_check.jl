include("common.jl")
using Test

original_path=joinpath(ROOT,"results","certificate_d7.json")
replay_path=joinpath(ROOT,"results","replay_certificate_d7.json")
original=JSON3.read(read(original_path,String),Dict{String,Any})
replay=JSON3.read(read(replay_path,String),Dict{String,Any})

# Timing and provenance are expected to differ. Every mathematical output,
# including all block error bounds, must reproduce exactly in this environment.
function mathematical(x)
    x isa AbstractDict && return Dict(k=>mathematical(v) for (k,v) in x if k ∉ ("metadata","seconds"))
    x isa AbstractVector && return mathematical.(x)
    x
end
@testset "Fresh full spectral replay" begin
    @test mathematical(original)==mathematical(replay)
    @test replay["bound"]==certified_gap_bound(parserat(replay["epsilon_s"]),parserat(replay["epsilon_t"]),parserat(replay["beta0"])).bound
    for (path,hash) in replay["metadata"]["source_sha256"]
        @test sha256_file(joinpath(ROOT,path))==hash
    end
    for (key,cells) in (("sectors_n",3),("sectors_2n",6))
        @test sum(s["dimension"]*s["multiplicity"] for s in replay[key])==9^cells
    end
end
meta=provenance(Dict("test"=>"replay_check.jl");certified=true)
meta["source_sha256"]["replay_check.jl"]=sha256_file(@__FILE__)
write_result(joinpath(ROOT,"results","replay_comparison.json"),Dict(
    "metadata"=>meta,"status"=>"passed","all_mathematical_outputs_identical"=>true,
    "original_sha256"=>sha256_file(original_path),"replay_sha256"=>sha256_file(replay_path),
    "original_seconds"=>original["seconds"],"replay_seconds"=>replay["seconds"],
    "bound"=>replay["bound"],"epsilon_s"=>replay["epsilon_s"],"epsilon_t"=>replay["epsilon_t"]))
