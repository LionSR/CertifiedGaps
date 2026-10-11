# One command for the fast tests, independent spectral replay, and provenance checks.
start=time()
include("test/runtests.jl")
fast_seconds=time()-start
paths=sort(filter(p->startswith(basename(p),"g") && endswith(p,".json"),
    readdir(joinpath(C.ROOT,"results");join=true)))
@assert length(paths)==30
ENV["DEEP_CHECK"]="1"
start=time()
for path in paths
    @assert check_result(path)
    println("Replayed ",basename(path)); flush(stdout)
end
deep_seconds=time()-start
using SHA
archive=JSON3.read(read(joinpath(C.ROOT,"results/original_generation_sources.json"),String),Dict{String,Any})
hash_count=0
for path in paths
    d=JSON3.read(read(path,String),Dict{String,Any})
    @assert d["metadata"]["certified"]
    for (name,hash) in d["metadata"]["source_sha256"]
        current=joinpath(C.ROOT,name)
        live=isfile(current) && sha256_file(current)==hash
        historical=haskey(archive["sources"],name) &&
            bytes2hex(sha256(archive["sources"][name]))==hash
        @assert live || historical
        global hash_count+=1
    end
end
write_result(joinpath(C.ROOT,"results","validation.json"),Dict(
    "status"=>"passed", "parameter_certificates"=>length(paths),
    "fast_test_seconds_including_compilation"=>fast_seconds,
    "deep_replay_seconds"=>deep_seconds,"source_hash_comparisons"=>hash_count,
    "blas"=>string(C.BLAS.get_config()),"gradual_underflow"=>!get_zero_subnormals(),
    "certificate_sha256"=>Dict(basename(p)=>sha256_file(p) for p in paths),
    "metadata"=>C.metadata(Dict("kind"=>"independent complete replay");certified=true)))
println("Fast checks: ",fast_seconds," s; deep replay: ",deep_seconds," s")
