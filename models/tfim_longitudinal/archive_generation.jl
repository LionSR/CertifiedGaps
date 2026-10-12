# Preserve exact historical generation sources after the upstream namespace rename.
using CertifiedGaps.Results, JSON3, SHA
const ROOT=@__DIR__
const REPO=normpath(joinpath(ROOT,"../.."))
const ORIGINAL="e67c0925f6d069d2dccf34bfbb3dbfc9971c6d0e"
archive=joinpath(ROOT,"results/original_generation_sources.json")
# The portable archive also works when the original local commit is unavailable.
sources=isfile(archive) ? Dict{String,String}(
    JSON3.read(read(archive,String),Dict{String,Any})["sources"]) : Dict{String,String}()
for path in readdir(joinpath(ROOT,"results");join=true)
    startswith(basename(path),"g") && endswith(path,".json") || continue
    d=JSON3.read(read(path,String),Dict{String,Any})
    for (name,hash) in d["metadata"]["source_sha256"]
        relative=relpath(normpath(joinpath(ROOT,name)),REPO)
        source=get!(sources,name) do
            read(`git -C $REPO show $(ORIGINAL*":"*relative)`,String)
        end
        @assert bytes2hex(sha256(source))==hash
    end
end
write_result(joinpath(ROOT,"results/original_generation_sources.json"),Dict(
    "generation_commit"=>ORIGINAL,"sources"=>sources))
# Mathematical generation sources differ only by the upstream package rename.
checked=String[]
for (name,source) in sources
    (startswith(name,"../../julia/src/") || name in
        ("certify.jl","scan.jl","check.jl","test/runtests.jl")) || continue
    current=replace(name,"ThermalMPOGap"=>"CertifiedGaps")
    @assert read(joinpath(ROOT,current),String)==replace(source,
        "ThermalMPOGap"=>"CertifiedGaps")
    push!(checked,current)
end
write_result(joinpath(ROOT,"results/namespace_migration.json"),Dict(
    "status"=>"passed","generation_commit"=>ORIGINAL,"checked_sources"=>sort(checked),
    "current_hashes"=>Dict(p=>sha256_file(joinpath(ROOT,p)) for p in checked),
    "method"=>"Exact byte equality after replacing only the package namespace"))
println("Archived ",length(sources)," sources; namespace-only equality passed for ",
    length(checked)," mathematical sources")
