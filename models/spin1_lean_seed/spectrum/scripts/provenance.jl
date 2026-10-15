# Record SHA-256 of every source and stored result plus the generation environment.
# The generation-time run metadata (CPU, packages, base commit) is carried over from an existing
# results/provenance.json, because the data were produced once and only re-hashed afterwards.
include("../src/common.jl")
using Pkg
sources=String[]
for (dir,_,files) in walkdir(ROOT)
    any(d->occursin(joinpath(ROOT,d),dir),("runtime","results")) && continue
    for f in files
        endswith(f,".jl") || endswith(f,".lean") || endswith(f,".toml") || continue
        push!(sources,relpath(joinpath(dir,f),ROOT))
    end
end
path=joinpath(DAT,"provenance.json")
previous=isfile(path) ? JSON3.read(read(path,String),Dict{String,Any}) : Dict{String,Any}()
records=metadata(sources)
records["source_hashes_scope"]=string("Current layout. The metadata inside each spectrum_*.json ",
    "and summary.json/analysis.json records hashes of the flat layout at generation time.")
records["run_metadata"]=get(previous,"run_metadata",Dict("CPU"=>Sys.CPU_NAME,
    "visible_cpu_threads"=>Sys.CPU_THREADS,"BLAS_threads"=>BLAS.get_num_threads(),
    "packages"=>Dict(p.name=>string(p.version) for p in values(Pkg.dependencies())
                     if !isnothing(p.version))))
haskey(previous,"development_note") &&
    (records["development_note"]=previous["development_note"])
records["hashed_on"]=Dict("CPU"=>Sys.CPU_NAME,"os"=>string(Sys.KERNEL))
records["input_hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(ROOT,f)))) for f in
    ("../../spin1/src/thermal.jl","../../spin1/results/source1.json",
     "../../spin1/results/source2.json"))
records["data_hashes"]=Dict(basename(p)=>bytes2hex(sha256(read(p)))
    for p in readdir(DAT;join=true)
    if (endswith(p,".csv") || endswith(p,".json") || endswith(p,".txt")) &&
       basename(p)!="provenance.json" && !startswith(basename(p),"filter_s") &&
       basename(p) ∉ ("eigenvalues.csv","cost_sectors.csv","dense_validation.csv"))
writejson(path,records)
