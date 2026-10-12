# Preserve full dimension audits and beta-grid partition sums, rather than millions of
# unverified individual levels. The rigorous eigenvalue enclosures are never compacted.
include("ladder.jl")
using JSON3
for name in readdir(joinpath(@__DIR__,"results"))
    startswith(name,"spectrum_J") && endswith(name,".json") || continue
    path=joinpath(@__DIR__,"results",name)
    data=JSON3.read(read(path,String),Dict{String,Any})
    haskey(data,"energies") || continue
    es=Float64.(pop!(data,"energies"))
    data["lowest_32_numerical_levels"]=es[1:min(32,length(es))]
    data["logZ_grid"]=[Dict("beta"=>string(Rational{Int}(beta)),
        "logZ"=>log_partition(es,beta)) for beta in 0.25:0.25:32]
    data["storage_note"]="Full spectrum computed; retained logZ grid and low levels. Rerun scan.jl to regenerate all levels."
    write_result(path,data)
end
