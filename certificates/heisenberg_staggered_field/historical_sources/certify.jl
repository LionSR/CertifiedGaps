include("common.jl")
h=isempty(ARGS) ? 1//1 : rat(ARGS[1]); maxL=length(ARGS)>1 ? parse(Int,ARGS[2]) : 16
scan=JSON3.read(read(joinpath(ROOT,"results","scan_$(tag(h))_$(maxL).json"),String),Dict{String,Any})
best=scan["best"];@assert best!==nothing
n=best["n0"];beta=Rational{BigInt}(best["beta"]);exact_rp(h)
start=time();a=all_spectra(n,h;verified=true);b=all_spectra(2n,h;verified=true)
seeds=enclosed_seeds(a.records,b.records,beta)
bound=certified_gap_bound(seeds.s,seeds.t,beta);@assert bound.bound!==nothing
write_result(joinpath(ROOT,"results","certificate_$(tag(h)).json"),Dict(
    "metadata"=>provenance(Dict("h"=>string(h),"n0"=>n);certified=true),
    "h"=>string(h),"n0"=>n,"sites_per_cell"=>1,"beta"=>string(beta),
    "epsilon_s"=>string(seeds.s),"epsilon_t"=>string(seeds.t),
    "bound"=>bound.bound,"logs"=>seeds.logs,"seconds"=>time()-start,
    "small_sectors"=>a.records,"large_sectors"=>b.records,
    "reference_gap"=>b.spec.E[2]-b.spec.E[1],"ratio"=>bound.bound/(b.spec.E[2]-b.spec.E[1])))
println("CERTIFIED h=$h n=$n beta=$beta bound=$(bound.bound) seconds=$(time()-start)")
