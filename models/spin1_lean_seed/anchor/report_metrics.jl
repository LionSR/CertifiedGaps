include("common.jl")
using .AnchorStudy, Printf
using CertifiedGaps.Results
const A=AnchorStudy
const R=joinpath(@__DIR__,"results")

data=read_exact_json(joinpath(R,"analysis.json"))
println("SUMMARY")
for x in data["summary"]
    println("n=$(x["n"]) g=$(x["twist"]) E0=$(x["E0"]) losses="*
        "$(x["ground_loss"]),$(x["count_ground_loss"]) counts=$(x["totals"])")
end
println("EXACT ANCHORS")
cert=read_exact_json(joinpath(@__DIR__,"certificate.json"))
for p in cert["parameters"]
    p["kind"]=="count_optimized" || continue
    println("n=$(p["n"]) u=$(p["u"]) v=$(p["v"]) c=$(p["c"])")
    for x in p["sector_anchors"][1:5]
        println("  m=$(x["Sz"]) k=$(x["k"]) rho=$(x["rho"])")
    end
    println("  binned counts=",[x["k"] for x in p["binned_tail"]["sector_anchors"]])
    println("  binned tail=",p["binned_tail"]["tail_upper"])
end
println("ALTERNATIVES")
for n in (8,10,12)
    alts=filter(x->A.ival(x["n"])==n,data["alternative_partitions"])
    sort!(alts;by=x->A.ival(x["k0"]))
    println(first(alts))
end
println("ANDERSON SCALARS (non-certifying)")
for l in 2:6
    e=first(A.block_spectra(l,0.,0.)[0])
    println("cluster=$l min=$e per_bond=$(e/(l-1)) max_Sz_dim=$(length(A.basis(l,0)))")
end
println("HOMOTOPY")
for f in sort(readdir(R))
    startswith(f,"homotopy_n") && endswith(f,".json") || continue
    x=read_exact_json(joinpath(R,f));ss=x["stages"]
    ks=[A.ival(s["k_in"]) for s in ss]
    margin=minimum(A.num(s["margin"]) for s in ss)
    products=sum(3*k*(k+1)÷2 for k in ks)
    println("$f complete=$(x["complete"]) stages=$(length(ss)) attempts=$(x["attempts"]) "*
        "k=$(x["initial_k"])->$(x["final_k"]) sumK=$(sum(ks)) sumK2=$(sum(abs2,ks)) "*
        "pair_products=$products seconds=$(x["seconds"]) min_margin=$margin")
end
