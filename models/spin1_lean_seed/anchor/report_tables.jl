include("common.jl")
include("exact.jl")
using .AnchorStudy, Printf
using CertifiedGaps.Results
const A=AnchorStudy
const R=joinpath(@__DIR__,"results")
fmt(x,d=6)=@sprintf("%.*f",d,A.num(x))
function table(io,heads,rows)
    println(io,"| ",join(heads," | ")," |")
    println(io,"|",join(fill("---",length(heads)),"|"),"|")
    for row in rows
        println(io,"| ",join(row," | ")," |")
    end
    println(io)
end

function generate()
    data=read_exact_json(joinpath(R,"analysis.json"))
    cert=read_exact_json(joinpath(@__DIR__,"certificate.json"))
    hom=Any[]
    for f in sort(readdir(R))
        startswith(f,"homotopy_n") && endswith(f,"_binned.json") || continue
        x=read_exact_json(joinpath(R,f));x["file"]=f;push!(hom,x)
    end
    sort!(hom;by=x->(A.ival(x["n"]),x["twist"],A.ival(x["Sz"])))
    @assert length(hom)==26 "expected all nonzero-count sectors"
    @assert all(x["complete"] for x in hom) "continuation incomplete"
    path=joinpath(@__DIR__,"REPORT.md")
    prefix=first(split(read(path,String),"<!-- GENERATED TABLES -->"))
    open(path,"w") do io
        print(io,prefix);println(io,"<!-- GENERATED TABLES -->\n")
        println(io,"## Numerical and exact-check results\n")
        println(io,"Ring energies and losses in the next table are NON-CERTIFYING. "*
            "K counts are for the uniform-cutoff comparison.\n")
        rows=Any[]
        for x in data["summary"]
            sector=only(filter(s->s["n"]==x["n"] && s["twist"]==x["twist"] &&
                A.ival(s["Sz"])==0,data["sectors"]))
            ks=sector["k_H0"]
            push!(rows,[x["n"],x["twist"],fmt(x["E0"]),fmt(x["ground_loss"]),
                fmt(x["count_ground_loss"]),sector["k_H_tail"],ks["bare"],ks["ground"],ks["count"]])
        end
        table(io,["n","g","E0(H)","ground-opt loss","count-opt loss",
            "physical k (m=0)","bare K0","ground-opt K0","count-opt K0"],rows)

        println(io,"Exact comparison counts and next-level lower anchors (both twists). "*
            "Decimals here display the rational values in certificate.json. "*
            "All unlisted sectors have K=0 at both cutoffs.\n")
        rows=Any[]
        for p in cert["parameters"]
            p["kind"]=="count_optimized" || continue
            for (x,y) in zip(p["sector_anchors"],p["binned_tail"]["sector_anchors"])
                A.ival(x["k"])+A.ival(y["k"])==0 && continue
                push!(rows,[p["n"],x["Sz"],x["k"],fmt(Float64(rat(x["rho"]))),
                    y["k"],fmt(Float64(rat(y["rho"])))])
            end
        end
        table(io,["n","Sz","uniform K","uniform next rho","binned K","binned next rho"],rows)

        println(io,"Exact binned-tail arithmetic totals; the analytic Taylor and spectral "*
            "inequalities remain to be formalized in Lean.\n")
        rows=Any[]
        for p in cert["parameters"]
            haskey(p,"binned_tail") || continue
            b=p["binned_tail"]
            push!(rows,[p["n"],fmt(Float64(rat(b["target"])),3),
                @sprintf("%.12g",Float64(rat(b["tail_upper"]))),"0.000015"])
        end
        table(io,["n","target rho","tail upper","budget"],rows)

        println(io,"First five ordered losses in Sz=0 for the count-optimized comparison "*
            "(NON-CERTIFYING). Full per-sector losses and both definitions of k are in "*
            "[results/sectors.csv](results/sectors.csv).\n")
        rows=[[x["n"],x["twist"],fmt.(x["low_level_losses"]["count"])...]
            for x in data["sectors"] if A.ival(x["Sz"])==0]
        table(io,["n","g","loss 1","loss 2","loss 3","loss 4","loss 5"],rows)

        println(io,"Best sampled multiblock comparison at the uniform cutoff "*
            "(NON-CERTIFYING selection).\n")
        rows=Any[]
        for n in (8,10,12)
            xs=filter(x->A.ival(x["n"])==n,data["alternative_partitions"])
            sort!(xs;by=x->A.ival(x["k0"]))
            x=first(xs)
            push!(rows,[n,join(x["lengths"],"+"),x["family"],x["k0"],x["ktotal"]])
        end
        table(io,["n","partition","parameters","Sz=0 K","total K"],rows)

        println(io,"Binned-tail continuation, all initially nonzero sectors "*
            "(**NON-CERTIFYING**). Final K is the retained count, which may exceed the "*
            "actual number below the target. Timings include numerical eigensolves, "*
            "compilation, and contention; they are not Lean timings. Negative Sz shares "*
            "the same spectrum. Each linked file lists k and rho at every stage.\n")
        rows=Any[]
        csv=open(joinpath(R,"homotopy_summary.csv"),"w")
        println(csv,"n,twist,Sz,stages,attempts,initial_k,final_k,physical_k_binned,seconds,min_margin,sum_K,sum_K2,pair_products")
        for x in hom
            ss=x["stages"];ks=[A.ival(s["k_in"]) for s in ss]
            margin=minimum(A.num(s["margin"]) for s in ss)
            ref=read_exact_json(joinpath(R,"n$(x["n"])_g$(x["twist"]).json"))
            s=only(filter(s->s["Sz"]==x["Sz"],ref["rows"]))
            phys=count(<(A.num(x["target"])),A.num.(s["values"]))
            @assert A.ival(x["final_k"])>=phys
            push!(rows,[x["n"],x["twist"],x["Sz"],length(ss),
                "$(x["initial_k"]) → $(x["final_k"])",phys,fmt(x["seconds"],1),
                @sprintf("%.2g",margin),"[stages](results/$(x["file"]))"])
            vals=[x["n"],x["twist"],x["Sz"],length(ss),x["attempts"],
                x["initial_k"],x["final_k"],phys,x["seconds"],margin,sum(ks),
                sum(abs2,ks),sum(3*k*(k+1)÷2 for k in ks)]
            println(csv,join(vals,","))
        end
        close(csv)
        table(io,["n","g","Sz","stages","initial → final K","physical k",
            "Julia seconds","min margin","data"],rows)

        println(io,"Naive three-Gram packed-dot work for completed binned runs. "*
            "Only nonnegative Sz is computed; negative Sz is reused. The final column "*
            "applies a hypothetical 2.5-second budget per dot, not a measured cost.\n")
        rows=Any[]
        for n in (8,10,12),g in ("I","P")
            hs=filter(x->A.ival(x["n"])==n && x["twist"]==g,hom)
            ks=[A.ival(s["k_in"]) for x in hs for s in x["stages"]]
            dots=sum(3*k*(k+1)÷2 for k in ks)
            push!(rows,[n,g,length(ks),sum(ks),dots,fmt(dots*2.5/3600,1)])
        end
        table(io,["n","g","sector-stages","sum K","pair products","budget hours"],rows)
    end
    write(path,rstrip(read(path,String))*"\n")
    println("Wrote REPORT.md and results/homotopy_summary.csv from 26 completed runs.")
end
generate()
