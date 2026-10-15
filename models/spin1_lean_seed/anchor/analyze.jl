include("common.jl")
using .AnchorStudy, LinearAlgebra, Printf
using CertifiedGaps.Results
const A=AnchorStudy

function run_analysis()
    ground=read_exact_json(joinpath(@__DIR__,"results","optimized.json"))["rows"]
    countopt=read_exact_json(joinpath(@__DIR__,"results","count_optimized.json"))["rows"]
    rows=Any[]; summary=Any[]; alternatives=Any[]
    for n in (8,10,12)
        p=only(filter(x->A.ival(x["n"])==n,countopt))
        q=only(filter(x->A.ival(x["n"])==n,ground))
        u,v,c=A.num.((p["u"],p["v"],p["c"]))
        specs=Dict("bare"=>A.comparison_spectra([n÷2,n÷2],0.,0.;c=-2.),
            "ground"=>A.comparison_spectra([n÷2,n÷2],A.num(q["u"]),A.num(q["v"]);
                c=A.num(q["c"])),
            "count"=>A.comparison_spectra([n÷2,n÷2],u,v;c))
        cut=A.tail_cutoff(n)
        partitions=n==8 ? [[4,2,2],[3,3,2],[2,2,2,2]] :
            n==10 ? [[4,3,3],[4,4,2],[3,3,2,2],[2,2,2,2,2]] :
            [[4,4,4],[3,3,3,3],[6,3,3],[6,4,2],[2,2,2,2,2,2]]
        for lengths in partitions, (label,x,y) in (("bare",0.,0.),("boundary",u,v))
            sp=A.comparison_spectra(lengths,x,y)
            ks=[count(<(cut),sp[m]) for m in 0:n]
            push!(alternatives,Dict("n"=>n,"lengths"=>lengths,"family"=>label,
                "k0"=>ks[1],"ktotal"=>ks[1]+2sum(ks[2:end]),"ground"=>sp[0][1]))
        end
        for g in ("I","P")
            r=read_exact_json(joinpath(@__DIR__,"results","n$(n)_g$g.json"))
            e0=A.num(r["E0"])
            totals=Dict("tail"=>0,"two"=>0,"bare"=>0,"ground"=>0,"count"=>0)
            z=0.
            for row in r["rows"]
                m=A.ival(row["Sz"]); d=A.ival(row["dimension"])
                e=A.num.(row["values"]); mult=m==0 ? 1 : 2
                kp=count(<(cut),e); k2=count(<(e0+2),e)
                counts=Dict(label=>count(<(cut),sp[m]) for (label,sp) in specs)
                K=counts["count"]
                rho=K<d ? specs["count"][m][K+1] : nothing
                # What anchor would be available using only the physically relevant count?
                rho_small=specs["count"][m][kp+1]
                losses=Dict(label=>e[1:min(5,length(e))]-sp[m][1:min(5,length(e))]
                    for (label,sp) in specs)
                totals["tail"]+=mult*kp; totals["two"]+=mult*k2
                for label in keys(counts)
                    totals[label]+=mult*counts[label]
                end
                z+=mult*sum(exp.(-A.BETA.*(e.+n*A.SHIFT)))
                push!(rows,Dict("n"=>n,"twist"=>g,"Sz"=>m,"dimension"=>d,
                    "k_H_tail"=>kp,"k_H_E0plus2"=>k2,"k_H0"=>counts,
                    "rho_at_physical_k"=>rho_small,"rho_at_comparison_k"=>rho,
                    "H_count_below_rho"=>count(<(rho),e),"H_min"=>first(e),
                    "low_level_losses"=>losses,"H0_ground"=>specs["count"][m][1]))
            end
            push!(summary,Dict("n"=>n,"twist"=>g,"E0"=>e0,"cutoff"=>cut,
                "cutoff_minus_E0"=>cut-e0,"totals"=>totals,"low_trace"=>z,
                "ground_loss"=>e0-specs["ground"][0][1],
                "count_ground_loss"=>e0-specs["count"][0][1]))
        end
    end
    write_result(joinpath(@__DIR__,"results","analysis.json"),Dict(
        "metadata"=>A.metadata(Dict("meaning"=>"all floating spectra non-certifying")),
        "summary"=>summary,"sectors"=>rows,"alternative_partitions"=>alternatives))
    open(joinpath(@__DIR__,"results","sectors.csv"),"w") do io
        println(io,"n,twist,Sz,dimension,k_H_tail,k_H_E0plus2,k_bare,k_ground,k_count,rho_small,rho_count,H_count_below_rho,loss1,loss2,loss3,loss4,loss5")
        for r in rows
            vals=Any[r[k] for k in ("n","twist","Sz","dimension","k_H_tail","k_H_E0plus2")]
            append!(vals,[r["k_H0"][k] for k in ("bare","ground","count")])
            append!(vals,[r["rho_at_physical_k"],r["rho_at_comparison_k"],r["H_count_below_rho"]])
            append!(vals,r["low_level_losses"]["count"])
            println(io,join(vals,","))
        end
    end
    for r in summary
        println(r)
    end
end
run_analysis()
