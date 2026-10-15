include("../src/common.jl")
using .AnchorStudy, LinearAlgebra
using CertifiedGaps.Results
const A=AnchorStudy

function counts(l,u,v,cut)
    c=A.local_min(u,v)
    es=A.block_spectra(l,u,v)
    ks=Int[]
    for m in 0:2l
        k=0
        for q in -l:l
            abs(m-q)<=l || continue
            e,f=es[abs(q)],es[abs(m-q)]
            for x in e
                k += searchsortedfirst(f,cut-2c-x)-1
            end
        end
        push!(ks,k)
    end
    return (k0=ks[1],ktotal=ks[1]+2sum(ks[2:end]),c=c,ks=ks)
end

function run_counts()
    rows=Any[]
    for l in 4:6
        cut=A.tail_cutoff(2l)
        best=(score=(typemax(Int),typemax(Int)),u=0.,v=0.)
        allrows=Any[]
        function consider(u,v)
            r=counts(l,u,v,cut)
            score=(r.k0,r.ktotal)
            push!(allrows,Dict("u"=>u,"v"=>v,"k0"=>r.k0,"ktotal"=>r.ktotal))
            if score<best.score
                best=(score=score,u=u,v=v)
            end
        end
        for u in -1.2:0.1:0.2, v in -0.3:0.05:0.3
            consider(u,v)
        end
        for step in (0.02,0.004,0.0008)
            u0,v0=best.u,best.v
            for du in -4:4,dv in -4:4
                consider(u0+step*du,v0+step*dv)
            end
        end
        u,v=round(best.u;digits=4),round(best.v;digits=4)
        r=counts(l,u,v,cut)
        c=floor((r.c-1e-8)*1e7)/1e7
        push!(rows,Dict("n"=>2l,"u"=>u,"v"=>v,"c"=>c,
            "kind"=>"count_optimized","counts"=>r.ks,"k0"=>r.k0,
            "ktotal"=>r.ktotal,"cutoff"=>cut,"samples"=>length(allrows)))
        write_result(joinpath(A.RESULTS,"scan_n$(2l).json"),allrows)
        println("COUNT n=$(2l) u=$u v=$v c=$c k0=$(r.k0) total=$(r.ktotal)")
        flush(stdout)
    end
    write_result(joinpath(A.RESULTS,"count_optimized.json"),Dict(
        "metadata"=>A.metadata(Dict("objective"=>"lexicographic (Sz=0 count, total count)",
            "global_optimality"=>false)),"rows"=>rows))
end
run_counts()
