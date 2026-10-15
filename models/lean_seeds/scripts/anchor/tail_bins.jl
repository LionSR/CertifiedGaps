# NON-CERTIFYING search for a cheaper tail bound using the comparison density of states.
include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.AnchorStudy
using CertifiedGaps.Results
const A = LeanSeeds.AnchorStudy

function binned_tail(spec,n,rho;step=1/8,bins=32)
    total=0.
    for (m,es) in spec
        mult=m==0 ? 1 : 2
        last=0
        for j in 0:bins-1
            k=searchsortedlast(es,rho+(j+1)*step)
            total+=mult*(k-last)*exp(-A.BETA*(rho+j*step+n*A.SHIFT))
            last=k
        end
        total+=mult*(length(es)-last)*exp(-A.BETA*(rho+bins*step+n*A.SHIFT))
    end
    return total
end

function run_bins()
    ps=read_exact_json(joinpath(A.RESULTS,"count_optimized.json"))["rows"]
    rows=Any[]
    for p in ps
        n=A.ival(p["n"]);u,v,c=A.num.((p["u"],p["v"],p["c"]))
        spec=A.comparison_spectra([n÷2,n÷2],u,v;c)
        lo,hi=A.tail_cutoff(n)-2,A.tail_cutoff(n)
        for _ in 1:50
            mid=(lo+hi)/2
            if binned_tail(spec,n,mid)<=A.TOLERANCE/2
                hi=mid
            else
                lo=mid
            end
        end
        rho=ceil(hi*1000)/1000
        ks=[count(<(rho),spec[m]) for m in 0:n]
        push!(rows,Dict("n"=>n,"rho"=>rho,"k"=>ks,"step"=>1/8,"bins"=>32,
            "tail_upper_float"=>binned_tail(spec,n,rho),"k_total"=>ks[1]+2sum(ks[2:end])))
        println("BINNED n=$n rho=$rho counts=$ks tail=$(binned_tail(spec,n,rho))")
    end
    write_result(joinpath(A.RESULTS,"binned_tails.json"),Dict(
        "metadata"=>A.metadata(Dict("tail_budget"=>A.TOLERANCE/2)),"rows"=>rows))
end
run_bins()
