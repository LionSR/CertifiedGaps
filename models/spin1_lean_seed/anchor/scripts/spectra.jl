include("../src/common.jl")
using .AnchorStudy, LinearAlgebra, Serialization
using CertifiedGaps.Results
const A = AnchorStudy

function run_spectra(ns)
    for n in ns, twist in (1,-1)
        name="n$(n)_g$(twist==1 ? "I" : "P")"
        path=joinpath(A.RESULTS,name*".json")
        rows=Any[]
        e0=Inf
        for m in 0:n
            t=time()
            H=A.ring(n,m;twist)
            cutoff=m==0 ? A.tail_cutoff(n)+0.3 : max(A.tail_cutoff(n),e0+2)+0.08
            spec=A.low_spectrum(H,cutoff)
            m==0 && (e0=first(spec.values))
            push!(rows,Dict("Sz"=>m,"dimension"=>size(H,1),"values"=>spec.values,
                "residual"=>spec.residual,"covered_cutoff"=>spec.complete,
                "seconds"=>time()-t))
            println("$name m=$m d=$(size(H,1)) E=$(first(spec.values)) "*
                "count2=$(count(<(e0+2),spec.values)) "*
                "counttail=$(count(<(A.tail_cutoff(n)),spec.values)) seconds=$(time()-t)")
            flush(stdout)
            write_result(path,Dict("metadata"=>A.metadata(Dict("n"=>n,"twist"=>twist)),
                "n"=>n,"twist"=>twist,"E0"=>e0,"tail_cutoff"=>A.tail_cutoff(n),
                "rows"=>rows))
        end
    end
end
run_spectra(isempty(ARGS) ? [8,10,12] : parse.(Int,ARGS))
