include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.SeedDesign
using CertifiedGaps.ExactArithmetic, CertifiedGaps.Results
using LinearAlgebra, JSON3, SHA, Dates, TOML
function main()
    records=Any[]
    for source in (1,2), nmax in (8,10,12), mode in ("Z","JN","ION"),
        eps in (3e-5,1e-4,3e-4,1e-3,3e-3,1e-2), odd in (false,true), cap in (false,true)
        source==2 && (nmax==12 || odd) && continue
        source==1 && mode=="ION" && continue
        !cap && odd && continue
        table=stored_table(source)
        input=prepare(table,nmax,eps,mode;allowodd=odd,masscap=cap,beta=source==1 ? 12.25 : 10.5)
        isnothing(input) && continue
        rows=scan_seeds(input)
        record=Dict("source"=>source,"nmax"=>nmax,"mode"=>mode,"eps"=>eps,
            "odd"=>odd,"capped"=>cap,"caps"=>[f.cap for f in input.filters],
            "filters"=>input.filters,"masses"=>input.masses,"viable"=>length(rows),
            "best"=>isempty(rows) ? nothing : first(rows),
            "best_per_trial"=>[first(filter(x->x.trial==n,rows))
                for n in unique(x.trial for x in rows)])
        push!(records,record)
        if eps in (3e-5,1e-3) && cap
            println((source,nmax,mode,eps,odd)," caps=",record["caps"],
                " best=",isempty(rows) ? nothing :
                    (rows[1].gap,rows[1].p,rows[1].q,rows[1].k,rows[1].r))
            flush(stdout)
        end
    end
    write_compact(joinpath(RAW,"existing_search.json"),Dict(
        "metadata"=>metadata("existing inputs"),"records"=>records))
end
main()
