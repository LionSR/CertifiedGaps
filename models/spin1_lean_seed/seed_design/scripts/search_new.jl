include(joinpath(@__DIR__, "..", "src", "common.jl"))

function load_spectra(nmax)
    Dict((n,g)=>JSON3.read(read(joinpath(CACHE,"spectrum-$n-$g.json"),String),
        Dict{String,Any})
        for n in 4:nmax for g in 0:2)
end

function thermal_table(spectra,beta)
    table=Dict{Tuple{Int,Int},Float64}(); tails=Dict{Tuple{Int,Int},Float64}()
    for ((n,g),spec) in spectra
        z=0.0; tail=0.0
        for s in spec["sectors"]
            es=s["eigenvalues"]
            z+=s["multiplicity"]*sum(exp(-beta*(e+n*SHIFT)) for e in es)
            # Heuristic ONLY: Lanczos has no certified counting/completeness anchor.
            tail+=s["multiplicity"]*(s["dim"]-length(es))*exp(-beta*(last(es)+n*SHIFT))
        end
        table[(n,g)]=z;tails[(n,g)]=tail
    end
    table,tails
end

function main()
    nmax=isempty(ARGS) ? 8 : parse(Int,first(ARGS))
    spectra=load_spectra(nmax)
    checks=Any[]
    for source in (1,2)
        beta=source==1 ? 12.25 : 10.5
        t,tails=thermal_table(spectra,beta)
        for (key,v) in stored_table(source)
            haskey(t,key) || continue
            push!(checks,Dict("source"=>source,"n"=>key[1],"g"=>key[2],
                "center_difference"=>t[key]-Float64(v),"tail_proxy"=>tails[key]))
        end
    end
    records=Any[]
    for beta in sort(unique(vcat(collect(4.0:0.5:20.0),[10.5,12.25]))),
        mode in ("Z","JN","ION"), eps in (3e-5,1e-4,1e-3), odd in (false,true)
        table,tails=thermal_table(spectra,beta)
        maximum(values(tails))<eps/100 || continue
        input=prepare(table,nmax,eps,mode;allowodd=odd,beta)
        rows=scan_seeds(input;trials=(48,60,72,84,96,108,120,144,168),multipliers=1.0:0.25:6.0)
        rec=Dict("nmax"=>nmax,"mode"=>mode,"eps"=>eps,"odd"=>odd,
            "base_beta"=>beta,"caps"=>[f.cap for f in input.filters],
            "filters"=>input.filters,"masses"=>input.masses,"viable"=>length(rows),
            "best"=>isempty(rows) ? nothing : first(rows),
            "best_per_trial"=>[first(filter(x->x.trial==n,rows))
                for n in unique(x.trial for x in rows)],
            "best_existing_trial"=>let rr=filter(r->r.trial in (72,120),rows)
                isempty(rr) ? nothing : first(rr)
            end)
        push!(records,rec)
        if !isempty(rows)
            println((nmax,beta,mode,eps,odd)," best=",
                (rows[1].gap,rows[1].p,rows[1].q,rows[1].k,rows[1].r,rows[1].trial))
            flush(stdout)
        end
    end
    write_compact(joinpath(RAW,"new_search_$nmax.json"),Dict(
        "metadata"=>metadata("new spectra; all non-certifying"),
        "checks_against_stored_centers"=>checks,"records"=>records))
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
