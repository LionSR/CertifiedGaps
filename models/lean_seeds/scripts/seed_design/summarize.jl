include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.SeedDesign
using CertifiedGaps.ExactArithmetic, CertifiedGaps.Results
using LinearAlgebra, JSON3, SHA, Dates, TOML

loadjson(path)=JSON3.read(read(path,String),Dict{String,Any})

function used_pairs(filters,mode,nmax;exact=false)
    supports=[Set([nmax]) for _ in filters]
    for (i,f) in enumerate(filters)
        isnothing(f) && continue
        P=exact ? parse_rational.(f["P"]) : f["P"]
        union!(supports[i],[j-1 for (j,c) in enumerate(P) if !iszero(c)])
    end
    pairs=Set{Tuple{Int,Int}}()
    for (i,ns) in enumerate(supports)
        twists=mode=="ION" && i<3 ? (0:2) : mode=="Z" ? (0:0) : (0:1)
        for n in ns,g in twists
            push!(pairs,(n,g))
        end
    end
    sort(collect(pairs))
end

"Work forecasts only: no eigenvalue counts or spectral anchors are certified."
function work(pairs,beta,eps;a=SHIFT)
    rows=Any[]
    dimension_sum=proxy=matvec=gram=complex_proxy=0.0
    for (n,g) in pairs
        path=joinpath(CACHE,"spectrum-$n-$g.json")
        isfile(path) || error("missing spectra for work proxy: $path")
        spec=loadjson(path)
        cutoff=log(2*3.0^n/eps)/beta-n*a
        K=0;mv=0;gr=0;censored=false;persector=Any[]
        for s in spec["sectors"]
            k=count(e->e<=cutoff,s["eigenvalues"])
            censored |= k==length(s["eigenvalues"]) && k<s["dim"]
            K+=s["multiplicity"]*k
            mv+=s["multiplicity"]*s["dim"]*k
            gr+=s["multiplicity"]*s["dim"]*k^2
            push!(persector,Dict("m"=>s["m"],"dimension"=>s["dim"],"low_states"=>k))
        end
        @assert !censored "increase nev for n=$n g=$g before using workload estimate"
        dimension_sum+=3^n;proxy+=3^n*K;matvec+=mv;gram+=gr
        complex_proxy+=(g==2 ? 2 : 1)*3^n*K
        push!(rows,Dict("n"=>n,"g"=>g,"cutoff"=>cutoff,"K"=>K,"censored"=>censored,
            "sectors"=>persector))
    end
    Dict("dimension_sum"=>dimension_sum,"state_proxy"=>proxy,"complex_state_proxy"=>complex_proxy,
        "matvec_volume"=>matvec,"gram_volume"=>gram,"traces"=>length(pairs),"rows"=>rows)
end

function main_summary()
    certificate=loadjson(joinpath(OUT,"certificate.json"))
    design_rows=Any[]
    for d in certificate["designs"]
        input=d["input"];seed=d["seed"];chain=d["chain"]
        pairs=used_pairs(input["filters"],input["mode"],input["nmax"];exact=true)
        W=work(pairs,Float64(parse_rational(input["base_beta"])),
            Float64(parse_rational(input["eps"])))
        push!(design_rows,Dict("name"=>d["name"],"nmax"=>input["nmax"],"mode"=>input["mode"],
            "eps"=>input["eps"],"base_beta"=>input["base_beta"],"beta"=>chain["beta"],
            "trial"=>seed["trial"],"n0"=>seed["n0"],"p"=>Float64(parse_rational(seed["p"])),
            "q"=>Float64(parse_rational(seed["q"])),"gap"=>chain["gap_display"],"work"=>W,
            "caps"=>[isnothing(f) ? nothing : Float64(parse_rational(f["cap"]))
                for f in input["filters"]]))
    end
    baseline=work([(n,g) for n in 4:12 for g in 0:1],12.25,3e-5)
    # A frontier over recorded grid representatives; never a global optimality claim.
    candidates=Any[]
    for (filename,sourcekind) in (("existing_search.json","stored tables"),
        ("new_search_8.json","Float64 spectra"),("new_search_10.json","Float64 spectra"),
        ("refinement.json","stored tables"))
        for row in loadjson(joinpath(RAW,filename))["records"]
            reps=get(row,"best_per_trial",Any[])
            isempty(reps) && continue
            source=get(row,"source",2)
            nmax=get(row,"nmax",10);mode=get(row,"mode","ION")
            beta=get(row,"base_beta",source==1 ? 12.25 : 10.5)
            pairs=used_pairs(row["filters"],mode,nmax)
            W=work(pairs,beta,row["eps"])
            for r in reps
                push!(candidates,Dict("nmax"=>nmax,"mode"=>mode,"eps"=>row["eps"],
                    "base_beta"=>beta,"beta"=>r["beta"],"n0"=>r["n0"],"trial"=>r["trial"],
                    "gap"=>r["gap"],"p"=>r["p"],"q"=>r["q"],"k"=>r["k"],"r"=>r["r"],
                    "trace_count"=>W["traces"],"dimension_sum"=>W["dimension_sum"],
                    "state_proxy"=>W["state_proxy"],"complex_state_proxy"=>W["complex_state_proxy"],
                    "input_kind"=>sourcekind,"certified"=>false))
            end
        end
    end
    dominates(a,b)=a["state_proxy"]<=b["state_proxy"] && a["eps"]>=b["eps"] &&
        a["gap"]>=b["gap"] && a["n0"]<=b["n0"] &&
        (a["state_proxy"]<b["state_proxy"] || a["eps"]>b["eps"] ||
            a["gap"]>b["gap"] || a["n0"]<b["n0"])
    front=unique([b for b in candidates if !any(a->dominates(a,b),candidates)])
    sort!(front;by=r->(r["state_proxy"],r["eps"],-r["gap"]))
    cols=["nmax","mode","eps","base_beta","beta","n0","trial","gap","p","q","k","r",
        "trace_count","dimension_sum","state_proxy","complex_state_proxy","input_kind","certified"]
    open(joinpath(OUT,"pareto.csv"),"w") do io
        println(io,join(cols,","))
        for r in front
            println(io,join((r[c] for c in cols),","))
        end
    end
    control=Any[]
    for n in 4:12,g in 0:(n==12 ? 1 : 2)
        path=joinpath(CACHE,"spectrum-$n-$g.json")
        isfile(path) || continue
        spec=loadjson(path)
        push!(control,Dict("n"=>n,"g"=>g,"seconds"=>sum(s["seconds"] for s in spec["sectors"]),
            "largest_residual"=>maximum(s["maximum_ritz_residual"] for s in spec["sectors"]),
            "spectrum_sha256"=>sha256_file(path),"certified"=>false))
    end
    write_compact(joinpath(OUT,"summary.json"),Dict(
        "metadata"=>metadata("work and Pareto summary; estimates non-certifying"),
        "designs"=>design_rows,"original_work"=>baseline,"pareto"=>front,
        "candidate_count"=>length(candidates),"spectrum_controls"=>control,
        "upstream_sha256"=>Dict(n=>sha256_file(joinpath(UPSTREAM,n))
            for n in ("initialization.tex","transfer.tex","thermal72.tex",
                "trials.tex","bootstrap.tex"))))
    for r in design_rows
        println(r["name"]," ",r["caps"]," work=",r["work"]["dimension_sum"],
            " proxy=",r["work"]["state_proxy"]," weighted=",r["work"]["complex_state_proxy"],
            " matvec=",r["work"]["matvec_volume"]," gram=",r["work"]["gram_volume"])
    end
    println("baseline proxy=",baseline["state_proxy"]," matvec=",baseline["matvec_volume"],
        " gram=",baseline["gram_volume"])
    println("Pareto rows=",length(front)," from ",length(candidates)," representatives")
end
main_summary()
