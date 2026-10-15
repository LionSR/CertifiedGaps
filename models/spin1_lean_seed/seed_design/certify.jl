include(joinpath(@__DIR__, "src", "exact.jl"))

function audit_stored(source,eps)
    @assert eps >= (source==1 ? frac(3,100000) : frac(1,100000))
    verification=read_exact_json(joinpath(SPIN1_RESULTS,"verification.json"))
    @assert verification["status"]=="passed"
    paths=["source$source","trials"]
    hashes=Dict{String,String}()
    for name in paths
        path=joinpath(SPIN1_RESULTS,"$name.json")
        h=sha256_file(path)
        @assert verification["components"][name]["sha256"]==h
        hashes[name]=h
    end
    raw=read_exact_json(joinpath(SPIN1_RESULTS,"source$source.json"))
    for row in raw["certificates"]
        c=parse_rational(row["center"])
        lo=parse_rational(row[source==1 ? "lower" : "true_trace_lower"])
        hi=parse_rational(row[source==1 ? "upper" : "true_trace_upper"])
        @assert c-eps<lo<=hi<c+eps
    end
    hashes
end

function proposals_for(source,eps,mode,nmax;odd=false)
    table=stored_table(source)
    input=prepare(table,nmax,Float64(eps),mode;allowodd=odd,beta=source==1 ? 12.25 : 10.5)
    [Dict("q"=>f.q,"shift"=>f.shift) for f in input.filters]
end

function main_certify()
    rows=Any[]
    for (name,source,eps,mode,nmax,trial,k,r) in (
        ("small_simple",2,frac(3,10000),"ION",10,120,30,frac(5,2)),
        ("small_relaxed",2,frac(4,10000),"ION",10,144,30,frac(25,8)),
        ("real_even",1,frac(1,1000),"JN",12,120,30,frac(2)),
    )
        started=time()
        inp=exact_input(source,mode,nmax,eps,proposals_for(source,eps,mode,nmax);
            uncappedO=mode=="ION")
        trial_record=verify_trial(trial)
        seed=exact_seed(inp,trial,k,r)
        beta=parse_rational(inp["base_beta"])*r
        chain=exact_chain(parse_rational(seed["p"]),parse_rational(seed["q"]),beta)
        ratio=frac(784)/beta
        @assert parse_rational(chain["r"])^numerator(ratio)>frac(20)^denominator(ratio)
        chain["dominates_original_target"]=true
        chain["target_power_ratio"]=string(ratio)
        push!(rows,Dict("name"=>name,"input"=>inp,"seed"=>seed,"chain"=>chain,
            "trial"=>trial_record,"input_sha256"=>audit_stored(source,eps),
            "seconds"=>time()-started))
        println(name,": p=",Float64(parse_rational(seed["p"])),
            " q=",Float64(parse_rational(seed["q"])),
            " gap > ",chain["gap_display"]," every even L >= ",trial÷2," seconds=",time()-started)
        flush(stdout)
    end
    write_result(joinpath(OUT,"certificate.json"),Dict(
        "metadata"=>metadata("three exact scalar designs";
        status="exact rational scalar checks passed; " *
            "physical premises and analytic lemmas not Lean-checked"),
        "designs"=>rows))
    compact=[Dict("name"=>row["name"],"input"=>row["input"],
        "trial"=>row["seed"]["trial"],"k"=>row["seed"]["k"],
        "temperature_multiplier"=>row["seed"]["r"],"D"=>row["seed"]["D"],
        "E"=>row["seed"]["E"],"Q"=>row["seed"]["Q"],
        "dominant_lower_steps"=>vcat([row["seed"]["initial_t"]],
            [r["new"] for r in row["seed"]["refinements"]]),
        "chain_rate"=>row["chain"]["r"],"chain_steps"=>length(row["chain"]["s"])-1,
        "chain_grid"=>"10^512") for row in rows]
    write_result(joinpath(OUT,"compact_witnesses.json"),Dict(
        "scope"=>"compact replay recipes; not standalone physical or Lean certificates",
        "designs"=>compact))
end
if abspath(PROGRAM_FILE) == @__FILE__
    main_certify()
end
