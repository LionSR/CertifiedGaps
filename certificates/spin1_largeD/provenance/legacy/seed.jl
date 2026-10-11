include("model.jl")
root=@__DIR__
function enclosed_logz(data,beta)
    @assert data["complete"] && data["rigorous"]
    energies=Float64.(data["energies"]); radii=Float64.(data["radii"])
    shift=energies[1]
    b=ThermalMPOGap.to_interval(beta)
    total=interval(0.0)
    for (e,r) in zip(energies,radii)
        energy=interval(e)-interval(-r,r)
        total+=exp(-b*(energy-interval(shift)))
    end
    -b*interval(shift)+log(total)
end
function make_seed(D,n,beta)
    a=JSON3.read(read(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$n.json"),String),Dict{String,Any})
    b=JSON3.read(read(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$(2n).json"),String),Dict{String,Any})
    finite_gaps=Dict{String,Any}()
    for spec in (a,b)
        radius=maximum(spec["radii"])
        gap=interval(spec["energies"][2])-interval(spec["energies"][1])+interval(-2radius,2radius)
        @assert inf(gap)>0 # The seed rings themselves have verified unique ground states.
        finite_gaps[string(spec["L"])]=[inf(gap),sup(gap)]
    end
    za=enclosed_logz(a,beta); zb=enclosed_logz(b,beta); zbb=enclosed_logz(b,2beta)
    es=Rational{BigInt}(sup(1-exp(zb-2za))); et=Rational{BigInt}(sup(1-exp(zbb-2zb)))
    @assert 0<=es<1 && 0<=et<1
    result=certified_gap_bound(es,et,beta)
    @assert result.bound!==nothing
    @assert iseven(n)
    data=Dict("D"=>string(D),"n0"=>n,"beta"=>string(beta),"eps_s"=>string(es),"eps_t"=>string(et),
        "bound"=>result.bound,"validity"=>"every even L >= $n cells = sites", "certified"=>true,
        "j_admissible"=>result.j_admissible,"bound_upper"=>result.upper,
        "log_rJ_upper"=>result.log_rJ,
        "finite_gap_intervals"=>finite_gaps,
        "logZ"=>[[inf(z),sup(z)] for z in (za,zb,zbb)],"ED_gap_L2n"=>b["energies"][2]-b["energies"][1],
        "source_hashes"=>source_hashes(root;dirs=("test",),files=("model.jl","seed.jl")))
    data["enclosure_hashes"]=file_hashes(root,
        ["results/enclosed_D$(numerator(D))_$(denominator(D))_L$l.json" for l in (n,2n)])
    data["core_hashes"]=source_hashes(joinpath(root,"../../julia");dirs=("src",),
        files=("Project.toml","Manifest.toml"))
    data["metadata"]=run_metadata(Dict("D"=>string(D),"n0"=>n,"beta"=>string(beta));
        source_sha256=data["source_hashes"],packages=(ThermalMPOGap,IntervalArithmetic),
        status="certified uniform gap under the repository two-purity theorem",certified=true)
    data["ratio_to_ED_gap"]=data["bound"]/data["ED_gap_L2n"]
    write_result(joinpath(root,"results","certificate_D$(numerator(D))_$(denominator(D)).json"),data)
    println(data)
    data
end
if abspath(PROGRAM_FILE)==@__FILE__
    make_seed(parse(Rational{Int},ARGS[1]),parse(Int,ARGS[2]),parse(Rational{Int},ARGS[3]))
end
