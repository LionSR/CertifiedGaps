# Assemble only fresh block files; historical files enter the final comparison only.
include("../seed.jl")
using SHA, Test
const auditroot=@__DIR__
readjson(path)=JSON3.read(read(path,String),Dict{String,Any})
function fresh_spectrum(D)
    tag="D$(numerator(D))_$(denominator(D))"
    folder=joinpath(auditroot,"results",tag)
    E=Float64[]; radii=Float64[]; records=Any[]; inputs=String[]
    manifest_expected=Dict{String,Tuple{Int,Int,Int}}()
    total=0
    expected_hashes=Dict(p=>bytes2hex(sha256(read(joinpath(auditroot,p))))
        for p in ("replay_worker.jl","../model.jl"))
    core=readchomp(`git rev-parse HEAD:julia`)
    for q in 0:12
        o=orbit_data(12,q)
        for parity in (q==0 ? (1,-1) : (1,)), k in 0:6
        m=length(active_orbits(o,k,parity))
        m==0 && continue
        filename=joinpath(folder,"L12_q$(q)_k$(k)_p$(parity).json")
        data=readjson(filename)
        @assert data["source_hashes"]==expected_hashes && data["core_tree"]==core
        @assert data["recomputed_from_model"] && data["D"]==string(D)
        @assert (data["L"],data["q"],data["k"],data["parity"])==(12,q,k,parity)
        @assert data["dimension"]==length(data["energies"])==m
        mult=(q==0 ? 1 : 2)*(k==0 || k==6 ? 1 : 2)
        @assert mult==data["multiplicity"]
        @assert 0<=data["orthogonality_upper"]<1 && isfinite(data["radius"]) && data["radius"]>0
        @assert isfinite(data["residual_upper"]) && data["residual_upper"]>=0
        @assert isfinite(data["construction_upper"]) && data["construction_upper"]>=0
        @assert issorted(data["energies"]) && all(isfinite,data["energies"])
        for _ in 1:mult
            append!(E,data["energies"]); append!(radii,fill(data["radius"],m))
        end
        total+=m*mult
        manifest_expected[basename(filename)]=(m,mult,q)
        push!(inputs,relpath(filename,auditroot))
        push!(records,filter(p->first(p)!="energies",data))
        end
    end
    @assert length(E)==length(radii)==total==3^12
    @assert length(records)==92
    perm=sortperm(E)
    a=readjson(joinpath(folder,"L6.json"))
    @assert a["source_hashes"]==expected_hashes && a["core_tree"]==core
    @assert a["recomputed_from_model"] && a["D"]==string(D) && a["L"]==6
    @assert a["complete"] && a["rigorous"]
    @assert a["dimension"]==length(a["energies"])==length(a["radii"])==3^6
    @assert issorted(a["energies"]) && all(isfinite,a["energies"])
    @assert all(r->isfinite(r)&&r>0,a["radii"])
    expected6=Dict{Tuple{Int,Int,Int},Tuple{Int,Int}}()
    for q in 0:6
        o=orbit_data(6,q)
        for parity in (q==0 ? (1,-1) : (1,)), k in 0:3
            m=length(active_orbits(o,k,parity)); m==0 && continue
            mult=(q==0 ? 1 : 2)*(k==0 || k==3 ? 1 : 2)
            expected6[(q,k,parity)]=(m,mult)
        end
    end
    @assert length(a["sectors"])==length(expected6)==29
    seen6=Set{Tuple{Int,Int,Int}}()
    for s in a["sectors"]
        key=(s["q"],s["k"],s["parity"])
        @assert !(key in seen6); push!(seen6,key)
        @assert (s["dimension"],s["multiplicity"])==expected6[key]
        @assert s["levels"]==s["dimension"] && 0<=s["orthogonality_upper"]<1
    end
    @assert sum(s["dimension"]*s["multiplicity"] for s in a["sectors"])==3^6
    push!(inputs,relpath(joinpath(folder,"L6.json"),auditroot))
    workers=[readjson(joinpath(folder,"worker_$group.json")) for group in ("heavy","mixed")]
    listed=String[]
    for (group,worker) in zip(("heavy","mixed"),workers)
        @assert worker["source_hashes"]==expected_hashes && worker["core_tree"]==core
        @assert worker["recomputed_from_model"] && worker["blas_threads"]==1
        @assert worker["D"]==string(D) && worker["group"]==group
        @assert worker["julia_version"]=="1.12.4"
        @assert worker["seconds"]<45*60
        expected_charges=group=="heavy" ? Set(vcat([1],collect(4:12))) : Set([0,2,3])
        for r in worker["records"]
            m,mult,q=manifest_expected[r["file"]]
            @assert (r["dimension"],r["multiplicity"])==(m,mult)
            @assert q in expected_charges && isfinite(r["seconds"]) && r["seconds"]>=0
        end
        append!(listed,[r["file"] for r in worker["records"]])
        push!(inputs,relpath(joinpath(folder,"worker_$group.json"),auditroot))
    end
    expected_files=[basename(p) for p in inputs if startswith(basename(p),"L12_")]
    @assert length(listed)==length(unique(listed))==92
    @assert Set(listed)==Set(expected_files)
    b=Dict("L"=>12,"energies"=>E[perm],"radii"=>radii[perm],"complete"=>true,
        "rigorous"=>true,"dimension"=>total,"sectors"=>records)
    a,b,inputs,workers
end
function assemble(D,beta,published;save=true)
    @assert VERSION==v"1.12.4" && pkgversion(IntervalArithmetic)==v"1.0.12"
    @assert beta>0 && published>0
    # The tree name alone does not detect uncommitted changes to loaded core files.
    repository=normpath(joinpath(auditroot,"../../.."))
    @assert success(`git -C $repository diff --quiet HEAD -- julia`)
    current_core=readchomp(`git -C $repository rev-parse HEAD:julia`)
    audited_core=readchomp(`git -C $repository rev-parse 22b2759e0790197b83e2571b465aa8dc69ebc0be:julia`)
    @assert current_core==audited_core
    a,b,inputs,workers=fresh_spectrum(D)
    za=enclosed_logz(a,beta); zb=enclosed_logz(b,beta); zbb=enclosed_logz(b,2beta)
    es=Rational{BigInt}(sup(1-exp(zb-2za)))
    et=Rational{BigInt}(sup(1-exp(zbb-2zb)))
    @assert 0<=es<1 && 0<=et<1
    result=certified_gap_bound(es,et,beta)
    @assert result.bound!==nothing && Rational{BigInt}(result.bound)>=published
    gaps=Dict{String,Any}()
    comparison=Dict{String,Any}()
    tag="D$(numerator(D))_$(denominator(D))"
    for spec in (a,b)
        rad=maximum(spec["radii"])
        gap=interval(spec["energies"][2])-interval(spec["energies"][1])+interval(-2rad,2rad)
        @assert inf(gap)>0
        gaps[string(spec["L"])]=[inf(gap),sup(gap)]
        # Historical values are deliberately read only after fresh certification.
        old=readjson(joinpath(auditroot,"../results/enclosed_$(tag)_L$(spec["L"]).json"))
        centers=Float64.(spec["energies"]); oldcenters=Float64.(old["energies"])
        @assert all(abs.(centers.-oldcenters).<=Float64.(spec["radii"]).+Float64.(old["radii"]))
        comparison[string(spec["L"])]=Dict("centers_bitwise_equal"=>centers==oldcenters,
            "radii_bitwise_equal"=>spec["radii"]==old["radii"],
            "maximum_center_difference"=>maximum(abs.(centers.-oldcenters)),
            "fresh_maximum_radius"=>maximum(spec["radii"]))
    end
    oldcert=readjson(joinpath(auditroot,"../results/certificate_$tag.json"))
    data=Dict("D"=>string(D),"n0"=>6,"beta"=>string(beta),"eps_s"=>string(es),"eps_t"=>string(et),
        "bound"=>result.bound,"published_bound_rational"=>string(published),
        "j_admissible"=>result.j_admissible,"validity"=>"every even L >= 6 one-site cells = sites",
        "certified"=>true,"recomputed_from_model"=>true,"logZ"=>[[inf(z),sup(z)] for z in (za,zb,zbb)],
        "finite_gap_intervals"=>gaps,"sector_count_L12"=>length(b["sectors"]),
        "dimension_L6"=>a["dimension"],"dimension_L12"=>b["dimension"],
        "comparison"=>comparison,"original_bound"=>oldcert["bound"],
        "original_rational_seeds_equal"=>string(es)==oldcert["eps_s"] && string(et)==oldcert["eps_t"],
        "worker_seconds"=>[w["seconds"] for w in workers],
        "input_hashes"=>file_hashes(auditroot,inputs),
        "source_hashes"=>file_hashes(auditroot,["replay_worker.jl","assemble_replay.jl","../model.jl","../seed.jl"]),
        "core_tree"=>readchomp(`git rev-parse HEAD:julia`),
        "core_hashes"=>source_hashes(joinpath(repository,"julia");dirs=("src",),
            files=("Project.toml","Manifest.toml")),
        "core_main_revision"=>"22b2759e0790197b83e2571b465aa8dc69ebc0be")
    data["metadata"]=run_metadata(Dict("D"=>string(D),"n0"=>6,"beta"=>string(beta));
        source_sha256=data["source_hashes"],packages=(CertifiedGaps,IntervalArithmetic),
        status="fresh full dense spectral replay under the repository theorem",certified=true)
    save && write_result(joinpath(auditroot,"results/certificate_$tag.json"),data)
    println("FRESH CERTIFICATE $tag bound=$(result.bound) original=$(oldcert["bound"]) seeds_equal=$(data["original_rational_seeds_equal"]) worker_seconds=$(data["worker_seconds"])")
    return data
end
if abspath(PROGRAM_FILE)==@__FILE__
    assemble(4//1,9//4,694694192//1000000000)
    assemble(1603//500,381//125,2213456//1000000000)
end
