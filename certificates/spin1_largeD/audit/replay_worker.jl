# Fresh dense replay: this driver never reads historical spectra or radii.
include("../model.jl")
using SHA
const auditroot=@__DIR__
D=parse(Rational{Int},ARGS[1]); group=ARGS[2]
@assert VERSION==v"1.12.4"
@assert group in ("heavy","mixed")
rp_identity(D)
tag="D$(numerator(D))_$(denominator(D))"
out=joinpath(auditroot,"results",tag); mkpath(out)
started=time()
hashes=Dict(p=>bytes2hex(sha256(read(joinpath(auditroot,p)))) for p in ("replay_worker.jl","../model.jl"))
core=readchomp(`git rev-parse HEAD:julia`)
if group=="mixed"
    data=compute_spectrum(6,D;rigorous=true)
    data["source_hashes"]=hashes; data["core_tree"]=core
    data["recomputed_from_model"]=true
    write_result(joinpath(out,"L6.json"),data)
end
charges=group=="heavy" ? vcat([1],collect(4:12)) : [0,2,3]
records=Any[]
for q in charges
    o=orbit_data(12,q)
    for parity in (q==0 ? (1,-1) : (1,)), k in 0:6
        tick=time()
        H,c=block(o,k,parity,D;verified=true); m=size(H,1); m==0 && continue
        mult=(q==0 ? 1 : 2)*(k==0 || k==6 ? 1 : 2)
        E,eta,norminfo=enclose_eigen(H,c,12,D)
        info=Dict{String,Any}(norminfo)
        merge!(info,Dict("D"=>string(D),"L"=>12,"q"=>q,"k"=>k,"parity"=>parity,
            "dimension"=>m,"multiplicity"=>mult,"energies"=>E,"radius"=>eta,
            "seconds"=>time()-tick,"source_hashes"=>hashes,"core_tree"=>core,
            "recomputed_from_model"=>true))
        file="L12_q$(q)_k$(k)_p$(parity).json"
        write_result(joinpath(out,file),info)
        push!(records,Dict("file"=>file,"dimension"=>m,"multiplicity"=>mult,"seconds"=>info["seconds"]))
        println("FRESH $tag q=$q k=$k parity=$parity dim=$m radius=$eta seconds=$(info["seconds"])"); flush(stdout)
        GC.gc()
    end
end
summary=Dict("D"=>string(D),"group"=>group,"seconds"=>time()-started,"records"=>records,
    "source_hashes"=>hashes,"core_tree"=>core,"julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
    "recomputed_from_model"=>true)
write_result(joinpath(out,"worker_$group.json"),summary)
println("DONE $tag $group seconds=$(summary["seconds"])")
