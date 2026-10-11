# Independent charge groups run concurrently; each BLAS operation still uses one thread.
include("model.jl")
root=@__DIR__
D=parse(Rational{Int},ARGS[1]); L=parse(Int,ARGS[2])
@assert VERSION==v"1.12.4"
rp_identity(D)
started=time()
function certify_charge(q,L,D)
    levels=Float64[]; radii=Float64[]; records=Any[]; total=0
    o=orbit_data(L,q)
    for parity in (q==0 ? (1,-1) : (1,)), k in 0:div(L,2)
        tick=time()
        H,c=block(o,k,parity,D;verified=true); m=size(H,1); m==0 && continue
        mult=(q==0 ? 1 : 2)*(k==0 || 2k==L ? 1 : 2)
        E,eta,info=enclose_eigen(H,c,L,D)
        merge!(info,Dict("q"=>q,"k"=>k,"parity"=>parity,"dimension"=>m,
            "multiplicity"=>mult,"levels"=>length(E),"seconds"=>time()-tick))
        push!(records,info); total+=mult*m
        for _ in 1:mult
            append!(levels,E); append!(radii,fill(eta,length(E)))
        end
        println("verified D=$D L=$L q=$q k=$k parity=$parity dim=$m eta=$eta seconds=$(time()-tick)")
        flush(stdout)
    end
    (;levels,radii,records,total)
end
answers=Vector{Any}(undef,L+1)
Threads.@threads :greedy for q in 0:L
    answers[q+1]=certify_charge(q,L,D)
end
levels=vcat([a.levels for a in answers]...); radii=vcat([a.radii for a in answers]...)
records=vcat([a.records for a in answers]...); total=sum(a.total for a in answers)
@assert total==length(levels)==3^L
perm=sortperm(levels)
data=Dict("D"=>string(D),"L"=>L,"energies"=>levels[perm],"radii"=>radii[perm],
    "dimension"=>total,"complete"=>true,"rigorous"=>true,"sectors"=>records,
    "seconds"=>time()-started,
    "source_hashes"=>source_hashes(root;dirs=("test",),files=("model.jl","certify_parallel.jl")))
data["metadata"]=run_metadata(Dict("D"=>string(D),"L"=>L);
    source_sha256=data["source_hashes"],packages=(ThermalMPOGap,IntervalArithmetic),
    status="complete dense eigenvalue enclosures",certified=true)
write_result(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$L.json"),data)
println("DONE D=$D L=$L seconds=$(data["seconds"])")
