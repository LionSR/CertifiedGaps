# Retry only the killed J=7/2 large ring, retaining every completed dense block.
# Run with --threads=2, BLAS threads remain one. No certified computation uses this driver.
include("ladder.jl")
include("sparse_scan_basis.jl")
using JSON3
function retry_large()
    L=10; J=7//2; start=time()
    meta=provenance(Dict("J"=>string(J),"L"=>L,"mode"=>"threaded block checkpoints"))
    parts=Vector{Any}(undef,L+1)
    Threads.@threads :greedy for m in 0:L
        o=orbit_data(L,m); levels=Float64[]; rows=Any[]
        for k in 0:fld(L,2), p in 0:1, q in 0:Int(o.flip)
            mult=(m==L ? 1 : 2)*(mod(2k,L)==0 ? 1 : 2)
            path=joinpath(@__DIR__,"results","checkpoints_J3p5","m$(m)_k$(k)_p$(p)_q$(q).json")
            if isfile(path)
                saved=JSON3.read(read(path,String)); ev=Float64.(saved.energies)
                row=saved.sector
            else
                tick=time(); H,_=block_matrix(o,k,p,q,J)
                d=size(H,1); d==0 && continue
                if eltype(H)==ComplexF64; H=realify!(H,o,k,p,q); end
                ev=eigvals(Hermitian(H))
                row=Dict("m"=>m,"k"=>k,"leg_parity"=>(-1)^p,
                    "spin_parity"=>(o.flip ? (-1)^q : 0),"dimension"=>d,
                    "multiplicity"=>mult,"seconds"=>time()-tick)
                write_result(path,Dict("metadata"=>meta,"sector"=>row,"energies"=>ev))
            end
            @assert length(ev)==length(first(block_basis(o,k,p,q)))
            push!(rows,row)
            for _ in 1:mult; append!(levels,ev); end
        end
        parts[m+1]=(;levels,rows)
        @printf("Checkpointed J=3.5 L=10 m=%d elapsed=%.1fs\n",m,time()-start)
        flush(stdout)
    end
    es=sort!(reduce(vcat,[x.levels for x in parts])); all_rows=reduce(vcat,[x.rows for x in parts])
    @assert length(es)==sum(r["dimension"]*r["multiplicity"] for r in all_rows)==4^L
    # Keep elapsed ED time from the completed checkpoint-producing run on replay.
    oldlog=joinpath(@__DIR__,"results","retry_J3p5_compute.log")
    mt=isfile(oldlog) ? match(r"m=10 elapsed=([0-9.]+)s",read(oldlog,String)) : nothing
    prior_seconds=mt===nothing ? 0.0 : parse(Float64,mt.captures[1])
    write_result(joinpath(@__DIR__,"results","spectrum_J3p5_L10.json"),Dict(
        "metadata"=>meta,"J"=>string(J),"L"=>L,"dimension"=>4^L,"sectors"=>all_rows,
        "seconds"=>prior_seconds+time()-start,"checkpoint_compute_seconds"=>prior_seconds,
        "reassembly_seconds"=>time()-start,"energies"=>es,"gap"=>es[2]-es[1]))
end
retry_large()
empty!(ARGS); push!(ARGS,"7//2")
include_string(Main,replace(read(joinpath(@__DIR__,"scan.jl"),String),
    "include(\"ladder.jl\")"=>""),joinpath(@__DIR__,"scan.jl"))
