using ThermalMPOGap, ThermalMPOGap.Results, ThermalMPOGap.Seeds
using ThermalMPOGap.ExactDiagonalization: PartialSpectrum
using LinearAlgebra, Printf
include("sectors.jl")
BLAS.set_num_threads(1)
const ROOT=@__DIR__
const THETAS=[0.0,-0.1,-0.3,-0.5,-0.6,-0.7]

function provenance(config)
    hashes=source_hashes(ROOT;dirs=("test",),files=("model.jl","sectors.jl","scan.jl"))
    merge!(hashes,Dict("../../julia/"*k=>v for (k,v) in source_hashes(joinpath(ROOT,"../../julia");dirs=("src",),files=("Project.toml","Manifest.toml"))))
    run_metadata(config; source_sha256=hashes,packages=(ThermalMPOGap,),status="numerical, not certified")
end

function scan_length(L; maxdim=7000)
    start=time(); energies=[Float64[] for _ in THETAS]; rows=Any[]; dimsum=0
    for M in 0:L
        o=orbits(L,M)
        for k in 0:L÷2
            d=count(p->mod(k*p,L)==0,o.periods)
            mult=(M==0 ? 1 : 2)*((k==0 || 2k==L) ? 1 : 2)
            dimsum+=mult*d
            d==0 && continue
            build=time(); blocks=reduced_components(L,M,k,o)
            build=time()-build
            for block in blocks
            A,B=block.A,block.B
            d=size(A,1)
            d<=maxdim || error("dense block $d exceeds configured maxdim $maxdim")
            vals=Vector{Vector{Float64}}(undef,length(THETAS)); seconds=zeros(length(THETAS))
            Threads.@threads for t in eachindex(THETAS)
                seconds[t]=@elapsed vals[t]=eigvals(Symmetric(Matrix(cos(THETAS[t])*A+sin(THETAS[t])*B)))
            end
            for t in eachindex(THETAS), _ in 1:mult
                append!(energies[t],vals[t])
            end
            push!(rows,(;M,k,symmetry=block.label,dimension=d,multiplicity=mult,build_seconds=build,eigensolve_seconds=seconds))
            @printf("L=%d M=%d k=%d dim=%d times=%s\n",L,M,k,d,string(round.(seconds;digits=2)));flush(stdout)
            end
        end
    end
    @assert dimsum==3^L
    @assert all(length(e)==3^L for e in energies)
    for e in energies; sort!(e); end
    data=(;metadata=provenance((;L,thetas=THETAS,maxdim)),L,dimension_sum=dimsum,
        seconds=time()-start,sectors=rows,spectra=[(;theta=THETAS[t],energies=energies[t]) for t in eachindex(THETAS)])
    write_result(joinpath(ROOT,"results","spectrum_L$(L).json"),data)
    println("completed L=$L in $(time()-start) s");flush(stdout)
end

if abspath(PROGRAM_FILE)==@__FILE__
    exact_rp()
    for L in parse.(Int,ARGS)
        scan_length(L)
    end
end
