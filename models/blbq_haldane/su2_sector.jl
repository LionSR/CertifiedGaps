using CertifiedGaps, CertifiedGaps.Results, LinearAlgebra, Serialization
include("su2.jl")
BLAS.set_num_threads(1)
const THETAS=[0.0,-0.1,-0.3,-0.5,-0.6,-0.7]
L,S,k=parse.(Int,ARGS)
started=time()
hashes=source_hashes(@__DIR__;dirs=(),
    files=("su2.jl","su2_sector.jl","sectors.jl","model.jl"))
output=joinpath(@__DIR__,"results/su2_L$(L)_S$(S)_k$(k).json")
cache=joinpath(@__DIR__,"runtime/su2_L$(L)_S$(S)_k$(k).bin")
if isfile(output)
    println("Already complete: $output")
    exit()
end
if isfile(cache)
    saved=deserialize(cache)
    @assert saved.source_sha256==hashes "cached matrices use different source"
    A,B,diagnostics=saved.A,saved.B,saved.diagnostics
    build_seconds=saved.build_seconds
else
    println("Constructing SU2 L=$L S=$S k=$k");flush(stdout)
    build_seconds=@elapsed A,B,diagnostics=su2_components(L,S,k)
    serialize(cache,(;A,B,diagnostics,build_seconds,source_sha256=hashes))
end
GC.gc()
println("Constructed $(size(A,1))-state block in $build_seconds seconds");flush(stdout)
levels=Vector{Vector{Float64}}(undef,6);times=zeros(6)
Threads.@threads for t in eachindex(THETAS)
    times[t]=@elapsed levels[t]=eigvals(Symmetric(cos(THETAS[t])*A+sin(THETAS[t])*B))
end
mult=(2S+1)*((k==0 || 2k==L) ? 1 : 2)
write_result(output,(;metadata=run_metadata((;L,S,k,theta=THETAS);
    source_sha256=hashes,packages=(CertifiedGaps,),
    status="complete dense SU2 sector; numerical, not certified"),
    L,S,k,dimension=size(A,1),multiplicity=mult,build_seconds,diagnostics,
    eigensolve_seconds=times,seconds=time()-started,
    spectra=[(;theta=THETAS[t],energies=levels[t]) for t in eachindex(THETAS)]))
println("Saved $output in $(time()-started) seconds")
rm(cache;force=true)
