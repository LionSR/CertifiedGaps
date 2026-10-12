using CertifiedGaps,CertifiedGaps.Results,JSON3,LinearAlgebra
BLAS.set_num_threads(1)
include("dimensions.jl")
const ROOT=@__DIR__;const THETAS=[0.0,-0.1,-0.3,-0.5,-0.6,-0.7]
L=parse(Int,only(ARGS));files=String[];rows=Any[]
levels=[Float64[] for _ in 1:6];weights=[Int[] for _ in 1:6]
for S in 0:L,k in 0:L÷2
    name="su2_L$(L)_S$(S)_k$(k).json";path=joinpath(ROOT,"results",name)
    @assert isfile(path) "Full partition function unavailable: missing $name"
    d=JSON3.read(read(path,String));push!(files,name)
    @assert (d.L,d.S,d.k)==(L,S,k)
    @assert d.dimension==dimension(L,S,k)-dimension(L,S+1,k)
    expected=(2S+1)*((k==0 || 2k==L) ? 1 : 2)
    @assert d.multiplicity==expected
    for t in 1:6
        @assert d.spectra[t].theta==THETAS[t]
        @assert length(d.spectra[t].energies)==d.dimension
        append!(levels[t],d.spectra[t].energies)
        append!(weights[t],fill(expected,d.dimension))
    end
    push!(rows,(;S,k,dimension=d.dimension,multiplicity=d.multiplicity,
        seconds=d.seconds,build_seconds=d.build_seconds,diagnostics=d.diagnostics))
end
@assert sum(r.dimension*r.multiplicity for r in rows)==3^L
for t in 1:6
    p=sortperm(levels[t]);levels[t]=levels[t][p];weights[t]=weights[t][p]
    @assert sum(weights[t])==3^L
end
write_result(joinpath(ROOT,"results/spectrum_L$(L).json"),(;
    metadata=run_metadata((;L,method="all SU2 highest-weight/momentum sectors, dense solves");
        source_sha256=source_hashes(ROOT;dirs=(),files=("collect_su2.jl","dimensions.jl")),
        packages=(CertifiedGaps,),status="complete weighted spectrum, numerical not certified"),
    input_hashes=file_hashes(joinpath(ROOT,"results"),files),
    L,dimension_sum=3^L,seconds=sum(r.seconds for r in rows),sectors=rows,
    spectra=[(;theta=THETAS[t],energies=levels[t],weights=weights[t]) for t in 1:6]))
println("Collected all sectors for L=$L; exact dimension sum $(3^L)")
