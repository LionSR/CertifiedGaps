using CertifiedGaps.Results, LinearAlgebra
include("su2.jl")
BLAS.set_num_threads(1)
L,S,k=parse.(Int,ARGS)
start=time()
A,B,diagnostics=su2_components(L,S,k)
build_seconds=time()-start
seconds=@elapsed energies=eigvals(Symmetric(cos(-0.3)*A+sin(-0.3)*B))
row=(;L,S,k,dimension=size(A,1),build_seconds,eigensolve_seconds=seconds,diagnostics,
    ground=isempty(energies) ? nothing : minimum(energies),total_seconds=time()-start)
write_result(joinpath(@__DIR__,"results","su2_benchmark_L$(L)_S$(S)_k$(k).json"),(;
    metadata=run_metadata((;L,S,k);source_sha256=source_hashes(@__DIR__;dirs=(),
        files=("su2.jl","su2_benchmark.jl","sectors.jl","model.jl")),
        status="numerical highest-weight construction benchmark, not certified"),row))
println(row)
