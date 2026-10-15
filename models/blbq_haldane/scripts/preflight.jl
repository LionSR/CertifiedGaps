# Exact integer sector dimensions (characters) and dense-memory projections for L = 10..16;
# writes results/preflight.json. The dimensions are exact; the GiB figures are projections.
using CertifiedGaps.Results, LinearAlgebra
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"src","dimensions.jl"))
rows=Any[]
for L in (10,12,14,16)
    dims=[dimension(L,M,k) for M in -L:L,k in 0:L-1]
    @assert sum(dims)==3^L
    highest=[dimension(L,S,k)-dimension(L,S+1,k) for S in 0:L,k in 0:L-1]
    @assert minimum(highest)>=0
    @assert sum((2S+1)*highest[S+1,k+1] for S in 0:L,k in 0:L-1)==3^L
    d=dimension(L,1,1) # neither flip nor reflection acts internally here
    push!(rows,(;L,total_dimension=3^L,dimension_sum=sum(dims),
        magnetization_zero=charge_count(L,0),momentum_zero_M0=dimension(L,0,0),
        M1_k1_dimension=d,single_dense_matrix_GiB=8d^2/2.0^30,
        SU2_highest_weight_max_dimension=maximum(highest),
        SU2_single_matrix_GiB=8maximum(highest)^2/2.0^30,
        SU2_basis_implemented=true,
        two_dense_matrices_GiB=16d^2/2.0^30,
        cubic_work_relative_L12=(d/dimension(12,1,1))^3))
end
write_result(joinpath(ROOT,"results/preflight.json"),(;
    metadata=run_metadata((;memory_limit_bytes=17179869184,cpu_quota=4);
        source_sha256=source_hashes(ROOT;dirs=(),
            files=("scripts/preflight.jl","src/dimensions.jl")),
        status="exact integer dimensions, floating resource projections"),rows))
println(rows)
