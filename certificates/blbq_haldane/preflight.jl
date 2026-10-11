using CertifiedGaps.Results, LinearAlgebra
BLAS.set_num_threads(1)
function charge_count(L,M)
    abs(M)>L && return 0
    coeff=[1]
    for _ in 1:L
        next=zeros(Int,length(coeff)+2)
        for i in eachindex(coeff), j in 0:2; next[i+j]+=coeff[i]; end
        coeff=next
    end
    coeff[M+L+1]
end
function mobius(n)
    sign=1
    for p in 2:n
        n%p==0 || continue
        n÷=p;sign=-sign
        n%p==0 && return 0
    end
    sign
end
function dimension(L,M,k)
    value=0
    for r in 1:L
        L%r==0 && M%r==0 || continue
        ramanujan=sum(d*mobius(r÷d) for d in 1:r if gcd(r,k)%d==0)
        value+=ramanujan*charge_count(L÷r,M÷r)
    end
    @assert value%L==0
    value÷L
end
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
        SU2_basis_implemented=false,
        two_dense_matrices_GiB=16d^2/2.0^30,
        cubic_work_relative_L12=(d/dimension(12,1,1))^3))
end
write_result(joinpath(@__DIR__,"results/preflight.json"),(;
    metadata=run_metadata((;memory_limit_bytes=17179869184,cpu_quota=4);
        source_sha256=source_hashes(@__DIR__;dirs=(),files=("preflight.jl",)),
        status="exact integer dimensions, floating resource projections"),rows))
println(rows)
