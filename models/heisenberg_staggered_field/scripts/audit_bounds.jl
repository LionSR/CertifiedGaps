# Independent interval accumulation checks of the residual/orthogonality bound.
# This audits small dense blocks, not a new parameter or seed search. Writes
# results/audit_bounds.json.
#     julia --project=models/heisenberg_staggered_field \
#         models/heisenberg_staggered_field/scripts/audit_bounds.jl              (≈ 25 s)
include(joinpath(dirname(@__DIR__), "src", "common.jl"))
using Test
rows=Any[]
@testset "spin-half normalization and independently accumulated residuals" begin
    X=[0//1 1//2;1//2 0//1]; J=[0//1 1//2;-1//2 0//1]
    Z=[1//2 0//1;0//1 -1//2]; H2=kron(X,X)-kron(J,J)+kron(Z,Z)
    @test X^2-J^2+Z^2 == 3//4*Matrix{Rational{Int}}(I,2,2)
    @test (H2-1//4*I)*(H2+3//4*I)==zeros(Rational{Int},4,4)
    @test tr(H2)==0
    for L in (4,6,8), k in 0:(L÷2)÷2
        o=orbits(L,L÷2); A,a=block(L,7//20,k,o;verified=true)
        F=eigen(Hermitian(Matrix(A))); b=enclose(A,F,a,L,7//20)
        Q=F.vectors; lambda=F.values; n=length(lambda)
        ci(z)=complex(iv(real(z)),iv(imag(z)))
        r2=iv(0); e2=iv(0)
        # Every product and sum below is an interval operation. No BLAS residual is reused.
        for i in 1:n,j in 1:n
            r=complex(iv(0),iv(0));e=complex(iv(0),iv(0))
            for t in 1:n
                r+=ci(A[i,t])*ci(Q[t,j])
                e+=conj(ci(Q[t,i]))*ci(Q[t,j])
            end
            r-=ci(Q[i,j])*iv(lambda[j])
            e-=complex(iv(i==j ? 1 : 0),iv(0))
            r2+=real(r)^2+imag(r)^2
            e2+=real(e)^2+imag(e)^2
        end
        ru=sup(sqrt(r2));eu=sup(sqrt(e2))
        @test ru<=b.residual_upper
        @test eu<=b.orthogonality_upper
        push!(rows,Dict("L"=>L,"k"=>k,"dimension"=>n,
            "direct_interval_residual_upper"=>ru,"gamma_residual_upper"=>b.residual_upper,
            "direct_interval_orthogonality_upper"=>eu,
            "gamma_orthogonality_upper"=>b.orthogonality_upper))
    end
end
write_result(joinpath(ROOT,"results","audit_bounds.json"),Dict(
    "description"=>"Independent scalar interval accumulation on small dense blocks; no search",
    "source_hashes"=>file_hashes(ROOT,["src/common.jl","scripts/audit_bounds.jl"]),
    "julia"=>string(VERSION),"rows"=>rows,"passed"=>true))
