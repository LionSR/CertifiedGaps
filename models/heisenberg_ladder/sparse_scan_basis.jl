# Memory-saving construction for additional numerical scans. Verified blocks continue
# to dispatch to the original interval implementation, without any method replacement.
function block_matrix(o::NamedTuple,k,p,q,J;verified=false)
    verified && return invoke(block_matrix,Tuple{Any,Any,Any,Any,Any},o,k,p,q,J;verified=true)
    sel,idx=block_basis(o,k,p,q); d=length(sel)
    ii=Int[]; jj=Int[]; vv=ComplexF64[]
    phases=[cispi(-e/o.L) for e in 0:2o.L-1]
    for (col,r) in enumerate(sel), (s,h) in transitions(o.reps[r],o.L,J)
        rr=o.rid[s+1]; row=idx[rr]; row==0 && continue
        e=exponent(k,p,q,o.tt[s+1],o.aa[s+1],o.bb[s+1],o.L)
        push!(ii,row); push!(jj,col)
        push!(vv,Float64(h)*sqrt(o.sizes[r]/o.sizes[rr])*phases[e+1])
    end
    H=sparse(ii,jj,vv,d,d)
    @assert isapprox(H,H';atol=1e-10,rtol=1e-12)
    H=sparse(Hermitian(H))
    mod(2k,o.L)==0 ? (Matrix(real(H)),0.0) : (H,0.0)
end

function realify!(H::SparseMatrixCSC,o,k,p,q)
    sel,idx=block_basis(o,k,p,q); n=length(sel); seen=falses(n)
    ii=Int[]; jj=Int[]; vv=ComplexF64[]
    for (i,r) in enumerate(sel)
        seen[i] && continue
        s=reflectrungs(o.reps[r],o.L); j=idx[o.rid[s+1]]
        @assert j!=0
        e=exponent(k,p,q,o.tt[s+1],o.aa[s+1],o.bb[s+1],o.L)
        t=cispi(-e/(2o.L)); seen[i]=true; seen[j]=true
        if i==j
            push!(ii,i); push!(jj,i); push!(vv,t)
        else
            append!(ii,(i,j,i,j)); append!(jj,(i,i,j,j))
            append!(vv,(t/sqrt(2),t/sqrt(2),im*t/sqrt(2),-im*t/sqrt(2)))
        end
    end
    U=sparse(ii,jj,vv,n,n); A=U'*H*U
    @assert maximum(abs,imag(A);init=0.0)<1e-9
    Matrix(real(A))
end
