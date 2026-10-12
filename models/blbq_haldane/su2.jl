# Numerical highest-weight reduction. No iterative eigenvalue solver is used here.
using LinearAlgebra, SparseArrays, Random
include("sectors.jl")

"Total raising operator S+ from magnetization M to M+1, at fixed momentum."
function raising_operator(L,M,k,o,op,U,Up)
    active=findall(p->mod(k*p,L)==0,o.periods)
    next=findall(p->mod(k*p,L)==0,op.periods)
    lookup=zeros(Int,length(op.reps));lookup[next]=1:length(next)
    rows=Int[];cols=Int[];vals=ComplexF64[]
    for (j,orb) in enumerate(active)
        s=o.reps[orb];p=o.periods[orb]
        for site in 0:L-1
            place=3^site
            s÷place%3==2 && continue
            t=s+place;other=op.owner[t+1];i=lookup[other]
            i==0 && continue
            value=sqrt(2p/op.periods[other])*cis(2pi*k/L*op.shift[t+1])
            push!(rows,i);push!(cols,j);push!(vals,value)
        end
    end
    D=Up'*sparse(rows,cols,vals,length(next),length(active))*U
    @assert norm(imag(D))<1e-9*max(1,norm(D))
    real(D)
end

"Choose root order to control amplification in every contiguous polynomial product."
function projector_order(L,S;trials=2000)
    roots=Float64[j*(j+1)-S*(S+1) for j in S+1:L]
    isempty(roots) && return roots,1.0
    function growth(order)
        worst=1.0
        for x in roots, first in eachindex(order)
            product=1.0
            for last in first:length(order)
                product*=abs(1-x/order[last]);worst=max(worst,product)
            end
        end
        worst
    end
    rng=MersenneTwister(100L+S)
    best=copy(roots);cost=growth(best)
    for _ in 1:trials
        candidate=shuffle(rng,roots);g=growth(candidate)
        if g<cost;best=candidate;cost=g;end
    end
    best,cost
end

"Independent column batches allow sparse multiplication with one BLAS thread."
function column_mul!(Y,A,X)
    p=size(X,2)
    Threads.@threads for chunk in 1:Threads.nthreads()
        lo=fld((chunk-1)*p,Threads.nthreads())+1
        hi=fld(chunk*p,Threads.nthreads())
        lo>hi && continue
        mul!(view(Y,:,lo:hi),A,view(X,:,lo:hi))
    end
    Y
end

"CholeskyQR2: two Gram/triangular normalizations, with independent BLAS-one batches."
function orthonormalize!(V)
    r=size(V,2);G=zeros(r,r);first_diagonal=Float64[]
    for pass in 1:2
        column_mul!(G,transpose(V),V)
        F=cholesky!(Symmetric(G))
        pass==1 && (first_diagonal=abs.(diag(F.U)))
        Threads.@threads for chunk in 1:Threads.nthreads()
            lo=fld((chunk-1)*size(V,1),Threads.nthreads())+1
            hi=fld(chunk*size(V,1),Threads.nthreads())
            lo>hi && continue
            rdiv!(view(V,lo:hi,:),F.U)
        end
    end
    V,first_diagonal
end

"Project onto ker S+ using its known SU(2) singular values, then orthonormalize."
function highest_basis(D,L,S,k;memory_limit=12*2.0^30)
    n=size(D,2);m=size(D,1);r=n-m
    @assert r>=0
    r==0 && return zeros(n,0),(;rank=0,seconds=0.0)
    # Three dense work arrays; QR and its thin Q use at most two such arrays.
    estimate=8.0r*(2n+m)+8.0r^2
    estimate<memory_limit || error("projector memory estimate $(estimate/2.0^30) GiB")
    start=time();rng=MersenneTwister(10000L+100S+k)
    V=randn(rng,n,r);Y=zeros(m,r);W=zeros(n,r)
    roots,amplification=projector_order(L,S)
    projection_seconds=@elapsed for lambda in roots
        column_mul!(Y,D,V);column_mul!(W,transpose(D),Y)
        @. V=V-W/lambda
    end
    Y=nothing;W=nothing;GC.gc()
    qr_seconds=@elapsed begin
        Q,diag_abs=orthonormalize!(V)
        # Numerical rank diagnostic only, not an interval proof.
        @assert minimum(diag_abs)>1e-8*maximum(diag_abs)
    end
    V=nothing;GC.gc()
    residual=norm(D*Q)/max(1,sqrt(r))
    gram=zeros(r,r);column_mul!(gram,transpose(Q),Q)
    orthogonality=norm(gram-I)/max(1,sqrt(r))
    @assert residual<1e-5
    @assert orthogonality<1e-9
    Q,(;rank=r,seconds=time()-start,projection_seconds,qr_seconds,
        residual,orthogonality,amplification,estimated_peak_bytes=estimate,
        min_qr_diagonal=minimum(diag_abs),max_qr_diagonal=maximum(diag_abs),
        orthonormalization="CholeskyQR2")
end

function su2_components(L,S,k;memory_limit=12*2.0^30)
    o=orbits(L,S);A,B,active,U=momentum_components(L,k,o)
    n=length(active)
    if S==L
        return Matrix(A),Matrix(B),(;rank=n,seconds=0.0)
    end
    op=orbits(L,S+1);_,_,next,Up=momentum_components(L,k,op)
    D=raising_operator(L,S,k,o,op,U,Up)
    Q,diagnostics=highest_basis(D,L,S,k;memory_limit)
    r=size(Q,2)
    if r==0
        return zeros(0,0),zeros(0,0),diagnostics
    end
    temp=zeros(n,r)
    Ar=zeros(r,r);Br=zeros(r,r)
    column_mul!(temp,A,Q);column_mul!(Ar,transpose(Q),temp)
    column_mul!(temp,B,Q);column_mul!(Br,transpose(Q),temp)
    @assert norm(Ar-Ar')<1e-9*max(1,norm(Ar))
    @assert norm(Br-Br')<1e-9*max(1,norm(Br))
    Matrix(Symmetric(Ar)),Matrix(Symmetric(Br)),diagnostics
end
