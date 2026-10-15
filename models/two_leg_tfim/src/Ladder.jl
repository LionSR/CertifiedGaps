module Ladder
using LinearAlgebra, SparseArrays, CertifiedGaps
using CertifiedGaps.Results
export rp_identity, blocks, block_matrix, full_spectrum, verify_eigen, logz_interval
BLAS.set_num_threads(1)

"Exact rational check of −h_bond = Σ C⊗C (C = Z⊗I, I⊗Z) and of a real symmetric on-site term."
function rp_identity(g::Rational)
    Z = Rational{BigInt}[1 0; 0 -1]; X = Rational{BigInt}[0 1; 1 0]
    I2 = Matrix{Rational{BigInt}}(I,2,2)
    C = [kron(Z,I2), kron(I2,Z)]
    # Independent four-spin construction, order (a,b,a′,b′).
    h = -kron(Z,I2,Z,I2)-kron(I2,Z,I2,Z)
    o = -kron(Z,Z)-g*(kron(X,I2)+kron(I2,X))
    @assert -h == sum(kron(c,c) for c in C)
    @assert all(isreal,o) && issymmetric(o)
    return true
end

function permute_bits(s,L,leg,half,refl)
    t=0
    for j in 0:L-1, a in 0:1
        k = refl ? mod(-j,L) : j
        half && (k=mod(k+L÷2,L))
        b = leg ? 1-a : a
        t |= ((s >> (2j+a)) & 1) << (2k+b)
    end
    t
end

"Complete real character blocks of parity, leg exchange, reflection, and half translation."
function blocks(L)
    even=iseven(L); ng=even ? 8 : 4; N=4^L
    # Group elements are commuting involutions; half translation is used only for even L.
    actions=[(isodd(a),even && isodd(a>>1),isodd(a>>(even ? 2 : 1))) for a in 0:ng-1]
    orbits=Vector{Vector{Int}}(); elements=Vector{Vector{Int}}(); seen=falses(N)
    for s in 0:N-1
        seen[s+1] && continue
        ys=[permute_bits(s,L,v...) for v in actions]
        seen[ys .+ 1].=true
        push!(orbits,sort(unique(ys))); push!(elements,ys)
    end
    out=[]
    for parity in 0:1, ch in 0:ng-1
        reps=Int[]; sizes=Int[]; lookup=zeros(Int,N); signs=zeros(Int8,N)
        for (orb,ys) in zip(orbits,elements)
            s=first(orb)
            count_ones(s)%2 == parity || continue
            cs=[isodd(count_ones(ch&a)) ? -1 : 1 for a in 0:ng-1]
            any(ys[k]==s && cs[k]!=1 for k in eachindex(ys)) && continue
            push!(reps,s); push!(sizes,length(orb)); j=length(reps)
            for k in eachindex(ys)
                lookup[ys[k]+1]=j; signs[ys[k]+1]=cs[k]
            end
        end
        push!(out,(;L,parity,ch,reps,sizes,lookup,signs))
    end
    @assert sum(length(b.reps) for b in out)==N
    out
end

"Integer coefficient matrix: exact block entry A[i,j]/2 * sqrt(size[j]/size[i])."
function block_matrix(b,g)
    L=b.L; m=length(b.reps); ii=Int[]; jj=Int[]; vv=Int[]
    masks=Int[]
    for j in 0:L-1
        push!(masks,(1<<(2j)) ⊻ (1<<(2j+1)))
        for a in 0:1
            push!(masks,(1<<(2j+a)) ⊻ (1<<(2mod(j+1,L)+a)))
        end
    end
    for (j,s) in enumerate(b.reps)
        push!(ii,j);push!(jj,j);push!(vv,Int(-2g*(2L-2count_ones(s))))
        for mask in masks
            t=s ⊻ mask; i=b.lookup[t+1]; i==0 && continue
            push!(ii,i); push!(jj,j); push!(vv,-2Int(b.signs[t+1]))
        end
    end
    A=sparse(ii,jj,vv,m,m); dropzeros!(A)
    I0,J0,V0=findnz(A)
    for (i,j,v) in zip(I0,J0,V0)
        @assert v*b.sizes[j] == A[j,i]*b.sizes[i]
    end
    H=zeros(m,m)
    for (i,j,v) in zip(I0,J0,V0)
        i<=j || continue
        H[i,j]=H[j,i]=(v/2)*sqrt(b.sizes[j]/b.sizes[i])
    end
    H,A
end

up(f)=setrounding(BigFloat,RoundUp) do; f(); end
down(f)=setrounding(BigFloat,RoundDown) do; f(); end

"Polar/Weyl eigenvalue enclosure with rigorous norm upper bounds on BLAS errors."
function verify_eigen(H,A,b,F)
    m=size(H,1); Q=F.vectors; lam=F.values
    # Exact algebraic matrix construction error: interval sqrt, then signed multiplication.
    rows=zeros(BigFloat,m)
    I0,J0,V0=findnz(A)
    for (i,j,v) in zip(I0,J0,V0)
        lo=down(()->sqrt(BigFloat(b.sizes[j])/BigFloat(b.sizes[i])))
        hi=up(()->sqrt(BigFloat(b.sizes[j])/BigFloat(b.sizes[i])))
        a=BigFloat(v)/2
        l=down(()->a*(v>=0 ? lo : hi)); u=up(()->a*(v>=0 ? hi : lo))
        err=up(()->max(abs(BigFloat(H[i,j])-l),abs(u-BigFloat(H[i,j]))))
        rows[i]=up(()->rows[i]+err)
    end
    dH=maximum(rows;init=BigFloat(0))
    R=H*Q-Q*Diagonal(lam); E=Q'Q-I
    # Frobenius norm of stored residual is bounded by m*maxabs. Products' error
    # Frobenius norm is gamma*(||H||F||Q||F+||Q||F||Lambda||2).
    # All scalar norms/bounds below are accumulated with MPFR upward rounding.
    f_upper(M)=up(()->sqrt(sum(x->BigFloat(x)^2,M)))
    qf=f_upper(Q); hf=f_upper(H); lf=BigFloat(maximum(abs,lam))
    gamma=up(()->BigFloat(2m+8)*big(2.0)^(-53)/(1-BigFloat(2m+8)*big(2.0)^(-53)))
    # Additive underflow allowance, valid also for subnormal intermediate results.
    tiny=up(()->BigFloat(16)*m^3*BigFloat(2.0)^(-1074))
    r=up(()->BigFloat(m)*BigFloat(maximum(abs,R))+gamma*(hf*qf+qf*lf)+tiny)
    e=up(()->BigFloat(m)*BigFloat(maximum(abs,E))+gamma*(qf*qf+sqrt(BigFloat(m)))+tiny)
    @assert e<1
    # U=Q(Q'Q)^(-1/2); ||U-Q|| <= e/(1+sqrt(1-e)).
    den=down(()->1+sqrt(1-e)); a=up(()->e/den)
    eta=up(()->dH+r+(hf+lf)*a)
    return (;radius=string(Rational{BigInt}(eta)), residual=string(r),orthogonality=string(e),
        construction=string(dH), dimension=m)
end

"""
All eigenvalues of the L-rung ladder at field g, block by block; with `verified=true` each
block also carries its eigenvalues and a whole-spectrum enclosure radius (`verify_eigen`).
"""
function full_spectrum(L,g; verified=false)
    es=Float64[]; sectors=[]; bs=blocks(L); t=time()
    for b in bs
        m=length(b.reps); m==0 && continue
        st=time(); H,A=block_matrix(b,g)
        if verified
            F=eigen(Symmetric(H)); ev=F.values
            check=verify_eigen(H,A,b,F)
        else
            ev=eigvals(Symmetric(H)); check=nothing
        end
        append!(es,ev)
        push!(sectors,(;parity=b.parity,character=b.ch,dimension=m,seconds=time()-st,
            verification=check, eigenvalues=verified ? ev : nothing))
        println("L=$L g=$g parity=$(b.parity) character=$(b.ch) dim=$m seconds=$(round(time()-st,digits=2))")
        flush(stdout)
    end
    @assert length(es)==4^L
    (;L,g=string(g),energies=sort!(es),sectors,seconds=time()-t,dimension=length(es))
end

"Directed MPFR interval for log Z, using all eigenvalue enclosures."
function logz_interval(spec,beta)
    beta=BigFloat(beta); lows=BigFloat[]; highs=BigFloat[]
    for sec in spec.sectors
        rad=parse(Rational{BigInt},sec.verification.radius)
        rr=up(()->BigFloat(numerator(rad))/BigFloat(denominator(rad)))
        for e in sec.eigenvalues
            push!(lows,down(()->BigFloat(e)-rr));push!(highs,up(()->BigFloat(e)+rr))
        end
    end
    lo=down(()->log(sum(e->exp(-beta*e),highs)))
    hi=up(()->log(sum(e->exp(-beta*e),lows)))
    (lo,hi)
end
end
