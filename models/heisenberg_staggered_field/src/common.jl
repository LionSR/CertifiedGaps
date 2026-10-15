using CertifiedGaps, LinearAlgebra, SparseArrays, IntervalArithmetic, JSON3
using CertifiedGaps.Results, CertifiedGaps.Seeds, CertifiedGaps.ExactDiagonalization
using KrylovKit, Random
BLAS.set_num_threads(1)
# The model directory (this file is in src/); results/ and hashed paths are relative to it.
const ROOT = dirname(@__DIR__)
const REPO = normpath(joinpath(ROOT,"../.."))
const FIELDS = [1//1,3//4,1//2,7//20,1//4]
rat(x::AbstractString) = (p=split(x,"//"); parse(BigInt,p[1])//parse(BigInt,p[2]))
tag(h) = replace(string(h),"//"=>"_")
iv(x) = CertifiedGaps.Recurrence.to_interval(x)
iv(x::Irrational) = interval(Float64,x)
"Run metadata with SHA-256 of this model's sources and of the core (`core/...`)."
function provenance(config; certified=false)
    hashes = source_hashes(ROOT; files=("certify.jl",))
    merge!(hashes, Dict("core/"*k=>v for (k,v) in source_hashes(joinpath(REPO,"julia"))))
    run_metadata(config;source_sha256=hashes,packages=(CertifiedGaps,IntervalArithmetic),
        status=certified ? "certified" : "numerical, not certified",certified)
end
"Exact rational check of the sublattice rotation, the RP identity −h = Σ C⊗C and o = −h Z."
function exact_rp(h)
    X = [0//1 1//2;1//2 0//1]; J=[0//1 1//2;-1//2 0//1]
    Z=[1//2 0//1;0//1 -1//2]; U=2J; Id=Matrix{Rational{Int}}(I,2,2)
    original=kron(X,X)-kron(J,J)+kron(Z,Z)
    rotated=-kron(X,X)-kron(J,J)-kron(Z,Z)
    @assert kron(U,Id)*original*kron(U',Id)==rotated
    @assert kron(Id,U)*original*kron(Id,U')==rotated
    @assert -rotated==sum(kron(C,C) for C in (X,J,Z))
    @assert U*Z*U' == -Z
    @assert issymmetric(-h*Z) && isreal(-h*Z)
    true
end
rot(s,L)=((s<<2)&((1<<L)-1)) | (s>>(L-2))
"Exact orbit partition under translation by two sites in a fixed particle sector."
function orbits(L,p)
    reps=Int[]; lens=Int[]; owner=zeros(Int32,1<<L); shift=zeros(Int8,1<<L)
    for s in 0:(1<<L)-1
        count_ones(s)==p && owner[s+1]==0 || continue
        push!(reps,s); t=s; ell=0
        while true
            owner[t+1]=length(reps);shift[t+1]=ell;ell+=1;t=rot(t,L)
            t==s && break
        end
        push!(lens,ell)
    end
    @assert sum(lens)==binomial(L,p)
    (;reps,lens,owner,shift)
end
"Momentum block, with optional interval enclosure of every exact matrix entry."
function block(L,h,k,o; verified=false)
    m=L÷2; ids=findall(l->(k*l)%m==0,o.lens); n=length(ids)
    loc=zeros(Int,length(o.reps));loc[ids]=1:n
    V=verified ? Complex{Interval{Float64}} : ComplexF64
    entries=Dict{Tuple{Int,Int},V}()
    phases=verified ? [complex(cos(iv(2*k*t)*iv(pi)/iv(m)),sin(iv(2*k*t)*iv(pi)/iv(m))) for t in 0:m-1] :
        [cispi(2*k*t/m) for t in 0:m-1]
    for (col,a) in enumerate(ids)
        s=o.reps[a]; diag=0//1
        for j in 0:L-1
            z=1-2*((s>>j)&1); zn=1-2*((s>>mod(j+1,L))&1)
            diag+=z*zn//4-h*(-1)^(j+1)*z/2
            if z!=zn
                t=s ⊻ (1<<j) ⊻ (1<<mod(j+1,L));b=o.owner[t+1];row=loc[b]
                row==0 && continue
                f=verified ? sqrt(iv(o.lens[a])/iv(o.lens[b]))/iv(2) : sqrt(o.lens[a]/o.lens[b])/2
                v=f*phases[Int(o.shift[t+1])+1]
                entries[(row,col)]=get(entries,(row,col),zero(v))+v
            end
        end
        v=verified ? complex(iv(diag),iv(0)) : complex(Float64(diag))
        entries[(col,col)]=get(entries,(col,col),zero(v))+v
    end
    ii=Int[];jj=Int[];vv=ComplexF64[]; err2=iv(0)
    for ((i,j),v) in entries
        i>j && continue
        z=verified ? complex(mid(real(v)),mid(imag(v))) : v
        i==j && (z=complex(real(z),0))
        push!(ii,i);push!(jj,j);push!(vv,z)
        if i!=j
            push!(ii,j);push!(jj,i);push!(vv,conj(z))
        end
        if verified
            er=max(abs(inf(real(v)-iv(real(z)))),abs(sup(real(v)-iv(real(z)))))
            ei=max(abs(inf(imag(v)-iv(imag(z)))),abs(sup(imag(v)-iv(imag(z)))))
            err2+=iv(i==j ? 1 : 2)*(iv(er)^2+iv(ei)^2)
        end
    end
    sparse(ii,jj,vv,n,n), verified ? sup(sqrt(err2)) : 0.0
end
# Outward Frobenius norm; every scalar multiply/add is rounded upward. nextfloat also
# covers underflow. No Float64 reduction is accepted as an inequality.
function frob(A)
    s=0.0
    for z in A
        a=nextfloat(real(z)*real(z));b=nextfloat(imag(z)*imag(z))
        s=nextfloat(s+nextfloat(a+b))
    end
    nextfloat(sqrt(s))
end
"Residual/polar/Weyl enclosure, with conservative IEEE complex matmul error bounds."
function enclose(H,F,entry_error,L,h)
    Q=F.vectors; lam=F.values;n=length(lam)
    R=H*Q-Q*Diagonal(lam)
    E=Q'*Q-I
    q=iv(frob(Q));hf=iv(frob(nonzeros(H)));lf=iv(frob(lam))
    u=iv(1//big(2)^53); g=iv(32n+32)*u/(iv(1)-iv(32n+32)*u)
    # Standard dot-product gamma bound, enlarged for complex products and subtraction.
    # Frobenius product inequality avoids an extra dense abs(A)*abs(B) multiplication.
    r=iv(frob(R))+g*(hf*q+q*lf)+iv(1//big(10)^280)
    e=iv(frob(E))+g*(q^2+sqrt(iv(n)))+iv(1//big(10)^280)
    @assert sup(e)<1
    d=e/(iv(1)+sqrt(iv(1)-e))
    hn=iv(L)*(iv(3)/iv(4)+iv(h)/iv(2))+iv(entry_error)
    eta=iv(entry_error)+r+(hn+iv(maximum(abs,lam)))*d
    (;radius=sup(eta),residual_upper=sup(r),orthogonality_upper=sup(e),
        entry_error,dimension=n)
end
"""
Spectrum of the L-site ring at field h from one representative per (S^z, momentum) equivalence
class. With `verified=true` every block is dense and its record carries the eigenvalues and the
`enclose` radius; otherwise blocks above `dense_max` are partial Krylov spectra (numerical).
"""
function all_spectra(L,h; verified=false,dense_max=2500,nev=80)
    start=time();levels=Float64[];records=Any[];N=0;cut=Inf
    for p in 0:L÷2
        o=orbits(L,p)
        for k in 0:(L÷2)÷2
            H,err=block(L,h,k,o;verified);n=size(H,1);n==0 && continue
            # Global spin flip followed by one-site translation maps p to L-p.
            # Complex conjugation maps k to m-k. Both are exact unitary/antiunitary maps.
            mult=(p==L÷2 ? 1 : 2)*(k==0 || 2k==L÷2 ? 1 : 2)
            t=time();complete=verified || n<=dense_max
            if complete
                if verified
                    F=eigen(Hermitian(Matrix(H))); vals=F.values
                    check=enclose(H,F,err,L,h)
                else
                    vals=eigvals(Hermitian(Matrix(H)))
                end
            else
                Random.seed!(L*10000+p*100+k)
                es,_,info=eigsolve(H,randn(ComplexF64,n),min(nev,n-2),:SR;
                    ishermitian=true,krylovdim=min(n,max(2nev+30,180)),tol=1e-11,maxiter=400)
                @assert info.converged>=min(nev,n-2)
                vals=sort(real.(es[1:min(nev,n-2)]));cut=min(cut,vals[end])
            end
            for _ in 1:mult; append!(levels,vals);end
            N+=n*mult
            rec=Dict{String,Any}("p"=>p,"k"=>k,"dimension"=>n,"multiplicity"=>mult,
                "dense"=>complete,"seconds"=>time()-t)
            if verified
                rec["eigenvalues"]=vals
                for (key,v) in pairs(check);rec[string(key)]=v;end
            end
            push!(records,rec)
            println("L=$L h=$h p=$p k=$k dim=$n mult=$mult dense=$complete seconds=$(round(time()-t;digits=2))");flush(stdout)
        end
    end
    @assert N==2^L
    isfinite(cut) && filter!(x->x<=cut,levels)
    spec=PartialSpectrum(sort!(levels),N,isinf(cut))
    (;spec,records,seconds=time()-start)
end
"Outward logZ from ALL enclosed levels; energy shift is an exact stored Float64."
function enclosed_logz(records,beta)
    shift=minimum(minimum(r["eigenvalues"]) for r in records)
    z=iv(0);b=iv(beta)
    for r in records
        eta=iv(-r["radius"],r["radius"])
        for e in r["eigenvalues"]
            z+=iv(r["multiplicity"])*exp(-b*(iv(e)-iv(shift)+eta))
        end
    end
    log(z)-b*iv(shift)
end
iv(a,b)=interval(a,b)
"Rational upper bounds on 1 − Pₛ(n, β) and 1 − P_τ(2n, β) from the block records of n and 2n."
function enclosed_seeds(a,b,beta)
    x=enclosed_logz(a,beta);y=enclosed_logz(b,beta);z=enclosed_logz(b,2beta)
    es=iv(1)-exp(y-iv(2)*x);et=iv(1)-exp(z-iv(2)*y)
    # Rationalize the binary upper endpoints exactly; no decimal conversion.
    s=Rational{BigInt}(sup(es));t=Rational{BigInt}(sup(et))
    @assert 0<=s<1 && 0<=t<1
    (;s,t,logs=[[inf(x),sup(x)],[inf(y),sup(y)],[inf(z),sup(z)]])
end
