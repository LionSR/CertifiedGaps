using CertifiedGaps, LinearAlgebra, SparseArrays, Printf
using CertifiedGaps.Results, CertifiedGaps.ExactDiagonalization, CertifiedGaps.Seeds
using IntervalArithmetic
BLAS.set_num_threads(1)
include("rp.jl")
const ROOT = normpath(joinpath(@__DIR__, "../.."))

"Record all executable sources and the unchanged core/environment manifests."
function provenance(config; certified=false)
    hashes = source_hashes(ROOT; dirs=("julia/src", "models/heisenberg_ladder"),
        files=("julia/Project.toml", "julia/Manifest.toml", "outputs/thermal-mpo-gap-scope.tex"))
    run_metadata(config; source_sha256=hashes, packages=(CertifiedGaps,IntervalArithmetic),
                 status=certified ? "certified" : "numerical, not certified", certified)
end

# Bits 2j and 2j+1 are the two spins of rung j; bit value one is spin down.
translate(s,L) = ((s << 2) & ((1 << (2L))-1)) | (s >> (2L-2))
function swaplegs(s,L)
    mask = sum(1 << (2j) for j in 0:L-1)
    ((s & mask) << 1) | ((s >> 1) & mask)
end

"Exact permutation orbits for C_L x leg exchange x (spin reversal at zero magnetization)."
function orbit_data(L,m)
    N = 1 << (2L)
    rid = zeros(Int,N); tt=zeros(Int,N); aa=zeros(Int,N); bb=zeros(Int,N)
    reps=Int[]; sizes=Int[]; stabs=Vector{NTuple{3,Int}}[]
    flip = m == L
    for s in 0:N-1
        count_ones(s) == m && rid[s+1] == 0 || continue
        push!(reps,s); r=length(reps); o=0; stab=NTuple{3,Int}[]
        for b in 0:Int(flip), a in 0:1
            v = b==0 ? s : xor(s,N-1)
            v = a==0 ? v : swaplegs(v,L)
            for t in 0:L-1
                v==s && push!(stab,(t,a,b))
                if rid[v+1]==0
                    rid[v+1]=r; tt[v+1]=t; aa[v+1]=a; bb[v+1]=b; o+=1
                end
                v=translate(v,L)
            end
        end
        push!(sizes,o); push!(stabs,stab)
    end
    @assert sum(sizes) == binomial(2L,m)
    (; L,m,reps,sizes,stabs,rid,tt,aa,bb,flip)
end
exponent(k,p,q,t,a,b,L) = mod(2k*t + L*(p*a+q*b),2L)
function block_basis(o,k,p,q)
    selected = findall(st -> all(exponent(k,p,q,t,a,b,o.L)==0 for (t,a,b) in st),o.stabs)
    index=zeros(Int,length(o.reps))
    for (i,r) in enumerate(selected); index[r]=i; end
    selected,index
end

"Exact Heisenberg transitions in the original spin basis (rational at stated couplings)."
function transitions(s,L,J)
    out=Dict{Int,Rational{Int}}(s=>0//1)
    for j in 0:L-1
        pairs=((2j,2j+1,J),(2j,2mod(j+1,L),1//1),
               (2j+1,2mod(j+1,L)+1,1//1))
        for (a,b,w) in pairs
            different = ((s>>a)&1) != ((s>>b)&1)
            out[s] += different ? -w/4 : w/4
            if different
                t=xor(s,(1<<a)|(1<<b)); out[t]=get(out,t,0//1)+w/2
            end
        end
    end
    out
end

"Dense symmetry block, optionally enclosed entrywise using outward interval arithmetic."
function block_matrix(o,k,p,q,J; verified=false)
    sel,idx=block_basis(o,k,p,q); d=length(sel)
    realblock = mod(2k,o.L)==0
    H=zeros(realblock ? Float64 : ComplexF64,d,d)
    if verified
        re=fill(interval(0.0),d,d); im=copy(re)
        phases=[(cospi(to_interval(e//o.L)), -sinpi(to_interval(e//o.L))) for e in 0:2o.L-1]
    else
        phases=[cispi(-e/o.L) for e in 0:2o.L-1]
    end
    for (col,r) in enumerate(sel)
        for (s,h) in transitions(o.reps[r],o.L,J)
            rr=o.rid[s+1]; row=idx[rr]; row==0 && continue
            e=exponent(k,p,q,o.tt[s+1],o.aa[s+1],o.bb[s+1],o.L)
            if verified
                amp=to_interval(h)*sqrt(to_interval(o.sizes[r]//o.sizes[rr]))
                re[row,col] += amp*phases[e+1][1]
                im[row,col] += amp*phases[e+1][2]
            else
                v=Float64(h)*sqrt(o.sizes[r]/o.sizes[rr])*phases[e+1]
                H[row,col] += realblock ? real(v) : v
            end
        end
    end
    radius=interval(0.0)
    if verified
        errors=fill(interval(0.0),d)
        for j in 1:d, i in 1:j
            a=mid(re[i,j]); b=realblock || i==j ? 0.0 : mid(im[i,j])
            H[i,j]=realblock ? a : complex(a,b); H[j,i]=conj(H[i,j])
            err=abs(re[i,j]-a)+abs(im[i,j]-b)
            # Exact block is Hermitian; use upper triangle for both halves.
            errors[i]+=err
            i!=j && (errors[j]+=err)
        end
        radius=interval(0.0,maximum(sup,errors;init=0.0))
    else
        @assert isapprox(H,H';atol=1e-10,rtol=1e-12)
        H=Matrix(Hermitian(H))
    end
    H,sup(radius)
end

"Safe Frobenius upper bound using interval accumulation of the exact floating-point entries."
function frob(A)
    s=interval(0.0)
    for z in A
        s+=interval(real(z))^2 + interval(imag(z))^2
    end
    sqrt(s)
end

"Polar-factor/Weyl enclosure, including explicit IEEE dot-product error allowances."
function verified_eigen(H,herror)
    n=size(H,1); f=eigen(Hermitian(H)); Q=f.vectors; lam=f.values
    R=H*Q-Q*Diagonal(lam); E=Q'*Q-I
    # 8n+32 bounds the real operations in complex dot products plus subtraction.
    u=interval(2.0)^(-53); g=(8n+32)*u/(1-(8n+32)*u)
    hn=frob(H); qn=frob(Q); ln=frob(lam)
    # Absolute underflow allowance, far larger than 8n^3 smallest subnormals.
    tiny=interval(2.0)^(-1000)
    r=frob(R)+g*(hn*qn+qn*ln)+tiny
    eta=frob(E)+g*(qn^2+sqrt(interval(n)))+tiny
    @assert sup(eta)<1
    dist=eta/(1+sqrt(1-eta))
    rho=r+(hn+maximum(abs,lam))*dist+interval(herror)
    audit=Dict("dimension"=>n,"matrix_error"=>string(Rational{BigInt}(herror)),
        "residual_upper"=>string(Rational{BigInt}(sup(r))),
        "orthogonality_upper"=>string(Rational{BigInt}(sup(eta))),
        "eigenvalue_radius"=>string(Rational{BigInt}(sup(rho))))
    lam,sup(rho),audit
end

"Full dense spectrum with exact symmetry multiplicities and dimension accounting."
function full_spectrum(L,J; verified=false)
    start=time(); energies=Float64[]; radii=Float64[]; rows=Any[]; dimension=0
    for m in 0:L
        o=orbit_data(L,m)
        for k in 0:fld(L,2), p in 0:1, q in 0:Int(o.flip)
            # Complex conjugate momenta and opposite magnetizations are isospectral.
            mult=(m==L ? 1 : 2)*(mod(2k,L)==0 ? 1 : 2)
            tick=time(); H,err=block_matrix(o,k,p,q,J;verified)
            d=size(H,1); d==0 && continue
            if verified
                ev,rho,audit=verified_eigen(H,err)
            else
                if eltype(H)==ComplexF64; H=realify!(H,o,k,p,q); end
                ev=eigvals(Hermitian(H)); rho=0.0; audit=Dict{String,Any}()
            end
            merge!(audit,Dict("m"=>m,"k"=>k,"leg_parity"=>(-1)^p,
                "spin_parity"=>o.flip ? (-1)^q : 0,"dimension"=>d,
                "multiplicity"=>mult,"seconds"=>time()-tick))
            push!(rows,audit)
            for _ in 1:mult; append!(energies,ev); append!(radii,fill(rho,d)); end
            dimension+=mult*d
        end
        @printf("J=%g L=%d m=%d levels=%d elapsed=%.1fs\n",Float64(J),L,m,dimension,time()-start)
        flush(stdout)
    end
    @assert dimension == 4^L == length(energies)
    perm=sortperm(energies); energies=energies[perm]; radii=radii[perm]
    (; energies,radii,rows,seconds=time()-start,dimension)
end

function reflectrungs(s,L)
    out=0
    for j in 0:L-1
        out |= ((s >> (2j)) & 3) << (2mod(-j,L))
    end
    out
end

"Real basis fixed by the antiunitary rung reflection times complex conjugation (scan only)."
function realify!(H,o,k,p,q)
    sel,idx=block_basis(o,k,p,q); n=length(sel); seen=falses(n)
    transformations=Tuple{Int,Int,ComplexF64}[]
    for (i,r) in enumerate(sel)
        seen[i] && continue
        s=reflectrungs(o.reps[r],o.L); j=idx[o.rid[s+1]]
        @assert j!=0
        e=exponent(k,p,q,o.tt[s+1],o.aa[s+1],o.bb[s+1],o.L)
        t=cispi(-e/(2o.L))
        seen[i]=true; seen[j]=true
        push!(transformations,(i,j,t))
    end
    for (i,j,t) in transformations
        if i==j
            for r in 1:n; H[r,i]*=t; end
        else
            a=t/sqrt(2)
            for r in 1:n
                x,y=H[r,i],H[r,j]
                H[r,i]=a*(x+y); H[r,j]=im*a*(x-y)
            end
        end
    end
    for (i,j,t) in transformations
        if i==j
            for c in 1:n; H[i,c]*=conj(t); end
        else
            a=conj(t)/sqrt(2)
            for c in 1:n
                x,y=H[i,c],H[j,c]
                H[i,c]=a*(x+y); H[j,c]=-im*a*(x-y)
            end
        end
    end
    @assert maximum(abs,imag(H);init=0.0)<1e-9
    real.(H)
end

"""
    read_result_text(path)

Read a result file, falling back to `path * ".gz"` (decompressed with the system `gzip`).
The 8-rung verified spectrum is stored gzipped to keep the repository small.
"""
function read_result_text(path)
    isfile(path) && return read(path, String)
    gz = path * ".gz"
    isfile(gz) || error("missing result file $path (or $gz)")
    return read(pipeline(`gzip -dc $gz`), String)
end

"Compress a large result file in place (deterministic header, no timestamp)."
compress_result(path) = run(`gzip -9 -n -f $path`)
