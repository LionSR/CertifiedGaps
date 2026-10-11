using LinearAlgebra, SparseArrays, CertifiedGaps, IntervalArithmetic
using CertifiedGaps.Results, CertifiedGaps.ExactDiagonalization, CertifiedGaps.Seeds
using KrylovKit, JSON3, Random
BLAS.set_num_threads(1)

"Exact Cartesian spin algebra and independent bond identity P-|Omega><Omega|."
function rp_identity(D::Rational)
    G = [[0 0 0; 0 0 1; 0 -1 0], [0 0 -1; 0 0 0; 1 0 0],
         [0 1 0; -1 0 0; 0 0 0]]
    S = [-im .* g for g in G]
    @assert all(g == -transpose(g) for g in G)
    @assert S[1]*S[2]-S[2]*S[1] == im*S[3]
    @assert sum(s*s for s in S) == 2I
    swap = zeros(Int,9,9)
    for a in 1:3, b in 1:3
        swap[3(a-1)+b,3(b-1)+a] = 1
    end
    omega = vec(Matrix{Int}(I,3,3))
    h = swap - omega*omega'
    @assert sum(kron(s,s) for s in S) == h
    @assert -h == sum(kron(g,g) for g in G)
    onsite = D * Diagonal([1,1,0])
    @assert onsite == D*S[3]^2 && isreal(onsite) && issymmetric(onsite)
    return true
end

rot(s,L) = (s%3)*3^(L-1)+s÷3
flip(s,L) = 3^L-1-s
charge(s,L) = sum(digits(s,base=3,pad=L)) - L

"Translation orbits, with spin inversion included inside charge zero."
function orbit_data(L,q)
    N=3^L
    owner=zeros(Int,N); shift=zeros(Int,N); flipped=zeros(Int,N)
    reps=Int[]; sizes=Int[]; periods=Int[]; flipshifts=Int[]
    for s in 0:N-1
        charge(s,L)==q && owner[s+1]==0 || continue
        push!(reps,s); id=length(reps)
        t=s; p=0
        while true
            owner[t+1]=id; shift[t+1]=p
            p+=1; t=rot(t,L); t==s && break
        end
        fs=-1; size=p
        if q==0
            f=flip(s,L)
            if owner[f+1]==id
                fs=shift[f+1]
            else
                t=f
                for u in 0:p-1
                    @assert owner[t+1]==0
                    owner[t+1]=id; shift[t+1]=u; flipped[t+1]=1
                    t=rot(t,L)
                end
                size=2p
            end
        end
        push!(sizes,size); push!(periods,p); push!(flipshifts,fs)
    end
    (;L,q,owner,shift,flipped,reps,sizes,periods,flipshifts)
end

function active_orbits(o,k,parity)
    [a for a in eachindex(o.reps) if mod(k*o.periods[a],o.L)==0 &&
      (o.q!=0 || o.flipshifts[a]<0 ||
       mod(2k*o.flipshifts[a],2o.L)==(parity==1 ? 0 : o.L))]
end

"Exact integer transitions in the normalized Sz basis, including two bonds at L=2."
function transitions(s,L,D)
    ds=digits(s,base=3,pad=L); out=Dict{Int,Rational{Int}}()
    diagonal=D*sum((x-1)^2 for x in ds)
    for j in 1:L
        l=mod1(j+1,L); a=ds[j]; b=ds[l]
        diagonal+=(a-1)*(b-1)
        if a<2 && b>0
            t=s+3^(j-1)-3^(l-1); out[t]=get(out,t,0//1)+1
        end
        if a>0 && b<2
            t=s-3^(j-1)+3^(l-1); out[t]=get(out,t,0//1)+1
        end
    end
    out[s]=get(out,s,0//1)+diagonal
    out
end

"Build a Hermitian momentum block and a rigorous operator norm construction error."
function block(o,k,parity,D; verified=false)
    active=active_orbits(o,k,parity); m=length(active)
    index=zeros(Int,length(o.reps)); index[active]=1:m
    rows=Int[]; cols=Int[]; vals=ComplexF64[]
    errors=zeros(Float64,m)
    piI=interval(Float64,pi)
    for (col,a) in enumerate(active)
        sums=Dict{Int,ComplexF64}()
        boxes=Dict{Int,Complex{Interval{Float64}}}()
        for (t,h) in transitions(o.reps[a],o.L,D)
            b=o.owner[t+1]; row=index[b]; row==0 && continue
            row>col && continue # Hermiticity supplies the lower triangle.
            u=o.shift[t+1]; sign=parity^o.flipped[t+1]
            v=Float64(h)*sqrt(o.sizes[a]/o.sizes[b])*cis(-2pi*k*u/o.L)*sign
            sums[row]=get(sums,row,0.0+0.0im)+v
            if verified
                theta=-2*piI*interval(k*u)/interval(o.L)
                scale=CertifiedGaps.to_interval(h)*sqrt(interval(o.sizes[a])/interval(o.sizes[b]))*sign
                z=complex(scale*cos(theta),scale*sin(theta))
                boxes[row]=get(boxes,row,complex(interval(0),interval(0)))+z
            end
        end
        for (row,v0) in sums
            v=row==col || k==0 || 2k==o.L ? complex(real(v0),0) : v0
            push!(rows,row); push!(cols,col); push!(vals,v)
            if row!=col
                push!(rows,col); push!(cols,row); push!(vals,conj(v))
            end
            if verified
                z=boxes[row]
                er=sup(abs(real(z)-real(v))+abs(imag(z)-imag(v)))
                errors[row]=sup(interval(errors[row])+interval(er))
                row==col || (errors[col]=sup(interval(errors[col])+interval(er)))
            end
        end
    end
    H=sparse(rows,cols,vals,m,m)
    (k==0 || 2k==o.L ? real(H) : H),maximum(errors;init=0.0)
end

"Dense eigensystem enclosure using polar decomposition, Weyl and BLAS error bounds."
function enclose_eigen(H,construction,L,D)
    A=Matrix(H); F=eigen(Hermitian(A)); Q=F.vectors; lam=F.values; m=length(lam)
    # Standard real-operation model for complex dot products, conservatively 8m+32.
    # All nonzero data are normal; the additive term covers gradual underflow anyway.
    u=interval(2.0^-53); gamma=(8m+32)*u/(1-(8m+32)*u)
    function fnormupper(X)
        # Each squared complex magnitude needs <=3 operations, followed by <=2m^2 sums.
        n=length(X); g=(4n+16)*u/(1-(4n+16)*u)
        sqrt(interval(sum(abs2,X))/(1-g)+interval(1e-290))
    end
    qf=fnormupper(Q)
    M=interval(L)*(CertifiedGaps.to_interval(abs(D))+2)+interval(construction)
    af=sqrt(interval(m))*M
    R=A*Q-Q*Diagonal(lam)
    r=fnormupper(R)+gamma*(af*qf+qf*maximum(abs,lam))+interval(1e-280)
    E=Q'*Q-Matrix{ComplexF64}(I,m,m)
    e=fnormupper(E)+gamma*(qf*qf+sqrt(interval(m)))+interval(1e-280)
    @assert sup(e)<1
    # Q=U P, ||Q-U|| <= e/(1+sqrt(1-e)). U is unitary.
    d=e/(1+sqrt(1-e))
    eta=sup(r+(M+maximum(abs,lam))*d+interval(construction))
    lam,eta,Dict("dimension"=>m,"residual_upper"=>sup(r),
        "orthogonality_upper"=>sup(e),"construction_upper"=>construction,"radius"=>eta)
end

function compute_spectrum(L,D; rigorous=false,low=false)
    started=time(); levels=Float64[]; radii=Float64[]; sectors=Any[]; total=0; cut=Inf
    for q in 0:L
        o=orbit_data(L,q)
        for parity in (q==0 ? (1,-1) : (1,)), k in 0:div(L,2)
            H,c=block(o,k,parity,D;verified=rigorous); m=size(H,1); m==0 && continue
            mult=(q==0 ? 1 : 2)*(k==0 || 2k==L ? 1 : 2)
            tick=time()
            if rigorous
                E,eta,info=enclose_eigen(H,c,L,D)
            elseif low && m>150
                nkeep=min(48,m-1)
                ev,_,conv=eigsolve(H,randn(eltype(H),m),nkeep,:SR;
                    ishermitian=true,krylovdim=max(120,2nkeep+10),tol=1e-12,maxiter=500)
                @assert conv.converged>=nkeep
                E=sort(real(ev[1:nkeep])); eta=0.0; info=Dict{String,Any}()
                cut=min(cut,E[end])
            else
                E=eigvals(Hermitian(Matrix(H))); eta=0.0; info=Dict{String,Any}()
            end
            merge!(info,Dict("q"=>q,"k"=>k,"parity"=>parity,"dimension"=>m,
                "multiplicity"=>mult,"levels"=>length(E),"seconds"=>time()-tick))
            push!(sectors,info); total+=m*mult
            for _ in 1:mult
                append!(levels,E); append!(radii,fill(eta,length(E)))
            end
            rigorous && println("verified D=$D L=$L q=$q k=$k parity=$parity dim=$m eta=$eta seconds=$(time()-tick)")
            flush(stdout)
        end
    end
    @assert total==3^L
    if isfinite(cut)
        keep=findall(<=(cut),levels); levels=levels[keep]; radii=radii[keep]
    end
    perm=sortperm(levels)
    Dict("D"=>string(D),"L"=>L,"energies"=>levels[perm],"radii"=>radii[perm],
        "dimension"=>total,"complete"=>length(levels)==total,"rigorous"=>rigorous,
        "sectors"=>sectors,"seconds"=>time()-started)
end
