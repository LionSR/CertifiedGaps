using LinearAlgebra, SparseArrays
include("model.jl")

"Translate the base-3 product state `s` of an `L`-site ring by one site."
rotate_state(s,L) = (s % 3^(L-1))*3 + s ÷ 3^(L-1)
"Reflect the base-3 product state `s` of an `L`-site ring (reverse its digits)."
function reflect_state(s,L)
    r=0
    for _ in 1:L
        r=3r+s%3; s÷=3
    end
    r
end
"Total magnetization M of the base-3 product state `s` (digit 0, 1, 2 means m = -1, 0, 1)."
function charge(s,L)
    m=-L
    for _ in 1:L
        m+=s%3; s÷=3
    end
    m
end

"Translation orbits, their periods, and a lookup (representative, translation exponent)."
function orbits(L,M)
    N=3^L
    owner=zeros(Int32,N); shift=zeros(Int8,N)
    reps=Int[]; periods=Int[]
    for s in 0:N-1
        owner[s+1]!=0 && continue
        charge(s,L)==M || continue
        push!(reps,s); idx=length(reps)
        t=s; p=0
        while true
            owner[t+1]=idx; shift[t+1]=p
            p+=1; t=rotate_state(t,L)
            t==s && break
        end
        push!(periods,p)
    end
    (;reps,periods,owner,shift)
end

"Sparse momentum Hamiltonian components; only transitions out of orbit representatives needed."
function momentum_components(L,k,o)
    active=findall(p -> mod(k*p,L)==0, o.periods)
    lookup=zeros(Int,length(o.reps)); lookup[active]=1:length(active)
    I=Int[]; J=Int[]; VA=ComplexF64[]; VB=ComplexF64[]
    X=sz_bond(); Y=X^2; place=[3^j for j in 0:L-1]
    q=2pi*k/L
    for (col,orb) in enumerate(active)
        s=o.reps[orb]; p=o.periods[orb]
        for j in 1:L
            l=mod1(j+1,L); a=s÷place[j]%3; b=s÷place[l]%3; old=3a+b+1
            for new in 1:9
                u=X[new,old]; v=Y[new,old]
                u==0 && v==0 && continue
                aa=(new-1)÷3; bb=(new-1)%3
                t=s+(aa-a)*place[j]+(bb-b)*place[l]
                roworb=o.owner[t+1]; row=lookup[roworb]
                row==0 && continue
                phase=cis(q*o.shift[t+1])*sqrt(p/o.periods[roworb])
                push!(I,row); push!(J,col); push!(VA,u*phase); push!(VB,v*phase)
            end
        end
    end
    d=length(active)
    A=sparse(I,J,VA,d,d); B=sparse(I,J,VB,d,d)
    @assert norm(A-A')<1e-9*max(1,norm(A))
    @assert norm(B-B')<1e-9*max(1,norm(B))
    # Reflection followed by conjugation is an antiunitary involution within k.
    # Its fixed vectors give a real Hamiltonian, avoiding complex dense eigensolves.
    UI=Int[]; UJ=Int[]; UV=ComplexF64[]; seen=falses(d); c=0
    for i in 1:d
        seen[i] && continue
        r=reflect_state(o.reps[active[i]],L)
        j=lookup[o.owner[r+1]]; @assert j!=0
        phase=cis(q*o.shift[r+1]/2)
        if i==j
            c+=1; push!(UI,i);push!(UJ,c);push!(UV,phase)
        else
            c+=1
            append!(UI,[i,j]);append!(UJ,[c,c]);append!(UV,[phase,phase]./sqrt(2))
            c+=1
            append!(UI,[i,j]);append!(UJ,[c,c]);append!(UV,[im*phase,-im*phase]./sqrt(2))
        end
        seen[i]=true;seen[j]=true
    end
    U=sparse(UI,UJ,UV,d,d)
    @assert norm(U'*U-spdiagm(ones(d))) < 1e-10
    Ar=U'*A*U; Br=U'*B*U
    @assert norm(imag(Ar))<1e-9*max(1,norm(Ar))
    @assert norm(imag(Br))<1e-9*max(1,norm(Br))
    return real(Ar),real(Br),active,U
end

"Split a sparse real symmetric involution into its +/- eigenspaces locally."
function involution_bases(S)
    d=size(S,1); @assert norm(S*S-spdiagm(ones(d)))<1e-8*max(1,sqrt(d))
    seen=falses(d); ii=[Int[],Int[]]; jj=[Int[],Int[]]; vv=[Float64[],Float64[]]; cols=[0,0]
    for seed in 1:d
        seen[seed] && continue
        component=[seed]; seen[seed]=true; pos=1
        while pos<=length(component)
            j=component[pos];pos+=1
            for t in nzrange(S,j)
                abs(nonzeros(S)[t])<1e-12 && continue
                i=rowvals(S)[t]
                if !seen[i];push!(component,i);seen[i]=true;end
            end
        end
        @assert length(component)<=16
        eig=eigen(Symmetric(Matrix(S[component,component])))
        for j in eachindex(eig.values)
            @assert abs(abs(eig.values[j])-1)<1e-9
            signidx=eig.values[j]>0 ? 1 : 2;cols[signidx]+=1
            for (loc,i) in enumerate(component)
                v=eig.vectors[loc,j]
                abs(v)<1e-14 && continue
                push!(ii[signidx],i);push!(jj[signidx],cols[signidx]);push!(vv[signidx],v)
            end
        end
    end
    [sparse(ii[s],jj[s],vv[s],d,cols[s]) for s in 1:2]
end

"""
    reduced_components(L, M, k, o)

Real symmetric bilinear (`A`) and biquadratic (`B`) components of every block of magnetization
`M` and momentum `k`, further split by spin flip (M = 0) and reflection (k = 0, π). The ring
Hamiltonian at angle θ on a block is `cos(θ) A + sin(θ) B`.
"""
function reduced_components(L,M,k,o)
    A,B,active,U=momentum_components(L,k,o)
    lookup=zeros(Int,length(o.reps));lookup[active]=1:length(active)
    syms=Pair{String,SparseMatrixCSC{Float64,Int}}[]
    operations=Pair{String,Function}[]
    M==0 && push!(operations,"flip" => (s->3^L-1-s))
    (k==0 || 2k==L) && push!(operations,"reflection" => (s->reflect_state(s,L)))
    d=length(active)
    for (name,op) in operations
        rows=Int[];vals=ComplexF64[]
        for orb in active
            s=op(o.reps[orb]);push!(rows,lookup[o.owner[s+1]])
            push!(vals,cis(2pi*k/L*o.shift[s+1]))
        end
        S=U'*sparse(rows,collect(1:d),vals,d,d)*U
        @assert norm(imag(S))<1e-8
        push!(syms,name=>real(S))
    end
    blocks=[(;A,B,label="",syms)]
    while !isempty(blocks[1].syms)
        next=typeof(blocks)()
        for block in blocks
            name,S=block.syms[1]
            for (i,V) in enumerate(involution_bases(S))
                size(V,2)==0 && continue
                @assert norm(block.A*V-V*(V'*block.A*V))<1e-8*max(1,norm(block.A))
                push!(next,(;A=V'*block.A*V,B=V'*block.B*V,
                    label=block.label*name*(i==1 ? "+" : "-"),
                    syms=[n=>(V'*s*V) for (n,s) in block.syms[2:end]]))
            end
        end
        blocks=next
    end
    @assert sum(size(b.A,1) for b in blocks)==d
    blocks
end
