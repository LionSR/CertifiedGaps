using ThermalMPOGap, LinearAlgebra, SparseArrays, JSON3
using ThermalMPOGap.Results, ThermalMPOGap.ExactDiagonalization, ThermalMPOGap.Seeds
using IntervalArithmetic
BLAS.set_num_threads(1)
const ROOT = @__DIR__
const REPO = normpath(joinpath(ROOT,"../.."))
const DELTAS = [9//10, 4//5, 7//10, 3//5, 1//2]
ratstr(x::Rational) = string(numerator(x),"//",denominator(x))
parserat(s) = (a=split(s,"//");parse(BigInt,a[1])//parse(BigInt,a[2]))
function provenance(config; certified=false)
    hashes = source_hashes(ROOT; dirs=("test",), files=("common.jl","scan.jl","certify.jl","rp.jl","verified_blocks.jl","check.jl","Project.toml"))
    merge!(hashes, Dict("../../"*k=>v for (k,v) in source_hashes(REPO;dirs=("julia/src",),files=("julia/Project.toml","julia/Manifest.toml"))))
    run_metadata(config;source_sha256=hashes,packages=(ThermalMPOGap,IntervalArithmetic,JSON3),status=certified ? "certified" : "numerical",certified)
end
"Exact spin-1 Sz-basis Hamiltonian, scaled to an integer matrix."
function integer_hamiltonian(cells, delta, states)
    L=2cells; den=denominator(delta); num=numerator(delta)
    index=Dict(s=>i for (i,s) in enumerate(states)); place=3 .^ (0:L-1)
    rows=Int[]; cols=Int[]; vals=Int[]
    for (col,s) in enumerate(states)
        diag=0
        for j in 1:L
            k=mod1(j+1,L); a=(s÷place[j])%3; b=(s÷place[k])%3
            w=den+(iseven(j) ? num : -num)
            diag+=w*(a-1)*(b-1)
            for e in (-1,1)
                if 0<=a+e<=2 && 0<=b-e<=2
                    t=s+e*(place[j]-place[k]); row=index[t]
                    push!(rows,row);push!(cols,col);push!(vals,w)
                end
            end
        end
        push!(rows,col);push!(cols,col);push!(vals,diag)
    end
    A=sparse(rows,cols,vals,length(states),length(states))
    @assert A==transpose(A)
    return A,den
end
function charge_groups(cells)
    L=2cells; groups=[Int[] for _ in 0:L]
    for s in 0:3^L-1
        x=s; m=-L
        for j in 1:L; m+=x%3; x÷=3; end
        m>=0 && push!(groups[m+1],s)
    end
    @assert length(groups[1])+2sum(length,groups[2:end])==9^cells
    groups
end
translate(s,L) = (s÷9)+(s%9)*3^(L-2)
"Build normalized character projections of cell translation and (at M=0) spin flip."
function orbit_basis(states,cells,k,parity)
    L=2cells; flipmax=3^L-1; idx=Dict(s=>i for (i,s) in enumerate(states))
    visited=Set{Int}(); rows=Int[];cols=Int[];vals=ComplexF64[]; col=0
    for s in states
        s in visited && continue
        orbit=Dict{Int,ComplexF64}(); allowed=true
        for f in 0:(parity==0 ? 0 : 1)
            t=f==0 ? s : flipmax-s
            for r in 0:cells-1
                phase=cispi(-2*k*r/cells)*(f==0 ? 1 : parity)
                if haskey(orbit,t) && abs(orbit[t]-phase)>1e-8; allowed=false; end
                orbit[t]=phase;push!(visited,t);t=translate(t,L)
            end
        end
        allowed || continue
        col+=1; normfactor=sqrt(length(orbit))
        for (t,v) in orbit
            push!(rows,idx[t]);push!(cols,col);push!(vals,v/normfactor)
        end
    end
    sparse(rows,cols,vals,length(states),col)
end
"Full dense block spectra; charge and conjugate momentum equivalents carry multiplicity."
function dense_spectrum(cells,delta; vectors=false)
    levels=Float64[]; records=Any[]; start=time(); groups=charge_groups(cells)
    for m in 0:2cells
        states=groups[m+1]; A,den=integer_hamiltonian(cells,delta,states)
        for k in 0:fld(cells,2), p in (m==0 ? (-1,1) : (0,)), reflection in (k==0 || 2k==cells ? (-1,1) : (1,))
            B=real_orbit_basis(states,cells,k,p,reflection); dim=size(B,2); dim==0 && continue
            H=Matrix(B'*(A*B))/den
            @assert maximum(abs,H-H')<1e-10
            mult=(m==0 ? 1 : 2)*(k==0 || 2k==cells ? 1 : 2)
            t=@elapsed E=eigvals(Hermitian(k==0 || 2k==cells ? real(H) : H))
            for _ in 1:mult;append!(levels,E);end
            push!(records,(charge=m,momentum=k,spin_flip=p,reflection=reflection,dimension=dim,multiplicity=mult,seconds=t))
            println("cells=$cells delta=$delta M=$m k=$k p=$p reflection=$reflection dim=$dim seconds=$(round(t,digits=2))");flush(stdout)
        end
    end
    @assert sum(r.dimension*r.multiplicity for r in records)==9^cells==length(levels)
    sort!(levels)
    return PartialSpectrum(levels,9^cells,true),records,time()-start
end
reflect_state(s,L)=foldl((a,_) -> (a[1]*3+a[2]%3,a[2]÷3),1:L;init=(0,s))[1]
"Real dihedral character basis; generic momenta paired with exact multiplicity two."
function real_orbit_basis(states,cells,k,parity,reflection)
    L=2cells;flipmax=3^L-1;idx=Dict(s=>i for (i,s) in enumerate(states))
    visited=Set{Int}();rows=Int[];cols=Int[];vals=Float64[];col=0
    generic=!(k==0 || 2k==cells)
    factor=(generic ? 2 : 1)/(cells*2*(parity==0 ? 1 : 2))
    for s in states
        s in visited && continue
        orbit=Set{Int}()
        for f in 0:(parity==0 ? 0 : 1),h in 0:1
            t=f==0 ? s : flipmax-s;t=h==0 ? t : reflect_state(t,L)
            for r in 0:cells-1;push!(orbit,t);t=translate(t,L);end
        end
        union!(visited,orbit);os=sort!(collect(orbit)); oi=Dict(t=>j for (j,t) in enumerate(os))
        P=zeros(length(os),length(os))
        for (j,v) in enumerate(os), f in 0:(parity==0 ? 0 : 1), h in 0:1
            t=f==0 ? v : flipmax-v;t=h==0 ? t : reflect_state(t,L)
            for r in 0:cells-1
                P[oi[t],j]+=factor*cospi(2*k*r/cells)*(f==0 ? 1 : parity)*(h==0 ? 1 : reflection)
                t=translate(t,L)
            end
        end
        F=eigen(Symmetric(P))
        @assert maximum(min.(abs.(F.values),abs.(F.values.-1)))<1e-10
        for a in findall(>(0.5),F.values)
            col+=1
            for (i,t) in enumerate(os)
                push!(rows,idx[t]);push!(cols,col);push!(vals,F.vectors[i,a])
            end
        end
    end
    sparse(rows,cols,vals,length(states),col)
end
