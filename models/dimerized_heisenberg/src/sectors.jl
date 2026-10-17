# Symmetry sectors of the dimerized ring for the Float64 diagnostics: exact orbit enumeration
# under cell translations and, at half filling, global spin flip. `dense_spectrum` is the plain
# Float64 solver behind the numerical scans; its eigenvalues are NOT certified (the certificates
# use the core verified spectra).
using CertifiedGaps, LinearAlgebra
using CertifiedGaps.ExactDiagonalization
BLAS.set_num_threads(1)

"Translate a configuration of `L` spins by one two-spin cell."
rotate_cell(s,L) = ((s << 2) & ((1<<L)-1)) | (s >> (L-2))

"""
    group_orbits(L, q)

Orbits of the configurations of `L` spins with `q` down spins under cell translations (and global
spin flip when `2q == L`), with each orbit's images and stabiliser; returns `(orbits, useflip)`.
"""
function group_orbits(L, q)
    n=L÷2; mask=(1<<L)-1
    states=[s for s in 0:mask if count_ones(s)==q]
    seen=falses(1<<L); orbits=Any[]
    useflip=(2q==L)
    for rep in states
        seen[rep+1] && continue
        images=Dict{Int,Tuple{Int,Int}}(); stabilizers=Tuple{Int,Int}[]
        for f in 0:Int(useflip)
            s=f==0 ? rep : xor(rep,mask)
            for t in 0:n-1
                if haskey(images,s)
                    t0,f0=images[s]
                    push!(stabilizers,(t-t0,f-f0))
                else
                    images[s]=(t,f); seen[s+1]=true
                end
                s=rotate_cell(s,L)
            end
        end
        push!(orbits,(;rep,images,stabilizers))
    end
    @assert sum(length(o.images) for o in orbits)==binomial(L,q)
    return orbits,useflip
end

"""
    dense_spectrum(delta, L)

Complete Float64 spectrum of the `L`-spin ring by dense LAPACK solves in every Sz, cell-momentum and
spin-flip sector, with per-sector records. Numerical only: no eigenvalue enclosure.
"""
function dense_spectrum(delta,L)
    n=L÷2; levels=Float64[]; records=Any[]; dimension_sum=0
    for q in 0:L÷2
        orbits,useflip=group_orbits(L,q)
        mult=useflip ? 1 : 2
        qdim=0
        for k in 0:n-1, parity in (useflip ? (1,-1) : (1,))
            odd=parity==-1 ? 1 : 0
            allowed=[o for o in orbits if all(mod(2*k*t+n*odd*f,2*n)==0 for (t,f) in o.stabilizers)]
            dim=length(allowed); qdim+=dim
            dim==0 && continue
            labels=zeros(Int,1<<L); coeff=zeros(ComplexF64,1<<L)
            for (i,o) in enumerate(allowed), (s,(t,f)) in o.images
                labels[s+1]=i
                coeff[s+1]=cispi(2*k*t/n)*parity^f/sqrt(length(o.images))
            end
            H=zeros(ComplexF64,dim,dim)
            for (j,o) in enumerate(allowed)
                s=o.rep; scale=sqrt(length(o.images))
                diag=0.0
                for site in 0:L-1
                    next=mod(site+1,L)
                    J=1+(iseven(site) ? delta : -delta)
                    opposite=((s>>site)&1)!=((s>>next)&1)
                    diag+=J*(opposite ? -0.25 : 0.25)
                    if opposite
                        target=xor(s,(1<<site)|(1<<next)); i=labels[target+1]
                        i==0 || (H[i,j]+=J/2*scale*conj(coeff[target+1]))
                    end
                end
                H[j,j]+=diag
            end
            defect=maximum(abs,H-H')
            @assert defect < 1e-10
            duration=@elapsed E=eigvals(Hermitian(H)) # Dense LAPACK only.
            for _ in 1:mult; append!(levels,E); end
            dimension_sum+=mult*dim
            push!(records,Dict("down_spins"=>q,"cell_momentum"=>k,
                "spin_flip_parity"=>(useflip ? parity : nothing),"dimension"=>dim,
                "multiplicity"=>mult,"hermiticity_defect"=>defect,"seconds"=>duration))
        end
        @assert qdim==binomial(L,q)
    end
    @assert dimension_sum==2^L
    @assert length(levels)==2^L
    return PartialSpectrum(sort!(levels),2^L,true),records
end

