# Verified spectral enclosures; all scalar error bounds are exact rationals.
module VerifiedSpectra
using LinearAlgebra, IntervalArithmetic
using CertifiedGaps.Recurrence: to_interval
using CertifiedGaps.Results
include("sectors.jl") # Exact orbit/stabilizer enumeration.
export verified_ring, scalar_witness, enclosed_logZ, rational_upper, rational_lower
const QQ=Rational{BigInt}
const CI=Complex{Interval{Float64}}
const UNIT=big(1)//big(2)^53
const TINY=big(8)//big(2)^1022 # Covers gradual underflow or FTZ/DAZ on intermediates.
qexact(x::Float64)=QQ(x)
rational_upper(x)=qexact(sup(x))
rational_lower(x)=qexact(inf(x))
maxcomponent(A)=max(maximum(abs,real.(A);init=0.0),maximum(abs,imag.(A);init=0.0))

"Exact forward-error and polar/Weyl witness from observed Float64 component maxima."
function scalar_witness(dim,L,hmax,qmax,lmax,rmax,emax,herror)
    # Four real GEMMs and ordinary component arithmetic. The depth of each
    # dependency path is below 4*dim+16 elementary IEEE operations.
    steps=4*dim+16
    @assert steps*UNIT<1
    gamma=(steps*UNIT)/(1-steps*UNIT)
    underflow=(steps*TINY)/(1-steps*UNIT)
    # First certify ||Q'Q-I||<1 using only component maxima. This implies every
    # column has Euclidean norm <sqrt(2). Cauchy-Schwarz then bounds each absolute
    # Gram dot sum by 2, and each residual dot sum by 2||H0||. This second pass
    # removes the artificial dimension factor in the dot-product absolute sums.
    coarse_eerr=gamma*(2*dim*qmax^2+1)+underflow
    coarse_ebound=2*dim*(emax+coarse_eerr)
    @assert coarse_ebound<1
    hnorm=big(3)*L//4+herror
    rerr=gamma*(2*hnorm+qmax*lmax)+underflow
    eerr=3gamma+underflow
    # ||A||_2 <= ||A||_F <= 2*dim*max(|Re Aij|,|Im Aij|).
    rbound=2*dim*(rmax+rerr)
    ebound=2*dim*(emax+eerr)
    @assert ebound<1
    # Polar factor U: ||Q-U|| <= ||Q*Q-I||, ||Q||<2.
    # ||Htrue|| <= 3L/4 because all AF weights are positive and sum to L.
    radius=rbound+2*herror+(big(3)*L//4+lmax)*ebound
    return (;steps,gamma,underflow,rerr,eerr,rbound,ebound,radius,coarse_ebound)
end

function phases8()
    z=interval(0.0); one=interval(1.0); root=sqrt(interval(2.0))/2
    [complex(one,z),complex(root,root),complex(z,one),complex(-root,root),
        complex(-one,z),complex(-root,-root),complex(z,-one),complex(root,-root)]
end

"Construct a binary64 Hermitian center and an outward enclosure of its matrix error."
function verified_block(delta::QQ,L,q,k,parity,orbits)
    n=L÷2; @assert 8%n==0
    odd=parity==-1 ? 1 : 0
    allowed=[o for o in orbits if all(mod(2*k*t+n*odd*f,2*n)==0 for (t,f) in o.stabilizers)]
    dim=length(allowed); dim==0 && return nothing
    labels=zeros(Int,1<<L); phasecodes=zeros(Int,1<<L)
    sizes=[length(o.images) for o in allowed]
    for (i,o) in enumerate(allowed), (s,(t,f)) in o.images
        labels[s+1]=i
        phasecodes[s+1]=mod((8÷n)*k*t+4*odd*f,8)
    end
    phases=phases8()
    sizeset=unique(sizes)
    ratios=Dict((a,b)=>sqrt(interval(a)/interval(b)) for a in sizeset for b in sizeset)
    weights=[to_interval(1+delta),to_interval(1-delta)]
    H=zeros(ComplexF64,dim,dim); maxerr=0.0
    for (j,o) in enumerate(allowed)
        s=o.rep; entries=Dict{Int,CI}(); diag=interval(0.0)
        for site in 0:L-1
            next=mod(site+1,L); J=weights[iseven(site) ? 1 : 2]
            opposite=((s>>site)&1)!=((s>>next)&1)
            diag+=J*(opposite ? -1 : 1)/4
            if opposite
                target=xor(s,(1<<site)|(1<<next)); i=labels[target+1]
                if 0<i<=j
                    value=(J/2)*ratios[(sizes[j],sizes[i])]*conj(phases[phasecodes[target+1]+1])
                    entries[i]=get(entries,i,complex(interval(0.0),interval(0.0)))+value
                end
            end
        end
        entries[j]=get(entries,j,complex(interval(0.0),interval(0.0)))+diag
        for (i,value) in entries
            re=(inf(real(value))+sup(real(value)))/2
            imv=i==j ? 0.0 : (inf(imag(value))+sup(imag(value)))/2
            @assert isfinite(re)&&isfinite(imv)
            if i==j; @assert inf(imag(value))<=0<=sup(imag(value)); end
            H[i,j]=complex(re,imv); H[j,i]=conj(H[i,j])
            err=max(sup(abs(real(value)-interval(re))),sup(abs(imag(value)-interval(imv))))
            maxerr=max(maxerr,err)
        end
    end
    return H,2*dim*qexact(maxerr)
end

"Residual/orthogonality centers from four explicit real GEMMs, then rigorous bounds."
function enclose_block(H,herror,L)
    dim=size(H,1)
    result=eigen(Hermitian(H)) # Dense complete spectrum and eigenvectors.
    vals=result.values; vectors=result.vectors
    Hr=real.(H); Hi=imag.(H); Qr=real.(vectors); Qi=imag.(vectors)
    @assert all(A -> all(x -> iszero(x)||abs(x)>=floatmin(Float64),A),
        (Hr,Hi,Qr,Qi,vals))
    # No complex GEMM (and no fast-math); standard real dot-product bounds apply.
    Rr=Hr*Qr-Hi*Qi-Qr.*transpose(vals)
    Ri=Hr*Qi+Hi*Qr-Qi.*transpose(vals)
    Er=transpose(Qr)*Qr+transpose(Qi)*Qi-Matrix{Float64}(I,dim,dim)
    Ei=transpose(Qr)*Qi-transpose(Qi)*Qr
    hmax=qexact(max(maximum(abs,Hr),maximum(abs,Hi)))
    qmax=qexact(max(maximum(abs,Qr),maximum(abs,Qi)))
    lmax=qexact(maximum(abs,vals))
    rmax=qexact(max(maximum(abs,Rr),maximum(abs,Ri)))
    emax=qexact(max(maximum(abs,Er),maximum(abs,Ei)))
    @assert all(isfinite,vals) && qmax<=2 && hmax<=3L && lmax<=3L
    # No overflow in any product or partial sum; conservative bound far below 2^1023.
    @assert 4*dim*max(hmax,qmax,lmax,1)^2 < big(2)^100
    w=scalar_witness(dim,L,hmax,qmax,lmax,rmax,emax,herror)
    return Dict("dimension"=>dim,"h_component_max"=>string(hmax),
        "q_component_max"=>string(qmax),"lambda_max_abs"=>string(lmax),
        "raw_residual_component_max"=>string(rmax),"raw_orthogonality_component_max"=>string(emax),
        "matrix_error_upper"=>string(herror),"residual_norm_upper"=>string(w.rbound),
        "coarse_orthogonality_norm_upper"=>string(w.coarse_ebound),
        "orthogonality_norm_upper"=>string(w.ebound),"eigenvalue_radius"=>string(w.radius),
        "eigenvalues"=>string.(QQ.(vals))),vals,w.radius
end

"Verified complete ring spectrum using dense Sz/momentum/spin-flip blocks."
function verified_ring(delta::QQ,L)
    @assert 0<=delta<1
    @assert L in (4,8,16)
    config=sprint(show,BLAS.get_config())
    @assert occursin("openblas",lowercase(config))
    BLAS.set_num_threads(1)
    # Sanity check that the BLAS kernel does not flush subnormal results to zero.
    @assert only(fill(2.0^-600,1,1)*fill(2.0^-450,1,1))==2.0^-1050
    blocks=Any[]; centers=Float64[]; radii=QQ[]; n=L÷2; dimension=0
    seconds=@elapsed for q in 0:L÷2
        orbits,useflip=group_orbits(L,q); mult=useflip ? 1 : 2; qdim=0
        conjugates=Dict{Tuple{Int,Int},Tuple{Any,Any,Any}}()
        for k in 0:n-1, parity in (useflip ? (1,-1) : (1,))
            # Exact character conjugation: real H implies H_(n-k,p)=conj(H_(k,p)).
            # Compatible stabilizers and orbit ordering agree, so the entire spectrum
            # and its verified radius can be reused, including all multiplicities.
            if k>n÷2
                haskey(conjugates,(n-k,parity)) || continue
                prior,vals,radius=conjugates[(n-k,parity)]
                block=copy(prior); dim=block["dimension"]; runtime=0.0
                block["conjugate_momentum"]=n-k
            else
                built=verified_block(delta,L,q,k,parity,orbits)
                built===nothing && continue
                H,herror=built; dim=size(H,1)
                runtime=@elapsed block,vals,radius=enclose_block(H,herror,L)
                conjugates[(k,parity)]=(block,vals,radius)
            end
            qdim+=dim
            merge!(block,Dict("down_spins"=>q,"cell_momentum"=>k,
                "spin_flip_parity"=>(useflip ? parity : nothing),"multiplicity"=>mult,
                "seconds"=>runtime))
            push!(blocks,block); dimension+=mult*dim
            for _ in 1:mult
                append!(centers,vals); append!(radii,fill(radius,dim))
            end
            println("verified delta=$delta L=$L q=$q k=$k parity=$parity dim=$dim ",
                "radius=$(Float64(radius)) sec=$runtime"); flush(stdout)
        end
        @assert qdim==binomial(L,q)
    end
    @assert dimension==2^L==length(centers)
    return Dict("spins"=>L,"cells"=>L÷2,"dimension"=>dimension,
        "seconds"=>seconds,"blas_config"=>config,"blocks"=>blocks,
        "maximum_eigenvalue_radius"=>string(maximum(radii))),centers,radii
end

"Outward log partition enclosure from every paired eigenvalue interval."
function enclosed_logZ(centers,radii,beta::QQ)
    @assert length(centers)==length(radii)
    b=to_interval(beta)
    shift=minimum(centers) # Arbitrary exact binary64 energy shift, no claimed E0.
    total=interval(0.0)
    radius_cache=Dict(r=>sup(to_interval(r)) for r in unique(radii))
    for (center,radius) in zip(centers,radii)
        rr=radius_cache[radius]
        enclosure=interval(center)+interval(-rr,rr)
        total+=exp(-b*(enclosure-interval(shift)))
    end
    return -b*interval(shift)+log(total)
end
end
