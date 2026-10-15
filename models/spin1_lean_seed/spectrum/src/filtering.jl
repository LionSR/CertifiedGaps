include("common.jl")

"Floor dyadic multiply without overflowing Int128; checked against BigInt in tests."
@inline function mulfloor(x::Int128,c::Int128)
    h=x>>64
    r=x & ((Int128(1)<<64)-1)
    ((c*h)<<24)+((c*r)>>40)
end

"Three-term coefficients of the normalized Chebyshev filter damping [alpha, upper] relative to lo."
function cheb_coefficients(lo,alpha,upper,degree)
    center=(upper+alpha)/2; radius=(upper-alpha)/2
    x0=(lo-center)/radius
    @assert x0 < -1
    sigma=1/x0
    coeffs=[(sigma/radius,-sigma*center/radius,0.0)]
    for j in 2:degree
        next=1/(2x0-sigma)
        push!(coeffs,(2next/radius,-2next*center/radius,-next*sigma))
        sigma=next
    end
    coeffs
end

"Integer a+b*omega representation, omega=exp(2pi*i/3)."
function integer_hamiltonian(H)
    if eltype(H)<:Real
        return round.(Int64,H),nothing
    end
    A=round.(Int64,real(H).+imag(H)./sqrt(3))
    B=round.(Int64,2imag(H)./sqrt(3))
    A,B
end

"""
    fixed_filter(Hint, Q, coeffs, bits)

Apply the Chebyshev recurrence to the block `Q` in signed fixed point with `bits-16` fractional
bits (Q48 for 64, Q112 for 128) and 40-bit coefficients; returns the Float64 block and max lane.
"""
function fixed_filter(Hint,Q,coeffs,bits)
    A,B=Hint
    frac=bits-16
    scale=exp2(frac)
    if isnothing(B)
        u=round.(Int128,real(Q).*scale); v=nothing
    else
        u=round.(Int128,(real(Q).+imag(Q)./sqrt(3)).*scale)
        v=round.(Int128,2imag(Q)./sqrt(3).*scale)
    end
    oldu=zeros(Int128,size(u)); oldv=isnothing(v) ? nothing : zero(v)
    hu=similar(u); hv=isnothing(v) ? nothing : similar(v)
    tmp=isnothing(v) ? nothing : similar(v)
    newu=similar(u); newv=isnothing(v) ? nothing : similar(v)
    maxlane=Int128(0)
    for coeff in coeffs
        ca,cb,cc=round.(Int128,coeff.*exp2(40))
        mul!(hu,A,u)
        if !isnothing(B)
            mul!(tmp,B,v); hu.-=tmp
            mul!(hv,A,v); mul!(tmp,B,u); hv.+=tmp
            mul!(tmp,B,v); hv.-=tmp
        end
        @inbounds for i in eachindex(u)
            newu[i]=mulfloor(hu[i],ca)+mulfloor(u[i],cb)+mulfloor(oldu[i],cc)
            maxlane=max(maxlane,abs(newu[i]))
        end
        if !isnothing(v)
            @inbounds for i in eachindex(v)
                newv[i]=mulfloor(hv[i],ca)+mulfloor(v[i],cb)+mulfloor(oldv[i],cc)
                maxlane=max(maxlane,abs(newv[i]))
            end
            oldv,v,newv=v,newv,oldv
        end
        oldu,u,newu=u,newu,oldu
    end
    @assert maxlane < (Int128(1)<<(bits-2))
    R=Float64.(u)./scale
    !isnothing(v) && (R=R.+cis(2pi/3).*(Float64.(v)./scale))
    R,Float64(maxlane)/scale
end

"Float64 reference implementation of the same Chebyshev recurrence."
function float_filter(H,Q,coeffs)
    old=zero(Q); cur=copy(Q); next=similar(Q); hv=similar(Q)
    for (a,b,c) in coeffs
        mul!(hv,H,cur)
        @. next=a*hv+b*cur+c*old
        old,cur,next=cur,next,old
    end
    cur,maximum(abs,cur)
end

"Deterministic block of `p` distinct computational-basis states from a low-diagonal pool."
function product_seeds(H,p,seed)
    rng=MersenneTwister(seed)
    # Independent coordinate states: random tie-breaking among low-diagonal states.
    ids=sortperm(real(diag(H)).+rand(rng,size(H,1)))
    # Broaden the pool to avoid concentrating on one symmetry orbit.
    pool=ids[1:min(length(ids),max(8p,div(length(ids),3)))]
    chosen=shuffle(rng,pool)[1:p]
    Q=zeros(eltype(H),size(H,1),p)
    for (j,i) in enumerate(chosen); Q[i,j]=1; end
    Q,chosen
end

"""
    filter_sector(n, twist, m, target; degree=32, bits=0, maxcycles=16)

Run filter cycles (filter, QR, Rayleigh-Ritz) in one sector until the linear residual target or
a plateau; writes `filter_s<source>_<n>_<twist>_<m>_<degree>_<bits>.json` unless it exists.
"""
function filter_sector(n,twist,m,target;degree=32,bits=0,maxcycles=16)
    source=target.beta==12.25 ? 1 : 2
    path=joinpath(DAT,"filter_s$(source)_$(n)_$(twist)_$(m)_$(degree)_$(bits).json")
    isfile(path) && return
    H,_=hamiltonian(n,m,twist); d=size(H,1); k=Int(target.k)
    k==0 && return
    p=min(d,k+8)
    F=low_spectrum(H;initial=min(p+2,d-2),dense=d<=160)
    ev=F.values
    alpha=p==d ? n-0.1 : (ev[p]+ev[p+1])/2
    lo=first(ev)-0.05
    coeffs=cheb_coefficients(lo,alpha,Float64(n),degree)
    Hint=bits==0 ? nothing : integer_hamiltonian(H)
    Q,seeds=product_seeds(H,p,1000n+10m+findfirst(==(twist),TWISTS))
    rows=[]; started=time(); hvtotal=0
    for cycle in 1:maxcycles
        Y,maxlane=bits==0 ? float_filter(H,Q,coeffs) : fixed_filter(Hint,Q,coeffs,bits)
        hvtotal+=degree*p
        qrY=qr(Y)
        # Materialize only the thin Q, never the d-by-d orthogonal matrix.
        Q=Matrix(qrY.Q*Matrix{eltype(H)}(I,d,p))
        HQ=H*Q; hvtotal+=p
        rr=eigen(Hermitian(Q'*HQ))
        Q=Q*rr.vectors; HQ=HQ*rr.vectors
        R=HQ[:,1:k]-Q[:,1:k]*Diagonal(rr.values[1:k])
        residual=opnorm(R)
        eigerror=maximum(abs.(rr.values[1:k]-ev[1:k]))
        push!(rows,(;cycle,total_degree=cycle*degree,hv=hvtotal,residual,eigerror,
                    linear_pass=residual<=target.residual_linear,
                    quadratic_pass=residual<=target.residual_quadratic,
                    maxlane,seconds=time()-started))
        residual<=target.residual_linear && break
        # A measured rounding plateau is a failure, not an extrapolated success.
        if cycle>=8 && minimum(r.residual for r in rows[end-3:end]) >=
                        0.8minimum(r.residual for r in rows[end-7:end-4])
            break
        end
    end
    writejson(path,Dict("status"=>"NON-CERTIFYING; integer recurrence, Float64 QR and residuals",
        "source"=>source,"n"=>n,"twist"=>twist,"m"=>m,"degree"=>degree,"bits"=>bits,
        "fraction_bits"=>bits==0 ? 0 : bits-16,"coefficient_bits"=>40,
        "k"=>k,"p"=>p,"dimension"=>d,"seed_indices"=>seeds,
        "lo"=>lo,"alpha"=>alpha,"upper"=>n,"coefficients"=>coeffs,
        "target"=>target,"rows"=>rows,"maxcycles"=>maxcycles,
        "stopping_reason"=>last(rows).linear_pass ? "linear_target" : "plateau_or_cap"))
    println("filter n=$n g=$twist m=$m deg=$degree bits=$bits k=$k p=$p ",
            "cycles=$(length(rows)) r=$(last(rows).residual)")
    flush(stdout)
end

"Run `filter_sector` for every active (k > 0) sector of the two sources at the given sizes."
function run_filters(ns; degrees=(16,32,64),bitslist=(0,))
    analysis=JSON3.read(read(joinpath(DAT,"analysis.json"),String))
    for t in analysis.precision
        n=Int(t.n); twist=String(t.twist)
        n in ns || continue
        relevant=(t.beta==12.25 && twist!="C") || (t.beta==10.5 && n in (6,8,10))
        relevant && t.budget==3e-5 && t.k>0 || continue
        for degree in degrees,bits in bitslist
            filter_sector(n,twist,Int(t.m),t;degree,bits)
        end
    end
end

