# Exact finite-matrix checks, not Lean kernel certificates.
using LinearAlgebra

const QQ = Rational{BigInt}
rat(x::AbstractString) = parse(BigInt,split(x,"//")[1]) //
    parse(BigInt,split(x,"//")[2])
ratstring(x::Rational) = "$(numerator(x))//$(denominator(x))"

"Exact rational LDL check; a zero pivot is allowed only with a zero remaining row."
function exact_psd(A)
    M=QQ.(A)
    @assert M == transpose(M)
    pivots=QQ[]
    for k in axes(M,1)
        p=M[k,k]
        p >= 0 || return (ok=false,pivots=pivots)
        push!(pivots,p)
        if p == 0
            all(iszero,M[k,k+1:end]) || return (ok=false,pivots=pivots)
            continue
        end
        for j in k+1:size(M,1), i in j:size(M,1)
            M[i,j] -= M[i,k]*M[j,k]/p
            M[j,i] = M[i,j]
        end
    end
    return (ok=true,pivots=pivots)
end

"Certify an entire half-chain spectrum by integer reconstruction and polar decomposition."
function exact_spectrum_enclosure(A)
    A=QQ.(A)
    den=foldl(lcm,denominator.(A);init=big(1))
    Ai=Int128.(numerator.(A.*den))
    F=eigen(Symmetric(Float64.(A))) # untrusted proposal only
    s=Int128(1)<<32
    Q=round.(Int128,F.vectors*s)
    ds=round.(Int128,F.values*s)
    d=size(A,1)
    # Prove every product and every partial sum in the Int128 calculations fits.
    qmax=maximum(abs,Q); amax=maximum(abs,Ai); dmax=maximum(abs,ds)
    bound=big(d)*qmax^2*max(big(1),dmax)*den + big(s)^3*amax + big(s)^2
    @assert bound < typemax(Int128)
    G=Q*transpose(Q)
    G -= s^2*I
    δ=maximum(sum(abs,BigInt.(G);dims=2))//big(s)^2
    @assert δ < 1//2
    T=Q*Diagonal(ds)*transpose(Q)
    E=s^3*Ai-Int128(den)*T
    err=maximum(sum(abs,BigInt.(E);dims=2))//(big(s)^3*den)
    M=big(dmax)//big(s)
    radius=err+3M*δ
    vals=[big(v)//big(s) for v in ds]
    return (lower=vals.-radius,upper=vals.+radius,radius=radius,
        orthogonality=δ,dimension=d,integer_bound=bound,
        vector_bytes=sizeof(Q),values=vals)
end

"Rational upper bound on exp(-x), x>=0, from a positive Taylor lower bound on exp(x)."
function exp_negative_upper(x::Rational;degree=256)
    @assert x>=0
    term=QQ(1); total=term
    for k in 1:degree
        term*=x/k
        total+=term
    end
    den=big(10)^24
    return ceil(BigInt,den/total)//den
end
