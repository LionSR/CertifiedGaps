"""
    AnchorStudy

Study `anchor`: exact counting anchors for the spin-1 seed from local PSD comparisons. Ring
spectra, scans and the continuation are NON-CERTIFYING Float64; `certify_anchor` and
`check_anchor` recompute the exact records (local PSD pivots, half-chain enclosures, comparison
counts, binned tails) in Rational{BigInt} / bounded Int128 arithmetic.
"""
module AnchorStudy

using LinearAlgebra, SparseArrays, Arpack, Optim, Random, Serialization
using CertifiedGaps.Results
using ..Common: MODEL, results_dir, study_hashes, sz_basis, export_all
import Arpack, Optim

BLAS.set_num_threads(1)
# ROOT is the model directory; this study's results live in results/anchor/.
const ROOT = MODEL
const RESULTS = results_dir("anchor")
const BETA = 49 / 4
const SHIFT = 700741 / 500000
const TOLERANCE = 3e-5
num(x) = x isa AbstractString ? parse(Float64,x) : Float64(x)
ival(x) = x isa AbstractString ? parse(Int,x) : Int(x)

"Sz product basis, with digit 0,1,2 representing m=-1,0,1."
basis(n, m) = sz_basis(n, m)

"Integer spin-one bond. A pi rotation around z changes only the exchange sign."
function bond(n, states, i, j; twist=1)
    lookup = Dict(x => r for (r, x) in enumerate(states))
    rows, cols, vals = Int[], Int[], Int[]
    pi, pj = 3^(i - 1), 3^(j - 1)
    for (r, x) in enumerate(states)
        a, b = (x ÷ pi) % 3 - 1, (x ÷ pj) % 3 - 1
        if a * b != 0
            push!(rows, r); push!(cols, r); push!(vals, a * b)
        end
        for δ in (-1, 1)
            abs(a + δ) <= 1 && abs(b - δ) <= 1 || continue
            push!(rows, lookup[x + δ * (pi - pj)])
            push!(cols, r); push!(vals, twist)
        end
    end
    return sparse(rows, cols, vals, length(states), length(states))
end

"Open chain and two endpoint-bond correction matrices, all integer."
function block_matrices(n, m)
    b = basis(n, m)
    h = [bond(n, b, i, i + 1) for i in 1:n-1]
    return (H=Matrix(sum(h)), L=Matrix(first(h) + last(h)),
        Q=Matrix(first(h)^2 + last(h)^2), dim=length(b))
end

const BLOCKS = Dict{Tuple{Int,Int},Any}()
function block(n, m)
    return get!(BLOCKS, (n, abs(m))) do
        block_matrices(n, abs(m))
    end
end

"Four-site comparison remainder h23 - u(h12+h34) - v(h12²+h34²)."
function remainder_matrices(m)
    b = basis(4, m)
    h = [bond(4, b, i, i+1) for i in 1:3]
    return (H=Matrix(h[2]), L=Matrix(h[1]+h[3]), Q=Matrix(h[1]^2+h[3]^2))
end
const LOCAL = [remainder_matrices(m) for m in 0:4]

local_min(u, v) = minimum(eigmin(Symmetric(a.H-u*a.L-v*a.Q)) for a in LOCAL)

function block_spectra(n, u, v)
    return Dict(m => eigvals(Symmetric(block(n,m).H + u*block(n,m).L +
        v*block(n,m).Q)) for m in 0:n)
end

"Tensor-sum spectrum, retaining Sz and all multiplicities."
function comparison_spectra(lengths, u, v; c=local_min(u,v))
    current = Dict(0 => [length(lengths)*c])
    for n in lengths
        spec = block_spectra(n,u,v)
        next = Dict{Int,Vector{Float64}}()
        for (m, es) in current, q in -n:n
            out = get!(next,m+q,Float64[])
            for x in es, y in spec[abs(q)]
                push!(out,x+y)
            end
        end
        current = next
    end
    return Dict(m => sort(es) for (m,es) in current if m >= 0)
end

function ring(n,m; twist=1)
    b = basis(n,m)
    H = sum(bond(n,b,i,i+1) for i in 1:n-1)
    return H + bond(n,b,n,1;twist)
end

"The two-block comparison, embedded in the same sector as the ring."
function comparison_matrix(n,m,u,v,c)
    b=basis(n,m); l=n÷2
    h=[bond(n,b,i,i+1) for i in 1:n-1]
    H0=sum(h[i] for i in 1:n-1 if i!=l)
    ends=(1,l-1,l+1,n-1)
    return H0 + u*sum(h[i] for i in ends) + v*sum(h[i]^2 for i in ends) + 2c*I
end

"NON-CERTIFYING lowest spectrum, extended until it crosses the requested cutoff."
function low_spectrum(H, cutoff; minimum_count=12, vectors=false)
    d = size(H,1)
    if d <= 1300
        F = eigen(Symmetric(Matrix{Float64}(H)))
        k = min(d,max(minimum_count,searchsortedlast(F.values,cutoff)+4))
        r = maximum(norm(H*F.vectors[:,i]-F.values[i]*F.vectors[:,i]) for i in 1:k)
        return (values=F.values[1:k], vectors=vectors ? F.vectors[:,1:k] : nothing,
            residual=r, complete=true, iterations=0)
    end
    k = min(max(minimum_count,20),d-2)
    Random.seed!(20261014 + d)
    while true
        vals,V,nconv,niter,nmult,resid = eigs(Float64.(H); nev=k, which=:SR,
            ncv=min(d,max(2k+12,80)), tol=2e-11, maxiter=10000,
            v0=randn(d))
        p = sortperm(real(vals)); vals=real(vals[p]); V=real(V[:,p])
        @assert nconv == k
        r = maximum(norm(H*V[:,i]-vals[i]*V[:,i]) for i in 1:k)
        if vals[end-3] > cutoff || k == d-2
            return (values=vals,vectors=vectors ? V : nothing,residual=r,
                complete=vals[end-3]>cutoff,iterations=niter)
        end
        k = min(2k,d-2)
    end
end

function metadata(config; status="NON-CERTIFYING Float64 experiment")
    hashes=study_hashes("anchor")
    return run_metadata(config;source_sha256=hashes,packages=(Arpack,Optim),
        status,certified=false)
end
"Exact absolute shifted-trace cutoff using a single uniform tail energy."
tail_cutoff(n; budget=TOLERANCE) = log(3.0^n / budget) / BETA - n*SHIFT

# ---------------------------------------------------------------------------------------------
# Exact finite-matrix checks, not Lean kernel certificates.
# ---------------------------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------------------------
# Certificate writer and checker. Proposals may come from floating eigensolves; every accepted
# finite inequality is recomputed with Rational{BigInt} or bounded Int128 arithmetic.
# ---------------------------------------------------------------------------------------------

const TARGETS=Dict(8=>-9644246//big(10)^6,10=>-12267845//big(10)^6,
    12=>-14891444//big(10)^6)

"""
    exact_counts(blocks, n, c; rho)

Exact per-Sz count of two-half comparison levels below `rho`, from the stored half-chain
enclosures on a 1e-7 integer grid. Asserts that lower and upper endpoints give the same count,
and records the next comparison level as the counting anchor.
"""
function exact_counts(blocks,n,c;rho=TARGETS[n])
    scale=big(10)^7
    lower=Dict(parse(Int,x["Sz"])=>BigInt.(rat.(x["lower"]).*scale) for x in blocks)
    upper=Dict(parse(Int,x["Sz"])=>BigInt.(rat.(x["upper"]).*scale) for x in blocks)
    shift=BigInt(2c*scale); target=BigInt(rho*scale); l=n÷2
    records=Any[]
    for m in 0:n
        lows=BigInt[]; highs=BigInt[]
        for q in -l:l
            abs(m-q)<=l || continue
            append!(lows,[x+y+shift for x in lower[abs(q)],y in lower[abs(m-q)]])
            append!(highs,[x+y+shift for x in upper[abs(q)],y in upper[abs(m-q)]])
        end
        sort!(lows);sort!(highs)
        klo=count(<(target),highs);khi=count(<(target),lows)
        @assert klo==khi "enclosures do not resolve the count"
        @assert length(lows)==length(basis(n,m))
        @assert khi<length(lows)
        push!(records,Dict("Sz"=>string(m),"k"=>string(khi),
            "rho"=>ratstring(lows[khi+1]//scale),
            "target"=>ratstring(rho),"dimension"=>string(length(lows))))
    end
    return records
end

"""
    exact_binned_tail(blocks, n, c)

Exact binned upper bound on the comparison-spectrum thermal tail above the target `rho`
(32 bins of width 1/8 plus an overflow bin, rational `exp(-x)` upper bounds). Asserts the
total is below the 1.5e-5 budget.
"""
function exact_binned_tail(blocks,n,c)
    targets=Dict(8=>-9925//big(1000),10=>-12687//big(1000),12=>-15453//big(1000))
    rho=targets[n];scale=big(10)^7;l=n÷2
    lower=Dict(parse(Int,x["Sz"])=>BigInt.(rat.(x["lower"]).*scale) for x in blocks)
    masses=zeros(BigInt,33)
    for m in 0:n
        es=BigInt[]
        for q in -l:l
            abs(m-q)<=l || continue
            append!(es,[x+y+BigInt(2c*scale) for x in lower[abs(q)],y in lower[abs(m-q)]])
        end
        sort!(es);last=0;mult=m==0 ? 1 : 2
        for j in 0:31
            k=searchsortedlast(es,BigInt((rho+(j+1)//8)*scale))
            masses[j+1]+=mult*(k-last);last=k
        end
        masses[33]+=mult*(length(es)-last)
    end
    @assert sum(masses)==big(3)^n
    uppers=[exp_negative_upper(49//4*(rho+j//8+n*(700741//big(500000)))) for j in 0:32]
    tail=sum(masses.*uppers)
    @assert tail<3//big(200000)
    return Dict("target"=>ratstring(rho),"bin_width"=>"1//8","masses"=>string.(masses),
        "exponential_uppers"=>ratstring.(uppers),"tail_upper"=>ratstring(tail),
        "budget"=>"3//200000","taylor_degree"=>"256",
        "sector_anchors"=>exact_counts(blocks,n,c;rho))
end

"""
    check_local(u, v, c)

Exact LDL pivots proving the four-site remainder `h23 - B12 - B34 - c` is PSD in every Sz
sector, plus the pi-twist and spin-flip conjugation identities.
"""
function check_local(u,v,c)
    records=Any[]
    for m in 0:4
        X=LOCAL[m+1]
        R=QQ.(X.H)-u*X.L-v*X.Q-c*I
        fac=exact_psd(R)
        @assert fac.ok
        # The pi-twisted bond is integer-diagonal unitarily equivalent.
        b=basis(4,m)
        signs=[isodd(((x÷9)%3-1)+((x÷27)%3-1)) ? -1 : 1 for x in b]
        Hp=Matrix(bond(4,b,2,3;twist=-1))
        Rp=QQ.(Hp)-u*X.L-v*X.Q-c*I
        @assert Rp == Diagonal(signs)*R*Diagonal(signs)
        Xneg=remainder_matrices(-m)
        Rneg=QQ.(Xneg.H)-u*Xneg.L-v*Xneg.Q-c*I
        @assert Rneg==R[end:-1:1,end:-1:1]
        push!(records,Dict("Sz"=>string(m),"dimension"=>string(length(b)),
            "pivots"=>ratstring.(fac.pivots),"twist_identity"=>true))
    end
    return records
end

"Recompute every exact record and write `results/anchor/certificate.json`."
function certify_anchor()
    opt=read_exact_json(joinpath(RESULTS,"optimized.json"))["rows"]
    params=[Dict("n"=>r["n"],"u"=>r["u"],"v"=>r["v"],"c"=>r["c"],
        "kind"=>"ground_optimized") for r in opt]
    countpath=joinpath(RESULTS,"count_optimized.json")
    if isfile(countpath)
        append!(params,read_exact_json(countpath)["rows"])
    end
    records=Any[]
    for p in params
        n=parse(Int,string(p["n"]))
        u=round(BigInt,num(p["u"])*10000)//big(10000)
        v=round(BigInt,num(p["v"])*10000)//big(10000)
        c=round(BigInt,num(p["c"])*10^7)//big(10)^7
        localrecords=check_local(u,v,c)
        blocks=Any[]
        for m in 0:n÷2
            X=block(n÷2,m)
            mat=QQ.(X.H)+u*X.L+v*X.Q
            t=time(); cert=exact_spectrum_enclosure(mat); seconds=time()-t
            lower=[floor(BigInt,x*10^6)//big(10)^6 for x in cert.lower]
            upper=[ceil(BigInt,x*10^6)//big(10)^6 for x in cert.upper]
            push!(blocks,Dict("Sz"=>string(m),"dimension"=>string(cert.dimension),
                "lower"=>ratstring.(lower),"upper"=>ratstring.(upper),
                "radius"=>ratstring(cert.radius),"seconds"=>seconds,
                "orthogonality"=>ratstring(cert.orthogonality),
                "Int128_bound"=>string(cert.integer_bound)))
            println("EXACT n=$n $(p["kind"]) m=$m dim=$(cert.dimension) "*
                "radius=$(Float64(cert.radius)) seconds=$seconds")
            flush(stdout)
        end
        result=Dict("n"=>string(n),"kind"=>p["kind"],"u"=>ratstring(u),
            "v"=>ratstring(v),"c"=>ratstring(c),"local_psd"=>localrecords,
            "blocks"=>blocks,"sector_anchors"=>exact_counts(blocks,n,c))
        p["kind"]=="count_optimized" && (result["binned_tail"]=exact_binned_tail(blocks,n,c))
        push!(records,result)
    end
    write_result(joinpath(RESULTS,"certificate.json"),Dict(
        "metadata"=>metadata(Dict("scope"=>"local PSD and half-chain enclosures");
            status="EXACT Julia checks; not a Lean certificate or gap proof"),
        "parameters"=>records))
end
"Recompute every certificate record exactly (stored enclosures must contain fresh ones)."
function check_anchor(path=joinpath(RESULTS,"certificate.json"))
    data=read_exact_json(path)
    for p in data["parameters"]
        n=parse(Int,p["n"]);u=rat(p["u"]);v=rat(p["v"]);c=rat(p["c"])
        localrows=check_local(u,v,c)
        @assert localrows==p["local_psd"]
        for old in p["blocks"]
            m=parse(Int,old["Sz"])
            X=block(n÷2,m)
            mat=QQ.(X.H)+u*X.L+v*X.Q
            Xneg=block_matrices(n÷2,-m)
            @assert QQ.(Xneg.H)+u*Xneg.L+v*Xneg.Q==mat[end:-1:1,end:-1:1]
            cert=exact_spectrum_enclosure(mat)
            @assert all(rat.(old["lower"]).<=cert.lower)
            @assert all(cert.upper.<=rat.(old["upper"]))
        end
        @assert exact_counts(p["blocks"],n,c)==p["sector_anchors"]
        if haskey(p,"binned_tail")
            @assert exact_binned_tail(p["blocks"],n,c)==p["binned_tail"]
        end
        println("PASS exact n=$n $(p["kind"]): local PSD, twists, half spectra, counts")
    end
    println("PASS all finite checks. No Lean kernel certification is claimed.")
end

export_all(@__MODULE__)

end
