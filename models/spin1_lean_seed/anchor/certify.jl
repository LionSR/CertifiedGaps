include("src/common.jl")
include("src/exact.jl")
using .AnchorStudy, LinearAlgebra
using CertifiedGaps.Results
const A=AnchorStudy
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
        @assert length(lows)==length(A.basis(n,m))
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
        X=A.LOCAL[m+1]
        R=QQ.(X.H)-u*X.L-v*X.Q-c*I
        fac=exact_psd(R)
        @assert fac.ok
        # The pi-twisted bond is integer-diagonal unitarily equivalent.
        b=A.basis(4,m)
        signs=[isodd(((x÷9)%3-1)+((x÷27)%3-1)) ? -1 : 1 for x in b]
        Hp=Matrix(A.bond(4,b,2,3;twist=-1))
        Rp=QQ.(Hp)-u*X.L-v*X.Q-c*I
        @assert Rp == Diagonal(signs)*R*Diagonal(signs)
        Xneg=A.remainder_matrices(-m)
        Rneg=QQ.(Xneg.H)-u*Xneg.L-v*Xneg.Q-c*I
        @assert Rneg==R[end:-1:1,end:-1:1]
        push!(records,Dict("Sz"=>string(m),"dimension"=>string(length(b)),
            "pivots"=>ratstring.(fac.pivots),"twist_identity"=>true))
    end
    return records
end

"Recompute every exact record and write `results/certificate.json`."
function run_certify()
    opt=read_exact_json(joinpath(A.RESULTS,"optimized.json"))["rows"]
    params=[Dict("n"=>r["n"],"u"=>r["u"],"v"=>r["v"],"c"=>r["c"],
        "kind"=>"ground_optimized") for r in opt]
    countpath=joinpath(A.RESULTS,"count_optimized.json")
    if isfile(countpath)
        append!(params,read_exact_json(countpath)["rows"])
    end
    records=Any[]
    for p in params
        n=parse(Int,string(p["n"]))
        u=round(BigInt,A.num(p["u"])*10000)//big(10000)
        v=round(BigInt,A.num(p["v"])*10000)//big(10000)
        c=round(BigInt,A.num(p["c"])*10^7)//big(10)^7
        localrecords=check_local(u,v,c)
        blocks=Any[]
        for m in 0:n÷2
            X=A.block(n÷2,m)
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
    write_result(joinpath(A.RESULTS,"certificate.json"),Dict(
        "metadata"=>A.metadata(Dict("scope"=>"local PSD and half-chain enclosures");
            status="EXACT Julia checks; not a Lean certificate or gap proof"),
        "parameters"=>records))
end
abspath(PROGRAM_FILE)==(@__FILE__) && run_certify()
