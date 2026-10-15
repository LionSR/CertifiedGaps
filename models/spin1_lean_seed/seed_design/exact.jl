# Exact scalar verifier. The finite physical inputs and analytic lemmas remain hypotheses
# until separately proved in Lean; no floating-point comparison accepts an inequality here.
include("common.jl")

function exact_sector_data(source,mode,eps)
    table=stored_table(source)
    weights=mode=="ION" ? [frac.([1,3,8],12),frac.([1,3,-4],12),frac.([1,-1,0],4)] :
        [frac.([1,3],4),frac.([1,-1],4)]
    ns=sort(unique(first(k) for k in keys(table)))
    moments=[Dict(n=>sum(w[g+1]*table[(n,g)] for g in 0:length(w)-1)
        for n in ns) for w in weights]
    errors=[eps*sum(abs,w) for w in weights]
    moments,errors
end

function filter_valid(P,B,v)
    for sign in (-1,1)
        shifted=[sum(binomial(i,j)*P[i+1]*sign^i*v^(i-j) for i in j:length(P)-1)
            for j in 0:length(P)-1]
        first(shifted)>B && all(>=(0),shifted[2:end]) || return false
    end
    true
end

"Recompute the filter bound and check the full exterior rays with exact rationals."
function exact_filter(mom,eps,shift,q; digits=10)
    @assert iseven(shift) && shift>=4
    P=convsquare(q,shift)
    B=sum(c*mom[i-1] for (i,c) in enumerate(P) if c!=0)+eps*sum(abs,P)
    @assert B>0
    _,cap=first_true(v->filter_valid(P,B,v),frac(1,20),frac(2);digits)
    @assert filter_valid(P,B,cap)
    Dict("shift"=>shift,"q"=>string.(q),"P"=>string.(P),"B"=>string(B),"cap"=>string(cap))
end

function capped_upper(m,cap,r)
    @assert m>=0 && cap>0 && r>=1
    count=fld(m,cap)
    count*power_upper(cap,r)+power_upper(m-count*cap,r)
end

function exact_bounds(input,k)
    nmax=input["nmax"];r=frac(k,nmax)
    masses=parse_rational.(input["masses"])
    bs=RationalBig[]
    for (m,f) in zip(masses,input["filters"])
        push!(bs,isnothing(f) ? power_upper(m,r) :
            capped_upper(m,parse_rational(f["cap"])^nmax,r))
    end
    mult=input["mode"]=="ION" ? [1,2,3] : [1,3]
    sum(mult.*bs),input["mode"]=="ION" ? bs[1]+2bs[2] : bs[1]
end

function exact_residual(D,tauj,taui,t,s,mode)
    @assert 0<t<=taui<=tauj<=D && s>=1
    vertices=mode=="ION" ? ((frac(0),frac(0)),(taui-t,frac(0)),
        (taui-t,tauj-taui),(frac(0),tauj-t)) :
        ((frac(0),frac(0)),(tauj-t,frac(0)))
    maximum(power_upper(x,s)+2power_upper(y/2,s)+3power_upper((D-t-x-y)/3,s)
        for (x,y) in vertices)
end

function exact_seed(input,trial,k,r)
    @assert trial%4==0 && iseven(k) && input["nmax"]<=k<=trial÷2 && r>=1
    M,E=exact_bounds(input,k)
    # Compact upper witnesses avoid carrying thousand-digit powers into every later step.
    witness_grid=big(10)^12
    D,E=round_up(power_upper(M,r),witness_grid),round_up(power_upper(E,r),witness_grid)
    h=frac(trial,2k)
    tauj=(D+3E)/4
    taui=input["mode"]=="ION" ? (D+11E)/12 : tauj
    t=round_down(first(power_interval(inv(D),inv(2h-1))),witness_grid)-inv(frac(witness_grid))
    multiplicity=input["mode"]=="ION" ? 2 : 3
    @assert 0<t<taui && multiplicity*t>D
    @assert D*power_upper(t,2h-1)<1
    initial=t
    refinements=Any[]
    for _ in 1:35
        residual=exact_residual(D,tauj,taui,t,2h,input["mode"])
        @assert residual<1
        next=round_down(first(power_interval(1-residual,inv(2h))),witness_grid)-
            inv(frac(witness_grid))
        next<=t && break
        @assert next<taui && power_upper(next,2h)+residual<1
        push!(refinements,Dict("old"=>string(t),"new"=>string(next),"residual"=>string(residual)))
        t=next
    end
    res=exact_residual(D,tauj,taui,t,h,input["mode"])
    dom=first(power_interval(t,h))
    p=1-(1+res/dom)^(-2)
    Q,_=exact_bounds(input,trial)
    Q=round_up(Q,witness_grid)
    @assert Q>1
    ratio=round_up(power_upper(Q-1,r),big(10)^18)
    q=1-(1+ratio)^(-2)
    Dict("trial"=>trial,"n0"=>trial÷2,"k"=>k,"r"=>string(r),"h"=>string(h),
        "D"=>string(D),"E"=>string(E),"tauj"=>string(tauj),"taui"=>string(taui),
        "initial_t"=>string(initial),"refinements"=>refinements,"t"=>string(t),
        "Q"=>string(Q),"spatial_residual"=>string(res),"dominant_power"=>string(dom),
        "physical_ratio"=>string(ratio),"p"=>string(p),"q"=>string(q),
        "diagonal_threshold_passes"=>max(p,q)<frac(12731,100000))
end

zq(a,b)=1-(1-b)*(1-a)^2
fq(z)=z^2/(2*(1-z)^2)

"Exact replay of ChainCertificate.Valid, including interpolation and C=5 tail."
function chain_valid(c)
    s,t=parse_rational.(c["s"]),parse_rational.(c["t"])
    r,epsilon=parse_rational(c["r"]),parse_rational(c["epsilon"])
    J=length(s)-1
    length(t)==J+1 && r>0 && 0<epsilon<=frac(1,100) || return false
    s[end]<=epsilon && t[end]<=epsilon && r^(2^J)*(5epsilon)<=1 || return false
    all(x->0<=x<=1,vcat(s,t)) || return false
    for j in 1:J
        z=zq(s[j],t[j]);zt=zq(t[j],s[j+1])
        z<1 && zt<1 && fq(z)<=s[j+1] && fq(zt)<=t[j+1] || return false
        r^(2^(j-1))*z<1-z || return false
    end
    true
end

function exact_chain(p,q,beta;steps=9)
    grid=big(10)^512
    s,t=[p],[q]
    for _ in 1:steps
        z=zq(last(s),last(t));@assert 0<z<frac(1,2)
        p=round_up(fq(z),grid)
        zt=zq(last(t),p);@assert 0<zt<1
        q=round_up(fq(zt),grid)
        push!(s,p);push!(t,q)
    end
    epsilon=max(last(s),last(t));@assert epsilon<=frac(1,100)
    validrate(r)=r^(2^steps)*(5epsilon)<=1 &&
        all(r^(2^(j-1))*zq(s[j],t[j])<1-zq(s[j],t[j]) for j in 1:steps)
    @assert validrate(frac(1)) && !validrate(frac(10))
    r,_=first_true(x->!validrate(x),frac(1),frac(10);digits=8)
    @assert r>1
    c=Dict("s"=>string.(s),"t"=>string.(t),"r"=>string(r),"epsilon"=>string(epsilon))
    @assert chain_valid(c)
    lower=first(log_enclosure(r))/beta
    claim=round_down(lower,big(10)^9)-frac(1,big(10)^9)
    @assert 0<claim<lower
    c["beta"]=string(beta);c["gap_strict_lower"]=string(claim)
    c["gap_display"]=Float64(claim)
    c
end

function exact_input(source,mode,nmax,eps,proposals;uncappedO=false)
    moments,errors=exact_sector_data(source,mode,eps)
    fs=Any[]
    for (j,(m,e,p)) in enumerate(zip(moments,errors,proposals))
        if uncappedO && j==2
            push!(fs,nothing);continue
        end
        # Floats propose coefficients only. The following proof obligations are exact.
        q=[frac(round(Int,x*10^6),10^6) for x in p["q"]]
        push!(fs,exact_filter(m,e,Int(p["shift"]),q))
    end
    Dict("source"=>source,"mode"=>mode,"nmax"=>nmax,"eps"=>string(eps),"filters"=>fs,
        "masses"=>[string(m[nmax]+e) for (m,e) in zip(moments,errors)],
        "base_beta"=>source==1 ? "49//4" : "21//2")
end

function verify_trial(trial)
    if trial in (72,120)
        path=joinpath(RESULTS,"trials.json")
        row=only(x for x in read_exact_json(path)["certificates"] if parse(Int,x["length"])==trial)
    else
        path=joinpath(@__DIR__,"trial_scan.json")
        row=read_exact_json(path)["exact"]
        @assert parse(Int,row["length"])==trial
    end
    u,v=parse(BigInt,row["U"]),parse(BigInt,row["V"])
    @assert v>0 && 500000u+700741v<0
    Dict("path"=>relpath(path,@__DIR__),"sha256"=>sha256_file(path),"length"=>trial,
        "energy_per_bond"=>string(u//v),"positive_margin"=>string(-500000u-700741v),
        "scope"=>trial in (72,120) ? "stored exact trial output; no fresh contraction" :
            "fresh Julia integer contraction")
end
