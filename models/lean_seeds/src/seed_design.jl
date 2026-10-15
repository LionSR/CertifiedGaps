"""
    SeedDesign

Study `seed_design`: which finite thermal traces does the spin-1 seed need? Float64 search
utilities (NON-CERTIFYING) propose filters and seeds; the exact part replays three seed designs
in rational arithmetic down to `ChainCertificate.Valid` (`certify_seed_design`,
`check_seed_design`).
"""
module SeedDesign

using LinearAlgebra, JSON3, SHA, Dates, TOML
using CertifiedGaps.ExactArithmetic, CertifiedGaps.Results
using ..Common: Common, MODEL, REPO, results_dir, study_hashes, export_all
# Paths: ROOT is the repository, PROJECT the model directory. Generated data live under
# results/seed_design/ (raw search rows in its raw/, gitignored); spectra under cache/seed_design/
# and paper downloads under upstream/ (both gitignored).
const PROJECT = MODEL
const ROOT = REPO
const OUT = results_dir("seed_design")
const RAW = joinpath(OUT, "raw")
const CACHE = joinpath(MODEL, "cache", "seed_design")
const UPSTREAM = Common.UPSTREAM
const SPIN1_RESULTS = Common.SPIN1_RESULTS
const SHIFT = 700741 / 500000
const RATSHIFT = frac(700741, 500000)

const sqdef(x) = x * (2 + x) / (1 + x)^2

function stored_table(source)
    raw = read_exact_json(joinpath(SPIN1_RESULTS, "source$source.json"))
    Dict((parse(Int, row["n"]), parse(Int, row["twist"])) =>
        parse_rational(row["center"]) for row in raw["certificates"])
end

function sector_table(table, mode)
    weights = mode == "ION" ? [[1,3,8]./12, [1,3,-4]./12, [1,-1,0]./4] :
        mode == "JN" ? [[1,3]./4, [1,-1]./4] : [[1.0]]
    ns = sort(unique(first(k) for k in keys(table)))
    [Dict(n => sum(w[g+1] * Float64(table[(n,g)]) for g in 0:length(w)-1)
        for n in ns if all(haskey(table,(n,g)) for g in 0:length(w)-1)) for w in weights]
end

function convsquare(q, shift)
    out = zeros(eltype(q), 2length(q)-1+shift)
    for i in eachindex(q), j in eachindex(q)
        out[i+j-1+shift] += q[i]*q[j]
    end
    out
end
poly(P,x) = evalpoly(x,P)

function exterior(P,B,v; slack=1e-10)
    for sign in (-1,1)
        a = [sum(binomial(i,j)*P[i+1]*sign^i*v^(i-j) for i in j:length(P)-1)
            for j in 0:length(P)-1]
        a[1] > B || return false
        minimum(a[2:end]) >= -slack*maximum(abs,a) || return false
    end
    true
end

"Search squared filters; every selected filter must later pass exact exterior checks."
function best_filter(mom, eps, nmax; allowodd=true, capped=true, target=1.0)
    available = sort([n for n in keys(mom) if n <= nmax])
    nmax in available || return nothing
    best = (cap=(mom[nmax]+eps)^(1/nmax), shift=nmax, q=[1.0],
        P=vcat(zeros(nmax),[1.0]), B=mom[nmax]+eps)
    capped || return best
    for shift in 4:2:nmax-2, stride in (allowodd ? (1,2) : (2,))
        powers = collect(0:stride:div(nmax-shift,2))
        length(powers) >= 2 || continue
        all(shift+i+j in available for i in powers for j in powers) || continue
        M = [mom[shift+i+j] for i in powers, j in powers]
        for ridge in (0.0,eps,3eps,10eps,100eps), at0 in (1.005,1.0,0.95,0.9,0.875,0.85,0.8)
            at=target*at0
            v = at.^powers
            A = Symmetric(M + ridge*I)
            minimum(eigvals(A)) > 1e-13 || continue
            qq = A\v
            qq ./= maximum(abs,qq)
            q = zeros(maximum(powers)+1)
            q[powers.+1] = qq
            P = convsquare(q,shift)
            B = sum(c*mom[i-1] for (i,c) in enumerate(P) if c != 0) + eps*sum(abs,P)
            B > 0 || continue
            hi = best.cap
            exterior(P,B,hi) || continue
            lo = 0.05
            for _ in 1:42
                mid=(lo+hi)/2
                exterior(P,B,mid) ? (hi=mid) : (lo=mid)
            end
            best = (cap=hi, shift, q, P, B)
        end
    end
    best
end

function capped(m,h,r)
    r >= 1 || return Inf
    count=floor(Int,m/h)
    count*h^r + max(0,m-count*h)^r
end

function prepare(table, nmax, eps, mode; allowodd=true, masscap=true, shift=SHIFT, beta=12.25)
    # Input errors refer to the new shifted table; rescale centers accordingly.
    tab = Dict((n,g)=>Float64(v)*exp(-beta*n*(shift-SHIFT)) for ((n,g),v) in table)
    moments=sector_table(tab,mode)
    errors=eps .* (mode=="ION" ? [1,2/3,1/2] : mode=="JN" ? [1,1/2] : [1])
    fs=[best_filter(m,e,nmax;allowodd,capped=masscap,target=exp(-beta*(shift-SHIFT)))
        for (m,e) in zip(moments,errors)]
    any(isnothing,fs) && return nothing
    masses=[m[nmax]+e for (m,e) in zip(moments,errors)]
    multiplicities=mode=="ION" ? [1,2,3] : mode=="JN" ? [1,3] : [1]
    (nmax,eps,mode,filters=fs,masses,multiplicities,masscap,shift,beta)
end

function bounds(input,k)
    r=k/input.nmax
    b=[input.masscap ? capped(m,f.cap^input.nmax,r) : m^r
        for (m,f) in zip(input.masses,input.filters)]
    M=sum(input.multiplicities.*b)
    E=input.mode=="ION" ? b[1]+2b[2] : b[1]
    M,E
end

function residual(D,tauj,taui,t,s,mode)
    0<t<=tauj<=D || return Inf
    if mode=="Z"
        return (D-t)^s
    elseif mode=="JN"
        return max(3*((D-t)/3)^s,(tauj-t)^s+3*((D-tauj)/3)^s)
    end
    t<=taui<=tauj || return Inf
    maximum(x^s+2*(y/2)^s+3*((D-t-x-y)/3)^s
        for (x,y) in ((0.0,0.0),(taui-t,0.0),(taui-t,tauj-taui),(0.0,tauj-t)))
end

"Finite-prefix/tail estimate, NON-CERTIFYING; stop before Float64 underflows."
function gap_estimate(p,q,beta)
    g=Inf
    for _ in 1:20
        z=p*(2-p)+q*(1-p)^2
        0<z<0.5 || return 0.0
        g=min(g,(log1p(-z)-log(z))/beta)
        if max(p,q)<1e-40
            return min(g,-log(5max(p,q))/beta)
        end
        p=z^2/(2*(1-z)^2)
        zq=q*(2-q)+p*(1-q)^2
        0<zq<1 || return 0.0
        q=zq^2/(2*(1-zq)^2)
        beta*=2
    end
    0.0
end

function seed(input, trial, k, r; trial_floor=1.0)
    trial%4==0 && iseven(k) && input.nmax<=k<=trial/2 || return nothing
    M,E=bounds(input,k)
    D,E=M^r,E^r
    h=trial/(2k)
    tauj=input.mode=="Z" ? D : (D+3E)/4
    taui=input.mode=="ION" ? (D+11E)/12 : tauj
    spatial_floor=trial_floor^r
    t=(spatial_floor/D)^(1/(2h-1))*(1-1e-12)
    # J needs only the threefold exclusion; uniqueness is not required to remove one weight.
    multiplicity = input.mode == "ION" ? 2 : input.mode == "JN" ? 3 : Inf
    0<t<taui && multiplicity*t>D || return nothing
    R(s,t)=residual(D,tauj,taui,t,s,input.mode)
    for _ in 1:50
        v=R(2h,t)
        v<spatial_floor || return nothing
        next=(spatial_floor-v)^(1/(2h))*(1-1e-12)
        next<=t && break
        t=next
    end
    p=sqdef(R(h,t)/t^h)
    Q,_=bounds(input,trial)
    Q>=trial_floor || return nothing
    q=sqdef((Q/trial_floor-1)^r)
    gap=gap_estimate(p,q,input.beta*r)
    (p,q,gap,k,r,trial,n0=trial÷2,beta=input.beta*r,D,E,tauj,taui,t,Q,
        diagonal=max(p,q)<0.12731)
end

function scan_seeds(input; trials=(72,120), multipliers=1.0:0.125:5.0)
    rows=NamedTuple[]
    for trial in trials, k in input.nmax:2:trial÷2, r in multipliers
        row=seed(input,trial,k,r)
        isnothing(row) || row.gap<=0 || push!(rows,row)
    end
    sort!(rows;by=x->-x.gap)
end

"Provenance record: source and environment hashes plus package versions."
function metadata(config; status="non-certifying Float64 discovery")
    hashes=study_hashes("seed_design";extra=("../../julia/src/ExactArithmetic.jl",
        "../../julia/src/Results.jl","Project.toml","Manifest.toml"))
    manifest=TOML.parsefile(joinpath(PROJECT,"Manifest.toml"))
    result=run_metadata(config;source_sha256=hashes,packages=(JSON3,),status,certified=false)
    result["package_versions"]=Dict(name=>only(manifest["deps"][name])["version"]
        for name in ("JSON3","KrylovKit","Nemo","FLINT_jll"))
    result
end

"""
    write_compact(path, data)

Write a search or summary record as single-line JSON (via `JSON3.write`); these files are
regenerable floating-point outputs, so they are stored without pretty-printing.
"""
function write_compact(path, data)
    mkpath(dirname(abspath(path)))
    open(io -> (JSON3.write(io, data); println(io)), path, "w")
    path
end

# ---------------------------------------------------------------------------------------------
# Exact scalar verifier. The finite physical inputs and analytic lemmas remain hypotheses
# until separately proved in Lean; no floating-point comparison accepts an inequality here.
# ---------------------------------------------------------------------------------------------

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

"Check the stored exact MPS trial energy beats the shift: 500000 U + 700741 V < 0."
function verify_trial(trial)
    if trial in (72,120)
        path=joinpath(SPIN1_RESULTS,"trials.json")
        row=only(x for x in read_exact_json(path)["certificates"] if parse(Int,x["length"])==trial)
    else
        path=joinpath(OUT,"trial_scan.json")
        row=read_exact_json(path)["exact"]
        @assert parse(Int,row["length"])==trial
    end
    u,v=parse(BigInt,row["U"]),parse(BigInt,row["V"])
    @assert v>0 && 500000u+700741v<0
    Dict("path"=>relpath(path,PROJECT),"sha256"=>sha256_file(path),"length"=>trial,
        "energy_per_bond"=>string(u//v),"positive_margin"=>string(-500000u-700741v),
        "scope"=>trial in (72,120) ? "stored exact trial output; no fresh contraction" :
            "fresh Julia integer contraction")
end

# ---------------------------------------------------------------------------------------------
# Certificate writer and checker
# ---------------------------------------------------------------------------------------------
function audit_stored(source,eps)
    @assert eps >= (source==1 ? frac(3,100000) : frac(1,100000))
    verification=read_exact_json(joinpath(SPIN1_RESULTS,"verification.json"))
    @assert verification["status"]=="passed"
    paths=["source$source","trials"]
    hashes=Dict{String,String}()
    for name in paths
        path=joinpath(SPIN1_RESULTS,"$name.json")
        h=sha256_file(path)
        @assert verification["components"][name]["sha256"]==h
        hashes[name]=h
    end
    raw=read_exact_json(joinpath(SPIN1_RESULTS,"source$source.json"))
    for row in raw["certificates"]
        c=parse_rational(row["center"])
        lo=parse_rational(row[source==1 ? "lower" : "true_trace_lower"])
        hi=parse_rational(row[source==1 ? "upper" : "true_trace_upper"])
        @assert c-eps<lo<=hi<c+eps
    end
    hashes
end

function proposals_for(source,eps,mode,nmax;odd=false)
    table=stored_table(source)
    input=prepare(table,nmax,Float64(eps),mode;allowodd=odd,beta=source==1 ? 12.25 : 10.5)
    [Dict("q"=>f.q,"shift"=>f.shift) for f in input.filters]
end

"Recompute the three exact designs; write results/seed_design/{certificate,compact_witnesses}.json."
function certify_seed_design()
    rows=Any[]
    for (name,source,eps,mode,nmax,trial,k,r) in (
        ("small_simple",2,frac(3,10000),"ION",10,120,30,frac(5,2)),
        ("small_relaxed",2,frac(4,10000),"ION",10,144,30,frac(25,8)),
        ("real_even",1,frac(1,1000),"JN",12,120,30,frac(2)),
    )
        started=time()
        inp=exact_input(source,mode,nmax,eps,proposals_for(source,eps,mode,nmax);
            uncappedO=mode=="ION")
        trial_record=verify_trial(trial)
        seed=exact_seed(inp,trial,k,r)
        beta=parse_rational(inp["base_beta"])*r
        chain=exact_chain(parse_rational(seed["p"]),parse_rational(seed["q"]),beta)
        ratio=frac(784)/beta
        @assert parse_rational(chain["r"])^numerator(ratio)>frac(20)^denominator(ratio)
        chain["dominates_original_target"]=true
        chain["target_power_ratio"]=string(ratio)
        push!(rows,Dict("name"=>name,"input"=>inp,"seed"=>seed,"chain"=>chain,
            "trial"=>trial_record,"input_sha256"=>audit_stored(source,eps),
            "seconds"=>time()-started))
        println(name,": p=",Float64(parse_rational(seed["p"])),
            " q=",Float64(parse_rational(seed["q"])),
            " gap > ",chain["gap_display"]," every even L >= ",trial÷2," seconds=",time()-started)
        flush(stdout)
    end
    write_result(joinpath(OUT,"certificate.json"),Dict(
        "metadata"=>metadata("three exact scalar designs";
        status="exact rational scalar checks passed; " *
            "physical premises and analytic lemmas not Lean-checked"),
        "designs"=>rows))
    compact=[Dict("name"=>row["name"],"input"=>row["input"],
        "trial"=>row["seed"]["trial"],"k"=>row["seed"]["k"],
        "temperature_multiplier"=>row["seed"]["r"],"D"=>row["seed"]["D"],
        "E"=>row["seed"]["E"],"Q"=>row["seed"]["Q"],
        "dominant_lower_steps"=>vcat([row["seed"]["initial_t"]],
            [r["new"] for r in row["seed"]["refinements"]]),
        "chain_rate"=>row["chain"]["r"],"chain_steps"=>length(row["chain"]["s"])-1,
        "chain_grid"=>"10^512") for row in rows]
    write_result(joinpath(OUT,"compact_witnesses.json"),Dict(
        "scope"=>"compact replay recipes; not standalone physical or Lean certificates",
        "designs"=>compact))
end

"Replay every exact record of results/seed_design/certificate.json; throws on any mismatch."
function check_seed_design(path=joinpath(OUT,"certificate.json"))
    c=JSON3.read(read(path,String),Dict{String,Any})
    for row in c["designs"]
        inp=row["input"];seed=row["seed"];chain=row["chain"]
        eps=parse_rational(inp["eps"])
        @assert audit_stored(inp["source"],eps)==row["input_sha256"]
        moments,errors=exact_sector_data(inp["source"],inp["mode"],eps)
        @assert parse_rational.(inp["masses"])==[m[inp["nmax"]]+e for (m,e) in zip(moments,errors)]
        for (m,e,f) in zip(moments,errors,inp["filters"])
            isnothing(f) && continue
            q=parse_rational.(f["q"])
            @assert exact_filter(m,e,f["shift"],q)==f
        end
        @assert verify_trial(seed["trial"])==row["trial"]
        calculated=exact_seed(inp,seed["trial"],seed["k"],parse_rational(seed["r"]))
        @assert calculated==seed
        @assert parse_rational(chain["beta"])==
            parse_rational(inp["base_beta"])*parse_rational(seed["r"])
        @assert parse_rational(first(chain["s"]))==parse_rational(seed["p"])
        @assert parse_rational(first(chain["t"]))==parse_rational(seed["q"])
        @assert chain_valid(chain)
        ratio=frac(784)/parse_rational(chain["beta"])
        @assert parse_rational(chain["r"])^numerator(ratio)>frac(20)^denominator(ratio)
        @assert 0<parse_rational(chain["gap_strict_lower"])<
            first(log_enclosure(parse_rational(chain["r"])))/parse_rational(chain["beta"])
        println("PASS ",row["name"]," exact scalar inequalities and ChainCertificate.Valid")
    end
    true
end

export_all(@__MODULE__)

end
