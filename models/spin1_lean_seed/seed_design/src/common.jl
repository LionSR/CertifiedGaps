# Search utilities. Float64 functions in this file are NON-CERTIFYING.
using LinearAlgebra, JSON3, SHA, Dates, TOML
# Paths: ROOT is the repository, PROJECT this directory. Generated data live under results/
# (raw search rows in results/raw/, gitignored); spectra and downloads in cache/ and upstream/.
const PROJECT = normpath(joinpath(@__DIR__, ".."))
const ROOT = normpath(joinpath(PROJECT, "../../.."))
const OUT = joinpath(PROJECT, "results")
const RAW = joinpath(OUT, "raw")
const CACHE = joinpath(PROJECT, "cache")
const UPSTREAM = joinpath(PROJECT, "upstream")
include(joinpath(ROOT, "julia/src/ExactArithmetic.jl"))
include(joinpath(ROOT, "julia/src/Results.jl"))
using .ExactArithmetic, .Results
const SHIFT = 700741 / 500000
const RATSHIFT = frac(700741, 500000)
const SPIN1_RESULTS = joinpath(ROOT, "models/spin1/results")
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
    hashes=source_hashes(PROJECT;files=("certify.jl","check.jl"))
    for name in ("julia/src/ExactArithmetic.jl","julia/src/Results.jl")
        hashes["../../../$name"]=sha256_file(joinpath(ROOT,name))
    end
    for name in ("Project.toml","Manifest.toml")
        hashes[name]=sha256_file(joinpath(PROJECT,name))
    end
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
