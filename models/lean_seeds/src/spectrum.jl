"""
    Spectrum

Study `spectrum`: spin-1 seed spectra, filtering precision and Lean kernel cost model. All
floating-point results are NON-CERTIFYING; `check_hamiltonian` and `check_data` are consistency
checks of the stored numerical data, not certificate checks.
"""
module Spectrum

using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, Serialization, Test
using Arpack, JSON3
using ..Common: MODEL, SPIN1_RESULTS, results_dir, sz_basis, export_all
BLAS.set_num_threads(1)
# ROOT is the model directory; results go to results/spectrum/ unless SPIN1_SPECTRUM_RESULTS
# points elsewhere (used to regenerate into a scratch directory without touching stored data).
const ROOT = MODEL
const DAT = get(ENV, "SPIN1_SPECTRUM_RESULTS", results_dir("spectrum"))
const SHIFT = 700741 / 500000
const TWISTS = ("I", "P", "C")

"Computational basis with ternary digits s+1, sorted by integer encoding."
basis(n, m) = sz_basis(n, m)
"Sparse Hermitian spin-1 ring; the only twisted bond is (n,1)."
function hamiltonian(n, m, twist; periodic=true)
    states = basis(n,m)
    index = Dict(s=>i for (i,s) in enumerate(states))
    T = twist == "C" ? ComplexF64 : Float64
    rows, cols, vals = Int[], Int[], T[]
    powers = [3^(j-1) for j in 1:n]
    phase = twist == "I" ? one(T) : twist == "P" ? -one(T) : cis(2pi/3)
    for (i,word) in enumerate(states)
        spins = [mod(div(word,p),3)-1 for p in powers]
        diag = 0
        for j in 1:n
            !periodic && j==n && continue
            k = mod1(j+1,n)
            diag += spins[j]*spins[k]
            for delta in (-1,1)
                abs(spins[j]-delta) <= 1 && abs(spins[k]+delta) <= 1 || continue
                target = word + delta*(powers[k]-powers[j])
                coeff = j == n ? (delta == 1 ? phase : conj(phase)) : one(T)
                push!(rows,index[target]); push!(cols,i); push!(vals,coeff)
            end
        end
        push!(rows,i); push!(cols,i); push!(vals,T(diag))
    end
    H = sparse(rows,cols,vals,length(states),length(states))
    @assert ishermitian(H)
    H, states
end

"Write a vector of NamedTuples as CSV (header from the first row)."
function writecsv(path, rows)
    isempty(rows) && return
    names = collect(keys(first(rows)))
    open(path,"w") do io
        println(io,join(names,','))
        for row in rows
            println(io,join((string(getproperty(row,k)) for k in names),','))
        end
    end
end

"Write `data` as compact JSON atomically (temporary file, then rename)."
function writejson(path, data)
    # Serialize before opening a file, then rename so interrupted runs cannot
    # leave a cache entry that looks complete to a subsequent invocation.
    # Compact JSON: shortest round-trip floats, no pretty-printing (keeps results/ small).
    encoded=JSON3.write(data)
    open(path*".tmp","w") do io
        println(io,encoded)
    end
    mv(path*".tmp",path;force=true)
end

"Run metadata with SHA-256 of the given project-relative source files."
function metadata(files)
    Dict("status"=>"NON-CERTIFYING numerical experiment", "julia"=>string(VERSION),
         "utc"=>string(now(UTC)), "blas"=>string(BLAS.get_config()),
         "source_hashes"=>Dict(f=>bytes2hex(sha256(read(joinpath(ROOT,f)))) for f in files))
end

"Low Ritz spectrum with an explicit residual check and a guard above cutoff."
function low_spectrum(H; cutoff=Inf, initial=32, dense=false, seed=1,force_sparse=false)
    d = size(H,1)
    t = time()
    if dense || (d <= 160 && !force_sparse)
        F = eigen(Hermitian(Matrix(H)))
        res = maximum(norm(H*F.vectors[:,i]-F.values[i]*F.vectors[:,i]) for i in 1:d)
        return (; values=F.values, vectors=F.vectors, residual=res,
                 nmult=0, seconds=time()-t, method="dense_Float64")
    end
    k = min(initial,d-2)
    totalmult = 0
    while true
        rng = MersenneTwister(seed)
        v0 = randn(rng,eltype(H),d)
        vals,V,nconv,niter,nmult,_ = eigs(H;nev=k,ncv=min(d,max(4k+20,100)),
                which=:SR,tol=2e-12,maxiter=5000,v0=v0)
        totalmult += nmult
        @assert nconv == k
        perm = sortperm(real(vals)); vals = real(vals[perm]); V=V[:,perm]
        res = maximum(norm(H*V[:,i]-vals[i]*V[:,i]) for i in 1:k)
        @assert res < 1e-8
        if last(vals) > cutoff+1e-7 || k == d-2 || isinf(cutoff)
            return (;values=vals,vectors=V,residual=res,nmult=totalmult,
                    seconds=time()-t,method="ARPACK_Float64")
        end
        k = min(2k,d-2)
    end
end

# ---------------------------------------------------------------------------------------------
# Fixed-point Chebyshev filtering
# ---------------------------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------------------------
# Consistency checks of the stored numerical data (NOT a certificate checker)
# ---------------------------------------------------------------------------------------------

"""
    check_hamiltonian()

Recompute the Hamiltonians and compare with the stored spectra: basis sizes, traces, forced
ARPACK against the stored dense spectra through n = 8, and m <-> -m spectral equality.
"""
function check_hamiltonian()
    rows=[]
    @testset "Hamiltonian and multiplicities" begin
        for n in 4:12
            @test sum(length(basis(n,m)) for m in -n:n)==3^n
        end
        for n in 4:8, twist in (n in (6,8) ? TWISTS : TWISTS[1:2])
            ring=JSON3.read(read(joinpath(DAT,"spectrum_$(n)_$(twist).json"),String))
            for s in ring.sectors
                H,_=hamiltonian(n,Int(s.m),twist)
                @test tr(H)==sum(Float64.(s.values)) ||
                      abs(tr(H)-sum(Float64.(s.values)))<1e-9
                size(H,1)>40 || continue
                F=low_spectrum(H;initial=min(24,size(H,1)-2),seed=271+n,force_sparse=true)
                err=maximum(abs.(F.values.-Float64.(s.values[1:length(F.values)])))
                @test err<1e-8
                push!(rows,(;n,twist,m=Int(s.m),nev=length(F.values),error=err,
                            residual=F.residual))
            end
        end
        for twist in TWISTS
            Hp,_=hamiltonian(6,1,twist); Hm,_=hamiltonian(6,-1,twist)
            @test eigvals(Hermitian(Matrix(Hp))) ≈ eigvals(Hermitian(Matrix(Hm)))
        end
        @test first(eigvals(Hermitian(Matrix(first(hamiltonian(4,0,"I")))))) ≈ -6
    end
    writecsv(joinpath(DAT,"dense_validation.csv"),rows)
end

"""
    check_data()

Check the stored spectra (dimensions, residuals, guard levels, trace and second moment through
n = 8), the shifted traces against the centres of `models/spin1/results/source{1,2}.json`, and
that every adaptive filter schedule in summary.json met its target.
"""
function check_data()
    analysis=JSON3.read(read(joinpath(DAT,"analysis.json"),String))
    summary=JSON3.read(read(joinpath(DAT,"summary.json"),String))
    rows=[]
    @testset "Complete numerical experiment" begin
        @test length(analysis.traces)==42
        for n in 4:12, twist in (n in (6,8,10) ? TWISTS : TWISTS[1:2])
            f=JSON3.read(read(joinpath(DAT,"spectrum_$(n)_$(twist).json"),String))
            @test sum(s.dimension*(s.m==0 ? 1 : 2) for s in f.sectors)==3^n
            @test maximum(s.residual for s in f.sectors)<1e-8
            for s in f.sectors
                @test length(s.values)==s.dimension || last(s.values)>f.E0+3.15
            end
            if n<=8
                @test abs(sum(sum(s.values)*(s.m==0 ? 1 : 2) for s in f.sectors))<1e-8
                @test sum(sum(abs2,s.values)*(s.m==0 ? 1 : 2) for s in f.sectors) ≈
                      4n*3^(n-1)
            end
        end
        for source in (1,2)
            path=joinpath(SPIN1_RESULTS,"source$source.json")
            reference=JSON3.read(read(path,String))
            for c in reference.certificates
                twist=TWISTS[Int(c.twist)+1]; beta=source==1 ? 12.25 : 10.5
                t=only(filter(t->t.n==c.n && t.twist==twist && t.beta==beta,analysis.traces))
                terms=split(String(c.center),"//")
                center=parse(Float64,terms[1])/parse(Float64,terms[2])
                @test abs(t.Z_shifted-center)<3e-5
                push!(rows,(;source,n=Int(c.n),twist,center,
                            measured=t.Z_shifted,error=t.Z_shifted-center))
            end
            expected=source==1 ? 44 : 24
            for bits in (64,128),method in ("linear","quadratic")
                r=filter(r->r.source==source && r.bits==bits && r.degree==0 &&
                            r.method==method,summary.runs)
                @test length(r)==expected
                @test all(x->x.success,r)
                # A small residual alone cannot exclude a missed lower state;
                # independently compare the Ritz values with the ED/Lanczos list.
                @test maximum(x.eigerror for x in r)<1e-5
            end
        end
    end
    writecsv(joinpath(DAT,"reference_trace_checks.csv"),rows)
end

export_all(@__MODULE__)

end
