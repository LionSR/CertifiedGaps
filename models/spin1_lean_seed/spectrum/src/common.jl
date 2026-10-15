# All floating-point results in this project are NON-CERTIFYING.
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, Serialization
using Arpack, JSON3
BLAS.set_num_threads(1)
# ROOT is the project directory; results go to results/ unless SPIN1_SPECTRUM_RESULTS points
# elsewhere (used to regenerate into a scratch directory without touching the stored data).
const ROOT = normpath(joinpath(@__DIR__, ".."))
const DAT = get(ENV, "SPIN1_SPECTRUM_RESULTS", joinpath(ROOT, "results"))
const SHIFT = 700741 / 500000
const TWISTS = ("I", "P", "C")

"Computational basis with ternary digits s+1, sorted by integer encoding."
function basis(n, m)
    states = Int[]
    function visit(site, remaining, word)
        if site == n
            remaining == 0 && push!(states, word)
            return
        end
        for s in -1:1
            abs(remaining-s) <= n-site-1 || continue
            visit(site+1, remaining-s, word+(s+1)*3^site)
        end
    end
    visit(0,m,0)
    sort!(states)
end

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
