module AnchorStudy

using LinearAlgebra, SparseArrays, Arpack, Optim, Random, Serialization
using CertifiedGaps.Results
import Arpack, Optim

BLAS.set_num_threads(1)
const ROOT = @__DIR__
const BETA = 49 / 4
const SHIFT = 700741 / 500000
const TOLERANCE = 3e-5
num(x) = x isa AbstractString ? parse(Float64,x) : Float64(x)
ival(x) = x isa AbstractString ? parse(Int,x) : Int(x)

"Sz product basis, with digit 0,1,2 representing m=-1,0,1."
function basis(n, m)
    states = Int[]
    function visit(i, remaining, code)
        if i == n
            remaining == 0 && push!(states, code)
            return
        end
        for s in -1:1
            abs(remaining - s) <= n - i - 1 || continue
            visit(i + 1, remaining - s, code + (s + 1) * 3^i)
        end
    end
    visit(0, m, 0)
    sort!(states)
    return states
end

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
    hashes=source_hashes(ROOT;dirs=("test",),
        files=filter(f->endswith(f,".jl"),readdir(ROOT)))
    return run_metadata(config;source_sha256=hashes,packages=(Arpack,Optim),
        status,certified=false)
end

"Exact absolute shifted-trace cutoff using a single uniform tail energy."
tail_cutoff(n; budget=TOLERANCE) = log(3.0^n / budget) / BETA - n*SHIFT

end
