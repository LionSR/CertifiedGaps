include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.Spectrum
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, JSON3
# The elementary comparison cuts two bonds and bounds each by -2 I.
# This is a diagnostic for the naive comparison, not a test of tuned boundaries.
analysis=JSON3.read(read(joinpath(DAT,"analysis.json"),String))
cache=Dict((n,m)=>eigvals(Hermitian(Matrix(first(hamiltonian(n,m,"I";periodic=false)))))
           for n in 2:6 for m in -n:n)
rows=[]
for t in analysis.precision
    t.budget==3e-5 || continue
    n=Int(t.n); m=Int(t.m); twist=String(t.twist)
    relevant=(t.beta==12.25 && twist!="C") || (t.beta==10.5 && n in (6,8,10))
    relevant || continue
    a=n÷2; b=n-a
    vals=sort([x+y-4 for ma in -a:a if abs(m-ma)<=b
               for x in cache[(a,ma)] for y in cache[(b,m-ma)]])
    @assert length(vals)==t.dimension
    tr=only(filter(r->r.n==n && r.twist==twist && r.beta==t.beta,analysis.traces))
    rho=tr.E0+t.cutoff
    push!(rows,(;source=t.beta==12.25 ? 1 : 2,n,twist,m,dimension=Int(t.dimension),
                target_k=Int(t.k),rho,naive_anchor_count=count(<(rho),vals),
                naive_lower_lambda_kplus1=vals[Int(t.k)+1],
                deficit=rho-vals[Int(t.k)+1]))
end
writecsv(joinpath(DAT,"anchor_probe.csv"),rows)
