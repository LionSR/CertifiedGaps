# UNTRUSTED generator. Floating-point eigensolves only propose vectors. Every integer
# Gram/residual claim is independently re-evaluated by Lean's kernel in Generated.lean.
using LinearAlgebra, SparseArrays, JSON, SHA
BLAS.set_num_threads(1)
const SCALE = big(2)^52
const WIDTH = 56
const BIAS = big(2)^54
const OUT = @__DIR__
exactstr(x::Integer)=string(x,"//1")
function exactint(s::AbstractString)
    a=split(s,"//"); @assert length(a)==2 && a[2]=="1"
    parse(BigInt,a[1])
end

"Integer diagonal of 20 H for g=1/4, h=1/20, with periodic boundary conditions."
function diagonal(n,s)
    z(j)=1-2*((s>>j)&1)
    -20sum(z(j)*z(mod(j+1,n)) for j in 0:n-1)-sum(z(j) for j in 0:n-1)
end
"Apply the exact integer Hamiltonian, without BLAS or a floating-point model."
function apply_h(n,v)
    [diagonal(n,s)*v[s+1]-5sum(v[xor(s,1<<j)+1] for j in 0:n-1)
     for s in 0:(1<<n)-1]
end
"Pack signed fixed-point integers as biased base-2^56 digits."
function pack(v; width=WIDTH)
    @assert all(0 .<= v .+ BIAS .< big(2)^width)
    foldr((x,a)->(a<<width)+(x+BIAS),v;init=big(0))
end
function proposal(n)
    d=1<<n
    H=zeros(Float64,d,d)
    for s in 0:d-1
        H[s+1,s+1]=diagonal(n,s)
        for j in 0:n-1; H[xor(s,1<<j)+1,s+1]=-5; end
    end
    t=@elapsed F=eigen(Symmetric(H),1:2)
    V=[round.(BigInt,F.vectors[:,j].*Float64(SCALE)) for j in 1:2]
    lam=round.(BigInt,F.values.*Float64(SCALE))
    AV=apply_h.(n,V)
    R=[SCALE.*AV[j].-lam[j].*V[j] for j in 1:2]
    gram=[dot(V[i],V[j]) for i in 1:2,j in 1:2]
    res=[dot(r,r) for r in R]
    # All tests below are exact integers, even though their candidates came from LAPACK.
    @assert maximum(abs.(gram-SCALE^2*I)) * big(10)^12 <= SCALE^2
    @assert sum(res)*big(10)^18 <= SCALE^4 # residual Frobenius <= 1e-9 for 20 H
    anchor=-21n+82-5n # two domain walls cost 80, field flip costs 2
    @assert count(s->diagonal(n,s)<anchor+5n,0:d-1)==2
    @assert lam[1]+SCALE < lam[2] && lam[2]+SCALE<anchor*SCALE
    Dict("n"=>n,"packed"=>exactstr.(pack.(V)),"lambda"=>exactstr.(lam),
         "gram"=>exactstr.(vec(gram)),"residual_sq"=>exactstr.(res),
         "anchor_20H"=>exactstr(anchor),"generator_seconds"=>t,
         "numerical_E"=>F.values./20,"status"=>"integer claims; spectral interpretation pending")
end

function main()
    cert=Dict("scale"=>exactstr(SCALE),"width"=>WIDTH,"bias"=>exactstr(BIAS),
        "points"=>proposal.([6,12]),
        "source_hashes"=>Dict("certify.jl"=>bytes2hex(sha256(read(@__FILE__)))),
        "run_metadata"=>Dict("julia"=>string(VERSION),"blas_threads"=>1,
            "JSON"=>string(Base.pkgversion(JSON)), "blas"=>string(BLAS.get_config()),
            "generator"=>"untrusted Float64 LAPACK plus exact BigInt arithmetic",
            "status"=>"partial: not a seed proof"))
    open(joinpath(OUT,"certificate.json"),"w") do io; JSON.print(io,cert,2); end
    open(joinpath(OUT,"LeanSeedSmall","Generated.lean"),"w") do io
        println(io,"import LeanSeedSmall.Packed\n\nnamespace LeanSeedSmall\n")
        println(io,"set_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n")
        for p in cert["points"]
            n=p["n"]
            for i in 1:2
                println(io,"def compact$(n)_$(i) : Nat := $(exactint(p["packed"][i]))")
                println(io,"def v$(n)_$(i) : Nat := Packed.widen $n compact$(n)_$(i)")
                println(io,"def lam$(n)_$(i) : Int := $(exactint(p["lambda"][i]))")
            end
            println(io,"theorem check$n : Packed.check $n v$(n)_1 v$(n)_2 lam$(n)_1 lam$(n)_2 = true := by\n  decide +kernel")
            println(io,"theorem anchor$n : classicalCount $n ($(exactint(p["anchor_20H"]))+5*$n) = 2 := by\n  decide +kernel")
            println(io,"#print axioms check$n\n#print axioms anchor$n\n")
        end
        println(io,"theorem direct6 : checkVectors 6 compact6_1 compact6_2 lam6_1 lam6_2 = true := by\n  decide +kernel")
        println(io,"#print axioms direct6")
        println(io,"end LeanSeedSmall")
    end
    println("Generated exact integer proposals for n=6,12. No RPRing seed is claimed.")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
