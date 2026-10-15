"""
    IsingSeed

Study `ising`: kernel-checked small-seed arithmetic for the longitudinal Ising chain at
g = 1/4, h = 1/20. UNTRUSTED generator: Float64 eigensolves only propose vectors; every integer
Gram/residual claim is independently re-evaluated by Lean's kernel in
`lean/SeedKernel/Ising/Generated*.lean`. The exact Julia replay (`verify_certificate`) is a
diagnostic, not part of the trust base.
"""
module IsingSeed

using LinearAlgebra, SparseArrays, JSON, SHA
using ..Common: MODEL, LEAN, results_dir, export_all
BLAS.set_num_threads(1)
const SCALE = big(2)^52
const WIDTH = 56
const BIAS = big(2)^54
const OUT = results_dir("ising")
const LEANDIR = joinpath(LEAN, "SeedKernel", "Ising")
"The generator, whose SHA-256 the certificate records (path relative to MODEL)."
const GENERATOR = "src/ising.jl"
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


"Lean source of `Generated.lean`: the proposals as constants and their `decide +kernel` checks."
function generated_lean(cert)
    io=IOBuffer()
    println(io,"import SeedKernel.Ising.Packed\n\nnamespace SeedKernel.Ising\n")
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
    println(io,"end SeedKernel.Ising")
    String(take!(io))
end

"""
Lean source of `GeneratedChunks.lean`: separately checked direct scalar blocks (64 blocks of 64
states of the 12-site ring), which bound kernel memory.
"""
function chunks_lean(c)
    p=c["points"][2]; n=p["n"]; @assert n==12
    V=map(p["packed"]) do x
        u=exactint(x); [((u>>(WIDTH*s))&((big(1)<<WIDTH)-1))-BIAS for s in 0:4095]
    end
    lam=exactint.(p["lambda"])
    R=[SCALE.*apply_h(n,V[i]).-lam[i].*V[i] for i in 1:2]
    io=IOBuffer()
    println(io,"import SeedKernel.Ising.Generated\nimport SeedKernel.Ising.Chunked\n")
    println(io,"namespace SeedKernel.Ising\nset_option maxHeartbeats 0\nset_option maxRecDepth 1000000")
    for j in 0:63
        ix=64j+1:64j+64
        vals=[sum(V[1][ix].^2),sum(V[2][ix].^2),sum(V[1][ix].*V[2][ix]),
              sum(R[1][ix].^2),sum(R[2][ix].^2)]
        println(io,"theorem chunk$j : chunkStats 12 compact12_1 compact12_2 lam12_1 lam12_2 $(64j) 64 =\n  (",join(vals,", "),") := by decide +kernel")
    end
    println(io,"def allChunks12 : List Stats := [")
    println(io,join(["  chunkStats 12 compact12_1 compact12_2 lam12_1 lam12_2 $(64j) 64" for j in 0:63],",\n"),"]")
    println(io,"theorem direct12_chunked : goodStats (sumStats allChunks12) = true := by")
    println(io,"  unfold allChunks12\n  rw [",join(["chunk$j" for j in 0:63],", "),"]\n  decide +kernel")
    # Check the exact coverage separately, rather than leave a generator-controlled list trusted.
    println(io,"theorem coverage12 : ((List.range 64).flatMap (fun b => (List.range 64).map (fun i => 64*b+i))) = List.range 4096 := by decide +kernel")
    println(io,"#print axioms direct12_chunked\n#print axioms coverage12\nend SeedKernel.Ising")
    String(take!(io))
end

"""
    certify_ising()

Propose vectors for n = 6, 12 (Float64 LAPACK, untrusted), round them to exact integers, and
write `results/ising/certificate.json`, `lean/SeedKernel/Ising/Generated.lean` and
`GeneratedChunks.lean`. About 10 s. No RPRing seed is claimed.
"""
function certify_ising()
    cert=Dict("scale"=>exactstr(SCALE),"width"=>WIDTH,"bias"=>exactstr(BIAS),
        "points"=>proposal.([6,12]),
        "source_hashes"=>Dict(GENERATOR=>bytes2hex(sha256(read(joinpath(MODEL,GENERATOR))))),
        "run_metadata"=>Dict("julia"=>string(VERSION),"blas_threads"=>1,
            "JSON"=>string(Base.pkgversion(JSON)), "blas"=>string(BLAS.get_config()),
            "generator"=>"untrusted Float64 LAPACK plus exact BigInt arithmetic",
            "status"=>"partial: not a seed proof"))
    open(joinpath(OUT,"certificate.json"),"w") do io; JSON.print(io,cert); end
    # Round-trip through the file so that the Lean files are generated from what is stored.
    stored=JSON.parsefile(joinpath(OUT,"certificate.json"))
    @assert verify_certificate(stored)
    write(joinpath(LEANDIR,"Generated.lean"),generated_lean(stored))
    write(joinpath(LEANDIR,"GeneratedChunks.lean"),chunks_lean(stored))
    println("Generated exact integer proposals for n=6,12 and 64 direct chunks. ",
        "No RPRing seed is claimed.")
end

"Exact Julia replay of every integer claim of the certificate (a diagnostic)."
function verify_certificate(c)
    @assert exactint(c["scale"])==SCALE
    @assert c["width"]==WIDTH && exactint(c["bias"])==BIAS
    for p in c["points"]
        n=p["n"]; d=1<<n
        V=map(p["packed"]) do x
            packed=exactint(x)
            [((packed>>(WIDTH*s))&((big(1)<<WIDTH)-1))-BIAS for s in 0:d-1]
        end
        @assert exactstr.(pack.(V))==p["packed"] # rejects unconsumed high digits
        lam=exactint.(p["lambda"])
        G=[sum(V[i].*V[j]) for i in 1:2,j in 1:2]
        @assert exactstr.(vec(G))==p["gram"]
        R=[SCALE.*apply_h(n,V[i]).-lam[i].*V[i] for i in 1:2]
        rs=[sum(r.^2) for r in R]
        @assert exactstr.(rs)==p["residual_sq"]
        @assert maximum(abs.(G-SCALE^2*I))*big(10)^12 <= SCALE^2
        @assert sum(rs)*big(10)^18<=SCALE^4
        anchor=exactint(p["anchor_20H"])
        @assert anchor == -26n+82
        @assert count(s->diagonal(n,s)<anchor+5n,0:d-1)==2
        @assert lam[1]+SCALE<lam[2] && lam[2]+SCALE<anchor*SCALE
    end
    true
end

"""
    check_ising()

Exact replay of `results/ising/certificate.json`, the recorded generator hash, and byte
equality of `Generated.lean`/`GeneratedChunks.lean` with the text generated from the stored
certificate. The Lean kernel itself is run by `Common.build_seed_kernel`. Throws on mismatch.
"""
function check_ising()
    c=JSON.parsefile(joinpath(OUT,"certificate.json"))
    @assert verify_certificate(c)
    @assert c["source_hashes"]==Dict(GENERATOR=>bytes2hex(sha256(read(joinpath(MODEL,GENERATOR)))))
    generated_lean(c)==read(joinpath(LEANDIR,"Generated.lean"),String) ||
        error("lean/SeedKernel/Ising/Generated.lean differs from results/ising/certificate.json")
    chunks_lean(c)==read(joinpath(LEANDIR,"GeneratedChunks.lean"),String) ||
        error("lean/SeedKernel/Ising/GeneratedChunks.lean differs from the certificate")
    println("ising: exact Julia replay, generator hash and generated Lean files match ",
        "(this does not certify RPRing spectra)")
    true
end

export_all(@__MODULE__)

end
