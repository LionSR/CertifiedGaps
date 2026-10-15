# Generate separately checked direct scalar blocks to bound kernel memory.
include("check.jl")
c=JSON.parsefile(joinpath(OUT,"certificate.json")); @assert verify_certificate(c)
p=c["points"][2]; n=p["n"]; @assert n==12
V=map(p["packed"]) do x
    u=exactint(x); [((u>>(WIDTH*s))&((big(1)<<WIDTH)-1))-BIAS for s in 0:4095]
end
lam=exactint.(p["lambda"])
R=[SCALE.*apply_h(n,V[i]).-lam[i].*V[i] for i in 1:2]
open(joinpath(OUT,"LeanSeedSmall/GeneratedChunks.lean"),"w") do io
    println(io,"import LeanSeedSmall.Generated\nimport LeanSeedSmall.Chunked\n")
    println(io,"namespace LeanSeedSmall\nset_option maxHeartbeats 0\nset_option maxRecDepth 1000000")
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
    println(io,"#print axioms direct12_chunked\n#print axioms coverage12\nend LeanSeedSmall")
end
println("Generated 64 exact scalar chunks, covering all 4096 states.")
