# NON-CERTIFYING runtime model. The benchmark propositions themselves are kernel checked.
using JSON, Printf
const HERE=@__DIR__
lean=get(ENV,"LEAN_EXE",joinpath(HERE,".tools/lean-4.35.0-rc3-linux/bin/lean"))
mkpath(joinpath(HERE,"benchmarks"))
timings="--reuse" in ARGS ? JSON.parsefile(joinpath(HERE,"costs.json"))["benchmarks"] : []
for bits in ("--reuse" in ARGS ? Int[] : [262144,1179648,3706560,21251232])
    f=joinpath(HERE,"benchmarks","Mul$bits.lean")
    write(f,"""
    import Lean
    set_option maxHeartbeats 0
    set_option profiler true
    set_option profiler.threshold 1
    theorem product :
      ((1 <<< $bits : Nat)-1)*((1 <<< $bits)-3) =
      (1 <<< $(2bits))-(1 <<< $(bits+2))+3 := by decide +kernel
    #print axioms product
    """)
    log=joinpath(HERE,"logs","mul$bits.log")
    t=@elapsed open(log,"w") do io; run(pipeline(`$lean $f`,stdout=io,stderr=io)); end
    text=read(log,String)
    matches=collect(eachmatch(r"(?m)^type checking took ([\d.]+)(ms|s)",text))
    cpu=isempty(matches) ? missing : sum(parse(Float64,m[1])*(m[2]=="ms" ? .001 : 1) for m in matches)
    push!(timings,Dict("operand_bits"=>bits,"wall_seconds"=>t,"checking_seconds"=>cpu))
end
inv=JSON.parsefile(joinpath(HERE,"inventory.json"))["rows"]
function maxsector(q,n)
    q==9 && return 1107
    q==3 && return 73789
    q==4 && return 12870
    binomial(2n,n)
end
for r in inv
    d=r["dimension"]; k=r["levels_below_cutoff"]; name=r["model"]
    charge=!(startswith(name,"Ising"))
    m=charge ? maxsector(r["local_dimension"],r["n0"]) : d
    sectors=charge ? (r["local_dimension"] in (3,9) ?
        (r["local_dimension"]==9 ? 17 : 25) : 17) : 1
    bits=288m
    nearest=argmin(abs(log(t["operand_bits"]/bits)) for t in timings)
    t=timings[nearest]["checking_seconds"]
    # Interpolation is an assumption, not a guarantee; masks, generation and anchors excluded.
    per=t*(bits/timings[nearest]["operand_bits"])^1.15
    pairs=3k*(k+1)/2
    r["largest_unsymmetrized_or_Sz_block"]=m
    r["three_gram_products_worst"]=pairs
    r["three_gram_products_balanced"]=3*(k^2/sectors+k)/2
    r["gram_hours_balanced"]=r["three_gram_products_balanced"]*per/3600
    r["gram_hours_concentrated"]=pairs*per/3600
    r["two_vectors_per_level_binary_MiB"]=2*7*m*k/2.0^20 # forward+reverse, at largest block width
end
open(joinpath(HERE,"costs.json"),"w") do io
    JSON.print(io,Dict("status"=>"NON-CERTIFYING, optimistic arithmetic-only model",
        "benchmarks"=>timings,"rows"=>inv),2)
end
println(JSON.json(timings))
