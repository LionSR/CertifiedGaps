# Tiny independently kernel-checked operation timings; no native_decide.
include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.Spectrum
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, JSON3
lean=joinpath(ROOT,"runtime/lean-4.35.0-rc3-linux/bin/lean")
rows=[]
for bits in (1000000,4000000,10000000,20000000,68000000),
        kind in ("baseline","multiply","maskloop")
    predicate=kind=="multiply" ? "((a*b) >>> $bits) > 0" :
        kind=="baseline" ? "a > 0" : "work 100 a ((1 <<< $bits)-1) < (1 <<< $bits)"
    code="""
    import Lean
    set_option maxRecDepth 1000000
    set_option maxHeartbeats 0
    def a : Nat := ((1 <<< $bits)-1)/3
    def b : Nat := ((1 <<< $bits)-1)/5
    def work : Nat → Nat → Nat → Nat
      | 0, a, _ => a
      | n+1, a, mask => work n
          (((((a <<< 17) ||| (a >>> ($bits-17))) ^^^ 65537) &&& mask)) mask
    theorem bench : $predicate := by decide +kernel
    #print axioms bench
    """
    file=joinpath(ROOT,"runtime/bench.lean"); write(file,code)
    for repetition in 1:2
        started=time()
        output=read(`$lean $file`,String)
        elapsed=time()-started
        push!(rows,(;kind,bits,repetition,seconds=elapsed))
        println("$kind $bits $repetition $elapsed ",strip(output))
    end
end
for count in (10000,30000,100000)
    code="""
    import Lean
    set_option maxRecDepth 1000000
    set_option maxHeartbeats 0
    def loop : Nat → Nat → Nat
      | 0, a => a
      | n+1, a => loop n (a+1)
    theorem bench : loop $count 0 = $count := by decide +kernel
    #print axioms bench
    """
    file=joinpath(ROOT,"runtime/bench.lean"); write(file,code)
    started=time(); output=read(`$lean $file`,String); elapsed=time()-started
    push!(rows,(kind="small_loop",bits=count,repetition=1,seconds=elapsed))
    println("small_loop $count $elapsed ",strip(output))
end
writecsv(joinpath(DAT,"kernel_timings.csv"),rows)
