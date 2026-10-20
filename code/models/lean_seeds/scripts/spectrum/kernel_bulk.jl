include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.Spectrum
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, JSON3
lean=joinpath(ROOT,"runtime/lean-4.35.0-rc3-linux/bin/lean")
rows=[]
for (bits,reps) in ((4000000,100),(20000000,30),(68000000,10)),
        kind in ("baseline","multiply_bulk","mask_bulk")
    count=kind=="baseline" ? 0 : kind=="mask_bulk" ? 200 : reps
    expression=kind=="multiply_bulk" ? "((a*b) >>> $bits)" :
        "((((a <<< 17) ||| (a >>> ($bits-17))) ^^^ 65537) &&& b)"
    arg=kind=="multiply_bulk" ? "((1 <<< $bits)-1)/5" : "(1 <<< $bits)-1"
    write(joinpath(ROOT,"runtime/bulk.lean"),"""
    import Lean
    set_option maxRecDepth 1000000
    set_option maxHeartbeats 0
    def work : Nat → Nat → Nat → Nat
      | 0, a, _ => a
      | n+1, a, b => work n ($expression) b
    theorem bench : work $count (((1 <<< $bits)-1)/3) ($arg) > 0 := by decide +kernel
    #print axioms bench
    """)
    for repetition in 1:2
        started=time(); result=read(`$lean $(joinpath(ROOT,"runtime/bulk.lean"))`,String)
        elapsed=time()-started
        push!(rows,(;kind,bits,count,repetition,seconds=elapsed))
        println("$kind bits=$bits count=$count seconds=$elapsed ",strip(result)); flush(stdout)
    end
end
writecsv(joinpath(DAT,"kernel_bulk.csv"),rows)
