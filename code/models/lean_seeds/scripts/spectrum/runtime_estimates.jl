include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.Spectrum
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, JSON3
summary=JSON3.read(read(joinpath(DAT,"summary.json"),String))
rows=[]
for t in summary.totals
    t.degree==0 || continue
    traffic=parse(Float64,String(t.bit_pass_volume))
    low=(5e-11*traffic+0.015*t.gram_4M_equivalent+t.small_loop_seconds)/60
    high=(1e-10*traffic+0.035*t.gram_4M_equivalent+t.small_loop_seconds)/60
    supplied=(1e-10*traffic+t.gram_4M_equivalent+t.small_loop_seconds)/60
    push!(rows,(;source=Int(t.source),bits=Int(t.bits),method=String(t.method),
        layout=String(t.layout),sectors=Int(t.sectors),successes=Int(t.successes),
        hv=Int(t.hv),bit_pass_volume=String(t.bit_pass_volume),
        gram_mults=Int(t.gram_mults),gram_4M_equivalent=Float64(t.gram_4M_equivalent),
        measured_model_minutes_low=low,measured_model_minutes_high=high,
        supplied_1second_multiply_model_minutes=supplied))
end
writecsv(joinpath(DAT,"runtime_estimates.csv"),rows)
