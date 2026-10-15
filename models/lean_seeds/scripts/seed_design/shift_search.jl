include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.SeedDesign
using CertifiedGaps.ExactArithmetic, CertifiedGaps.Results
using LinearAlgebra, JSON3, SHA, Dates, TOML
function main_shift()
    records=Any[]
    # The original trial gives a known ground Boltzmann lower bound after rescaling.
    for a in (SHIFT-.02,SHIFT-.01,SHIFT-.005,SHIFT,SHIFT+1e-7),
        eps in (3e-4,1e-3,3e-3)
        input=prepare(stored_table(2),10,eps,"ION";allowodd=false,shift=a,beta=10.5)
        floor=exp(10.5*120*(SHIFT-a))
        row=seed(input,120,30,2.5;trial_floor=floor)
        push!(records,Dict("a"=>a,"eps"=>eps,"base_beta"=>10.5,"trial"=>120,
            "trial_floor"=>floor,"effective_eps_at_n10"=>eps*exp(10.5*10*(a-SHIFT)),
            "normalized_cap_I"=>input.filters[1].cap*exp(10.5*(a-SHIFT)),"seed"=>row))
    end
    write_compact(joinpath(OUT,"shift_search.json"),Dict(
        "metadata"=>metadata("shift scan; non-certifying"),"records"=>records))
    for r in records
        println((r["a"],r["eps"],r["effective_eps_at_n10"],
            isnothing(r["seed"]) ? nothing : r["seed"].gap))
    end
end
main_shift()
