# Chebyshev block-filter runs (NON-CERTIFYING). Writes one results/spectrum/filter_s*.json per run;
# existing files are treated as a cache, so move them aside to recompute.
#
#   julia --project=. scripts/spectrum/filters.jl float [n...]  # Float64, degrees 16/32/64
#   julia --project=. scripts/spectrum/filters.jl fixed [n...]  # Q48/Q112, degree 32
#   julia --project=. scripts/spectrum/filters.jl sweep         # Q48/Q112, deg. 16/64, n = 10, 12
include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.Spectrum
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, JSON3

if abspath(PROGRAM_FILE)==@__FILE__
    mode=popfirst!(ARGS)
    ns=isempty(ARGS) ? collect(4:12) : parse.(Int,ARGS)
    if mode=="float"
        run_filters(ns)
    elseif mode=="fixed"
        run_filters(ns;degrees=(32,),bitslist=(64,128))
    elseif mode=="sweep"
        run_filters([10,12];degrees=(16,64),bitslist=(64,128))
    else
        error("mode must be float, fixed or sweep")
    end
end
