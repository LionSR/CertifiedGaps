# Chebyshev block-filter runs (NON-CERTIFYING). Writes one results/filter_s*.json per run;
# existing files are treated as a cache, so move them aside to recompute.
#
#   julia --project=. scripts/filters.jl float [n...]  # Float64, degrees 16/32/64
#   julia --project=. scripts/filters.jl fixed [n...]  # Q48/Q112 integer recurrence, degree 32
#   julia --project=. scripts/filters.jl sweep         # Q48/Q112 at degrees 16/64 for n = 10, 12
include("../src/filtering.jl")

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
