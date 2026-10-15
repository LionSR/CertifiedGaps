# Produce the per-point certificates results/certificate_g<g>_h<h>.json.gz and results/seeds.json.
#
# For each grid point: diagonalise every ring of n and 2n sites (n ∈ 2:2:NMAX) densely in
# floating point, rank the seeds (n₀, β₀) numerically, then redo the best one with rigorous
# eigenvalue enclosures and interval partition sums and apply `certified_gap_bound`.
# Usage: julia --project=models/tfim_longitudinal models/tfim_longitudinal/certify.jl \
#     [point indices 1-30]
# (g-major order; NMAX defaults to 6). Takes about five minutes for all 30 points.
include(joinpath(@__DIR__, "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, CertifiedGaps, CertifiedGaps.Seeds, CertifiedGaps.Results
using Printf, IntervalArithmetic
const C=LongitudinalCertificate
nmax=parse(Int,get(ENV,"NMAX","6")); ns=collect(2:2:nmax); betas=collect(0.25:0.25:64)
points=[(g,h) for g in C.GS for h in C.HS]
!isempty(ARGS) && (points=points[parse.(Int,ARGS)])
for (g,h) in points
    start=time(); C.rp_check(g,h); spec=Dict{Int,Any}()
    tag="g$(numerator(g))_$(denominator(g))_h$(numerator(h))_$(denominator(h))"
    for L in sort(unique(vcat(ns,2ns)))
        spec[L]=C.ring_spectrum(L,g,h)
        @printf("g=%g h=%g L=%d dimmax=%d %.2fs gap=%.9f\n",g,h,L,
            maximum(spec[L].dimensions),spec[L].seconds,spec[L].E[2]-spec[L].E[1])
        flush(stdout)
    end
    # Numerical ranking only: the chosen seed is re-certified rigorously below.
    rows=scan_seeds((n,b)->C.defects(spec,n,b),ns,betas); best=best_seed(rows)
    summary=[Dict("n"=>n,
        "max_min_purity"=>maximum(min(1-r.εₛ,1-r.εₜ) for r in rows if r.n₀==n),
        "delta"=>spec[2n].E[1]-2spec[n].E[1]) for n in ns]
    data=Dict{String,Any}("g"=>C.str(g),"h"=>C.str(h),"nmax"=>nmax,
        "diagnostics"=>summary,"spectra"=>Dict(string(L)=>Dict("dimensions"=>s.dimensions,
        "seconds"=>s.seconds,"E0"=>s.E[1],"gap"=>s.E[2]-s.E[1]) for (L,s) in spec),
        "scan_seconds"=>time()-start,"beta_grid"=>["1//4","1//4","64//1"],
        "best"=>best===nothing ? nothing : Dict("n"=>best.n₀,"beta"=>C.str(best.β₀),
            "epsilon_s"=>best.εₛ,"epsilon_t"=>best.εₜ,"bound"=>best.bound))
    if best!==nothing
        n=best.n₀; b=C.exact(best.β₀); small=C.ring_spectrum(n,g,h;rigorous=true)
        large=C.ring_spectrum(2n,g,h;rigorous=true)
        seeds=C.seed_enclosure(small.blocks,large.blocks,b)
        bound=certified_gap_bound(seeds.es,seeds.et,b).bound
        data["certificate"]=Dict("n"=>n,"beta"=>C.str(b),"epsilon_s"=>C.str(seeds.es),
            "epsilon_t"=>C.str(seeds.et),"bound"=>bound,
            "bound_exact_binary"=>bound===nothing ? nothing : C.str(bound),
            "small"=>small.blocks,"large"=>large.blocks,
            "logZ_intervals"=>[[C.str(IntervalArithmetic.inf(x)),
                C.str(IntervalArithmetic.sup(x))] for x in seeds.logs],
            "seconds"=>small.seconds+large.seconds,
            "validity"=>"every even L >= $n cells = sites")
        @printf("CERTIFICATE g=%g h=%g n=%d beta=%g bound=%s\n",g,h,n,b,string(bound))
        flush(stdout)
    end
    data["seconds"]=time()-start
    certified=haskey(data,"certificate") && data["certificate"]["bound"]!==nothing
    data["metadata"]=C.metadata(Dict("g"=>C.str(g),"h"=>C.str(h),"nmax"=>nmax);certified)
    out=joinpath(C.ROOT,"results","certificate_"*tag*".json"); write_result(out,data)
    # Spectra make each file ~260 kB; store compressed (`gzip -dc` restores the JSON).
    run(`gzip -9 -n -f $out`)
end
# The seed table is rebuilt from every stored certificate, not only the points just run.
length(C.result_paths())==length(C.GS)*length(C.HS) && C.write_seed_table()
