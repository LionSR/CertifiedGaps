# Rebuilds results/scaling.json (grid maxima, best even seeds, energy scaling, extrapolations)
# from the complete spectra results/spectrum_L*.json[.gz] and the references. The spectra are
# not committed; regenerate them with scripts/scan.jl and the SU(2) pipeline first.
# Numerical only.
using CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
using JSON3, LinearAlgebra, Printf
const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"src","result_io.jl"))
BLAS.set_num_threads(1)
const THETAS=[0.0,-0.1,-0.3,-0.5,-0.6,-0.7]
ratstr(x)=string(ceil(BigInt,BigFloat(x)*big(10)^15)//big(10)^15)
"""
    forecast(ns, ds, gap, style)

First even n0 at which a dilute-triplet ansatz with assumed bulk `gap` and fitted mismatch
tail (`style`: "exponential" or "exponential_sqrt") gives a positive recurrence bound.
Linear least squares, explicitly a heuristic and never a seed enclosure.
"""
function forecast(ns,ds,gap,style)
    (length(ns)<2 || any(ds .<= 0) || gap<=0) && return nothing
    y=log.(ds) .+ (style=="exponential_sqrt" ? log.(ns)./2 : zeros(length(ns)))
    a,b=hcat(ones(length(ns)),Float64.(ns))\y
    b>=0 && return nothing
    delta(n)=exp(a+b*n)/(style=="exponential_sqrt" ? sqrt(n) : 1)
    # Dilute independent triplets: z_L=(1+3 exp(-beta gap))^L after ground-energy shift.
    # This makes Ps=exp(-beta delta(n)) and Pt=[(1+3x^2)/(1+3x)^2]^(2n).
    for n in 8:2:2000
        for bg in 2.0:0.125:30.0
            beta=bg/gap; x=exp(-bg)
            es=-expm1(-beta*delta(n))
            et=-expm1(2n*(log1p(3x*x)-2log1p(3x)))
            0<=es<1 && 0<=et<1 || continue
            bound=certified_gap_bound(es,et,beta).bound
            if bound!==nothing && bound>0
                return (;n,beta,epsilon_s=es,epsilon_t=et,conditional_bound=bound,
                    correlation_length=-1/b,log_amplitude=a,gap_assumed=gap,style,
                    fit_ns=ns,fit_deltas=ds)
            end
        end
    end
    nothing
end
"Rebuild results/scaling.json from the (uncommitted) spectra and the references."
function analyze()
    started=time()
    inputs=Dict{Int,Any}()
    for file in readdir(joinpath(ROOT,"results"))
        startswith(file,"spectrum_L") || continue
        d=read_result(joinpath(ROOT,"results",file);exact=false);inputs[d.L]=d
    end
    references=Dict{Int,Any}()
    for file in readdir(joinpath(ROOT,"results"))
        startswith(file,"reference_L") || continue
        d=JSON3.read(read(joinpath(ROOT,"results",file),String));references[d.L]=d
    end
    ns=sort([n for n in keys(inputs) if haskey(inputs,2n)])
    allrows=Any[];best_even=Any[];scaling=Any[];extrapolations=Any[]
    for (t,theta) in enumerate(THETAS)
        specs=Dict(L=>result_spectrum(d,t) for (L,d) in inputs)
        results=scan_seeds(spectral_defects(specs),ns,0.125:0.125:64.0)
        eligible=filter(r->iseven(r.n₀),results)
        @assert best_seed(eligible)===nothing "eligible numerical seed: implement rigorous recheck"
        for n in ns
            row=argmax(r->min(1-r.εₛ,1-r.εₜ),filter(r->r.n₀==n,results))
            push!(allrows,(;theta,n,beta=row.β₀,spatial=1-row.εₛ,thermal=1-row.εₜ,
                max_min=min(1-row.εₛ,1-row.εₜ),eligible=iseven(n),bound=row.bound,
                grid_boundary=row.β₀ in (0.125,64.0)))
        end
        best=argmax(r->min(1-r.εₛ,1-r.εₜ),eligible)
        push!(best_even,(;theta,n=best.n₀,beta=best.β₀,
            epsilon_s_rational=ratstr(best.εₛ),epsilon_t_rational=ratstr(best.εₜ),
            beta_rational=ratstr(best.β₀),bound=best.bound,
            seed_status="numerical decimals rounded upward, NOT eigenvalue enclosures"))
        grounds=Dict(L=>s.E[1] for (L,s) in specs)
        gaps=Dict(L=>spectrum_gap(s) for (L,s) in specs)
        for (L,d) in references
            haskey(inputs,L) && continue
            grounds[L]=d.summary[t].ground;gaps[L]=d.summary[t].gap
        end
        dns=sort([n for n in keys(grounds) if haskey(grounds,2n)])
        for n in dns
            push!(scaling,(;theta,n,delta=grounds[2n]-2grounds[n],ground_n=grounds[n],
                ground_2n=grounds[2n],gap_2n=gaps[2n],complete_spectrum_2n=haskey(inputs,2n)))
        end
        fitns=filter(n->iseven(n)&&n>=4,dns)
        ds=[grounds[2n]-2grounds[n] for n in fitns]
        glens=sort(filter(iseven,collect(keys(gaps))))
        glens=glens[max(1,end-2):end]
        estimates=Any[]
        for power in (1,2)
            design=hcat(ones(length(glens)),Float64.(glens).^(-power))
            intercept,slope=design\[gaps[L] for L in glens]
            forecasts=[forecast(fitns,ds,intercept,style)
                for style in ("exponential","exponential_sqrt")]
            push!(estimates,(;method="gap intercept vs L^(-$power)",gap=intercept,
                lengths=glens,finite_gaps=[gaps[L] for L in glens],forecasts))
        end
        push!(extrapolations,(;theta,estimates))
    end
    sources=source_hashes(ROOT;dirs=("src","test"),
        files=("scripts/analyze.jl","scripts/scan.jl","scripts/reference.jl"))
    inputs_used=filter(f->startswith(f,"spectrum_L")||startswith(f,"reference_L"),
        readdir(joinpath(ROOT,"results")))
    data=(;metadata=run_metadata((;thetas=THETAS,beta_grid="1//8:1//8:64//1",n_grid=ns);
            source_sha256=sources,status="numerical negative grid survey, not certified"),
        input_hashes=file_hashes(joinpath(ROOT,"results"),inputs_used),
        seconds=time()-started,certified_gap=nothing,validity_range=nothing,best_even,
        grid_maxima=allrows,
        energy_scaling=scaling,extrapolations,
        rp=(;identity="-X=A; -X^2=-A/2-B",coefficient_G="cos(theta)-sin(theta)/2",
            coefficient_Q="-sin(theta)",
            theta_exact=["0//1","-1//10","-3//10","-1//2","-3//5","-7//10"],
            identity_arithmetic="Rational{BigInt} coefficient comparison; " *
                "trigonometric parameters symbolic"))
    write_result(joinpath(ROOT,"results","scaling.json"),data)
    for row in allrows;println(row);end
    for row in extrapolations;println(row);end
end
analyze()
