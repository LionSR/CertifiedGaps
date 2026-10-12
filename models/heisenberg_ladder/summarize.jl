include("ladder.jl")
using JSON3
readj(name)=JSON3.read(read(joinpath(@__DIR__,"results",name),String))
function summarize()
    c=readj("certificate.json"); threshold=readj("threshold.json")
    scans=[readj("scan_J$(tag).json") for tag in ("4p0","3p5","3p0","2p5","2p0","1p5")]
    extrap=Any[]
    for s in scans
        s.best!==nothing && continue
        ds=Dict(Int(d.n0)=>Float64(d.delta) for d in s.diagnostics)
        xi=2/log(ds[2]/ds[4]); gap=Float64(s.rings[Symbol("10")].gap)
        prediction=nothing
        for n in 6:2:128
            delta=ds[4]*exp(-(n-4)/xi)
            defects=(nn,beta)->(-expm1(-beta*delta),-expm1(-6nn*exp(-beta*gap)))
            b=best_seed(scan_seeds(defects,(n,),0.25:0.25:32))
            if b!==nothing
                prediction=Dict("n0"=>n,"beta"=>b.β₀,"eps_s"=>b.εₛ,"eps_t"=>b.εₜ)
                break
            end
        end
        push!(extrap,Dict("J"=>s.J,"xi_from_even_deltas"=>xi,"ed_gap_10"=>gap,
            "prediction"=>prediction))
    end
    write_result(joinpath(@__DIR__,"results","extrapolation.json"),Dict(
        "metadata"=>provenance(Dict("model"=>"delta4 exp(-(n-4)/xi), P_t=exp(-6n exp(-beta gap10))")),
        "estimates"=>extrap))
    report=joinpath(@__DIR__,"REPORT.md"); body=read(report,String)
    marker="<!-- GENERATED RESULTS -->"
    body=first(split(body,marker))
    open(report,"w") do io
        print(io,body,marker,"\n\n")
        println(io,"**Certified at J_perp=4:** Δ_L ≥ **",c.bound,"** for **every even L ≥ ",
            c.n0," rungs**, i.e. N=2L ≥ ",2c.n0," sites with N divisible by 4.\n")
        println(io,"The best seed is (n0,beta)=(",c.n0,",",c.beta,
            "), with exact upper bounds eps_s=",c.eps_s," and eps_t=",c.eps_t,". ",
            "The rigorous run took ",round(c.seconds;digits=2)," s, including output. ",
            "Maximum eigenvalue radius: ",c.maximum_eigenvalue_radius,".\n")
        println(io,"| J_perp | Best eligible (n0,beta) | Numerical eps_s | Numerical eps_t | Numerical recurrence bound | ED gap (10 rungs) | Dense 10-rung seconds | Driver seconds |")
        println(io,"|---:|:---|---:|---:|---:|---:|---:|---:|")
        for s in scans
            b=s.best; gap=s.rings[Symbol("10")].gap
            if b===nothing
                @printf(io,"| %s | none | — | — | — | %.8f | %.2f | %.2f |\n",s.J,gap,s.rings[Symbol("10")].seconds,s.seconds)
            else
                @printf(io,"| %s | (%d, %.2f) | %.10f | %.10f | %.10f | %.8f | %.2f | %.2f |\n",
                    s.J,b.n0,b.beta,b.eps_s,b.eps_t,b.bound,gap,s.rings[Symbol("10")].seconds,s.seconds)
            end
        end
        gap10=scans[1].rings[Symbol("10")].gap
        println(io,"\nThe certified bound is ",round(100c.bound/gap10;digits=3),
            "% of the numerical 10-rung ED gap ",gap10,
            ". The 8-rung ED reference is ",c.ed_gap_2n,
            " (ratio ",c.bound_over_ed_gap,
            "). These finite-ring gaps are numerical references, not thermodynamic-limit estimates.\n")
        println(io,"All scans use beta=0.25:0.25:16 and n0=2:5, restricting valid seeds to even n0. ",
            "Every numerical scan still computes the full dense spectrum through 10 rungs.\n")
        println(io,"| Rungs L | Sites | Total dimension | Largest dense block | Numerical seconds at J=4 |")
        println(io,"|---:|---:|---:|---:|---:|")
        for L in (2,3,4,5,6,8,10)
            r=scans[1].rings[Symbol(string(L))]
            @printf(io,"| %d | %d | %d | %d | %.3f |\n",L,2L,r.dimension,r.maximum_block_dimension,r.seconds)
        end
        println(io,"\nExact sector dimensions and multiplicities are individually recorded in each spectrum JSON.\n")
        lo=Float64(parse(Rational{Int},threshold.last_failure)); hi=Float64(parse(Rational{Int},threshold.first_success))
        println(io,"### Numerical threshold\n")
        println(io,"Fixed-grid bisection finds a failure at J=",lo," and a success at J=",hi,
            ". Exact sampled endpoints are ",threshold.last_failure," and ",threshold.first_success,
            ". This locates an observed transition near ",round((lo+hi)/2;digits=4),
            " for n0∈{2,4} and this beta grid. It does not certify the endpoint, global ",
            "monotonicity, a continuum-beta optimum, or a physical gap threshold.\n")
        println(io,"### Failed seeds and size extrapolation\n")
        println(io,"| J_perp | n0 | max_beta min(P_s,P_t) | beta at maximum | delta(n0) | Eligible even seed? |")
        println(io,"|---:|---:|---:|---:|---:|:---|")
        for s in scans
            s.best!==nothing && continue
            for d in s.diagnostics
                @printf(io,"| %s | %d | %.9f | %.2f | %.10g | %s |\n",s.J,d.n0,
                    d.max_min_purity,d.beta_at_max_min,d.delta,d.eligible_even_seed ? "yes" : "no")
            end
        end
        println(io,"\nA deliberately tentative extrapolation fits delta(n)=delta(4) exp(-(n-4)/xi) ",
            "using only the even n=2,4 data, xi=2/log(delta(2)/delta(4)). It models ",
            "P_s≈exp(-beta delta(n)) and P_t≈exp(-6n exp(-beta gap_10)). The latter ",
            "is a dilute-triplet estimate. The first predicted even seed, maximizing the ",
            "recurrence bound on beta=0.25:0.25:32 at that n, is:\n")
        for e in extrap
            println(io,"* J=",e["J"],": xi≈",round(e["xi_from_even_deltas"];digits=3),
                ", predicted (n0,beta)=",e["prediction"]===nothing ? "not found through 128" :
                string((e["prediction"]["n0"],e["prediction"]["beta"])),".")
        end
        println(io,"\nThese are numerical hypotheses, not seed enclosures. Two even sizes cannot ",
            "validate exponential scaling, and the doubled leg bond at L=2 biases that fit. ",
            "Odd-n delta values come from frustrated rings and are excluded from the fit. ",
            "Actual even n0=6 testing would require a 12-rung ring beyond the requested cap.\n")
        println(io,"### Certified versus numerical\n\n",
            "* Certified: rational RP identity and real onsite term at every stated J; ",
            "the J=4 full-spectrum enclosures, rational seed defects, and uniform gap above.\n",
            "* Exact combinatorial accounting: all orbit stabilizers, sector dimensions, and multiplicities.\n",
            "* Numerical, not certified: parameter scans, other positive recurrence values, ",
            "threshold refinement, finite-ring reference gaps and ratios, failed-grid evidence, ",
            "and extrapolated future seeds. No claim of absence at all beta or all sizes.\n")
    end
end
summarize()
