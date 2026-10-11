include("ladder.jl")
function threshold()
    start=time(); trials=Any[]
    function trial(J)
        specs=Dict{Int,PartialSpectrum}()
        for L in (2,4,8)
            s=full_spectrum(L,J)
            specs[L]=PartialSpectrum(s.energies,4^L)
        end
        best=best_seed(scan_seeds(spectral_defects(specs),(2,4),0.25:0.25:16))
        push!(trials,Dict("J"=>string(J),"succeeds"=>best!==nothing,
            "best"=>best===nothing ? nothing : Dict("n0"=>best.n₀,"beta"=>best.β₀,
                "eps_s"=>best.εₛ,"eps_t"=>best.εₜ,"bound"=>best.bound)))
        best!==nothing
    end
    lo=3//2; hi=4//1
    @assert !trial(lo) && trial(hi)
    for _ in 1:16
        mid=(lo+hi)/2
        if trial(mid); hi=mid; else; lo=mid; end
    end
    write_result(joinpath(@__DIR__,"results","threshold.json"),Dict(
        "metadata"=>provenance(Dict("ns"=>[2,4],"beta_grid"=>"1//4:1//4:16//1")),
        "last_failure"=>string(lo),"first_success"=>string(hi),
        "interpretation"=>"Numerical local transition bracket; monotonicity is not proved.",
        "trials"=>trials,"seconds"=>time()-start))
    println("NUMERICAL BRACKET ",lo," ",hi)
end
threshold()
