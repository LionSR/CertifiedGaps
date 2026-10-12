include("model.jl")
root=@__DIR__
for arg in ARGS
    D=parse(Rational{Int},arg)
    spectra=Dict{Int,PartialSpectrum}()
    for L in (6,12)
        data=JSON3.read(read(joinpath(root,"results","scan_D$(numerator(D))_$(denominator(D))_L$L.json"),String),Dict{String,Any})
        spectra[L]=PartialSpectrum(Float64.(data["energies"]),Int(data["dimension"]),data["complete"])
    end
    results=scan_seeds(spectral_defects(spectra),[6],2.75:0.002:3.75)
    best=best_seed(results)
    data=Dict("D"=>string(D),"certified"=>false,
        "beta_grid"=>"2.75:0.002:3.75", "n0"=>6,
        "best"=>best===nothing ? nothing : Dict("n0"=>best.n₀,"beta"=>best.β₀,
            "beta_exact"=>string(round(Int,1000best.β₀)//1000),"eps_s"=>best.εₛ,
            "eps_t"=>best.εₜ,"bound"=>best.bound),
        "rows"=>[Dict("beta"=>r.β₀,"eps_s"=>r.εₛ,"eps_t"=>r.εₜ,"bound"=>r.bound) for r in results],
        "source_hashes"=>source_hashes(root;dirs=(),files=("model.jl","refine_beta.jl",)))
    write_result(joinpath(root,"results","beta_refinement_D$(numerator(D))_$(denominator(D)).json"),data)
    println("FINE D=$D best=$best"); flush(stdout)
end
