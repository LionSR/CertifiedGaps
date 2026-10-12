include("ladder.jl")
using JSON3
J=parse(Rational{Int},ARGS[1]); tag=replace(string(Float64(J)),"."=>"p")
spectra=Dict{Int,PartialSpectrum}(); records=Dict{String,Any}()
start=time()
for L in (2,3,4,5,6,8,10)
    path=joinpath(@__DIR__,"results","spectrum_J$(tag)_L$(L).json")
    if isfile(path) && haskey(JSON3.read(read(path,String)),:energies)
        data=JSON3.read(read(path,String)); es=Float64.(data.energies)
    else
        s=full_spectrum(L,J)
        es=s.energies
        data=Dict("metadata"=>provenance(Dict("J"=>string(J),"L"=>L)),
            "J"=>string(J),"L"=>L,"dimension"=>s.dimension,"sectors"=>s.rows,
            "seconds"=>s.seconds,"energies"=>es,"gap"=>es[2]-es[1])
        write_result(path,data)
    end
    spectra[L]=PartialSpectrum(es,4^L)
    records[string(L)]=Dict("ground_energy"=>es[1],"gap"=>es[2]-es[1],
        "dimension"=>data["dimension"],"seconds"=>data["seconds"],
        "maximum_block_dimension"=>maximum(Int(r["dimension"]) for r in data["sectors"]))
end
rows=scan_seeds(spectral_defects(spectra),2:5,0.25:0.25:16.0)
eligible=filter(r->iseven(r.n₀),rows)
seedrow(r)=r===nothing ? nothing : Dict("n0"=>r.n₀,"beta"=>r.β₀,
    "eps_s"=>r.εₛ,"eps_t"=>r.εₜ,"bound"=>r.bound,"eligible_even_seed"=>iseven(r.n₀))
diagnostics=Any[]
for n in 2:5
    rs=filter(r->r.n₀==n,rows); best=argmax(r->min(1-r.εₛ,1-r.εₜ),rs)
    push!(diagnostics,Dict("n0"=>n,"max_min_purity"=>min(1-best.εₛ,1-best.εₜ),
        "beta_at_max_min"=>best.β₀,"delta"=>spectra[2n].E[1]-2spectra[n].E[1],
        "best_bound_seed"=>seedrow(best_seed(rs)),"eligible_even_seed"=>iseven(n)))
end
write_result(joinpath(@__DIR__,"results","scan_J$(tag).json"),Dict(
    "metadata"=>provenance(Dict("J"=>string(J),"beta_grid"=>"1//4:1//4:16//1")),
    "J"=>string(J),"best"=>seedrow(best_seed(eligible)),"rings"=>records,
    "diagnostics"=>diagnostics,"grid"=>seedrow.(rows),"seconds"=>time()-start))
println("BEST ",best_seed(eligible))
