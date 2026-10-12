include("ladder.jl")
using JSON3
rows=Any[]
for tag in ("4p0","3p5","3p0","2p5","2p0","1p5")
    specs=Dict{Int,PartialSpectrum}()
    for L in (2,4,8)
        d=JSON3.read(read(joinpath(@__DIR__,"results","spectrum_J$(tag)_L$(L).json"),String))
        specs[L]=PartialSpectrum(Float64.(d.energies),4^L)
    end
    results=scan_seeds(spectral_defects(specs),(2,4),0.25:0.25:16)
    b=best_seed(results)
    row=Dict("J_tag"=>tag,"best"=>b===nothing ? nothing : Dict("n0"=>b.n₀,
        "beta"=>b.β₀,"eps_s"=>b.εₛ,"eps_t"=>b.εₜ,"bound"=>b.bound),
        "max_min_purity"=>[Dict("n0"=>n,"value"=>maximum(min(1-r.εₛ,1-r.εₜ)
            for r in results if r.n₀==n),"delta"=>specs[2n].E[1]-2specs[n].E[1]) for n in (2,4)])
    push!(rows,row); println(row)
end
write_result(joinpath(@__DIR__,"results","small_grid.json"),Dict(
    "metadata"=>provenance(Dict("ns"=>[2,4],"beta_grid"=>"1//4:1//4:16//1")),"parameters"=>rows))
