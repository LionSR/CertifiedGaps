include("src/Ladder.jl")
using .Ladder, ThermalMPOGap, ThermalMPOGap.Results, ThermalMPOGap.Seeds
using ThermalMPOGap.ExactDiagonalization
using LinearAlgebra
BLAS.set_num_threads(1)
g=parse(Rational{Int},get(ARGS,1,"2//1"))
tag=replace(string(g),"//"=>"_")
root=@__DIR__; rp_identity(g); specs=Dict{Int,PartialSpectrum}(); records=[]
for L in [2,3,4,6,8]
    rec=full_spectrum(L,g)
    push!(records,rec); specs[L]=PartialSpectrum(rec.energies,4^L)
    write_result(joinpath(root,"results","spectra_g$tag.json"),(;status="numerical, not certified",records,
        source_hashes=source_hashes(root;files=["scan.jl"])))
end
rs=scan_seeds(spectral_defects(specs),2:4,1/8:1/8:32)
best=best_seed(filter(r->iseven(r.n₀),rs))
summary=[(;n, max_min_purity=maximum(min(1-r.εₛ,1-r.εₜ) for r in rs if r.n₀==n),
    delta=specs[2n].E[1]-2specs[n].E[1]) for n in 2:4]
write_result(joinpath(root,"results","scan_g$tag.json"),(;g=string(g),best,summary,grid=rs,
    metadata=run_metadata((;ns=2:4,beta="1//8:1//8:32//1",max_rungs=8);
        source_sha256=source_hashes(root;files=["scan.jl"]),packages=(ThermalMPOGap,),
        status="numerical, not certified")))
println("BEST ",best); println(summary)
