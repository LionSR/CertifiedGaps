# Numerical seed scan (not certified): floating-point spectra of the rings of 2, 3, 4, 6, 8
# rungs at field g from the core's exact block reduction (`Ladder.ladder_energies`) and a β
# grid; writes results/spectra_g<g>.json (gitignored) and results/scan_g<g>.json.
#     julia --project=julia models/two_leg_tfim/scripts/scan.jl 3//1
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "src", "Ladder.jl"))
using .Ladder, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds
using CertifiedGaps.ExactDiagonalization
using LinearAlgebra
BLAS.set_num_threads(1)
g=parse(Rational{Int},get(ARGS,1,"2//1"))
tag=replace(string(g),"//"=>"_")
root=ROOT; specs=Dict{Int,PartialSpectrum}(); records=[]
for L in [2,3,4,6,8]
    seconds=@elapsed energies=ladder_energies(g,L)
    push!(records,(;L,g=string(g),energies,seconds)); specs[L]=PartialSpectrum(energies,4^L)
    println("L=$L g=$g seconds=$(round(seconds,digits=2))"); flush(stdout)
    write_result(joinpath(root,"results","spectra_g$tag.json"),
        (;status="numerical, not certified",records,
        source_hashes=source_hashes(root)))
end
rs=scan_seeds(spectral_defects(specs),2:4,1/8:1/8:32)
best=best_seed(filter(r->iseven(r.n₀),rs))
summary=[(;n, max_min_purity=maximum(min(1-r.εₛ,1-r.εₜ) for r in rs if r.n₀==n),
    delta=specs[2n].E[1]-2specs[n].E[1]) for n in 2:4]
write_result(joinpath(root,"results","scan_g$tag.json"),(;g=string(g),best,summary,grid=rs,
    metadata=run_metadata((;ns=2:4,beta="1//8:1//8:32//1",max_rungs=8);
        source_sha256=source_hashes(root),packages=(CertifiedGaps,),
        status="numerical, not certified")))
println("BEST ",best); println(summary)
