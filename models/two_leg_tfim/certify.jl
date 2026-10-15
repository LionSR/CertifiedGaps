# Certifies the best stored numerical seed at field g (default 3//1): verified dense spectra of
# the rings n0 and 2n0, outward log Z, rational seeds, recurrence (~11 min). Writes OUTPUT
# (default results/certificate_g<g>_certifiedgaps.json, which test/runtests.jl also checks;
# check.jl --full passes a temporary file instead).
#     julia --project=models/two_leg_tfim models/two_leg_tfim/certify.jl [g] [OUTPUT]
include(joinpath(@__DIR__, "src", "Ladder.jl"))
using .Ladder, CertifiedGaps, CertifiedGaps.Results, LinearAlgebra
BLAS.set_num_threads(1)
root=@__DIR__
g=parse(Rational{Int},get(ARGS,1,"3//1"))
tag=replace(string(g),"//"=>"_")
output=get(ARGS,2,joinpath(root,"results","certificate_g$(tag)_certifiedgaps.json"))
scan=read_exact_json(joinpath(root,"results","scan_g$tag.json"))
@assert scan["best"] !== nothing "No eligible numerical seed"
best=scan["best"]; n=parse(Int,best["n₀"])
beta=Rational{BigInt}(best["β₀"]); @assert iseven(n)
t=time(); specs=[]
for L in (n,2n)
    push!(specs,full_spectrum(L,g;verified=true))
end
(;es,et,logZ)=seed_defects(specs[1],specs[2],beta)
@assert 0<=es<1 && 0<=et<1
bound=certified_gap_bound(es,et,beta)
@assert bound.bound !== nothing
payload=(;g=string(g),n0=n,beta=string(beta),epsilon_s=string(es),epsilon_t=string(et),
    bound=bound.bound,recurrence=bound,seconds=time()-t,
    validity="even L >= $n rungs (2L spin sites, multiples of four, >= $(2n) sites)",
    logZ=[map(x->string(Rational{BigInt}(x)),v) for v in logZ],
    spectra=specs,
    metadata=run_metadata((;g=string(g),n0=n,beta=string(beta));
        source_sha256=merge(source_hashes(root;files=["certify.jl"]),
            Dict("../../julia/"*k=>v for (k,v) in source_hashes(joinpath(root,"../../julia")))),
        packages=(CertifiedGaps,),status="certified",certified=true))
write_result(output,payload)
println("CERTIFIED g=$g n0=$n beta=$beta bound=$(bound.bound), seconds=$(time()-t), output=$output")
