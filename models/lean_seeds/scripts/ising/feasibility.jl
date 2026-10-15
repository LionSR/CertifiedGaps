# NON-CERTIFYING: test the proposed two-eigenvalue plus classical-anchor sandwich.
# The exact integer residuals are real evidence, but the eigenvalue theorem is not in Lean yet.
using JSON, Printf
const HERE=joinpath(@__DIR__, "..", "..", "results", "ising")
c=JSON.parsefile(joinpath(HERE,"certificate.json"))
beta=97/4
eta=1e-10 # deliberately generous compared with exact normalized residuals
function zinterval(p,b)
    e=Float64.(p["numerical_E"])
    rho=Float64(parse(BigInt,split(p["anchor_20H"],"//")[1]))/20; d=2^p["n"]
    # Shift by e[1] for stable Float64 calculations; retain the log shift in output.
    lo=sum(exp.(-b.*(e.-e[1].+eta)))
    tail=(d-2)*exp(-b*(rho-e[1]))
    hi=sum(exp.(-b.*(e.-e[1].-eta)))+tail
    (lo=log(lo)-b*e[1], hi=log(hi)-b*e[1], tail=tail, cutoff=rho-e[1])
end
a=zinterval(c["points"][1],beta)
b=zinterval(c["points"][2],beta)
t=zinterval(c["points"][2],2beta)
es=-expm1(b.lo-2a.hi); et=-expm1(t.lo-2b.hi)
fits=a.hi < 2726798329802173/17592186044416 &&
     b.lo > 5453596500152007/17592186044416 &&
     b.hi < 2726798388197837/8796093022208 &&
     t.lo > 2726798250076001/4398046511104
out=Dict("status"=>"NON-CERTIFYING Float64 feasibility; assumes eigenvalue enclosure theorem",
    "energy_radius"=>eta,"small"=>a,"large"=>b,"double_beta"=>t,
    "eps_s_sandwich"=>es,"eps_t_sandwich"=>et,
    "eps_s_target"=>9064e-9,"eps_t_target"=>3141e-8,
    "fits_log_enclosures"=>fits,"passes"=>es<9064e-9 && et<3141e-8)
open(joinpath(HERE,"feasibility.json"),"w") do io; JSON.print(io,out); end
println(JSON.json(out))
