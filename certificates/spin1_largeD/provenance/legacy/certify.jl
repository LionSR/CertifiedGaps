include("model.jl")
root=@__DIR__
D=length(ARGS)>=1 ? parse(Rational{Int},ARGS[1]) : 4//1
L=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 12
@assert VERSION==v"1.12.4"
rp_identity(D)
data=compute_spectrum(L,D;rigorous=true)
@assert data["complete"] && data["rigorous"]
data["source_hashes"]=source_hashes(root;dirs=("test",),files=("model.jl","certify.jl"))
data["metadata"]=run_metadata(Dict("D"=>string(D),"L"=>L);
    source_sha256=data["source_hashes"],packages=(ThermalMPOGap,IntervalArithmetic),
    status="complete dense eigenvalue enclosures",certified=true)
write_result(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$L.json"),data)
println("DONE D=$D L=$L seconds=$(data["seconds"])")
