using CertifiedGaps.Results,JSON3,LinearAlgebra
BLAS.set_num_threads(1)
root=@__DIR__;name="su2_benchmark_L14_S2_k1.json"
raw=joinpath(root,"historical",name)
isfile(raw) || cp(joinpath(root,"results",name),raw)
mutable(x)=x isa JSON3.Object ? Dict{String,Any}(string(k)=>mutable(v) for (k,v) in pairs(x)) :
    x isa JSON3.Array ? [mutable(v) for v in x] : x
d=mutable(JSON3.read(read(raw,String)))
old=d["metadata"]["source_sha256"]["su2.jl"]
d["metadata"]["post_run_disk_su2_sha256"]=old
d["metadata"]["source_sha256"]["su2.jl"]=sha256_file(joinpath(root,"historical/su2_householder.jl"))
d["metadata"]["source_hash_note"]="The pilot loaded the Householder version before optimization edits. Original end-of-run disk metadata is preserved in historical/. The corrected loaded-source hash matches the independently recorded L12 pilot hash."
d["row"]["diagnostics"]["orthonormalization"]="HouseholderQR"
write_result(joinpath(root,"results",name),d)
