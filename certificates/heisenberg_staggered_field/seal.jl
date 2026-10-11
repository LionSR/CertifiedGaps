# Seal the saved deliverables after the independent audit; no computation/search.
using ThermalMPOGap.Results
const root=@__DIR__
const output=joinpath(root,"results","audit_manifest.json")
paths=sort([joinpath(folder,name) for (folder,_,names) in walkdir(root) for name in names
    if joinpath(folder,name)!=output && !endswith(name,".tmp")])
write_result(output,Dict("description"=>"SHA-256 of every saved deliverable except this manifest",
    "file_sha256"=>Dict(relpath(path,root)=>sha256_file(path) for path in paths)))
println("Sealed ",length(paths)," files")
