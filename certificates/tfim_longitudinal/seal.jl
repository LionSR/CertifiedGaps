# Fingerprint final audit deliverables; the manifest excludes itself to avoid a cycle.
using CertifiedGaps.Results
const SEAL_ROOT=@__DIR__
const OUTPUT=joinpath(SEAL_ROOT,"results","audit_manifest.json")
paths=String[]
for (dir,_,names) in walkdir(SEAL_ROOT), name in names
    path=joinpath(dir,name)
    path==OUTPUT && continue
    endswith(name,".tmp") && continue
    push!(paths,path)
end
for name in ("../../julia/Project.toml","../../julia/Manifest.toml",
    "../../outputs/thermal-mpo-gap-scope.tex","../../AGENTS.md")
    push!(paths,joinpath(SEAL_ROOT,name))
end
for (dir,_,names) in walkdir(joinpath(SEAL_ROOT,"../../julia/src")), name in names
    endswith(name,".jl") && push!(paths,joinpath(dir,name))
end
write_result(OUTPUT,Dict("status"=>"sealed after independent review",
    "sha256"=>Dict(relpath(p,SEAL_ROOT)=>sha256_file(p) for p in sort(paths))))
println("Sealed ",length(paths)," files")
