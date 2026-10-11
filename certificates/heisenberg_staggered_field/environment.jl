include("common.jl")
using Pkg
@assert VERSION==v"1.12.4"
certfiles=filter(f->startswith(f,"certificate_") && endswith(f,".json"),readdir(joinpath(ROOT,"results")))
for f in certfiles
    c=JSON3.read(read(joinpath(ROOT,"results",f),String),Dict{String,Any})
    for (name,hash) in c["metadata"]["source_sha256"]
        path=startswith(name,"core/") ? joinpath(REPO,"julia",name[6:end]) : joinpath(ROOT,name)
        @assert sha256_file(path)==hash "source changed since certificate: $name"
    end
end
write_result(joinpath(ROOT,"results","environment.json"),Dict(
    "julia"=>string(VERSION),"git_commit"=>readchomp(`git -C $REPO rev-parse HEAD`),
    "packages"=>Dict(x.name=>string(x.version) for x in values(Pkg.dependencies())),
    "project_hashes"=>file_hashes(joinpath(REPO,"julia"),["Project.toml","Manifest.toml"]),
    "source_hashes_match"=>true,"certificates_audited"=>certfiles,
    "blas"=>string(BLAS.get_config()),"blas_threads"=>BLAS.get_num_threads()))
println("Julia 1.12.4; all stored certificate source hashes match current sources.")
