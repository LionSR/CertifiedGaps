using CertifiedGaps, CertifiedGaps.Results, LinearAlgebra, Pkg
include("model.jl")
BLAS.set_num_threads(1)
@assert VERSION==v"1.12.4"
r=exact_rp()
root=@__DIR__
ratrows(M)=[[string(M[i,j]) for j in axes(M,2)] for i in axes(M,1)]
write_result(joinpath(root,"results/rp.json"),(;
    metadata=run_metadata((;claim="local Cartesian reflection positivity only");
        source_sha256=source_hashes(root;dirs=(),files=("audit.jl","model.jl")),
        status="exact local coefficient identities proved; no uniform gap certificate"),
    G=ratrows.(r.G),Q=[(;a,b,entries=ratrows(r.Q[a,b])) for a in 1:3 for b in 1:3],
    onsite=ratrows(r.onsite),coefficient_G="cos(theta)-sin(theta)/2",
    coefficient_Q="-sin(theta)",parameter_domain="-pi/4 < theta <= 0",
    exact_identity="X=-A; X^2=A/2+B",identity_verified=true,uniform_gap_certified=false))
versions=Dict(string(info.name)=>string(info.version) for info in values(Pkg.dependencies()) if info.version!==nothing)
resource_files=Dict(name=>strip(read(joinpath("/sys/fs/cgroup",name),String))
    for name in ("memory.max","memory.peak","memory.events","cpu.max")
    if isfile(joinpath("/sys/fs/cgroup",name)))
write_result(joinpath(root,"results/environment.json"),(;
    metadata=run_metadata((;memory_limit_bytes=17179869184,cpu_quota=4);
        source_sha256=source_hashes(root;dirs=("test",),files=[filter(f->endswith(f,".jl"),readdir(root));"run.sh";"run_su2_grid.sh";"Project.toml"]),
        status="environment ready; Julia version asserted and dependencies loaded"),
    packages=versions,core_environment_hashes=file_hashes(joinpath(root,"../../julia"),("Project.toml","Manifest.toml")),
    resource_observations=(;scope="entire environment lifetime, including pilots; not a per-sector peak",
        cgroup=resource_files),
    setup="juliaup installer HTTP 403; official julialang-s3 Julia 1.12.4 tarball succeeded; instantiate/precompile succeeded"))
println("Environment audit and exact RP result written")
