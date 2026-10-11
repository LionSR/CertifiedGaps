using Pkg, LinearAlgebra, CertifiedGaps, CertifiedGaps.Results
BLAS.set_num_threads(1)
root=@__DIR__; repo=normpath(joinpath(root,"../.."))
packages=Dict(string(p.name)=>string(p.version) for p in values(Pkg.dependencies()))
write_result(joinpath(root,"results","migration_2026-10-11_environment.json"),(;julia=string(VERSION),
    blas=string(BLAS.get_config()),blas_threads=BLAS.get_num_threads(),packages,
    source_hashes=file_hashes(repo,["julia/Project.toml","julia/Manifest.toml",
        "AGENTS.md","julia/README.md","outputs/thermal-mpo-gap-scope.tex"]),
    setup="Official Julia 1.12.4 tarball; /tmp/tfim-julia-depot; juliaup endpoint returned HTTP 403."))
