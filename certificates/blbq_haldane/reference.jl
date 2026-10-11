# Numerical reference only: these incomplete spectra are NEVER used in partition functions.
using CertifiedGaps, CertifiedGaps.Results, KrylovKit, LinearAlgebra, Random
include("sectors.jl")
BLAS.set_num_threads(1)
const THETAS=[0.0,-0.1,-0.3,-0.5,-0.6,-0.7]
function reference(L)
    start=time();o=orbits(L,0);rows=Any[]; levels=[Float64[] for _ in THETAS]
    for k in 0:L÷2
        blocks=reduced_components(L,0,k,o)
        for block in blocks
            d=size(block.A,1)
            for t in eachindex(THETAS)
                H=cos(THETAS[t])*block.A+sin(THETAS[t])*block.B
                Random.seed!(1000L+100k+t)
                elapsed=@elapsed begin
                    E,Q,info=eigsolve(H,randn(d),min(6,d),:SR;ishermitian=true,
                        krylovdim=min(d,80),tol=1e-11,maxiter=300)
                    @assert info.converged>=min(6,d)
                end
                # Count each eigenvalue once: magnetization zero contains one state per SU(2) multiplet.
                append!(levels[t],real(E[1:min(6,d)]))
                push!(rows,(;k,symmetry=block.label,theta=THETAS[t],dimension=d,
                    energies=real(E[1:min(6,d)]),residual_norms=info.normres[1:min(6,d)],seconds=elapsed))
                println("reference L=$L k=$k $(block.label) theta=$(THETAS[t]) E0=$(E[1]) seconds=$elapsed");flush(stdout)
            end
        end
    end
    for e in levels;sort!(e);end
    hashes=source_hashes(@__DIR__;dirs=(),files=("reference.jl","sectors.jl","model.jl"))
    data=(;metadata=run_metadata((;L,method="KrylovKit low multiplets, NEVER used for Z");
        source_sha256=hashes,status="numerical reference, not certified"),L,seconds=time()-start,sectors=rows,
        summary=[(;theta=THETAS[t],ground=levels[t][1],gap=levels[t][2]-levels[t][1]) for t in eachindex(THETAS)])
    write_result(joinpath(@__DIR__,"results","reference_L$(L).json"),data)
end
for L in parse.(Int,ARGS);reference(L);end
