# Low-energy M = 0 reference levels at L = 14, 16 (KrylovKit, six levels per symmetry block).
# Numerical reference only: these incomplete spectra are NEVER used in partition functions.
# Usage: models/blbq_haldane/run.sh models/blbq_haldane/scripts/reference.jl 14
using CertifiedGaps, CertifiedGaps.Results, KrylovKit, LinearAlgebra, Random
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "sectors.jl"))
BLAS.set_num_threads(1)
const THETAS = [0.0, -0.1, -0.3, -0.5, -0.6, -0.7]

"Lowest six levels of every M = 0 symmetry block of the length-`L` ring; write the summary."
function reference(L)
    start = time(); o = orbits(L, 0); rows = Any[]; levels = [Float64[] for _ in THETAS]
    for k in 0:L÷2, block in reduced_components(L, 0, k, o)
        d = size(block.A, 1); nev = min(6, d)
        for t in eachindex(THETAS)
            H = cos(THETAS[t]) * block.A + sin(THETAS[t]) * block.B
            Random.seed!(1000L + 100k + t)
            elapsed = @elapsed begin
                E, _, info = eigsolve(H, randn(d), nev, :SR; ishermitian=true,
                    krylovdim=min(d, 80), tol=1e-11, maxiter=300)
                @assert info.converged >= nev
            end
            # Count each level once: M = 0 contains one state per integer-spin multiplet.
            append!(levels[t], real(E[1:nev]))
            push!(rows, (; k, symmetry=block.label, theta=THETAS[t], dimension=d,
                energies=real(E[1:nev]), residual_norms=info.normres[1:nev], seconds=elapsed))
            println("reference L=$L k=$k $(block.label) theta=$(THETAS[t]) E0=$(E[1]) ",
                "seconds=$elapsed"); flush(stdout)
        end
    end
    foreach(sort!, levels)
    hashes = source_hashes(ROOT; dirs=(),
        files=("scripts/reference.jl", "src/sectors.jl", "src/model.jl"))
    metadata = run_metadata((; L, method="KrylovKit low multiplets, NEVER used for Z");
        source_sha256=hashes, status="numerical reference, not certified")
    summary = [(; theta=THETAS[t], ground=levels[t][1], gap=levels[t][2] - levels[t][1])
        for t in eachindex(THETAS)]
    write_result(joinpath(ROOT, "results", "reference_L$(L).json"),
        (; metadata, L, seconds=time() - start, sectors=rows, summary))
end

for L in parse.(Int, ARGS)
    reference(L)
end
