# Check every study of the Lean-seed project. Julia never decides the truth of a Lean statement:
# it replays the exact records, compares the generated Lean files with the stored witnesses, and
# runs the Lean kernel on lean/SeedKernel (`lake build SeedKernel` plus the axiom audit).
#
#     julia --project=. check.jl            # fast: exact replays, data checks, Lake build + audit,
#                                           # compiled evaluation of the integer hypotheses
#     julia --project=. check.jl --full     # also regenerate the spin-1 witnesses from scratch,
#                                           # rebuild SeedKernel from clean, run every benchmark
#
# Set LAKE_BIN if `lake` is not on PATH. Throws (nonzero exit) on any mismatch.
include(joinpath(@__DIR__, "src", "LeanSeeds.jl"))
using .LeanSeeds
using .LeanSeeds.Common: build_seed_kernel, lean_output, LAKE, LEAN

"Run one named check and report its wall time."
function step(f, name)
    println("== ", name)
    t = @elapsed f()
    println("== ", name, " passed in ", round(t; digits=1), " s")
end

function main(full)
    step("seed_design: exact designs and ChainCertificate.Valid") do
        LeanSeeds.SeedDesign.check_seed_design()
    end
    step("spectrum: Float64 data consistency (NOT a certificate)") do
        LeanSeeds.Spectrum.check_hamiltonian()
        LeanSeeds.Spectrum.check_data()
    end
    step("anchor: exact local PSD, enclosures, counts and tails") do
        LeanSeeds.AnchorStudy.check_anchor()
    end
    step("spin1_kernel: hashes and witness files") do
        LeanSeeds.Spin1Kernel.check_spin1_kernel(; full)
    end
    step("ising: exact replay and generated Lean files") do
        LeanSeeds.IsingSeed.check_ising()
    end
    step("lean: lake build SeedKernel and axiom audit") do
        seconds, names = build_seed_kernel(; fresh=full)
        println("lake build SeedKernel", full ? " (clean)" : "", " took ", round(seconds; digits=1),
            " s; ", length(names), " declarations use only standard axioms")
    end
    # Compiled code, not the kernel: shows that the integer hypotheses hold for the Lean
    # definitions as written, without proving them (SeedKernel/Eval/*.lean).
    step("lean: compiled evaluation of the integer hypotheses (SeedKernel/Eval, not a proof)") do
        for file in sort(readdir(joinpath(LEAN, "SeedKernel", "Eval")))
            t = @elapsed lean_output(`$LAKE env lean $(joinpath("SeedKernel", "Eval", file))`)
            println("Eval/$file passed in $(round(t; digits=1)) s")
        end
    end
    full || return
    step("lean: kernel benchmarks (SeedKernel/Bench)") do
        for file in sort(readdir(joinpath(LEAN, "SeedKernel", "Bench")))
            t = @elapsed lean_output(`$LAKE env lean $(joinpath("SeedKernel", "Bench", file))`)
            println("Bench/$file checked in $(round(t; digits=1)) s")
        end
    end
end

main("--full" in ARGS)
