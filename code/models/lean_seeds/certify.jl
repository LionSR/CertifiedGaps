# Regenerate the exact records of the Lean-seed studies and the Lean witness files under
# lean/SeedKernel. Nothing written here is a gap certificate (see README "What is rigorous").
#
#     julia --project=. certify.jl [study ...]     # default: every study below, in order
#
#   seed_design   results/seed_design/{certificate,compact_witnesses}.json        ~10 s
#   anchor        results/anchor/certificate.json                                  ~10 s
#   spin1_kernel  lean/SeedKernel/Spin1/H{6,8,10}.lean, results/spin1_kernel/
#                 {exact_witnesses.tsv,witnesses.json}                            ~2 min
#   ising         results/ising/certificate.json, lean/SeedKernel/Ising/
#                 Generated{,Chunks}.lean                                          ~15 s
#
# The spectrum study is NON-CERTIFYING Float64 data; its drivers are under scripts/spectrum.
# Run spin1_kernel last when other sources changed: it records the hashes of the shared drivers.
include(joinpath(@__DIR__, "src", "LeanSeeds.jl"))
using .LeanSeeds

const TASKS = [
    "seed_design" => LeanSeeds.SeedDesign.certify_seed_design,
    "anchor" => LeanSeeds.AnchorStudy.certify_anchor,
    "ising" => LeanSeeds.IsingSeed.certify_ising,
    "spin1_kernel" => LeanSeeds.Spin1Kernel.certify_spin1_kernel,
]

function main(args)
    known = first.(TASKS)
    unknown = setdiff(args, known)
    isempty(unknown) || error("unknown study $(unknown); expected any of $(known)")
    for (name, task) in TASKS
        isempty(args) || name in args || continue
        println("== certify ", name)
        t = @elapsed task()
        println("== ", name, " done in ", round(t; digits=1), " s")
    end
end

main(ARGS)
