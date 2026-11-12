# Checker for the spin-1 small-seed certificate.
#
#     julia --project=code/julia code/models/spin1_small_seed/check.jl          # exact parts
#     julia --project=code/julia code/models/spin1_small_seed/check.jl --full   # + all 18 traces
#
# Verifies the model-local source hashes, re-audits the frozen spin-1 trace inputs against
# their recorded SHA-256 hashes and the 10⁻³ tolerance, re-proves the reflection-positivity
# identity exactly, replays the whole exact rational chain from the stored certificate, and
# requires every stored field to match the recomputation. `--full` additionally audits all 18
# source-1 trace rows. Exits with an error on any mismatch.
include(joinpath(@__DIR__, "src", "Spin1SmallSeed.jl"))
using .Spin1SmallSeed

Spin1SmallSeed.check(full = "--full" in ARGS)
