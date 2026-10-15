# Verify results/certificate.json; throws on any mismatch.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder \
#         models/heisenberg_ladder/check.jl            # fast, seconds: no spectra
#     ... check.jl --full                               # from scratch, about 1.5 minutes
#
# The fast check re-derives the hashes, exact RP, every grid decision from its stored rational
# defects, the continuation and the README claims. `--full` also recomputes every verified
# spectrum, the defect and log Z enclosures and the all-β exclusions.
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
"--full" in ARGS ? HeisenbergLadder.check() : HeisenbergLadder.check_stored()
