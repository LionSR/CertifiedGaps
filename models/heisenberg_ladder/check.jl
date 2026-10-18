# Verify results/suite.json and the core certificates; throws on any mismatch.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=julia \
#         models/heisenberg_ladder/check.jl            # fast, seconds: no spectra
#     ... check.jl --full                               # recomputes every verified spectrum
#
# The fast check re-derives the hashes, the exact frame and RP identities, every grid decision
# from its stored rational defects, each core certificate against its point and the Lean seed
# with the exact rational chain (`check_point(...; spectra = false)`), the continuations, the
# coverage and witnesses of the all-β exclusions and the README claims. `--full` also
# recomputes every verified spectrum, every grid defect enclosure, each certificate's seed
# (`check_point` with the fresh spectra) and every exclusion bound.
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
HeisenbergLadder.check(; full = "--full" in ARGS)
