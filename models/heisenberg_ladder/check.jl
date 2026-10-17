# Verify results/suite.json and the core certificates; throws on any mismatch.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=julia \
#         models/heisenberg_ladder/check.jl            # fast, seconds: no spectra
#     ... check.jl --full                               # recomputes every verified spectrum
#
# The fast check re-derives the hashes, the exact frame and RP identities, every grid decision
# from its stored rational defects, the exact rational chain of each core certificate
# (`check_certificate(cert; spectra = false)`), the continuations and the README claims.
# `--full` also recomputes every verified spectrum, every grid defect enclosure, each
# certificate's seed (`check_certificate` with the fresh spectra) and the all-β exclusions.
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
HeisenbergLadder.check(; full = "--full" in ARGS)
