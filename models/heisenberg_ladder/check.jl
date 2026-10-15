# Verify results/certificate.json from scratch: recomputes every spectrum and decision and
# throws on any mismatch (about 5 minutes).
#
#     OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder \
#         models/heisenberg_ladder/check.jl
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
HeisenbergLadder.check()
