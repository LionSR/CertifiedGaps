# Produce results/certificate.json from fresh verified spectra (about 6 minutes).
#
#     OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder \
#         models/heisenberg_ladder/certify.jl
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
HeisenbergLadder.produce()
