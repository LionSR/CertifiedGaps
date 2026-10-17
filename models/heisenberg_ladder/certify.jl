# Produce one core gap certificate per certified coupling (results/certificate_J<tag>.json)
# and the suite record results/suite.json from fresh verified spectra.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=julia \
#         models/heisenberg_ladder/certify.jl
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
HeisenbergLadder.produce()
