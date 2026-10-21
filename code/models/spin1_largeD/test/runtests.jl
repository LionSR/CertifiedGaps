# Unit tests, then the full fast checker.
#
#   OPENBLAS_NUM_THREADS=1 julia --project=julia \
#       models/spin1_largeD/test/runtests.jl
include("frame.jl")      # exact Cartesian/Sᶻ frame, RP data, spectra against dense ED, L ≤ 6
include("transfer.jl")   # transfer-squaring and Weyl inequalities; mutation rejection
include("../check.jl")
