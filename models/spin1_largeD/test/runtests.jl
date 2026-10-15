# Unit tests, then the full checker (about one minute in total).
#
#   OPENBLAS_NUM_THREADS=1 julia --project=models/spin1_largeD \
#       models/spin1_largeD/test/runtests.jl
include("small.jl")      # dense blocks against the full product-space Hamiltonian, L <= 6
include("symmetry.jl")   # exact symmetry covariance and orbit counting, L <= 12
include("reviewer.jl")   # Cartesian/Sz basis bridge and BigFloat residual diagnostics
include("transfer.jl")   # transfer-squaring, tail and Weyl inequalities; mutation rejection
include("../check.jl")
