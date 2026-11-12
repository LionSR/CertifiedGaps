# Generates the spin-1 small-seed certificate: results/certificate.json.
#
#     julia --project=code/julia code/models/spin1_small_seed/certify.jl
#
# Recomputes every exact scalar step from the hash-audited frozen spin-1 trace inputs and
# writes the certificate (exact rational defects, doubling chain, and gap bound) into
# results/ (a symlink to data/models/spin1_small_seed/).
include(joinpath(@__DIR__, "src", "Spin1SmallSeed.jl"))
using .Spin1SmallSeed

Spin1SmallSeed.certify()
