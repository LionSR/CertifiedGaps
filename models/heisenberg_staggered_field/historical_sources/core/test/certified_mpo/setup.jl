# Shared loader for the CertifiedMPO tests; each test file also runs on its own:
#   julia --project=julia julia/test/certified_mpo/mpo_norms.jl
using ThermalMPOGap
const CM = ThermalMPOGap.CertifiedMPO
