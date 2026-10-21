import SeedKernel.Spin1.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace SeedKernel.Spin1
-- Arithmetic witnesses, NOT a theorem about ordered physical eigenvalues.
theorem h6_ritz_witness : ritzCheck 262 6 24 8495447 1000000 274967 1000000 = true := by
  decide +kernel
theorem h6_lane_bound : 729 * ((24 : Nat)^50) < 2^262 := by
  decide +kernel
#print axioms h6_ritz_witness
end SeedKernel.Spin1
