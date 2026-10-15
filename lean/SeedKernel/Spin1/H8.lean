import SeedKernel.Spin1.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace SeedKernel.Spin1
-- Arithmetic witnesses, NOT a theorem about ordered physical eigenvalues.
theorem h8_ritz_witness : ritzCheck 345 8 32 11179504 1000000 279904 1000000 = true := by
  decide +kernel
theorem h8_lane_bound : 6561 * ((32 : Nat)^66) < 2^345 := by
  decide +kernel
#print axioms h8_ritz_witness
end SeedKernel.Spin1
