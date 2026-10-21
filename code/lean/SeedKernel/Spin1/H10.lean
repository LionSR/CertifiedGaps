import SeedKernel.Spin1.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace SeedKernel.Spin1
-- Arithmetic witnesses, NOT a theorem about ordered physical eigenvalues.
theorem h10_ritz_witness : ritzCheck 222 10 16 13548710 1000000 632355 1000000 = true := by
  decide +kernel
theorem h10_lane_bound : 59049 * ((40 : Nat)^34) < 2^222 := by
  decide +kernel
#print axioms h10_ritz_witness
end SeedKernel.Spin1
