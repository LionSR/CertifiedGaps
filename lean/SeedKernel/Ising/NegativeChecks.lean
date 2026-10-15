import SeedKernel.Ising.Generated

namespace SeedKernel.Ising
set_option maxHeartbeats 0
set_option maxRecDepth 1000000

/-- A false eigenvalue proposal is rejected by the direct kernel computation. -/
theorem wrong_eigenvalue_rejected :
    checkVectors 6 compact6_1 compact6_2 0 lam6_2 = false := by decide +kernel

/-- A malformed vector proposal is rejected, even with the correct eigenvalue numbers. -/
theorem wrong_vector_rejected :
    Packed.check 6 0 v6_2 lam6_1 lam6_2 = false := by decide +kernel

#print axioms wrong_eigenvalue_rejected
#print axioms wrong_vector_rejected
end SeedKernel.Ising
