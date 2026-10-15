import Lean

set_option maxHeartbeats 0
set_option maxRecDepth 1000000

namespace LeanSeedSmall.Benchmark

-- Dense operands, not multiplication by a power of two. This is a primitive throughput
-- measurement; it is not a proof of packed-dot-product decoding or Hamiltonian semantics.
theorem mul_262144 :
    ((1 <<< 262144 : Nat) - 1) * ((1 <<< 262144) - 1) =
      (1 <<< 524288) - (1 <<< 262145) + 1 := by decide +kernel

theorem mul_4194304 :
    ((1 <<< 4194304 : Nat) - 1) * ((1 <<< 4194304) - 1) =
      (1 <<< 8388608) - (1 <<< 4194305) + 1 := by decide +kernel

#print axioms mul_262144
#print axioms mul_4194304

end LeanSeedSmall.Benchmark
