import Lean
set_option maxHeartbeats 0
set_option profiler true
set_option profiler.threshold 1
theorem product :
  ((1 <<< 262144 : Nat)-1)*((1 <<< 262144)-3) =
  (1 <<< 524288)-(1 <<< 262146)+3 := by decide +kernel
#print axioms product
