import Lean
set_option maxHeartbeats 0
set_option profiler true
set_option profiler.threshold 1
theorem product :
  ((1 <<< 3706560 : Nat)-1)*((1 <<< 3706560)-3) =
  (1 <<< 7413120)-(1 <<< 3706562)+3 := by decide +kernel
#print axioms product
