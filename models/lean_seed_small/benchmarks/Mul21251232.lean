import Lean
set_option maxHeartbeats 0
set_option profiler true
set_option profiler.threshold 1
theorem product :
  ((1 <<< 21251232 : Nat)-1)*((1 <<< 21251232)-3) =
  (1 <<< 42502464)-(1 <<< 21251234)+3 := by decide +kernel
#print axioms product
