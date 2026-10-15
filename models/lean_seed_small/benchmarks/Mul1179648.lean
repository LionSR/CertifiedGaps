import Lean
set_option maxHeartbeats 0
set_option profiler true
set_option profiler.threshold 1
theorem product :
  ((1 <<< 1179648 : Nat)-1)*((1 <<< 1179648)-3) =
  (1 <<< 2359296)-(1 <<< 1179650)+3 := by decide +kernel
#print axioms product
