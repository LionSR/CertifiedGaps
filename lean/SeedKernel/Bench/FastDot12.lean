import SeedKernel.Spin1.Optimized
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedKernel SeedKernel.Spin1
example : dotFast 128 (3^12) (repunit 128 (3^12)) (repunit 128 (3^12)) = 3^12 := by
  decide +kernel
