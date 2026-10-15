import SeedKernel.Spin1.Optimized
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedKernel SeedKernel.Spin1
example : dot 128 (3^12) (repunit 128 (3^12))
    (filterPrepared 128 (prepare 128 12) (repunit 128 (3^12))) = 17*12*3^10 := by
  decide +kernel
