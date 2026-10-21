import SeedKernel.Spin1.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedKernel SeedKernel.Spin1
-- Dense 128-bit lanes; one matvec, reversal, and dot, with a closed exact answer.
example : dot 128 (3^8) (repunit 128 (3^8))
    (filterApply 128 8 (repunit 128 (3^8))) = 17*8*3^(8-2) := by
  decide +kernel
