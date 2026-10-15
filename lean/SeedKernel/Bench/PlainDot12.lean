import SeedKernel.Spin1.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedKernel SeedKernel.Spin1
-- Dot with its already-reversed operand available: one multiplication, no butterfly.
example : lane 128 (3^12-1) (repunit 128 (3^12) * repunit 128 (3^12)) = 3^12 := by
  decide +kernel
