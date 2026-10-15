import SeedPrototype.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedPrototype
-- Dense 128-bit lanes; one matvec, reversal, and dot, with a closed exact answer.
example : dot 128 (3^12) (repunit 128 (3^12))
    (filterApply 128 12 (repunit 128 (3^12))) = 17*12*3^(12-2) := by
  decide +kernel
