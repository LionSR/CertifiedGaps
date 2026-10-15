import SeedPrototype.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedPrototype
-- Dense 128-bit lanes; one matvec, reversal, and dot, with a closed exact answer.
example : dot 128 (3^6) (repunit 128 (3^6))
    (filterApply 128 6 (repunit 128 (3^6))) = 17*6*3^(6-2) := by
  decide +kernel
