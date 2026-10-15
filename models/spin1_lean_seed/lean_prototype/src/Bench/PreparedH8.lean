import SeedPrototype.Optimized
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedPrototype
-- End-to-end equality of the optimized and original arithmetic moments.
example : momentsPrepared 345 8 32 = moments 345 8 32 := by decide +kernel
