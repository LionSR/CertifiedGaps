import SeedKernel.Spin1.Optimized
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedKernel SeedKernel.Spin1
-- End-to-end equality of the optimized and original arithmetic moments.
example : momentsPrepared 345 8 32 = moments 345 8 32 := by decide +kernel
