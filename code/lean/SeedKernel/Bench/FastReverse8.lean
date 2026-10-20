import SeedKernel.Spin1.Optimized
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedKernel SeedKernel.Spin1
-- A non-palindromic regenerated trial, not merely a uniform-vector reversal test.
example : reverseFast 128 (3^8) (trial 128 8 3) =
    reverseLanes 128 (3^8) (trial 128 8 3) := by decide +kernel
