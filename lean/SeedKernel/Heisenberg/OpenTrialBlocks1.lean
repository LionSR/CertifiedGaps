/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.OpenTrialCompute

/-!
# Kernel evaluation of the open-chain block iterations, part 1

One `decide +kernel` per block runs the row iteration of `OpenTrialCompute` (240 states, three
packed vector-matrix products per state) and compares its outputs with `openOutputs`. This
file checks the blocks `0`–`3` (spins `0` and `1`, with the largest block, `m = 32`). The two
parts are independent, so that a parallel build runs them side by side.
-/

namespace SeedKernel.OpenTrial

set_option maxRecDepth 100000

theorem openBlockCheck_0 : OpenBlockCheck 0 := by decide +kernel
theorem openBlockCheck_1 : OpenBlockCheck 1 := by decide +kernel
theorem openBlockCheck_2 : OpenBlockCheck 2 := by decide +kernel
theorem openBlockCheck_3 : OpenBlockCheck 3 := by decide +kernel

end SeedKernel.OpenTrial
