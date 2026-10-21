/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.OpenTrialCompute

/-!
# Kernel evaluation of the open-chain block iterations, part 2

One `decide +kernel` per block runs the row iteration of `OpenTrialCompute` (240 states, three
packed vector-matrix products per state) and compares its outputs with `openOutputs`. This
file checks the blocks `4`–`10` (spins `2` to `5`). The two parts are independent, so that a
parallel build runs them side by side.
-/

namespace SeedKernel.OpenTrial

set_option maxRecDepth 100000

theorem openBlockCheck_4 : OpenBlockCheck 4 := by decide +kernel
theorem openBlockCheck_5 : OpenBlockCheck 5 := by decide +kernel
theorem openBlockCheck_6 : OpenBlockCheck 6 := by decide +kernel
theorem openBlockCheck_7 : OpenBlockCheck 7 := by decide +kernel
theorem openBlockCheck_8 : OpenBlockCheck 8 := by decide +kernel
theorem openBlockCheck_9 : OpenBlockCheck 9 := by decide +kernel
theorem openBlockCheck_10 : OpenBlockCheck 10 := by decide +kernel

end SeedKernel.OpenTrial
