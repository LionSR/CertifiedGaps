/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialChecksInsertion1

/-!
# Kernel check of `O W = W C`, blocks `4 ≤ b < 5`

The packed comparison for the columns of the blocks `4 ≤ b < 5` on the radix vector of
width `71` (see `TrialChecks` and `O_mul_Wm`). The import of `TrialChecksInsertion1` only orders
the kernel runs, so that a parallel build does not hold several of them in memory at once.
-/

namespace SeedKernel.Trial

set_option maxRecDepth 100000

/-- `O W = W C` on the `130` columns of the blocks `4 ≤ b < 5`, packed rows (slots of
`130 · 71` bits). -/
theorem checkO2 : OPack 9230 (YG 71 4 5) = ZPack 9230 (ZG 71 4 5 oB) := by
  decide +kernel

end SeedKernel.Trial
