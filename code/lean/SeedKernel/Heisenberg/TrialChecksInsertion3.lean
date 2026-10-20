/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialChecksInsertion2

/-!
# Kernel check of `O W = W C`, blocks `5 ≤ b < 7`

The packed comparison for the columns of the blocks `5 ≤ b < 7` on the radix vector of
width `71` (see `TrialChecks` and `O_mul_Wm`). The import of `TrialChecksInsertion2` only orders
the kernel runs, so that a parallel build does not hold several of them in memory at once.
-/

namespace SeedKernel.Trial

set_option maxRecDepth 100000

/-- `O W = W C` on the `180` columns of the blocks `5 ≤ b < 7`, packed rows (slots of
`180 · 71` bits). -/
theorem checkO3 : OPack 12780 (YG 71 5 7) = ZPack 12780 (ZG 71 5 7 oB) := by
  decide +kernel

end SeedKernel.Trial
