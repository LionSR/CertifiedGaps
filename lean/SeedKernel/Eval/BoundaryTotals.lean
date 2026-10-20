/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryData

/-!
# Compiled evaluation of the boundary-field integers (not a proof)

`SeedKernel/Heisenberg/BoundaryTraces.lean` proves the nine entries `k_{n,g}` with `n ≤ 4` in the
kernel and keeps the six with `n = 5, 6` as the hypothesis `BoundaryTotalsLarge`, because the
kernel would need about `10` GB and half an hour per `n = 5` entry.

The `#guard` commands below evaluate all fifteen entries of the Lean definition `quotient` with
compiled code (about `25` s in all) and fail if one differs from the printed `paperK`. This
trusts the compiler and the runtime's big-integer arithmetic, like `native_decide`, so it
proves nothing and no declaration depends on it. It shows that `BoundaryTotalsLarge` holds for
the definitions as written: the remaining hypothesis is a fact about the Lean code, not a hope.

Run with `lake env lean SeedKernel/Eval/BoundaryTotals.lean` in `lean/` (opt-in library
`SeedKernelEval`; `models/lean_seeds/check.jl` runs it).
-/

namespace CertifiedGaps.SpinOneHeisenberg.BoundaryThermal

-- The six entries of `BoundaryTotalsLarge`.
#guard quotient 5 1 = some (paperK 5 1)
#guard quotient 5 (-1) = some (paperK 5 (-1))
#guard quotient 5 0 = some (paperK 5 0)
#guard quotient 6 1 = some (paperK 6 1)
#guard quotient 6 (-1) = some (paperK 6 (-1))
#guard quotient 6 0 = some (paperK 6 0)

-- The whole table, including the nine entries the kernel proves.
#guard lengths.all fun n ↦ fields.all fun g ↦ quotient n g == some (paperK n g)

-- A wrong printed integer is detected.
#guard quotient 6 0 ≠ some (paperK 6 0 + 1)

end CertifiedGaps.SpinOneHeisenberg.BoundaryThermal
