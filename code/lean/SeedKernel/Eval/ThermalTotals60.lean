/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalPacked60

/-!
# Compiled evaluation of the `b = 21/2` totals (not a proof)

`HornerTotals60` (`Heisenberg/ThermalTable60.lean`) says that the integer recursion of
`thermal120.tex` reproduces the nine printed totals `t_{n,g}`, `n = 6, 8, 10`, `g = 1, P, C`.
`Heisenberg/ThermalPacked60.lean` proves `hornerTotals60_of_packed`: if the packed evaluator
`packedTotal60 n g` returns the printed `paperTotal60 n g` at the nine points, then
`HornerTotals60` holds.

The `#guard` commands below run `packedTotal60` with the evaluator of `#eval` for the six entries
with `n = 6, 8` and fail if one differs from the printed value. This trusts the compiler, the
interpreter and the runtime's big-integer arithmetic, like `native_decide`, so it proves nothing
and no declaration depends on it. It shows that these six instances of `HornerTotals60` hold for
the definitions as written. Measured (laptop): about `10` min for the file and `3.1` GB of
memory, most of it for the twist `C` at `n = 8`. The entries with `n = 10` would cost about as
much as `n = 10` of the `b = 49/4` table (`ThermalTotals.lean`, about `20` min and `25` GB per
entry), and twice that for the twist `C`, whose real form has twice the rows; they are not
evaluated here.

Run with `lake env lean code/lean/SeedKernel/Eval/ThermalTotals60.lean` from the repository root
(opt-in library `SeedKernelEval`).
-/

namespace CertifiedGaps.SpinOneHeisenberg.ThermalPacked

-- The twists `1` (`g = 0`), `P` (`g = 1`) and `C` (`g = 2`) at `n = 6, 8`.
#guard packedTotal60 6 0 = some (paperTotal60 6 0)
#guard packedTotal60 6 1 = some (paperTotal60 6 1)
#guard packedTotal60 6 2 = some (paperTotal60 6 2)
#guard packedTotal60 8 0 = some (paperTotal60 8 0)
#guard packedTotal60 8 1 = some (paperTotal60 8 1)
#guard packedTotal60 8 2 = some (paperTotal60 8 2)

-- A wrong printed integer is detected.
#guard packedTotal60 6 2 ≠ some (paperTotal60 6 2 + 1)

end CertifiedGaps.SpinOneHeisenberg.ThermalPacked
