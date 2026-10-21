/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalPacked60

/-!
# Compiled evaluation of the `b = 21/2` totals for the twists `1` and `P` (not a proof)

`HornerTotals60` (`Heisenberg/ThermalTable60.lean`) says that the integer recursion of
`thermal120.tex` reproduces the nine printed totals `t_{n,g}`, `n = 6, 8, 10`, `g = 1, P, C`.
`Heisenberg/ThermalPacked60.lean` proves `hornerTotal60_eq_of_packed`: for the twists `1` and
`P`, whenever the packed evaluator `packedHornerTotal60 n g` returns `some W`, then
`hornerTotal60 n g = W`.

The `#guard` commands below run `packedHornerTotal60` with the evaluator of `#eval` for the four
entries with `n = 6, 8` and fail if one differs from the printed `paperTotal60`. This trusts the
compiler, the interpreter and the runtime's big-integer arithmetic, like `native_decide`, so it
proves nothing and no declaration depends on it. It shows that these four instances of
`HornerTotals60` hold for the definitions as written. Measured (laptop, a build running in
parallel): about `3` min for the file, almost all of it at `n = 8`, and `2.9` GB of memory. The
entries with `n = 10` would cost about as
much as `n = 10` of the `b = 49/4` table (`ThermalTotals.lean`, about `20` min and `25` GB per
entry), and the `C` twist is not covered by the packed evaluator, so they are not evaluated here.

Run with `lake env lean code/lean/SeedKernel/Eval/ThermalTotals60.lean` from the repository root
(opt-in library `SeedKernelEval`).
-/

namespace CertifiedGaps.SpinOneHeisenberg.ThermalPacked

-- The twists `1` (`g = 0`) and `P` (`g = 1`) at `n = 6, 8`. The index `g.castSucc` of
-- `paperTotal60` is the same twist.
#guard packedHornerTotal60 6 0 = some (paperTotal60 6 0)
#guard packedHornerTotal60 6 1 = some (paperTotal60 6 1)
#guard packedHornerTotal60 8 0 = some (paperTotal60 8 0)
#guard packedHornerTotal60 8 1 = some (paperTotal60 8 1)

-- A wrong printed integer is detected.
#guard packedHornerTotal60 6 1 ≠ some (paperTotal60 6 1 + 1)

end CertifiedGaps.SpinOneHeisenberg.ThermalPacked
