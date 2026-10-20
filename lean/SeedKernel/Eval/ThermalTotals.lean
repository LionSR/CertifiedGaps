/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalPacked

/-!
# Compiled evaluation of the periodic integer totals (not a proof)

`HornerTotals` (`Heisenberg/ThermalTable.lean`) says that the integer recursion of
`thermal72.tex` reproduces the eighteen printed totals `W_{n,g}`, `n = 4, …, 12`.
`Heisenberg/ThermalPacked.lean` proves `hornerTotal_eq_of_packed`: whenever the packed evaluator
`packedHornerTotal n g` returns `some W`, then `hornerTotal n g = W`.

The `#guard` commands below run `packedHornerTotal` with the evaluator of `#eval` for the twelve
entries with `n ≤ 9` and fail if one differs from the printed `paperTotal`. This trusts the
compiler, the interpreter and the runtime's big-integer arithmetic, like `native_decide`, so it
proves nothing and no declaration depends on it. It shows that these twelve instances of
`HornerTotals` hold for the definitions as written. Measured times (laptop, other jobs
running): `n ≤ 6` under `2` s per entry, `n = 7` about `5` s, `n = 8` about `22` s, `n = 9`
about `150` s per entry; about `5.5` min for the file. The cost grows about eightfold per site
(`3ⁿ` rows times about `3ⁿ/(4n)` representative columns: `630` at `n = 9`, `1652` at
`n = 10`), and the memory with it (about `3` GB at `n = 9`). Extrapolated, `n = 10` would take
about `20` min and `25` GB per entry, `n = 11` about `3` h and `n = 12` about a day, so they
are not evaluated here.

Run with `lake env lean lean/SeedKernel/Eval/ThermalTotals.lean` from the repository root
(opt-in library `SeedKernelEval`).
-/

namespace CertifiedGaps.SpinOneHeisenberg.ThermalPacked

-- The four entries that `SeedKernel/Heisenberg/ThermalTraces.lean` also proves in the kernel.
#guard packedHornerTotal 4 0 = some (paperTotal 4 0)
#guard packedHornerTotal 4 1 = some (paperTotal 4 1)

-- The entries with `5 ≤ n ≤ 9`.
#guard packedHornerTotal 5 0 = some (paperTotal 5 0)
#guard packedHornerTotal 5 1 = some (paperTotal 5 1)
#guard packedHornerTotal 6 0 = some (paperTotal 6 0)
#guard packedHornerTotal 6 1 = some (paperTotal 6 1)
#guard packedHornerTotal 7 0 = some (paperTotal 7 0)
#guard packedHornerTotal 7 1 = some (paperTotal 7 1)
#guard packedHornerTotal 8 0 = some (paperTotal 8 0)
#guard packedHornerTotal 8 1 = some (paperTotal 8 1)
#guard packedHornerTotal 9 0 = some (paperTotal 9 0)
#guard packedHornerTotal 9 1 = some (paperTotal 9 1)

-- A wrong printed integer is detected.
#guard packedHornerTotal 5 1 ≠ some (paperTotal 5 1 + 1)

end CertifiedGaps.SpinOneHeisenberg.ThermalPacked
