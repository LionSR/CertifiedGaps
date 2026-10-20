/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalPacked

/-!
# The periodic integer totals with `n = 4` in the kernel

`HornerTotals` (`Heisenberg/ThermalTable.lean`) is the hypothesis that the integer recursion of
`thermal72.tex` reproduces the eighteen printed totals `W_{n,g}`, `n = 4, …, 12`.
`Heisenberg/ThermalPacked.lean` computes the same totals with the packed evaluator
`packedHornerTotal` and proves `hornerTotal_eq_of_packed`. This file checks the two entries with
`n = 4` by `decide +kernel` (each a closed run of the packed recursion on the `13`
representative columns of the `81` words), so that `HornerTotals` reduces to the sixteen totals
with `n ≥ 5` (`HornerTotalsLarge`, `hornerTotals_of_large`).

Kernel time (laptop, other jobs running): about `11` s per entry, the whole file about `30` s,
with a peak of about `7` GB (of which about `3` GB is the imported library). The entries with
`n ≥ 5` are left as hypotheses: the kernel keeps every intermediate packed row until the end of
a declaration, and `n = 5` has three times the rows and `22` columns instead of `13`.
`SeedKernel/Eval/ThermalTotals.lean` evaluates the entries with `n ≤ 9` with compiled code
(not a proof).
-/

namespace CertifiedGaps.SpinOneHeisenberg

open ThermalPacked

theorem packedHornerTotal_four_zero : packedHornerTotal 4 0 = some (paperTotal 4 0) := by
  decide +kernel

theorem packedHornerTotal_four_one : packedHornerTotal 4 1 = some (paperTotal 4 1) := by
  decide +kernel

/-- The two totals with `n = 4`, from the kernel-checked packed runs. -/
theorem hornerTotal_four (g : Fin 2) : hornerTotal 4 g = paperTotal 4 g := by
  fin_cases g
  · exact hornerTotal_eq_of_packed (by norm_num) packedHornerTotal_four_zero
  · exact hornerTotal_eq_of_packed (by norm_num) packedHornerTotal_four_one

/-- **The remaining hypothesis**: the sixteen printed totals with `n = 5, …, 12`. -/
def HornerTotalsLarge : Prop :=
  ∀ n ∈ [5, 6, 7, 8, 9, 10, 11, 12], ∀ g : Fin 2, hornerTotal n g = paperTotal n g

/-- The two kernel-checked totals and the sixteen remaining ones give all eighteen. -/
theorem hornerTotals_of_large (h : HornerTotalsLarge) : HornerTotals := by
  intro n hn g
  simp only [tableLengths, List.mem_cons, List.not_mem_nil, or_false] at hn
  rcases hn with rfl | hn
  · exact hornerTotal_four g
  · exact h n (by simp only [List.mem_cons, List.not_mem_nil, or_false]; omega) g

end CertifiedGaps.SpinOneHeisenberg
