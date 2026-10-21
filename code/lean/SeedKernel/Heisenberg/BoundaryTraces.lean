/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryThermal

/-!
# The boundary-field thermal table in the kernel

`CertifiedGaps.Heisenberg.BoundaryThermal` proves the fifteen rational enclosures of the
short-chain traces `Z_{n,g}` (`BoundaryTraceTable`) from `BoundaryTotals`: the packed integer
recurrence of the thermal appendix reproduces the paper's printed integers `k_{n,g}`. This file
checks the nine entries with `n ≤ 4` by `decide +kernel` (each a closed run of the recurrence on
`3ⁿ` columns), so that the table and its two consequences need only the six integers with
`n = 5, 6` as hypotheses (`BoundaryTotalsLarge`).

Kernel times (laptop, other jobs running): `n = 2` about 0.5 s, `n = 3` about 2 to 3 s, `n = 4`
about 10 to 12 s per entry; the whole file about 45 s and 6 GB of memory. The entries with
`n = 5` are left as hypotheses: one trial run of `k_{5,1}` needed about 10 GB and half an hour,
because the kernel keeps every intermediate packed row; `n = 6` would need several times that.
-/

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryThermal

theorem quotient_two_one : quotient 2 1 = some (paperK 2 1) := by decide +kernel
theorem quotient_two_neg_one : quotient 2 (-1) = some (paperK 2 (-1)) := by decide +kernel
theorem quotient_two_zero : quotient 2 0 = some (paperK 2 0) := by decide +kernel
theorem quotient_three_one : quotient 3 1 = some (paperK 3 1) := by decide +kernel
theorem quotient_three_neg_one : quotient 3 (-1) = some (paperK 3 (-1)) := by decide +kernel
theorem quotient_three_zero : quotient 3 0 = some (paperK 3 0) := by decide +kernel
theorem quotient_four_one : quotient 4 1 = some (paperK 4 1) := by decide +kernel
theorem quotient_four_neg_one : quotient 4 (-1) = some (paperK 4 (-1)) := by decide +kernel
theorem quotient_four_zero : quotient 4 0 = some (paperK 4 0) := by decide +kernel

/-- **The remaining hypothesis**: the six printed integers with `n = 5, 6`. -/
def BoundaryTotalsLarge : Prop :=
  ∀ n ∈ [5, 6], ∀ g ∈ fields, quotient n g = some (paperK n g)

/-- The nine kernel-checked entries and the six remaining ones give all fifteen. -/
theorem boundaryTotals_of_large (h : BoundaryTotalsLarge) : BoundaryTotals := by
  intro n hn g hg
  simp only [lengths, List.mem_cons, List.not_mem_nil, or_false] at hn
  have hg' := hg
  simp only [fields, List.mem_cons, List.not_mem_nil, or_false] at hg'
  rcases hn with rfl | rfl | rfl | rfl | rfl
  · rcases hg' with rfl | rfl | rfl
    exacts [quotient_two_one, quotient_two_neg_one, quotient_two_zero]
  · rcases hg' with rfl | rfl | rfl
    exacts [quotient_three_one, quotient_three_neg_one, quotient_three_zero]
  · rcases hg' with rfl | rfl | rfl
    exacts [quotient_four_one, quotient_four_neg_one, quotient_four_zero]
  · exact h 5 (by simp) g hg
  · exact h 6 (by simp) g hg

/-- **The table of the fifteen traces from the six integers with `n = 5, 6`.** -/
theorem boundaryTraceTable_of_large (h : BoundaryTotalsLarge) :
    BoundaryTraces.BoundaryTraceTable :=
  boundaryTraceTable_of_totals (boundaryTotals_of_large h)

/-- `Y₆(b₀) < .027` (the input `BoundaryInit.ThermalData.sixth`) from the six integers. -/
theorem sixth_of_large (h : BoundaryTotalsLarge) :
    BoundaryInit.shiftFactor (49 / 4) ^ 6 * oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) <
      27 / 1000 :=
  BoundaryTraces.sixth_of_table (boundaryTraceTable_of_large h)

/-- The leading-weight polynomial `< .003745` (eq. (appth:polynomial)) from the six integers. -/
theorem moment_polynomial_lt_of_large (h : BoundaryTotalsLarge) :
    BoundaryTraces.moment 6 - 28 / 100 * BoundaryTraces.moment 5 -
      3404 / 10000 * BoundaryTraces.moment 4 + 504 / 10000 * BoundaryTraces.moment 3 +
        324 / 10000 * BoundaryTraces.moment 2 < 3745 / 10 ^ 6 :=
  BoundaryTraces.moment_polynomial_lt_of_table (boundaryTraceTable_of_large h)

end BoundaryThermal

end CertifiedGaps.SpinOneHeisenberg
