/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialCompute

/-!
# Kernel checks of the intertwiner identities: data, `T W = W B`, `V W = c 1`

Each theorem here is one `decide +kernel` on the generated data:

* well-formedness of the sparse codes of `W` and `V` (`WIdxOk`, `VIdxOk`) and their entry
  bounds (`WAbsOk`, `VAbsOk`), the recipe bound `|Aₛ(i, j)| ≤ 640` (`ABoundOk`), and the offsets
  of the blocks;
* the packed comparison for `T W = W B` on the radix vector of width `42` over all columns
  (`checkT`) and the comparison for `V W = c 1` on the radix vector of width `39` (`checkV`).

The comparisons for `O W = W C` (width `71`) are in `TrialChecksInsertion1`–
`TrialChecksInsertion4`, one file per group of columns, to bound the memory of each kernel
run. The widths exceed the a priori entry bounds used in `T_mul_Wm`, `O_mul_Wm`, `Vm_mul_Wm`.
-/

namespace SeedKernel.Trial

open CertifiedGaps.SpinOneHeisenberg

set_option maxRecDepth 100000

/-- The Lean recipe gives the matrices written by the Julia generator (which asserts them equal
to the upstream `matrix-data.json`). -/
theorem trialRows_eq_code :
    trialRows 0 = shiftL trialMatrixCode0.1 (decodeSquare 26 trialMatrixCode0) ∧
      trialRows 1 = shiftL trialMatrixCode1.1 (decodeSquare 26 trialMatrixCode1) ∧
      trialRows 2 = shiftL trialMatrixCode2.1 (decodeSquare 26 trialMatrixCode2) := by
  decide +kernel

theorem wIdxOk : WIdxOk := by decide +kernel
theorem vIdxOk : VIdxOk := by decide +kernel
theorem wAbsOk : WAbsOk := by decide +kernel
theorem vAbsOk : VAbsOk := by decide +kernel
theorem aBoundOk : ABoundOk := by decide +kernel

theorem tB_offset : ∀ b : Fin 11, (tB b).1 ≤ 2 ^ 19 := by decide +kernel
theorem oB_offset : ∀ b : Fin 11, (oB b).1 ≤ 2 ^ 37 := by decide +kernel

/-- `T W = W B` on the radix vector of width `42`, packed rows (slots of `676 · 42` bits). -/
theorem checkT : TPack 28392 (YG 42 0 11) = ZPack 28392 (ZG 42 0 11 tB) := by decide +kernel

/-- `V W = c 1` on the radix vector of width `39`. -/
theorem checkV : VCheck 39 := by decide +kernel

end SeedKernel.Trial
