/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum

/-!
# Real spin one-half factors

In the basis `(↑, ↓)` the matrices `Sˣ`, `iSʸ` and `Sᶻ` are real: `Sˣ` and `Sᶻ` are symmetric
and `iSʸ` is antisymmetric. After a sublattice rotation `e^{iπSʸ}` on every other site, a
Heisenberg bond `S · S` becomes `-SˣSˣ + SʸSʸ - SᶻSᶻ`, that is
`-(Sˣ ⊗ Sˣ + iSʸ ⊗ iSʸ + Sᶻ ⊗ Sᶻ)`: minus a sum of squares of real tensor products.
-/

open Matrix

namespace CertifiedGaps.SpinHalf

/-- The real one-site factors `Sˣ, iSʸ, Sᶻ` of spin one-half in the basis `(↑, ↓)`. -/
noncomputable def P : Fin 3 → Matrix (Fin 2) (Fin 2) ℝ :=
  ![!![0, 1 / 2; 1 / 2, 0], !![0, 1 / 2; -1 / 2, 0], !![1 / 2, 0; 0, -1 / 2]]

/-- `Sˣ` and `Sᶻ` are symmetric, `iSʸ` is antisymmetric. -/
theorem P_symm (k : Fin 3) : (P k)ᵀ = P k ∨ (P k)ᵀ = -P k := by
  fin_cases k
  · left; ext i j; fin_cases i <;> fin_cases j <;> norm_num [P]
  · right; ext i j; fin_cases i <;> fin_cases j <;> norm_num [P]
  · left; ext i j; fin_cases i <;> fin_cases j <;> norm_num [P]

end CertifiedGaps.SpinHalf
