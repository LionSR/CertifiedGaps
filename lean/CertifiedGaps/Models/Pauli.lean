/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.FinCases

/-!
# The real Pauli matrices `X` and `Z`

The Ising models need only `X` and `Z`, which are real and symmetric, so their bond factors and
on-site terms are real without any basis change.
-/

open Matrix

namespace CertifiedGaps

/-- The Pauli matrix `Z`. -/
def pauliZ : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, -1]

/-- The Pauli matrix `X`. -/
def pauliX : Matrix (Fin 2) (Fin 2) ℝ := !![0, 1; 1, 0]

/-- `Z` is symmetric. -/
theorem pauliZ_transpose : pauliZᵀ = pauliZ := by
  ext i j; fin_cases i <;> fin_cases j <;> rfl

/-- `X` is symmetric. -/
theorem pauliX_transpose : pauliXᵀ = pauliX := by
  ext i j; fin_cases i <;> fin_cases j <;> rfl

end CertifiedGaps
