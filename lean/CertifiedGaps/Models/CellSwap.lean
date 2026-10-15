/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.Basic.Real.Basic

/-!
# The swap of a two-site cell

For a cell of two sites with local space `n`, the swap `F (a, b) = (b, a)` is a real symmetric
involution with `F (A ⊗ B) F = B ⊗ A`. It is the internal reflection of a dimerized cell: a bond
from the second site of one cell to the first site of the next is `-∑ (1 ⊗ C) ⊗ F (1 ⊗ C) F`.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps

variable {n : Type} [Fintype n] [DecidableEq n]

/-- The swap of the two sites of a cell. -/
def cellSwap : Matrix (n × n) (n × n) ℝ := of fun x y ↦ if x = y.swap then 1 else 0

/-- `F (A ⊗ B) F = B ⊗ A`. -/
theorem cellSwap_mul_kronecker_mul_cellSwap (A B : Matrix n n ℝ) :
    cellSwap * (A ⊗ₖ B) * cellSwap = B ⊗ₖ A := by
  ext ⟨a, b⟩ ⟨c, e⟩
  simp only [mul_apply, cellSwap, of_apply, kroneckerMap_apply, Fintype.sum_prod_type,
    Prod.swap_prod_mk, Prod.mk.injEq, ite_mul, one_mul, zero_mul]
  simp [Finset.sum_ite_eq, Finset.sum_ite_eq', ite_and, mul_comm]

/-- The swap is an involution. -/
theorem cellSwap_mul_cellSwap : (cellSwap : Matrix (n × n) (n × n) ℝ) * cellSwap = 1 := by
  have h := cellSwap_mul_kronecker_mul_cellSwap (n := n) 1 1
  rwa [one_kronecker_one, Matrix.mul_one] at h

/-- The swap commutes with an exchange `c ∑ₖ Aₖ ⊗ Aₖ` within the cell. -/
theorem cellSwap_conj_smul_sum_kronecker_self {κ : Type*} [Fintype κ] (c : ℝ)
    (A : κ → Matrix n n ℝ) :
    cellSwap * (c • ∑ k, A k ⊗ₖ A k) * cellSwap = c • ∑ k, A k ⊗ₖ A k := by
  rw [Matrix.mul_smul, Matrix.smul_mul, Finset.mul_sum, Finset.sum_mul]
  simp_rw [cellSwap_mul_kronecker_mul_cellSwap]

omit [Fintype n] in
/-- The swap is symmetric. -/
theorem cellSwap_transpose : (cellSwap : Matrix (n × n) (n × n) ℝ)ᵀ = cellSwap := by
  ext ⟨a, b⟩ ⟨c, e⟩
  simp only [transpose_apply, cellSwap, of_apply, Prod.swap_prod_mk, Prod.mk.injEq]
  by_cases h₁ : a = e <;> by_cases h₂ : b = c
  · subst h₁ h₂; simp
  · simp [h₂, Ne.symm h₂]
  · simp [h₁, Ne.symm h₁]
  · simp [h₁, Ne.symm h₁]

end CertifiedGaps
