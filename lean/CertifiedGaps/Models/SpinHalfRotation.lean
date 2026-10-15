/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.SpinHalf
import Mathlib.Analysis.Complex.Basic
import Mathlib.LinearAlgebra.Matrix.Hermitian

/-!
# The rotation `e^{iπSʸ}` of spin one-half

The complex spin matrices `S_α` (`spinHalf`), the real factors `P_α = Sˣ, iSʸ, Sᶻ` as complex
matrices (`Pc α = c_α S_α` with `c = (1, i, 1)`), and the rotation `rot = e^{iπSʸ}`, which maps
`P_α ↦ s_α P_α` with `s = (-1, 1, -1)`. Since `-s_α c_α² = 1`, a rotated bond
`-(s_α Pc α) ⊗ Pc α` is the exchange term `S_α ⊗ S_α`.
-/

open Matrix

namespace CertifiedGaps.SpinHalf

/-- The spin one-half matrices `Sˣ, Sʸ, Sᶻ`. -/
noncomputable def spinHalf : Fin 3 → Matrix (Fin 2) (Fin 2) ℂ :=
  ![!![0, 1 / 2; 1 / 2, 0], !![0, -Complex.I / 2; Complex.I / 2, 0], !![1 / 2, 0; 0, -1 / 2]]

/-- The spin matrices are Hermitian. -/
theorem spinHalf_isHermitian (α : Fin 3) : (spinHalf α).IsHermitian := by
  fin_cases α <;> ext i j <;> fin_cases i <;> fin_cases j <;>
    simp [spinHalf, conjTranspose_apply]

/-- The real factors `Sˣ, iSʸ, Sᶻ` as complex matrices. -/
noncomputable def Pc (α : Fin 3) : Matrix (Fin 2) (Fin 2) ℂ := (P α).map (↑)

/-- The phases `c = (1, i, 1)` with `Pc α = c_α S_α`. -/
noncomputable def c : Fin 3 → ℂ := ![1, Complex.I, 1]

/-- The sign of `P_α` under the rotation. -/
def s : Fin 3 → ℂ := ![-1, 1, -1]

/-- `Pc α = c_α S_α`. -/
theorem Pc_eq (α : Fin 3) : Pc α = c α • spinHalf α := by
  fin_cases α <;> ext i j <;> fin_cases i <;> fin_cases j <;>
    simp [Pc, P, c, spinHalf] <;> ring_nf <;> simp [Complex.ext_iff] <;> norm_num

/-- `-s_α c_α² = 1`: a rotated bond `-(s_α Pc α) ⊗ Pc α` is `S_α ⊗ S_α`. -/
theorem neg_s_mul_c_sq (α : Fin 3) : -(s α * c α ^ 2) = 1 := by
  fin_cases α <;> simp [s, c]

/-- The rotation `e^{iπSʸ}`. -/
def rot : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; -1, 0]

/-- Conjugation by the rotation in closed form. -/
theorem rot_conj (A : Matrix (Fin 2) (Fin 2) ℂ) :
    rot * A * rotᴴ = !![A 1 1, -A 1 0; -A 0 1, A 0 0] := by
  ext a b; fin_cases a <;> fin_cases b <;>
    simp [rot, mul_apply, conjTranspose_apply, Fin.sum_univ_two, vecMul, dotProduct]

/-- The rotation is unitary. -/
theorem rot_mul_conjTranspose : rot * rotᴴ = 1 := by
  ext a b; fin_cases a <;> fin_cases b <;> simp [rot, mul_apply, conjTranspose_apply]

/-- The rotation maps `P_α ↦ s_α P_α`. -/
theorem rot_conj_Pc (α : Fin 3) : rot * Pc α * rotᴴ = s α • Pc α := by
  rw [rot_conj]
  fin_cases α <;> ext i j <;> fin_cases i <;> fin_cases j <;> simp [Pc, P, s] <;> norm_num

/-- The rotation maps `S_α ↦ s_α S_α`. -/
theorem rot_conj_spinHalf (α : Fin 3) : rot * spinHalf α * rotᴴ = s α • spinHalf α := by
  have h := rot_conj_Pc α
  rw [Pc_eq, Matrix.mul_smul, Matrix.smul_mul, smul_comm] at h
  exact smul_right_injective _ (by fin_cases α <;> simp [c]) h

end CertifiedGaps.SpinHalf
