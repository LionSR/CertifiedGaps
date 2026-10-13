/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Model

/-!
# Sitewise unitary conjugation

Several models are reflection positive only in a rotated frame, for example after
`e^{iπSʸ}` on every other site. A product of one-site unitaries `U = ⨂ᵢ aᵢ` conjugates two-site
and one-site operators sitewise (`conj_bondOp_kronecker`, `conj_siteOp`), and unitary
conjugation preserves the spectrum with multiplicity (`eigenvalues_eq_of_conj`). So a gap
proved for the rotated Hamiltonian holds for the physical one.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps

variable {d : Type} [Fintype d] [DecidableEq d] {L : ℕ}

omit [DecidableEq d] in
/-- Conjugating a tensor product by a tensor product acts sitewise. -/
theorem piTensor_conj (a X : Fin L → Matrix d d ℂ) :
    piTensor a * piTensor X * piTensor (fun i ↦ (a i)ᴴ) =
      piTensor fun i ↦ a i * X i * (a i)ᴴ := by
  rw [piTensor_mul, piTensor_mul]

/-- A product of one-site unitaries is unitary. -/
theorem piTensor_conjTranspose_mul_self (a : Fin L → Matrix d d ℂ)
    (ha : ∀ i, (a i)ᴴ * a i = 1) : piTensor (fun i ↦ (a i)ᴴ) * piTensor a = 1 := by
  rw [piTensor_mul]
  simp only [ha, piTensor_one]

/-- Sitewise conjugation of a two-site product operator. -/
theorem conj_bondOp_kronecker (a : Fin L → Matrix d d ℂ) (ha : ∀ i, a i * (a i)ᴴ = 1)
    {j k : Fin L} (hjk : j ≠ k) (A B : Matrix d d ℂ) :
    piTensor a * bondOp j k (A ⊗ₖ B) * piTensor (fun i ↦ (a i)ᴴ) =
      bondOp j k ((a j * A * (a j)ᴴ) ⊗ₖ (a k * B * (a k)ᴴ)) := by
  rw [bondOp_kronecker hjk, bondOp_kronecker hjk, piTensor_conj]
  congr 1
  ext1 i
  by_cases h₁ : i = j
  · subst h₁; simp
  · by_cases h₂ : i = k
    · subst h₂; simp [h₁]
    · simp [h₁, h₂, ha i]

/-- Sitewise conjugation of a one-site operator. -/
theorem conj_siteOp (a : Fin L → Matrix d d ℂ) (ha : ∀ i, a i * (a i)ᴴ = 1) (j : Fin L)
    (A : Matrix d d ℂ) :
    piTensor a * siteOp j A * piTensor (fun i ↦ (a i)ᴴ) = siteOp j (a j * A * (a j)ᴴ) := by
  rw [siteOp, siteOp, piTensor_conj]
  congr 1
  ext1 i
  by_cases h : i = j
  · subst h; simp
  · simp [h, ha i]

/-- Unitary conjugation preserves the eigenvalues with multiplicity. -/
theorem eigenvalues_eq_of_conj {n : Type*} [Fintype n] [DecidableEq n] {A B U : Matrix n n ℂ}
    (hA : A.IsHermitian) (hB : B.IsHermitian) (hU : Uᴴ * U = 1) (h : B = U * A * Uᴴ) :
    hB.eigenvalues = hA.eigenvalues := by
  rw [IsHermitian.eigenvalues_eq_eigenvalues_iff, h, Matrix.mul_assoc, charpoly_mul_comm,
    Matrix.mul_assoc, hU, Matrix.mul_one]

end CertifiedGaps
