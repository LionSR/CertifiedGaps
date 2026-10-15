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
and one-site operators sitewise (`conj_bondOp_kronecker`, `conj_siteOp`), hence the ring
Hamiltonian bond by bond and site by site (`RPRing.conj_hamiltonian`), and unitary conjugation
preserves the spectrum with multiplicity (`eigenvalues_eq_of_conj`). So a gap proved for the
rotated Hamiltonian holds for the physical one (`RPRing.eigenvalues_eq_spectrum`).
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps

variable {d : Type} [Fintype d] [DecidableEq d] {L : ℕ}

/-- Conjugating a two-site product by a two-site product acts factorwise. -/
theorem kronecker_conj {m n : Type*} [Fintype m] [Fintype n] (X A : Matrix m m ℂ)
    (Y B : Matrix n n ℂ) :
    (X ⊗ₖ Y) * (A ⊗ₖ B) * (X ⊗ₖ Y)ᴴ = (X * A * Xᴴ) ⊗ₖ (Y * B * Yᴴ) := by
  rw [conjTranspose_kronecker, ← mul_kronecker_mul, ← mul_kronecker_mul]

omit [DecidableEq d] in
/-- Conjugating a tensor product by a tensor product acts sitewise. -/
theorem piTensor_conj (a X : Fin L → Matrix d d ℂ) :
    piTensor a * piTensor X * piTensor (fun i ↦ (a i)ᴴ) =
      piTensor fun i ↦ a i * X i * (a i)ᴴ := by
  rw [piTensor_mul, piTensor_mul]

/-- A product of one-site unitaries is unitary. -/
theorem piTensor_conjTranspose_mul_self (a : Fin L → Matrix d d ℂ)
    (ha : ∀ i, a i * (a i)ᴴ = 1) : piTensor (fun i ↦ (a i)ᴴ) * piTensor a = 1 := by
  rw [piTensor_mul]
  simp only [fun i ↦ mul_eq_one_comm.mp (ha i), piTensor_one]

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

namespace RPRing

variable {ι : Type} [Fintype ι] (R : RPRing d ι)

/-- Conjugating the ring Hamiltonian by a product of one-site unitaries `a` acts bond by bond
and site by site. (A ring of one site has a self-bond, which is excluded.) -/
theorem conj_hamiltonian (a : Fin L → Matrix d d ℂ) (ha : ∀ i, a i * (a i)ᴴ = 1) (hL : L ≠ 1) :
    piTensor a * R.hamiltonian L * piTensor (fun i ↦ (a i)ᴴ) =
      ∑ j, ∑ α, bondOp j (finRotate L j)
          (-((a j * R.Cc α * (a j)ᴴ) ⊗ₖ (a (finRotate L j) * R.Dc α * (a (finRotate L j))ᴴ))) +
        ∑ j, siteOp j (a j * R.oc * (a j)ᴴ) := by
  rw [hamiltonian, Matrix.mul_add, Matrix.add_mul, Finset.mul_sum, Finset.sum_mul,
    Finset.mul_sum, Finset.sum_mul]
  congr 1
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    have hL2 : 2 ≤ L := by have := j.isLt; omega
    rw [bond, bondMatrix, bondOp_neg, bondOp_sum, Matrix.mul_neg, Matrix.neg_mul,
      Finset.mul_sum, Finset.sum_mul, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun α _ ↦ ?_
    rw [conj_bondOp_kronecker a ha (ne_finRotate hL2 j), bondOp_neg]
  · exact Finset.sum_congr rfl fun j _ ↦ conj_siteOp a ha j _

/-- A Hermitian matrix that is the ring Hamiltonian conjugated by a product of one-site
unitaries has the spectrum of the ring, with multiplicity. -/
theorem eigenvalues_eq_spectrum (a : Fin L → Matrix d d ℂ) (ha : ∀ i, a i * (a i)ᴴ = 1)
    {H : Matrix (Fin L → d) (Fin L → d) ℂ} (hH : H.IsHermitian)
    (h : H = piTensor a * R.hamiltonian L * piTensor fun i ↦ (a i)ᴴ) :
    hH.eigenvalues = R.spectrum L := by
  refine eigenvalues_eq_of_conj _ _ (U := piTensor a) ?_ ?_
  · rw [conjTranspose_piTensor]; exact piTensor_conjTranspose_mul_self a ha
  · rw [conjTranspose_piTensor]; exact h

end RPRing

end CertifiedGaps
