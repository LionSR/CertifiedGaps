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

The same conjugation changes the local basis of a ring `∑ⱼ Bⱼ,ⱼ₊₁`, including a ring whose
closing bond carries a twist `Bc`: a one-site unitary `V` on every site turns every bond `B` into
`(V ⊗ V)ᴴ B (V ⊗ V)` (`conj_twistedRingHamiltonian`, `conj_ringHamiltonian`) and leaves the
Gibbs trace unchanged (`gibbsTrace_twistedRingHamiltonian_conj`,
`gibbsTrace_ringHamiltonian_conj`). This lets a model state its reflection-positivity
factorization in one basis and compute its seeds in another.
-/

open Matrix NormedSpace
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

omit [Fintype d] in
/-- A matrix unit on a pair of sites is a tensor product of matrix units. -/
theorem single_eq_kronecker (p q : d × d) (x : ℂ) :
    single p q x = single p.1 q.1 x ⊗ₖ single p.2 q.2 1 := by
  rw [single_kronecker_single, mul_one]

/-- Sitewise conjugation of any two-site operator on distinct sites: the product `⨂ᵢ aᵢ` acts
on the bond `j, k` as `aⱼ ⊗ aₖ`. Every two-site operator is a sum of tensor products of matrix
units, which reduces this to `conj_bondOp_kronecker`. -/
theorem conj_bondOp (a : Fin L → Matrix d d ℂ) (ha : ∀ i, a i * (a i)ᴴ = 1) {j k : Fin L}
    (hjk : j ≠ k) (B : Matrix (d × d) (d × d) ℂ) :
    piTensor a * bondOp j k B * piTensor (fun i ↦ (a i)ᴴ) =
      bondOp j k ((a j ⊗ₖ a k) * B * (a j ⊗ₖ a k)ᴴ) := by
  conv_lhs => rw [matrix_eq_sum_single B]
  conv_rhs => rw [matrix_eq_sum_single B]
  simp only [bondOp_sum, Finset.mul_sum, Finset.sum_mul, single_eq_kronecker,
    conj_bondOp_kronecker a ha hjk, kronecker_conj]

/-- The product `V^{⊗L}` of copies of one unitary is unitary. -/
theorem piTensor_const_conjTranspose_mul_self (V : Matrix d d ℂ) (hV : Vᴴ * V = 1) :
    (piTensor fun _ : Fin L ↦ V)ᴴ * piTensor (fun _ : Fin L ↦ V) = 1 := by
  rw [conjTranspose_piTensor, piTensor_mul, hV, piTensor_one]

/-- `conj_bondOp` for the same unitary `V` on every site, in the change-of-basis form
`Uᴴ X U` (the columns of `V` are the new basis vectors). -/
theorem conj_bondOp_const (V : Matrix d d ℂ) (hV : Vᴴ * V = 1) {j k : Fin L} (hjk : j ≠ k)
    (B : Matrix (d × d) (d × d) ℂ) :
    (piTensor fun _ : Fin L ↦ V)ᴴ * bondOp j k B * piTensor (fun _ : Fin L ↦ V) =
      bondOp j k ((V ⊗ₖ V)ᴴ * B * (V ⊗ₖ V)) := by
  have h := conj_bondOp (fun _ : Fin L ↦ Vᴴ) (fun _ ↦ by rw [conjTranspose_conjTranspose, hV])
    hjk B
  simp only [conjTranspose_conjTranspose] at h
  rw [conjTranspose_piTensor, h, conjTranspose_kronecker, conjTranspose_kronecker,
    conjTranspose_conjTranspose]

/-! ### Changing the local basis of a ring -/

/-- **Sitewise change of basis of a ring Hamiltonian**, with or without a twisted closing bond:
the same unitary `V` on every site turns each bond `B` into `(V ⊗ V)ᴴ B (V ⊗ V)`. A ring of one
site is excluded because its bond joins the site to itself. -/
theorem conj_twistedRingHamiltonian (V : Matrix d d ℂ) (hV : Vᴴ * V = 1) (hL : L ≠ 1)
    {B Bc B' Bc' : Matrix (d × d) (d × d) ℂ} (hB : (V ⊗ₖ V)ᴴ * B * (V ⊗ₖ V) = B')
    (hBc : (V ⊗ₖ V)ᴴ * Bc * (V ⊗ₖ V) = Bc') :
    (piTensor fun _ : Fin L ↦ V)ᴴ * twistedRingHamiltonian B Bc L *
        piTensor (fun _ : Fin L ↦ V) = twistedRingHamiltonian B' Bc' L := by
  rw [twistedRingHamiltonian, twistedRingHamiltonian, Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  have hL2 : 2 ≤ L := by have := j.isLt; omega
  rw [conj_bondOp_const V hV (ne_finRotate hL2 j)]
  split_ifs
  · rw [hBc]
  · rw [hB]

/-- Sitewise change of basis of a translation-invariant ring: the bond `B` becomes
`(V ⊗ V)ᴴ B (V ⊗ V)`. -/
theorem conj_ringHamiltonian (V : Matrix d d ℂ) (hV : Vᴴ * V = 1) (hL : L ≠ 1)
    (B : Matrix (d × d) (d × d) ℂ) :
    (piTensor fun _ : Fin L ↦ V)ᴴ * ringHamiltonian B L * piTensor (fun _ : Fin L ↦ V) =
      ringHamiltonian ((V ⊗ₖ V)ᴴ * B * (V ⊗ₖ V)) L := by
  simp only [← twistedRingHamiltonian_self]
  exact conj_twistedRingHamiltonian V hV hL rfl rfl

/-- A sitewise change of basis leaves the Gibbs trace of a (twisted) ring unchanged: rings whose
bonds are conjugate by `V ⊗ V` have the same partition function. -/
theorem gibbsTrace_twistedRingHamiltonian_conj (V : Matrix d d ℂ) (hV : Vᴴ * V = 1)
    (hL : L ≠ 1) {B Bc B' Bc' : Matrix (d × d) (d × d) ℂ}
    (hB : (V ⊗ₖ V)ᴴ * B * (V ⊗ₖ V) = B') (hBc : (V ⊗ₖ V)ᴴ * Bc * (V ⊗ₖ V) = Bc') (β : ℝ) :
    gibbsTrace (twistedRingHamiltonian B' Bc' L) β =
      gibbsTrace (twistedRingHamiltonian B Bc L) β := by
  rw [← conj_twistedRingHamiltonian V hV hL hB hBc,
    gibbsTrace_conj _ _ (piTensor_const_conjTranspose_mul_self V hV)]

/-- A sitewise change of basis of the bond leaves the Gibbs trace of a translation-invariant
ring unchanged. -/
theorem gibbsTrace_ringHamiltonian_conj (V : Matrix d d ℂ) (hV : Vᴴ * V = 1) (hL : L ≠ 1)
    {B B' : Matrix (d × d) (d × d) ℂ} (hB : (V ⊗ₖ V)ᴴ * B * (V ⊗ₖ V) = B') (β : ℝ) :
    gibbsTrace (ringHamiltonian B' L) β = gibbsTrace (ringHamiltonian B L) β := by
  simp only [← twistedRingHamiltonian_self]
  exact gibbsTrace_twistedRingHamiltonian_conj V hV hL hB hB β

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
