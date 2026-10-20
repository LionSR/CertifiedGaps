/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Hamiltonian

/-!
# The variational lower bound on a Gibbs trace

A trial vector with Rayleigh quotient below `E` bounds the Gibbs trace of any Hermitian matrix
from below by `e^{-βE}`. In the eigenbasis of `H` the Rayleigh quotient is an average of
eigenvalues, so it bounds the lowest one; the Gibbs trace is a sum of positive terms
`e^{-β Eₘ}`, so it exceeds the term of that eigenvalue.

## Main results

* `exists_eigenvalue_lt_of_re_dotProduct_lt`: the variational principle.
* `exp_lt_gibbsTrace_of_eigenvalue_lt`: one low eigenvalue bounds the Gibbs trace.
* `exp_lt_gibbsTrace_of_re_dotProduct_lt`: the two combined.
-/

open Matrix NormedSpace

namespace CertifiedGaps

variable {n : Type*} [Fintype n] [DecidableEq n] {H : Matrix n n ℂ}

/-- **Variational principle.** If the Rayleigh quotient of a vector `ψ` lies below `E`, some
eigenvalue of the Hermitian matrix `H` lies below `E`. -/
theorem exists_eigenvalue_lt_of_re_dotProduct_lt (hH : H.IsHermitian) (ψ : n → ℂ) {E : ℝ}
    (h : (star ψ ⬝ᵥ (H *ᵥ ψ)).re < E * (star ψ ⬝ᵥ ψ).re) : ∃ i, hH.eigenvalues i < E := by
  set U : Matrix n n ℂ := (hH.eigenvectorUnitary : Matrix n n ℂ)
  have hU2 : U * Uᴴ = 1 := Unitary.coe_mul_star_self _
  have hspec : H = U * diagonal (fun i ↦ (hH.eigenvalues i : ℂ)) * Uᴴ := by
    conv_lhs => rw [hH.spectral_theorem, Unitary.conjStarAlgAut_apply]
    rfl
  -- In the eigenbasis `y = Uᴴ ψ` both quadratic forms become diagonal.
  set y := Uᴴ *ᵥ ψ
  have hquad : ∀ M : Matrix n n ℂ, star ψ ⬝ᵥ ((U * M * Uᴴ) *ᵥ ψ) = star y ⬝ᵥ (M *ᵥ y) := by
    intro M
    rw [← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec, star_mulVec,
      conjTranspose_conjTranspose]
  have hnorm := hquad 1
  rw [Matrix.mul_one, hU2, one_mulVec, one_mulVec] at hnorm
  rw [hspec, hquad, hnorm] at h
  simp only [mulVec_diagonal, dotProduct, Pi.star_apply, Complex.re_sum, Finset.mul_sum] at h
  by_contra hcon
  push Not at hcon
  refine absurd h (not_lt.2 (Finset.sum_le_sum fun i _ ↦ ?_))
  have hre : (star (y i) * ((hH.eigenvalues i : ℂ) * y i)).re =
      hH.eigenvalues i * (star (y i) * y i).re := by
    rw [mul_left_comm, Complex.re_ofReal_mul]
  have hpos : 0 ≤ (star (y i) * y i).re := by
    rw [Complex.star_def, ← Complex.normSq_eq_conj_mul_self, Complex.ofReal_re]
    exact Complex.normSq_nonneg _
  rw [hre]
  exact mul_le_mul_of_nonneg_right (hcon i) hpos

/-- An eigenvalue below `E` bounds the Gibbs trace from below: `Re Tr e^{-βH} > e^{-βE}` for
`β > 0`, since every term of the eigenvalue sum is positive. -/
theorem exp_lt_gibbsTrace_of_eigenvalue_lt (hH : H.IsHermitian) {E β : ℝ} (hβ : 0 < β)
    {i : n} (hi : hH.eigenvalues i < E) :
    Real.exp (-(β * E)) < gibbsTrace H β := by
  rw [gibbsTrace_eq_sum hH]
  calc Real.exp (-(β * E)) < Real.exp (-(β * hH.eigenvalues i)) :=
        Real.exp_lt_exp.2 (by nlinarith)
    _ ≤ ∑ k, Real.exp (-(β * hH.eigenvalues k)) :=
        Finset.single_le_sum (f := fun k ↦ Real.exp (-(β * hH.eigenvalues k)))
          (fun k _ ↦ (Real.exp_pos _).le) (Finset.mem_univ i)

/-- **Variational bound on the Gibbs trace.** A trial vector `ψ` with Rayleigh quotient below
`E` gives `Re Tr e^{-βH} > e^{-βE}` for every `β > 0`. -/
theorem exp_lt_gibbsTrace_of_re_dotProduct_lt (hH : H.IsHermitian) (ψ : n → ℂ) {E β : ℝ}
    (h : (star ψ ⬝ᵥ (H *ᵥ ψ)).re < E * (star ψ ⬝ᵥ ψ).re) (hβ : 0 < β) :
    Real.exp (-(β * E)) < gibbsTrace H β :=
  let ⟨_, hi⟩ := exists_eigenvalue_lt_of_re_dotProduct_lt hH ψ h
  exp_lt_gibbsTrace_of_eigenvalue_lt hH hβ hi

end CertifiedGaps
