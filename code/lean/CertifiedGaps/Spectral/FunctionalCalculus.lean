/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Matrix.HermitianFunctionalCalculus
import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.Analysis.SpecialFunctions.Exponential
import CertifiedGaps.Spectral.Basic

/-!
# Real functions of a Hermitian matrix

The function `f(A) = U diag(f(λᵢ)) Uᴴ` of a Hermitian matrix `A`, computed in its eigenbasis
`U = eigU hA`, is Mathlib's `Matrix.IsHermitian.cfc` (`hA.cfc f`), which `IsHermitian.cfc_eq`
identifies with the continuous functional calculus `cfc f A`. The generic `cfc` API asks for
continuity side conditions and a C⋆-algebra structure on matrices; the lemmas here work with
`hA.cfc` directly, for arbitrary `f : ℝ → ℝ`, and are all proved from `Uᴴ U = U Uᴴ = 1`.

## Main results

* `hermFun_eq`: `hA.cfc f = U diag(f(λᵢ)) Uᴴ` with `U = eigU hA`.
* The algebra of `f ↦ hA.cfc f`: `hermFun_id` (the spectral theorem), `hermFun_one`,
  `hermFun_add`, `hermFun_sub`, `hermFun_smul`, `hermFun_sum`, `hermFun_mul`, `hermFun_pow`;
  and, over `ℂ`, `exp_smul_hermFun`: `exp (c f(A)) = (e^{c f})(A)`.
* `nsq_conj_diagonal_mulVec_le`: `‖U diag(g) Uᴴ v‖² ≤ η² ‖v‖²` for unitary `U` and `|gᵢ| ≤ η`,
  and its case `nsq_hermFun_mulVec_le`: `‖f(A) v‖ ≤ max |f(λᵢ)| ‖v‖`.
-/

open Matrix NormedSpace

namespace CertifiedGaps.Spectral

variable {𝕜 : Type*} [RCLike 𝕜] {n : Type*} [Fintype n] [DecidableEq n]

/-- The reusable norm estimate: if `U` is unitary and `|gᵢ| ≤ η` for every `i`, then
`‖U diag(g) Uᴴ v‖² ≤ η² ‖v‖²`. Both `U` and `Uᴴ` preserve the norm, and the diagonal
scales each coordinate by at most `η`. -/
theorem nsq_conj_diagonal_mulVec_le {U : Matrix n n 𝕜} (hU : Uᴴ * U = 1) (hU' : U * Uᴴ = 1)
    (g : n → ℝ) {η : ℝ} (hg : ∀ i, |g i| ≤ η) (v : n → 𝕜) :
    nsq ((U * diagonal (fun i ↦ (g i : 𝕜)) * Uᴴ) *ᵥ v) ≤ η ^ 2 * nsq v := by
  rw [← mulVec_mulVec, ← mulVec_mulVec, nsq_mulVec_of_unitary hU,
    ← nsq_conjTranspose_mulVec hU' v, nsq_eq_sum, nsq_eq_sum, Finset.mul_sum]
  refine Finset.sum_le_sum fun i _ ↦ ?_
  rw [mulVec_diagonal, norm_mul, mul_pow, RCLike.norm_ofReal]
  exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (abs_nonneg _) (hg i) 2) (by positivity)

variable {A : Matrix n n 𝕜} (hA : A.IsHermitian)

/-- Mathlib's `hA.cfc f` is `U diag(f(λᵢ)) Uᴴ` in the eigenbasis `U = eigU hA`. -/
theorem hermFun_eq (f : ℝ → ℝ) :
    hA.cfc f = eigU hA * diagonal (fun i ↦ (f (hA.eigenvalues i) : 𝕜)) * (eigU hA)ᴴ := by
  rw [IsHermitian.cfc, Unitary.conjStarAlgAut_apply]
  rfl

/-- Functional-calculus norm bound: `‖f(A) v‖² ≤ η² ‖v‖²` whenever `|f| ≤ η` on the
eigenvalues of `A`. -/
theorem nsq_hermFun_mulVec_le (f : ℝ → ℝ) {η : ℝ} (hf : ∀ i, |f (hA.eigenvalues i)| ≤ η)
    (v : n → 𝕜) : nsq (hA.cfc f *ᵥ v) ≤ η ^ 2 * nsq v := by
  rw [hermFun_eq]
  exact nsq_conj_diagonal_mulVec_le (eigU_star_mul hA) (eigU_mul_star hA) _ hf v

/-- The spectral theorem: the identity function gives back `A`. -/
theorem hermFun_id : hA.cfc id = A := by
  have hspec : A = eigU hA * diagonal (fun i ↦ (hA.eigenvalues i : 𝕜)) * (eigU hA)ᴴ := by
    conv_lhs => rw [hA.spectral_theorem, Unitary.conjStarAlgAut_apply]
    rfl
  rw [hermFun_eq]
  exact hspec.symm

/-- The constant function `1` gives the identity matrix. -/
theorem hermFun_one : hA.cfc (fun _ ↦ 1) = 1 := by
  rw [hermFun_eq]
  simp only [RCLike.ofReal_one]
  rw [diagonal_one, Matrix.mul_one, eigU_mul_star]

/-- `f(A) - g(A) = (f - g)(A)`. -/
theorem hermFun_sub (f g : ℝ → ℝ) :
    hA.cfc f - hA.cfc g = hA.cfc (fun x ↦ f x - g x) := by
  simp only [hermFun_eq, RCLike.ofReal_sub, ← diagonal_sub, Matrix.mul_sub, Matrix.sub_mul]

/-- `f(A) + g(A) = (f + g)(A)`. -/
theorem hermFun_add (f g : ℝ → ℝ) :
    hA.cfc f + hA.cfc g = hA.cfc (fun x ↦ f x + g x) := by
  simp only [hermFun_eq, RCLike.ofReal_add, ← diagonal_add, Matrix.mul_add, Matrix.add_mul]

/-- `c f(A) = (c f)(A)` for real `c`. -/
theorem hermFun_smul (c : ℝ) (f : ℝ → ℝ) :
    (c : 𝕜) • hA.cfc f = hA.cfc (fun x ↦ c * f x) := by
  simp only [hermFun_eq, RCLike.ofReal_mul]
  rw [← smul_mul_assoc, ← mul_smul_comm, ← diagonal_smul]
  rfl

/-- Finite sums of functions of `A` are functions of `A`. -/
theorem hermFun_sum {ι : Type*} (s : Finset ι) (f : ι → ℝ → ℝ) :
    ∑ j ∈ s, hA.cfc (f j) = hA.cfc (fun x ↦ ∑ j ∈ s, f j x) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [hermFun_eq]
  | insert a s ha ih => simp only [Finset.sum_insert ha, ih, hermFun_add]

/-- `f(A) g(A) = (f g)(A)`, because `Uᴴ U = 1`. -/
theorem hermFun_mul (f g : ℝ → ℝ) :
    hA.cfc f * hA.cfc g = hA.cfc (fun x ↦ f x * g x) := by
  simp only [hermFun_eq, RCLike.ofReal_mul]
  rw [← diagonal_mul_diagonal]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc (eigU hA)ᴴ, eigU_star_mul, Matrix.one_mul]

/-- Powers of `A` are the power functions of its eigenvalues. -/
theorem hermFun_pow (j : ℕ) : A ^ j = hA.cfc (fun x ↦ x ^ j) := by
  induction j with
  | zero => simp only [pow_zero, hermFun_one]
  | succ j ih =>
    conv_lhs => rw [pow_succ, ih]
    conv_lhs => rhs; rw [← hermFun_id hA]
    rw [hermFun_mul]
    simp only [id, pow_succ]

omit [RCLike 𝕜] in
/-- The matrix exponential of a real multiple of a function of a Hermitian matrix is the
exponential of that function. -/
theorem exp_smul_hermFun {A : Matrix n n ℂ} (hA : A.IsHermitian) (c : ℝ) (f : ℝ → ℝ) :
    exp ((c : ℂ) • hA.cfc f) = hA.cfc (fun x ↦ Real.exp (c * f x)) := by
  set U := eigU hA
  have hU : Uᴴ * U = 1 := eigU_star_mul hA
  have hUinv : U⁻¹ = Uᴴ := inv_eq_left_inv hU
  have hUunit : IsUnit U := (isUnit_iff_isUnit_det U).mpr (isUnit_det_of_left_inverse hU)
  have h : (c : ℂ) • hA.cfc f =
      U * diagonal (fun i ↦ ((c * f (hA.eigenvalues i) : ℝ) : ℂ)) * U⁻¹ := by
    rw [hUinv]
    exact (hermFun_smul hA c f).trans (hermFun_eq hA _)
  have hd : exp (fun i ↦ ((c * f (hA.eigenvalues i) : ℝ) : ℂ)) =
      fun i ↦ ((Real.exp (c * f (hA.eigenvalues i)) : ℝ) : ℂ) := by
    funext i
    rw [Pi.coe_exp, ← Complex.exp_eq_exp_ℂ, Complex.ofReal_exp]
  rw [h, exp_conj _ _ hUunit, exp_diagonal, hUinv, hd, hermFun_eq]
  rfl

end CertifiedGaps.Spectral
