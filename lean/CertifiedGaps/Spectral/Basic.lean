/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas

/-!
# Quadratic forms, eigenbases and eigenvalue counts

Shared vocabulary for the eigenvalue-counting layer (`CertifiedGaps.Spectral.*`). Everything is
stated for matrices over a field `𝕜` with `[RCLike 𝕜]`, so it applies verbatim to real symmetric
(`𝕜 = ℝ`) and complex Hermitian (`𝕜 = ℂ`) matrices on a finite index type.

* `qform A x = Re ⟪x, A x⟫` and `nsq x = ‖x‖² = Re ⟪x, x⟫`, both real and computed with the
  plain dot product of `n → 𝕜` (no `EuclideanSpace` coercions);
* `eigU hA`: the unitary of eigenvectors, with `qform_eigU`: in the eigenbasis the form is
  `∑ λᵢ ‖yᵢ‖²`;
* `countBelow s a t = #{i ∈ s : a i < t}` and `eigCount hA t = #{i : λᵢ < t}`, the eigenvalue
  counting function. Sorted-eigenvalue statements (Weyl monotonicity, Courant–Fischer) are
  phrased through counting functions, which avoids choosing a sorting of `hA.eigenvalues`.
-/

open Matrix

namespace CertifiedGaps.Spectral

variable {𝕜 : Type*} [RCLike 𝕜] {n : Type*} [Fintype n]

/-! ### Quadratic form and squared norm -/

/-- The real quadratic form `Re ⟪x, A x⟫` of a square matrix. -/
noncomputable def qform (A : Matrix n n 𝕜) (x : n → 𝕜) : ℝ := RCLike.re (star x ⬝ᵥ (A *ᵥ x))

/-- The squared Euclidean norm `Re ⟪x, x⟫ = ∑ᵢ ‖xᵢ‖²` of a vector. -/
noncomputable def nsq (x : n → 𝕜) : ℝ := RCLike.re (star x ⬝ᵥ x)

theorem nsq_eq_sum (x : n → 𝕜) : nsq x = ∑ i, ‖x i‖ ^ 2 := by
  simp only [nsq, dotProduct, Pi.star_apply, map_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [RCLike.star_def, RCLike.conj_mul, ← RCLike.ofReal_pow, RCLike.ofReal_re]

theorem nsq_nonneg (x : n → 𝕜) : 0 ≤ nsq x := by
  rw [nsq_eq_sum]; positivity

theorem nsq_eq_zero_iff {x : n → 𝕜} : nsq x = 0 ↔ x = 0 := by
  rw [nsq_eq_sum, Finset.sum_eq_zero_iff_of_nonneg (fun i _ ↦ by positivity)]
  simp [funext_iff]

theorem nsq_pos {x : n → 𝕜} (hx : x ≠ 0) : 0 < nsq x :=
  (nsq_nonneg x).lt_of_ne (Ne.symm (mt nsq_eq_zero_iff.1 hx))

@[simp] theorem nsq_zero : nsq (0 : n → 𝕜) = 0 := nsq_eq_zero_iff.2 rfl

theorem qform_one [DecidableEq n] (x : n → 𝕜) : qform (1 : Matrix n n 𝕜) x = nsq x := by
  simp [qform, nsq]

theorem qform_add_left (A B : Matrix n n 𝕜) (x : n → 𝕜) :
    qform (A + B) x = qform A x + qform B x := by
  simp [qform, add_mulVec, dotProduct_add]

theorem qform_sub_left (A B : Matrix n n 𝕜) (x : n → 𝕜) :
    qform (A - B) x = qform A x - qform B x := by
  simp [qform, sub_mulVec, dotProduct_sub]

theorem qform_smul_left (r : ℝ) (A : Matrix n n 𝕜) (x : n → 𝕜) :
    qform ((r : 𝕜) • A) x = r * qform A x := by
  simp [qform, smul_mulVec, dotProduct_smul, RCLike.smul_re]

/-- The form of `A - t • 1` is the form of `A` shifted by `t ‖x‖²`. -/
theorem qform_sub_smul_one [DecidableEq n] (A : Matrix n n 𝕜) (t : ℝ) (x : n → 𝕜) :
    qform (A - (t : 𝕜) • 1) x = qform A x - t * nsq x := by
  rw [qform_sub_left, qform_smul_left, qform_one]

theorem qform_zero_right (A : Matrix n n 𝕜) : qform A 0 = 0 := by simp [qform]

/-- `⟪M y, M z⟫ = ⟪y, Mᴴ M z⟫`. -/
theorem star_mulVec_dotProduct_mulVec {m : Type*} [Fintype m] (M : Matrix n m 𝕜)
    (y z : m → 𝕜) : star (M *ᵥ y) ⬝ᵥ (M *ᵥ z) = star y ⬝ᵥ ((Mᴴ * M) *ᵥ z) := by
  rw [star_mulVec, ← dotProduct_mulVec, mulVec_mulVec]

/-- Compression: the form of `A` at `M y` is the form of `Mᴴ A M` at `y`. -/
theorem qform_mulVec {m : Type*} [Fintype m] (A : Matrix n n 𝕜) (M : Matrix n m 𝕜)
    (y : m → 𝕜) : qform A (M *ᵥ y) = qform (Mᴴ * A * M) y := by
  simp only [qform]
  rw [star_mulVec, ← dotProduct_mulVec, mulVec_mulVec, mulVec_mulVec, Matrix.mul_assoc]

/-- The squared norm of `M y` is the form of the Gram matrix `Mᴴ M` at `y`. -/
theorem nsq_mulVec {m : Type*} [Fintype m] (M : Matrix n m 𝕜) (y : m → 𝕜) :
    nsq (M *ᵥ y) = qform (Mᴴ * M) y := by
  simp only [nsq, qform, star_mulVec_dotProduct_mulVec]

/-- For a Hermitian matrix the cross term of the form is symmetric. -/
theorem qform_add (A : Matrix n n 𝕜) (hA : A.IsHermitian) (x y : n → 𝕜) :
    qform A (x + y) = qform A x + 2 * RCLike.re (star y ⬝ᵥ (A *ᵥ x)) + qform A y := by
  have hsym : RCLike.re (star x ⬝ᵥ (A *ᵥ y)) = RCLike.re (star y ⬝ᵥ (A *ᵥ x)) := by
    rw [star_dotProduct, star_mulVec, hA.eq, ← dotProduct_mulVec, RCLike.star_def,
      RCLike.conj_re]
  simp only [qform, mulVec_add, dotProduct_add, add_dotProduct, star_add, map_add]
  rw [hsym]; ring

/-! ### The eigenbasis -/

variable [DecidableEq n] {A : Matrix n n 𝕜} (hA : A.IsHermitian)

/-- The unitary matrix whose columns are the eigenvectors `hA.eigenvectorBasis`. -/
noncomputable abbrev eigU : Matrix n n 𝕜 := (hA.eigenvectorUnitary : Matrix n n 𝕜)

theorem eigU_star_mul : (eigU hA)ᴴ * eigU hA = 1 := Unitary.coe_star_mul_self _

theorem eigU_mul_star : eigU hA * (eigU hA)ᴴ = 1 := Unitary.coe_mul_star_self _

/-- `Uᴴ A U` is the diagonal matrix of eigenvalues. -/
theorem eigU_conj :
    (eigU hA)ᴴ * A * eigU hA = diagonal (fun i ↦ (hA.eigenvalues i : 𝕜)) := by
  have hspec : A = eigU hA * diagonal (fun i ↦ (hA.eigenvalues i : 𝕜)) * (eigU hA)ᴴ := by
    conv_lhs => rw [hA.spectral_theorem, Unitary.conjStarAlgAut_apply]
    rfl
  generalize hU : eigU hA = U at hspec
  have h1 : Uᴴ * U = 1 := hU ▸ eigU_star_mul hA
  conv_lhs => rw [hspec]
  simp only [Matrix.mul_assoc]
  rw [h1, Matrix.mul_one, ← Matrix.mul_assoc, h1, Matrix.one_mul]

/-- A unitary preserves the squared norm. -/
theorem nsq_eigU (y : n → 𝕜) : nsq (eigU hA *ᵥ y) = nsq y := by
  rw [nsq_mulVec, eigU_star_mul, qform_one]

/-- In the eigenbasis the quadratic form is `∑ λᵢ ‖yᵢ‖²`. -/
theorem qform_eigU (y : n → 𝕜) :
    qform A (eigU hA *ᵥ y) = ∑ i, hA.eigenvalues i * ‖y i‖ ^ 2 := by
  rw [qform_mulVec, eigU_conj]
  simp only [qform, mulVec_diagonal, dotProduct, Pi.star_apply, map_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [mul_left_comm, RCLike.star_def, RCLike.conj_mul, ← RCLike.ofReal_pow,
    ← RCLike.ofReal_mul, RCLike.ofReal_re]

/-- Every vector is `U y` for `y = Uᴴ x`. -/
theorem eigU_mulVec_star (x : n → 𝕜) : eigU hA *ᵥ ((eigU hA)ᴴ *ᵥ x) = x := by
  rw [mulVec_mulVec, eigU_mul_star, one_mulVec]

/-! ### Counting functions -/

/-- The number `#{i ∈ s : a i < t}` of values strictly below `t`. -/
noncomputable def countBelow {ι : Type*} (s : Finset ι) (a : ι → ℝ) (t : ℝ) : ℕ :=
  (s.filter fun i ↦ a i < t).card

/-- The eigenvalue counting function `#{i : λᵢ(A) < t}` (with multiplicity). -/
noncomputable def eigCount (t : ℝ) : ℕ := countBelow Finset.univ hA.eigenvalues t

theorem eigCount_le_card (t : ℝ) : eigCount hA t ≤ Fintype.card n :=
  Finset.card_filter_le _ _

end CertifiedGaps.Spectral
