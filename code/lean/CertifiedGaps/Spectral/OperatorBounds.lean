/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Spectral.Basic

/-!
# Operator inequalities through quadratic forms

Two-sided bounds `a ≤ A ≤ b` on a matrix, stated through the quadratic form of
`Spectral/Basic.lean`: `a ≤ A` means `∀ x, a * nsq x ≤ qform A x`, and `|A| ≤ c` means
`∀ x, |qform A x| ≤ c * nsq x`. Everything holds over any `RCLike` field `𝕜`.

## Main results

* Sums of bounds (`le_qform_sum`, `qform_sum_le`) and invariance of bounds under unitary
  conjugation `M A Mᴴ` (`le_qform_conj`, `qform_conj_le`, `abs_qform_conj_le`).
* *Weighted Schur test* (`abs_qform_le_of_schur`). If `‖Aᵢⱼ‖ ≤ Kᵢⱼ` with `K` symmetric and
  `∑ⱼ Kᵢⱼ qⱼ ≤ c qᵢ` for positive weights `q`, then `|⟪x, A x⟫| ≤ c ‖x‖²`. This is the bound
  `‖D⁻¹ A D‖_∞ ≤ c`, proved directly on the quadratic form by
  `2 |xᵢ| |xⱼ| ≤ (qⱼ/qᵢ) |xᵢ|² + (qᵢ/qⱼ) |xⱼ|²` (`two_mul_le_weighted`), so no similarity
  transform is needed.
* Bounds from eigenvalues (`qform_le_of_eigenvalues_le`, `le_qform_of_le_eigenvalues`) and back
  (`abs_eigenvalues_le_one_of_qform`: a contraction has its spectrum in `[-1, 1]`).
* *Annihilating polynomial* (`eigenvalues_mem_of_prod_eq_zero`). If a Hermitian `A` satisfies
  `∏ᵣ (A - r) = 0` over a list of real roots, every eigenvalue of `A` is one of the roots, so `A`
  lies between the least and the largest root (`qform_le_of_prod_eq_zero`,
  `le_qform_of_prod_eq_zero`). A model checks the polynomial identity by finite computation.
* `exists_eq_conjTranspose_mul_self`: a Hermitian matrix with a nonnegative form is a Gram
  matrix `Cᴴ C`.
-/

open Matrix

namespace CertifiedGaps.Spectral

variable {𝕜 : Type*} [RCLike 𝕜] {n : Type*} [Fintype n]

/-! ### Sums and unitary conjugates -/

/-- Lower bounds add up: `aᵢ ≤ Aᵢ` for all `i ∈ s` gives `∑ aᵢ ≤ ∑ Aᵢ`. -/
theorem le_qform_sum {ι : Type*} (s : Finset ι) {A : ι → Matrix n n 𝕜} {a : ι → ℝ}
    (h : ∀ i ∈ s, ∀ x, a i * nsq x ≤ qform (A i) x) (x : n → 𝕜) :
    (∑ i ∈ s, a i) * nsq x ≤ qform (∑ i ∈ s, A i) x := by
  rw [qform_sum, Finset.sum_mul]
  exact Finset.sum_le_sum fun i hi ↦ h i hi x

/-- Upper bounds add up: `Aᵢ ≤ bᵢ` for all `i ∈ s` gives `∑ Aᵢ ≤ ∑ bᵢ`. -/
theorem qform_sum_le {ι : Type*} (s : Finset ι) {A : ι → Matrix n n 𝕜} {b : ι → ℝ}
    (h : ∀ i ∈ s, ∀ x, qform (A i) x ≤ b i * nsq x) (x : n → 𝕜) :
    qform (∑ i ∈ s, A i) x ≤ (∑ i ∈ s, b i) * nsq x := by
  rw [qform_sum, Finset.sum_mul]
  exact Finset.sum_le_sum fun i hi ↦ h i hi x

/-- A lower bound `a ≤ A` holds for every unitary conjugate `M A Mᴴ`. -/
theorem le_qform_conj [DecidableEq n] {M A : Matrix n n 𝕜} (hM : M * Mᴴ = 1) {a : ℝ}
    (h : ∀ y, a * nsq y ≤ qform A y) (x : n → 𝕜) : a * nsq x ≤ qform (M * A * Mᴴ) x := by
  rw [qform_conj, ← nsq_conjTranspose_mulVec hM x]
  exact h _

/-- An upper bound `A ≤ b` holds for every unitary conjugate `M A Mᴴ`. -/
theorem qform_conj_le [DecidableEq n] {M A : Matrix n n 𝕜} (hM : M * Mᴴ = 1) {b : ℝ}
    (h : ∀ y, qform A y ≤ b * nsq y) (x : n → 𝕜) : qform (M * A * Mᴴ) x ≤ b * nsq x := by
  rw [qform_conj, ← nsq_conjTranspose_mulVec hM x]
  exact h _

/-- A norm bound `|⟪y, A y⟫| ≤ c ‖y‖²` holds for every unitary conjugate `M A Mᴴ`. -/
theorem abs_qform_conj_le [DecidableEq n] {M A : Matrix n n 𝕜} (hM : M * Mᴴ = 1) {c : ℝ}
    (h : ∀ y, |qform A y| ≤ c * nsq y) (x : n → 𝕜) : |qform (M * A * Mᴴ) x| ≤ c * nsq x := by
  rw [qform_conj, ← nsq_conjTranspose_mulVec hM x]
  exact h _

/-! ### The weighted Schur test -/

/-- The weighted arithmetic–geometric mean inequality behind the Schur test:
`2 u v ≤ (q'/q) u² + (q/q') v²` for positive weights `q, q'`. -/
theorem two_mul_le_weighted {q q' : ℝ} (hq : 0 < q) (hq' : 0 < q') (u v : ℝ) :
    2 * (u * v) ≤ q' / q * u ^ 2 + q / q' * v ^ 2 := by
  rw [div_mul_eq_mul_div, div_mul_eq_mul_div, div_add_div _ _ hq.ne' hq'.ne',
    le_div_iff₀ (by positivity)]
  nlinarith [sq_nonneg (q' * u - q * v)]

/-- **Weighted Schur test.** Let `K` be a symmetric majorant of the entries of `A`,
`‖Aᵢⱼ‖ ≤ Kᵢⱼ = Kⱼᵢ`, and let `q > 0` be weights with `∑ⱼ Kᵢⱼ qⱼ ≤ c qᵢ` for every row `i`.
Then `|Re ⟪x, A x⟫| ≤ c ‖x‖²`.

For a Hermitian `A` one may take `Kᵢⱼ = ‖Aᵢⱼ‖`; the majorant form also allows a sum of several
Hermitian terms, majorized by the sum of the norms of their entries. -/
theorem abs_qform_le_of_schur {A : Matrix n n 𝕜} {K : n → n → ℝ} {q : n → ℝ} {c : ℝ}
    (hAK : ∀ i j, ‖A i j‖ ≤ K i j) (hK : ∀ i j, K i j = K j i) (hq : ∀ i, 0 < q i)
    (hrow : ∀ i, ∑ j, K i j * q j ≤ c * q i) (x : n → 𝕜) : |qform A x| ≤ c * nsq x := by
  -- Step 1: the triangle inequality bounds the form by `∑ᵢⱼ Kᵢⱼ |xᵢ| |xⱼ|`.
  have h1 : |qform A x| ≤ ∑ i, ∑ j, K i j * (‖x i‖ * ‖x j‖) := by
    calc |qform A x| ≤ ‖star x ⬝ᵥ (A *ᵥ x)‖ := RCLike.abs_re_le_norm _
      _ = ‖∑ i, ∑ j, star (x i) * (A i j * x j)‖ := by
          simp only [dotProduct, mulVec, Pi.star_apply, Finset.mul_sum]
      _ ≤ ∑ i, ∑ j, ‖star (x i) * (A i j * x j)‖ :=
          (norm_sum_le _ _).trans (Finset.sum_le_sum fun i _ ↦ norm_sum_le _ _)
      _ ≤ ∑ i, ∑ j, K i j * (‖x i‖ * ‖x j‖) := by
          refine Finset.sum_le_sum fun i _ ↦ Finset.sum_le_sum fun j _ ↦ ?_
          rw [norm_mul, norm_mul, norm_star]
          calc ‖x i‖ * (‖A i j‖ * ‖x j‖) = ‖A i j‖ * (‖x i‖ * ‖x j‖) := by ring
            _ ≤ K i j * (‖x i‖ * ‖x j‖) :=
                mul_le_mul_of_nonneg_right (hAK i j) (by positivity)
  -- Step 2: split each product by the weighted mean inequality; by the symmetry of `K` the two
  -- halves are equal after exchanging the summation indices.
  set S := ∑ i, ∑ j, K i j * (q j / q i) * ‖x i‖ ^ 2 with hS
  have hswap : ∑ i, ∑ j, K j i * (q i / q j) * ‖x j‖ ^ 2 = S := Finset.sum_comm
  have h2 : ∑ i, ∑ j, K i j * (‖x i‖ * ‖x j‖) ≤ S := by
    have hpt : ∀ i j, K i j * (‖x i‖ * ‖x j‖) ≤
        (K i j * (q j / q i) * ‖x i‖ ^ 2 + K j i * (q i / q j) * ‖x j‖ ^ 2) / 2 := by
      intro i j
      have hK0 : 0 ≤ K i j := (norm_nonneg _).trans (hAK i j)
      have := mul_le_mul_of_nonneg_left (two_mul_le_weighted (hq i) (hq j) ‖x i‖ ‖x j‖) hK0
      rw [← hK i j]
      linarith
    calc ∑ i, ∑ j, K i j * (‖x i‖ * ‖x j‖)
        ≤ ∑ i, ∑ j, (K i j * (q j / q i) * ‖x i‖ ^ 2 + K j i * (q i / q j) * ‖x j‖ ^ 2) / 2 :=
          Finset.sum_le_sum fun i _ ↦ Finset.sum_le_sum fun j _ ↦ hpt i j
      _ = (S + ∑ i, ∑ j, K j i * (q i / q j) * ‖x j‖ ^ 2) / 2 := by
          rw [hS, ← Finset.sum_add_distrib, Finset.sum_div]
          refine Finset.sum_congr rfl fun i _ ↦ ?_
          rw [← Finset.sum_add_distrib, Finset.sum_div]
      _ = S := by rw [hswap]; ring
  -- Step 3: the row-sum condition.
  have h3 : S ≤ c * nsq x := by
    rw [nsq_eq_sum, Finset.mul_sum]
    refine Finset.sum_le_sum fun i _ ↦ ?_
    have hqi := hq i
    calc ∑ j, K i j * (q j / q i) * ‖x i‖ ^ 2 = ‖x i‖ ^ 2 / q i * ∑ j, K i j * q j := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun j _ ↦ ?_
          field_simp
      _ ≤ ‖x i‖ ^ 2 / q i * (c * q i) := mul_le_mul_of_nonneg_left (hrow i) (by positivity)
      _ = c * ‖x i‖ ^ 2 := by field_simp
  linarith

/-! ### Eigenvalues and bounds -/

section Eigen

variable [DecidableEq n] {A : Matrix n n 𝕜} (hA : A.IsHermitian)
include hA

/-- An upper bound on all eigenvalues is an upper bound on the quadratic form. -/
theorem qform_le_of_eigenvalues_le {b : ℝ} (h : ∀ i, hA.eigenvalues i ≤ b) (x : n → 𝕜) :
    qform A x ≤ b * nsq x := by
  rw [← eigU_mulVec_star hA x, qform_eigU, nsq_eigU, nsq_eq_sum, Finset.mul_sum]
  exact Finset.sum_le_sum fun i _ ↦ mul_le_mul_of_nonneg_right (h i) (sq_nonneg _)

/-- A lower bound on all eigenvalues is a lower bound on the quadratic form. -/
theorem le_qform_of_le_eigenvalues {a : ℝ} (h : ∀ i, a ≤ hA.eigenvalues i) (x : n → 𝕜) :
    a * nsq x ≤ qform A x := by
  rw [← eigU_mulVec_star hA x, qform_eigU, nsq_eigU, nsq_eq_sum, Finset.mul_sum]
  exact Finset.sum_le_sum fun i _ ↦ mul_le_mul_of_nonneg_right (h i) (sq_nonneg _)

/-- On a Hermitian matrix, `|⟨x, A x⟩| ≤ ‖x‖²` for every `x` puts the eigenvalues in
`[-1, 1]`: evaluate the form on the eigenvector `U eᵢ`. -/
theorem abs_eigenvalues_le_one_of_qform (h : ∀ x, |qform A x| ≤ nsq x) (i : n) :
    |hA.eigenvalues i| ≤ 1 := by
  have hx := h (eigU hA *ᵥ Pi.single i 1)
  rw [qform_eigU, nsq_eigU, nsq_eq_sum] at hx
  simpa [Pi.single_apply, apply_ite] using hx

/-- **Eigenvalues are roots of an annihilating polynomial.** If `∏ᵣ (A - r) = 0` for a list of
real numbers `r`, every eigenvalue of the Hermitian matrix `A` is in the list. -/
theorem eigenvalues_mem_of_prod_eq_zero {rs : List ℝ}
    (h : (rs.map fun r : ℝ ↦ A - (r : 𝕜) • (1 : Matrix n n 𝕜)).prod = 0) (i : n) :
    hA.eigenvalues i ∈ rs := by
  set v := eigU hA *ᵥ Pi.single i (1 : 𝕜) with hv_def
  have hv : v ≠ 0 := by
    intro h0
    have hback : (eigU hA)ᴴ *ᵥ v = Pi.single i 1 := by
      rw [hv_def, mulVec_mulVec, eigU_star_mul, one_mulVec]
    rw [h0, mulVec_zero] at hback
    simpa using congrFun hback i
  -- On the eigenvector `v`, each factor `A - r` acts as the scalar `λᵢ - r`.
  have hprod : ∀ l : List ℝ, (l.map fun r : ℝ ↦ A - (r : 𝕜) • (1 : Matrix n n 𝕜)).prod *ᵥ v =
      (l.map fun r ↦ ((hA.eigenvalues i - r : ℝ) : 𝕜)).prod • v := by
    intro l
    induction l with
    | nil => simp
    | cons r l ih =>
      rw [List.map_cons, List.prod_cons, ← mulVec_mulVec, ih, mulVec_smul, sub_mulVec,
        smul_mulVec, one_mulVec, hv_def, mulVec_eigU_single, ← hv_def, ← sub_smul,
        smul_smul, List.map_cons, List.prod_cons, mul_comm]
      push_cast
      ring_nf
  have hzero := hprod rs
  rw [h, zero_mulVec, eq_comm, smul_eq_zero] at hzero
  rcases hzero with hz | hz
  · obtain ⟨r, hr, hr0⟩ := List.mem_map.1 (List.prod_eq_zero_iff.1 hz)
    rw [RCLike.ofReal_eq_zero, sub_eq_zero] at hr0
    rwa [hr0]
  · exact absurd hz hv

/-- **Upper bound from an annihilating polynomial.** If `∏ᵣ (A - r) = 0` and every root is at
most `b`, then `A ≤ b`. -/
theorem qform_le_of_prod_eq_zero {rs : List ℝ}
    (h : (rs.map fun r : ℝ ↦ A - (r : 𝕜) • (1 : Matrix n n 𝕜)).prod = 0) {b : ℝ}
    (hb : ∀ r ∈ rs, r ≤ b) (x : n → 𝕜) : qform A x ≤ b * nsq x :=
  qform_le_of_eigenvalues_le hA (fun i ↦ hb _ (eigenvalues_mem_of_prod_eq_zero hA h i)) x

/-- **Lower bound from an annihilating polynomial.** If `∏ᵣ (A - r) = 0` and every root is at
least `a`, then `a ≤ A`. -/
theorem le_qform_of_prod_eq_zero {rs : List ℝ}
    (h : (rs.map fun r : ℝ ↦ A - (r : 𝕜) • (1 : Matrix n n 𝕜)).prod = 0) {a : ℝ}
    (ha : ∀ r ∈ rs, a ≤ r) (x : n → 𝕜) : a * nsq x ≤ qform A x :=
  le_qform_of_le_eigenvalues hA (fun i ↦ ha _ (eigenvalues_mem_of_prod_eq_zero hA h i)) x

end Eigen

/-! ### Gram factorization -/

/-- A Hermitian matrix with a nonnegative quadratic form is a Gram matrix `Cᴴ C`, with
`C = diag(√λ) Uᴴ` in its eigenbasis. -/
theorem exists_eq_conjTranspose_mul_self {A : Matrix n n 𝕜} (hA : A.IsHermitian)
    (h : ∀ x, 0 ≤ qform A x) : ∃ C : Matrix n n 𝕜, A = Cᴴ * C := by
  classical
  have hev : ∀ i, 0 ≤ hA.eigenvalues i := by
    intro i
    have := h (eigU hA *ᵥ Pi.single i 1)
    rw [qform_eigU, Finset.sum_eq_single i (fun k _ hk ↦ by simp [hk]) (by simp)] at this
    simpa using this
  refine ⟨diagonal (fun i ↦ ((√(hA.eigenvalues i) : ℝ) : 𝕜)) * (eigU hA)ᴴ, ?_⟩
  have hsq : (diagonal fun i ↦ ((√(hA.eigenvalues i) : ℝ) : 𝕜))ᴴ *
      diagonal (fun i ↦ ((√(hA.eigenvalues i) : ℝ) : 𝕜)) =
      diagonal (fun i ↦ (hA.eigenvalues i : 𝕜)) := by
    rw [diagonal_conjTranspose, diagonal_mul_diagonal]
    congr 1
    funext i
    rw [Pi.star_apply, RCLike.star_def, RCLike.conj_ofReal, ← RCLike.ofReal_mul,
      Real.mul_self_sqrt (hev i)]
  rw [conjTranspose_mul, conjTranspose_conjTranspose, Matrix.mul_assoc,
    ← Matrix.mul_assoc _ _ (eigU hA)ᴴ, hsq, ← eigU_conj hA]
  simp only [← Matrix.mul_assoc]
  rw [eigU_mul_star, Matrix.one_mul, Matrix.mul_assoc, eigU_mul_star, Matrix.mul_one]

end CertifiedGaps.Spectral
