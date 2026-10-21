/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Spectral.Inertia

/-!
# Counting eigenvalues through a Schur-type compression

Setting. `A` is Hermitian (in applications `A = H - t`), `W : n × m` spans a "low" subspace, and
`A ≥ m₀ > 0` on `S = ker Wᴴ` (the orthogonal complement of the columns of `W`):
`m₀ ‖u‖² ≤ Re ⟪u, A u⟫` whenever `Wᴴ u = 0`. A witness `Z : n × k` (for example `W` plus a
correction) and arbitrary coefficients `B : m × k` define the residual `R = A Z - W B`.

**Compression bound** (`qform_schur_le`). For all `c` and all `u ∈ S`,
`Re ⟪Z c + u, A (Z c + u)⟫ ≥ Re ⟪c, (Zᴴ A Z - δ) c⟫` whenever `‖R‖_F² ≤ m₀ δ`.

Derivation. Expanding, `q(Zc + u) = q(Zc) + 2 Re ⟪u, A Z c⟫ + q(u)`. Since `Wᴴ u = 0`, the cross
term is `⟪u, R c⟫`, and by Cauchy–Schwarz twice `(Re ⟪u, R c⟫)² ≤ ‖u‖² ‖R‖_F² ‖c‖²`. Completing
the square, `2s + m₀ ‖u‖² ≥ -r²/m₀` whenever `s² ≤ r² ‖u‖²` (`qform_add_ge_of_gap`).

**Count** (`eigCount_le_of_schur`, `eigCount_le_of_schur_LDL`). If moreover `range Z + S` is the
whole space (for square `Wᴴ Z` it suffices that `Wᴴ Z` is invertible,
`range_sup_ker_eq_top_of_isUnit`), then the number of eigenvalues of `H` below `t` is at most the
codimension of any subspace on which `Zᴴ (H - t) Z - δ` is positive semidefinite; with an exact
rational `L D Lᵀ` certificate of a matrix `T ≤ Zᴴ (H - t) Z - δ`, it is at most the number of
negative pivots.
-/

open Matrix Module

namespace CertifiedGaps.Spectral

variable {𝕜 : Type*} [RCLike 𝕜] {n m k : Type*} [Fintype n] [Fintype m] [Fintype k]

/-- The squared Frobenius norm `∑ᵢⱼ ‖Rᵢⱼ‖²`. -/
noncomputable def frobSq (R : Matrix n k 𝕜) : ℝ := ∑ i, ∑ j, ‖R i j‖ ^ 2

/-- Cauchy–Schwarz for the real part of the dot product. -/
theorem re_dotProduct_sq_le (u r : n → 𝕜) :
    RCLike.re (star u ⬝ᵥ r) ^ 2 ≤ nsq u * nsq r := by
  have h1 : |RCLike.re (star u ⬝ᵥ r)| ≤ ∑ i, ‖u i‖ * ‖r i‖ := by
    refine (RCLike.abs_re_le_norm _).trans ((norm_sum_le _ _).trans (le_of_eq ?_))
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    rw [Pi.star_apply, norm_mul, norm_star]
  rw [nsq_eq_sum, nsq_eq_sum, ← sq_abs]
  exact (pow_le_pow_left₀ (abs_nonneg _) h1 2).trans (Finset.sum_mul_sq_le_sq_mul_sq _ _ _)

/-- `‖R c‖² ≤ ‖R‖_F² ‖c‖²`. -/
theorem nsq_mulVec_le (R : Matrix n k 𝕜) (c : k → 𝕜) : nsq (R *ᵥ c) ≤ frobSq R * nsq c := by
  rw [nsq_eq_sum, nsq_eq_sum, frobSq, Finset.sum_mul]
  refine Finset.sum_le_sum fun i _ ↦ ?_
  have h1 : ‖(R *ᵥ c) i‖ ≤ ∑ j, ‖R i j‖ * ‖c j‖ :=
    (norm_sum_le _ _).trans (le_of_eq (Finset.sum_congr rfl fun j _ ↦ norm_mul _ _))
  exact (pow_le_pow_left₀ (norm_nonneg _) h1 2).trans (Finset.sum_mul_sq_le_sq_mul_sq _ _ _)

/-- **Completing the square.** If `A ≥ m₀ > 0` at `u` and the cross term satisfies
`(Re ⟪u, A y⟫)² ≤ r² ‖u‖²`, then `Re ⟪y + u, A (y + u)⟫ ≥ Re ⟪y, A y⟫ - r²/m₀`. -/
theorem qform_add_ge_of_gap {A : Matrix n n 𝕜} (hA : A.IsHermitian) {m₀ r2 : ℝ} (hm₀ : 0 < m₀)
    (hr2 : 0 ≤ r2) {y u : n → 𝕜} (hu : m₀ * nsq u ≤ qform A u)
    (hcross : RCLike.re (star u ⬝ᵥ (A *ᵥ y)) ^ 2 ≤ r2 * nsq u) :
    qform A y - r2 / m₀ ≤ qform A (y + u) := by
  rw [qform_add A hA]
  set s := RCLike.re (star u ⬝ᵥ (A *ᵥ y))
  set N := nsq u
  have hN : 0 ≤ N := nsq_nonneg u
  -- Key inequality `m₀ (2s + m₀ N) + r2 ≥ 0`, from `(m₀ N + s)² ≥ 0` and `s² ≤ r2 N`.
  have key : 0 ≤ m₀ * (2 * s + m₀ * N) + r2 := by
    rcases hN.eq_or_lt with h0 | hpos
    · rw [← h0, mul_zero] at hcross
      have : s = 0 := by nlinarith
      rw [this, ← h0]; simpa using hr2
    · have h := sq_nonneg (m₀ * N + s)
      by_contra hneg
      push Not at hneg
      nlinarith [mul_neg_of_pos_of_neg hpos hneg]
  have : -(r2 / m₀) ≤ 2 * s + m₀ * N := by
    rw [neg_le, le_div_iff₀ hm₀]
    linarith
  linarith

variable [DecidableEq n]

omit [DecidableEq n] in
/-- **Compression bound.** If `A` is Hermitian, `A ≥ m₀ > 0` on `ker Wᴴ`, and the residual
`R = A Z - W B` satisfies `‖R‖_F² ≤ m₀ δ`, then for every `c` and every `u` with `Wᴴ u = 0`,
`Re ⟪c, (Zᴴ A Z - δ) c⟫ ≤ Re ⟪Z c + u, A (Z c + u)⟫`. -/
theorem qform_schur_le {A : Matrix n n 𝕜} (hA : A.IsHermitian) (W : Matrix n m 𝕜)
    (Z : Matrix n k 𝕜) (B : Matrix m k 𝕜) {m₀ δ : ℝ} (hm₀ : 0 < m₀)
    (hgap : ∀ u, Wᴴ *ᵥ u = 0 → m₀ * nsq u ≤ qform A u)
    (hδ : frobSq (A * Z - W * B) ≤ m₀ * δ) [DecidableEq k] (c : k → 𝕜) {u : n → 𝕜}
    (hu : Wᴴ *ᵥ u = 0) :
    qform (Zᴴ * A * Z - (δ : 𝕜) • 1) c ≤ qform A (Z *ᵥ c + u) := by
  set R := A * Z - W * B
  have hfrob : 0 ≤ frobSq R := Finset.sum_nonneg fun i _ ↦ Finset.sum_nonneg fun j _ ↦
    by positivity
  -- The `W B` part of `A Z c` is invisible to `u ∈ ker Wᴴ`.
  have hcross_eq : star u ⬝ᵥ (A *ᵥ (Z *ᵥ c)) = star u ⬝ᵥ (R *ᵥ c) := by
    have hW : star u ⬝ᵥ (W *ᵥ (B *ᵥ c)) = 0 := by
      rw [dotProduct_mulVec, ← conjTranspose_conjTranspose W, ← star_mulVec, hu, star_zero,
        zero_dotProduct]
    simp only [R, sub_mulVec, dotProduct_sub, ← mulVec_mulVec, hW, sub_zero]
  have hcross : RCLike.re (star u ⬝ᵥ (A *ᵥ (Z *ᵥ c))) ^ 2 ≤ (frobSq R * nsq c) * nsq u := by
    rw [hcross_eq]
    refine (re_dotProduct_sq_le u (R *ᵥ c)).trans ?_
    rw [mul_comm (frobSq R * nsq c)]
    exact mul_le_mul_of_nonneg_left (nsq_mulVec_le R c) (nsq_nonneg u)
  have h := qform_add_ge_of_gap hA hm₀ (mul_nonneg hfrob (nsq_nonneg c)) (hgap u hu) hcross
  have hδc : frobSq R * nsq c / m₀ ≤ δ * nsq c := by
    rw [div_le_iff₀ hm₀]
    nlinarith [nsq_nonneg c]
  rw [qform_sub_smul_one, ← qform_mulVec]
  linarith

omit [DecidableEq n] in
/-- If `Wᴴ Z` is invertible then every vector is `Z c + u` with `Wᴴ u = 0`. -/
theorem range_sup_ker_eq_top_of_isUnit [DecidableEq k] (W Z : Matrix n k 𝕜)
    (h : IsUnit (Wᴴ * Z).det) :
    LinearMap.range Z.mulVecLin ⊔ LinearMap.ker Wᴴ.mulVecLin = ⊤ := by
  refine eq_top_iff.2 fun x _ ↦ Submodule.mem_sup.2 ?_
  set G := (Wᴴ * Z)⁻¹
  refine ⟨Z *ᵥ (G *ᵥ (Wᴴ *ᵥ x)), ⟨_, rfl⟩, x - Z *ᵥ (G *ᵥ (Wᴴ *ᵥ x)), ?_, by abel⟩
  rw [LinearMap.mem_ker, Matrix.mulVecLin_apply, mulVec_sub, mulVec_mulVec, mulVec_mulVec,
    mul_nonsing_inv _ h, one_mulVec, sub_self]

/-- **Count through a Schur compression.** Let `H` be Hermitian, `t` a cut, and suppose
`H - t ≥ m₀ > 0` on `ker Wᴴ`, `‖(H - t) Z - W B‖_F² ≤ m₀ δ`, and `range Z + ker Wᴴ` is the
whole space. If `Zᴴ (H - t) Z - δ` is positive semidefinite on `P`, then
`#{i : λᵢ(H) < t} + dim P ≤ k`. -/
theorem eigCount_le_of_schur {H : Matrix n n 𝕜} (hH : H.IsHermitian) (t : ℝ)
    (W : Matrix n m 𝕜) (Z : Matrix n k 𝕜) (B : Matrix m k 𝕜) {m₀ δ : ℝ} (hm₀ : 0 < m₀)
    (hgap : ∀ u, Wᴴ *ᵥ u = 0 → (t + m₀) * nsq u ≤ qform H u)
    (hδ : frobSq ((H - (t : 𝕜) • 1) * Z - W * B) ≤ m₀ * δ)
    (hspan : LinearMap.range Z.mulVecLin ⊔ LinearMap.ker Wᴴ.mulVecLin = ⊤) [DecidableEq k]
    (P : Submodule 𝕜 (k → 𝕜))
    (hP : ∀ c ∈ P, 0 ≤ qform (Zᴴ * (H - (t : 𝕜) • 1) * Z - (δ : 𝕜) • 1) c) :
    eigCount hH t + finrank 𝕜 P ≤ Fintype.card k := by
  have hA : (H - (t : 𝕜) • 1).IsHermitian := by
    refine hH.sub ?_
    rw [IsHermitian, conjTranspose_smul, conjTranspose_one, RCLike.star_def, RCLike.conj_ofReal]
  refine eigCount_le_of_compression hH Z _ hspan P fun c hc u hu ↦ ?_
  have hgap' : ∀ u, Wᴴ *ᵥ u = 0 → m₀ * nsq u ≤ qform (H - (t : 𝕜) • 1) u := fun u hu ↦ by
    rw [qform_sub_smul_one]; linarith [hgap u hu]
  have h := qform_schur_le hA W Z B hm₀ hgap' hδ c (LinearMap.mem_ker.1 hu)
  have hPc := hP c hc
  simp only [qform_sub_smul_one] at h hPc
  linarith

/-- **Count through a Schur compression, rational certificate.** As `eigCount_le_of_schur`,
with an exact `L D Lᵀ` certificate for a rational matrix `T` whose form lies below that of
`Zᴴ (H - t) Z - δ`: then `#{i : λᵢ(H) < t} ≤ C.negCount`. -/
theorem eigCount_le_of_schur_LDL {K : ℕ} {H : Matrix n n 𝕜} (hH : H.IsHermitian) (t : ℝ)
    (W : Matrix n m 𝕜) (Z : Matrix n (Fin K) 𝕜) (B : Matrix m (Fin K) 𝕜) {m₀ δ : ℝ}
    (hm₀ : 0 < m₀) (hgap : ∀ u, Wᴴ *ᵥ u = 0 → (t + m₀) * nsq u ≤ qform H u)
    (hδ : frobSq ((H - (t : 𝕜) • 1) * Z - W * B) ≤ m₀ * δ)
    (hspan : LinearMap.range Z.mulVecLin ⊔ LinearMap.ker Wᴴ.mulVecLin = ⊤)
    (C : LDLCert K) (T : Matrix (Fin K) (Fin K) ℚ) (hC : C.Factors T)
    (hT : ∀ c, qform (T.map fun q : ℚ ↦ (q : 𝕜)) c ≤
      qform (Zᴴ * (H - (t : 𝕜) • 1) * Z - (δ : 𝕜) • 1) c) :
    eigCount hH t ≤ C.negCount := by
  obtain ⟨P, hPdim, hP⟩ := C.exists_nonneg_subspace (𝕜 := 𝕜) hC
  have h := eigCount_le_of_schur hH t W Z B hm₀ hgap hδ hspan P fun c hc ↦
    (hP c hc).trans (hT c)
  rw [Fintype.card_fin] at h
  omega

end CertifiedGaps.Spectral
