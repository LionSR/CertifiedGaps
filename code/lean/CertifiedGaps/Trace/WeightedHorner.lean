/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Trace.Horner
import CertifiedGaps.Spectral.FunctionalCalculus
import CertifiedGaps.Spectral.OperatorBounds

/-!
# Integer Horner columns in a weighted norm

`Trace/Horner.lean` bounds the rounding error of the integer recursion
`X ← ⌊D X / s⌋ + cⱼ e_r` when `R = D / s` is a Euclidean contraction. The thermal appendix of
OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026) uses an integer matrix
`D = s J A J⁻¹` that is similar to a contraction `A` by a positive diagonal matrix `J`, and is
itself not a Euclidean contraction. It is a contraction in the weighted norm
`‖v‖_ω² = ∑ᵢ ωᵢ vᵢ²` with `ωᵢ = Jᵢ⁻²`, because `‖v‖_ω = ‖J⁻¹ v‖`. This file proves the rounding
bound in that norm.

* `IsWeightedContraction D s ω`: `‖D v‖_ω ≤ s ‖v‖_ω` for every real `v`.
* `isWeightedContraction_of_similar`: `D = s J A J⁻¹` with `A` a Hermitian contraction is a
  weighted contraction with `ω = J⁻²`.
* `hornerFloor_error_weighted`: with `R = D / s` a weighted contraction and `|cⱼ - aⱼ| ≤ 1`,
  `‖X - (∑_{j ≤ J} aⱼ Rʲ) e_r‖_ω ≤ (J + 1)(δ + γ)` whenever `√(∑ᵢ ωᵢ) ≤ δ` and `√ω_r ≤ γ`.

Each step adds the floor error, of weighted norm at most `√(∑ ωᵢ)`, and the coefficient error
`(cⱼ - aⱼ) e_r`, of weighted norm at most `√ω_r`; the weighted contraction does not enlarge the
error already present. With `ω ≡ 1` this is `hornerFloor_error`.
-/

open Matrix

namespace CertifiedGaps.Trace

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- `R = D / s` is a contraction in the weighted norm `‖v‖_ω² = ∑ᵢ ωᵢ vᵢ²`, stated without
division: `‖D v‖_ω² ≤ s² ‖v‖_ω²` for every real `v`. -/
def IsWeightedContraction (D : Matrix ι ι ℤ) (s : ℕ) (ω : ι → ℝ) : Prop :=
  ∀ v : ι → ℝ, ∑ i, ω i * ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i ^ 2 ≤ (s : ℝ) ^ 2 * ∑ i, ω i * v i ^ 2

/-! ### The weighted norm -/

/-- The weighted norm `√(∑ ωᵢ vᵢ²)`, routed through `EuclideanSpace` (coordinates `√ωᵢ vᵢ`) to
reuse its triangle inequality. -/
private noncomputable abbrev lw (ω : ι → ℝ) (v : ι → ℝ) : ℝ :=
  ‖WithLp.toLp 2 fun i ↦ Real.sqrt (ω i) * v i‖

omit [DecidableEq ι] in
private theorem lw_eq {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) (v : ι → ℝ) :
    lw ω v = Real.sqrt (∑ i, ω i * v i ^ 2) := by
  simp only [lw, EuclideanSpace.norm_eq, Real.norm_eq_abs, sq_abs, mul_pow, Real.sq_sqrt (hω _)]

omit [DecidableEq ι] in
private theorem lw_add (ω : ι → ℝ) (u v : ι → ℝ) : lw ω (u + v) ≤ lw ω u + lw ω v := by
  have h : (WithLp.toLp 2 fun i ↦ Real.sqrt (ω i) * (u + v) i : EuclideanSpace ℝ ι) =
      WithLp.toLp 2 (fun i ↦ Real.sqrt (ω i) * u i) +
        WithLp.toLp 2 (fun i ↦ Real.sqrt (ω i) * v i) := by
    rw [← WithLp.toLp_add]; congr 1; funext i; simp [mul_add]
  simp only [lw, h]
  exact norm_add_le _ _

omit [DecidableEq ι] in
private theorem lw_le_of_abs_le_one {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) (u : ι → ℝ)
    (hu : ∀ i, |u i| ≤ 1) : lw ω u ≤ Real.sqrt (∑ i, ω i) := by
  rw [lw_eq hω]
  refine Real.sqrt_le_sqrt (Finset.sum_le_sum fun i _ ↦ ?_)
  exact mul_le_of_le_one_right (hω i) ((sq_le_one_iff_abs_le_one _).2 (hu i))

private theorem lw_single {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) (r : ι) (x : ℝ) :
    lw ω (Pi.single r x) = Real.sqrt (ω r) * |x| := by
  rw [lw_eq hω, ← Real.sqrt_sq_eq_abs, ← Real.sqrt_mul (hω r)]
  congr 1
  rw [Finset.sum_eq_single r (fun i _ hi ↦ by simp [hi]) (by simp)]
  simp

omit [DecidableEq ι] in
private theorem lw_scaled_le {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) {D : Matrix ι ι ℤ} {s : ℕ}
    (hs : 0 < s) (hD : IsWeightedContraction D s ω) (v : ι → ℝ) :
    lw ω (scaledMatrix D s *ᵥ v) ≤ lw ω v := by
  rw [lw_eq hω, lw_eq hω]
  refine Real.sqrt_le_sqrt ?_
  have hs' : (0 : ℝ) < s := Nat.cast_pos.2 hs
  have h := hD v
  simp only [scaledMatrix, smul_mulVec, Pi.smul_apply, smul_eq_mul]
  have e : ∀ i, ω i * ((s : ℝ)⁻¹ * ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i) ^ 2 =
      ((s : ℝ) ^ 2)⁻¹ * (ω i * ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i ^ 2) := fun i ↦ by ring
  rw [Finset.sum_congr rfl fun i _ ↦ e i, ← Finset.mul_sum]
  calc ((s : ℝ) ^ 2)⁻¹ * ∑ i, ω i * ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i ^ 2
      ≤ ((s : ℝ) ^ 2)⁻¹ * ((s : ℝ) ^ 2 * ∑ i, ω i * v i ^ 2) := by gcongr
    _ = ∑ i, ω i * v i ^ 2 := by field_simp

/-! ### One step and the accumulated error -/

omit [DecidableEq ι] in
/-- `R` applied to an integer vector is the exact integer product divided by `s`. -/
private theorem scaledMatrix_mulVec_intCast' (D : Matrix ι ι ℤ) (s : ℕ) (X : ι → ℤ) (i : ι) :
    (scaledMatrix D s *ᵥ fun j ↦ (X j : ℝ)) i = ((D *ᵥ X) i : ℝ) / s := by
  simp [scaledMatrix, mulVec, dotProduct, div_eq_inv_mul, Finset.mul_sum, mul_assoc]

/-- One floored step moves the weighted error by at most `√(∑ ωᵢ)` (vector floor) plus
`√ω_r` (coefficient error), on top of the propagated error, which `R` does not enlarge. -/
private theorem step_error_weighted {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) {D : Matrix ι ι ℤ}
    {s : ℕ} (hs : 0 < s) (hD : IsWeightedContraction D s ω) (r : ι) {c : ℤ} {a : ℝ}
    (hca : |(c : ℝ) - a| ≤ 1) (X : ι → ℤ) (Y : ι → ℝ) :
    lw ω (fun i ↦ (hornerStep D s r c X i : ℝ) - (scaledMatrix D s *ᵥ Y + Pi.single r a : ι → ℝ) i)
      ≤ lw ω (fun i ↦ (X i : ℝ) - Y i) + Real.sqrt (∑ i, ω i) + Real.sqrt (ω r) := by
  obtain ⟨u, hu⟩ : ∃ u : ι → ℝ,
      u = fun i ↦ (((D *ᵥ X) i / (s : ℤ) : ℤ) : ℝ) - ((D *ᵥ X) i : ℝ) / s := ⟨_, rfl⟩
  have hdecomp : (fun i ↦ (hornerStep D s r c X i : ℝ) - (scaledMatrix D s *ᵥ Y +
      Pi.single r a : ι → ℝ) i) = u + scaledMatrix D s *ᵥ (fun i ↦ (X i : ℝ) - Y i) +
      Pi.single r ((c : ℝ) - a) := by
    funext i
    have hsub : (fun i ↦ (X i : ℝ) - Y i) = (fun i ↦ (X i : ℝ)) - Y := rfl
    rw [hsub, mulVec_sub, hu]
    simp only [hornerStep, Pi.add_apply, Pi.sub_apply, scaledMatrix_mulVec_intCast',
      Pi.single_apply]
    push_cast
    split_ifs <;> ring
  have hu1 : ∀ i, |u i| ≤ 1 := by
    intro i
    rw [hu]; dsimp only; rw [intDiv_eq_floor]
    have h1 := Int.floor_le (((D *ᵥ X) i : ℝ) / s)
    have h2 := Int.lt_floor_add_one (((D *ᵥ X) i : ℝ) / s)
    rw [abs_le]; constructor <;> linarith
  have h1 := lw_le_of_abs_le_one hω u hu1
  have h2 := lw_scaled_le hω hs hD (fun i ↦ (X i : ℝ) - Y i)
  have h3 := lw_single hω r ((c : ℝ) - a)
  have h3' : Real.sqrt (ω r) * |(c : ℝ) - a| ≤ Real.sqrt (ω r) :=
    mul_le_of_le_one_right (Real.sqrt_nonneg _) hca
  have h4 := lw_add ω (u + scaledMatrix D s *ᵥ (fun i ↦ (X i : ℝ) - Y i))
    (Pi.single r ((c : ℝ) - a))
  have h5 := lw_add ω u (scaledMatrix D s *ᵥ (fun i ↦ (X i : ℝ) - Y i))
  rw [hdecomp]
  linarith

/-! ### Main error bounds -/

/-- **Weighted rounding error of the integer Horner recursion (coefficient lists).** If
`R = D / s` is a weighted contraction and every integer coefficient is within `1` of its real
target, the floored recursion stays within `(number of steps) · (δ + γ)` of the ideal recursion
in the weighted norm, for any `δ ≥ √(∑ ωᵢ)` and `γ ≥ √ω_r`. -/
theorem hornerList_error_weighted {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) {D : Matrix ι ι ℤ} {s : ℕ}
    (hs : 0 < s) (hD : IsWeightedContraction D s ω) {δ γ : ℝ}
    (hδ : Real.sqrt (∑ i, ω i) ≤ δ) (r : ι) (hγ : Real.sqrt (ω r) ≤ γ) {c : List ℤ}
    {a : List ℝ} (hca : List.Forall₂ (fun (cj : ℤ) (aj : ℝ) ↦ |(cj : ℝ) - aj| ≤ 1) c a) :
    Real.sqrt (∑ i, ω i *
        ((hornerList D s r c i : ℝ) - idealList (scaledMatrix D s) r a i) ^ 2) ≤
      c.length * (δ + γ) := by
  rw [← lw_eq hω]
  induction hca with
  | nil =>
    rw [show (fun i ↦ (hornerList D s r [] i : ℝ) - idealList (scaledMatrix D s) r [] i) = 0
      from funext fun i ↦ by simp]
    simp only [lw, Pi.zero_apply, mul_zero, List.length_nil, Nat.cast_zero, zero_mul]
    rw [show (fun _ : ι ↦ (0 : ℝ)) = 0 from rfl, WithLp.toLp_zero, norm_zero]
  | @cons cj aj c a hj _ ih =>
    rw [hornerList_cons, idealList_cons, List.length_cons]
    refine (step_error_weighted hω hs hD r hj _ _).trans ?_
    push_cast
    linarith

/-- **Weighted rounding error of `hornerFloor`.** With `R = D / s` a weighted contraction,
`|cⱼ - aⱼ| ≤ 1` for `j ≤ J`, `δ ≥ √(∑ ωᵢ)` and `γ ≥ √ω_r`, the integer column is within
`(J + 1)(δ + γ)` of `(∑_{j ≤ J} aⱼ Rʲ) e_r` in the weighted norm. -/
theorem hornerFloor_error_weighted {ω : ι → ℝ} (hω : ∀ i, 0 ≤ ω i) {D : Matrix ι ι ℤ}
    {s : ℕ} (hs : 0 < s) (hD : IsWeightedContraction D s ω) {δ γ : ℝ}
    (hδ : Real.sqrt (∑ i, ω i) ≤ δ) (r : ι) (hγ : Real.sqrt (ω r) ≤ γ) {c : ℕ → ℤ}
    {a : ℕ → ℝ} (J : ℕ) (hca : ∀ j ≤ J, |(c j : ℝ) - a j| ≤ 1) :
    Real.sqrt (∑ i, ω i * ((hornerFloor D s c r J i : ℝ) -
        ((∑ j ∈ Finset.range (J + 1), a j • scaledMatrix D s ^ j) *ᵥ Pi.single r 1) i) ^ 2) ≤
      (J + 1) * (δ + γ) := by
  have h2 : List.Forall₂ (fun (cj : ℤ) (aj : ℝ) ↦ |(cj : ℝ) - aj| ≤ 1)
      ((List.range (J + 1)).map c) ((List.range (J + 1)).map a) := by
    rw [List.forall₂_map_left_iff, List.forall₂_map_right_iff, List.forall₂_same]
    intro j hj
    exact hca j (Nat.lt_succ_iff.1 (List.mem_range.1 hj))
  have := hornerList_error_weighted hω hs hD hδ r hγ h2
  rw [idealList_range_map, List.length_map, List.length_range] at this
  simpa [hornerFloor] using this

/-! ### Weighted contractions from a diagonal similarity -/

omit [DecidableEq ι] in
/-- A Hermitian matrix with `|⟨x, A x⟩| ≤ ‖x‖²` does not enlarge norms: `‖A y‖ ≤ ‖y‖`. -/
theorem nsq_mulVec_le_of_qform {A : Matrix ι ι ℂ} (hA : A.IsHermitian)
    (hcon : ∀ x, |Spectral.qform A x| ≤ Spectral.nsq x) (y : ι → ℂ) :
    Spectral.nsq (A *ᵥ y) ≤ Spectral.nsq y := by
  classical
  have h := Spectral.nsq_hermFun_mulVec_le hA id (η := 1)
    (fun i ↦ by simpa using Spectral.abs_eigenvalues_le_one_of_qform hA hcon i) y
  rwa [Spectral.hermFun_id, one_pow, one_mul] at h

/-- **Weighted contraction from a diagonal similarity.** If `D = s J A J⁻¹` over `ℂ`, with `J`
the positive diagonal matrix of `j` and `A` a Hermitian contraction, then `R = D / s` is a
contraction in the norm with weights `ωᵢ = jᵢ⁻²`, since `‖v‖_ω = ‖J⁻¹ v‖`. -/
theorem isWeightedContraction_of_similar {A : Matrix ι ι ℂ} (hA : A.IsHermitian)
    (hcon : ∀ x, |Spectral.qform A x| ≤ Spectral.nsq x) {j : ι → ℝ} (hj : ∀ i, 0 < j i)
    {D : Matrix ι ι ℤ} {s : ℕ}
    (hD : D.map ((↑) : ℤ → ℂ) = (s : ℂ) • (diagonal (fun i ↦ (j i : ℂ)) * A *
      diagonal (fun i ↦ ((j i)⁻¹ : ℂ)))) :
    IsWeightedContraction D s fun i ↦ (j i)⁻¹ ^ 2 := by
  intro v
  -- `y = J⁻¹ v` and `J⁻¹ D v = s A y`.
  set y : ι → ℂ := fun i ↦ ((j i)⁻¹ * v i : ℝ)
  have hDv : ∀ i, ((((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i : ℝ) : ℂ) = (s : ℂ) * (j i : ℂ) *
      (A *ᵥ y) i := by
    intro i
    have h1 : ((((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i : ℝ) : ℂ) =
        ((D.map ((↑) : ℤ → ℂ)) *ᵥ fun k ↦ (v k : ℂ)) i := by
      simp [mulVec, dotProduct]
    rw [h1, hD]
    simp only [smul_mulVec, Pi.smul_apply, smul_eq_mul, ← mulVec_mulVec, mulVec_diagonal]
    have hy : (diagonal fun i ↦ ((j i : ℂ))⁻¹) *ᵥ (fun k ↦ (v k : ℂ)) = y := by
      funext k
      rw [mulVec_diagonal]
      simp only [y]
      push_cast
      ring
    rw [hy]
    ring
  have hlhs : ∑ i, (j i)⁻¹ ^ 2 * ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i ^ 2 =
      (s : ℝ) ^ 2 * Spectral.nsq (A *ᵥ y) := by
    rw [Spectral.nsq_eq_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    have h := congrArg (fun z : ℂ ↦ ‖z‖ ^ 2) (hDv i)
    simp only [norm_mul, Complex.norm_real, Complex.norm_natCast, Real.norm_eq_abs, mul_pow,
      sq_abs] at h
    have hji : (j i) ≠ 0 := (hj i).ne'
    rw [h]; field_simp
  have hrhs : ∑ i, (j i)⁻¹ ^ 2 * v i ^ 2 = Spectral.nsq y := by
    rw [Spectral.nsq_eq_sum]
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    simp only [y, Complex.norm_real, Real.norm_eq_abs, sq_abs, mul_pow]
  rw [hlhs, hrhs]
  exact mul_le_mul_of_nonneg_left (nsq_mulVec_le_of_qform hA hcon y) (by positivity)

end CertifiedGaps.Trace
