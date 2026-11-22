/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.PurityCriterion

/-!
# Three-direction spacetime bootstrap

The indices of `CubeData.Z` count doublings in two spatial directions and imaginary time.
The hypotheses are scalar moment inequalities. They do not assume a spectral gap.
`WindowInput` states the spatial moment interpolation inequality connecting these data to
the thermal purities of the square tori. Model-specific transfer and seed proofs are separate.
-/

open Real

namespace CertifiedGaps.Spacetime

/-- Purity for doubling the first spatial period. -/
noncomputable def px (Z : ℕ → ℕ → ℕ → ℝ) (x y t : ℕ) := Z (x + 1) y t / Z x y t ^ 2
/-- Purity for doubling the second spatial period. -/
noncomputable def py (Z : ℕ → ℕ → ℕ → ℝ) (x y t : ℕ) := Z x (y + 1) t / Z x y t ^ 2
/-- Purity for doubling inverse temperature. -/
noncomputable def pt (Z : ℕ → ℕ → ℕ → ℝ) (x y t : ℕ) := Z x y (t + 1) / Z x y t ^ 2

/-- Positivity, upper bounds, and the single-direction squaring inequalities. -/
structure CubeData where
  Z : ℕ → ℕ → ℕ → ℝ
  positive : ∀ x y t, 0 < Z x y t
  x_le : ∀ x y t, px Z x y t ≤ 1
  y_le : ∀ x y t, py Z x y t ≤ 1
  t_le : ∀ x y t, pt Z x y t ≤ 1
  x_square : ∀ x y t, 1 - px Z (x + 1) y t ≤ squaringDefect (1 - px Z x y t)
  y_square : ∀ x y t, 1 - py Z x (y + 1) t ≤ squaringDefect (1 - py Z x y t)
  t_square : ∀ x y t, 1 - pt Z x y (t + 1) ≤ squaringDefect (1 - pt Z x y t)

theorem px_pos (D : CubeData) (x y t : ℕ) : 0 < px D.Z x y t :=
  div_pos (D.positive _ _ _) (sq_pos_of_pos (D.positive _ _ _))

theorem py_pos (D : CubeData) (x y t : ℕ) : 0 < py D.Z x y t :=
  div_pos (D.positive _ _ _) (sq_pos_of_pos (D.positive _ _ _))

theorem pt_pos (D : CubeData) (x y t : ℕ) : 0 < pt D.Z x y t :=
  div_pos (D.positive _ _ _) (sq_pos_of_pos (D.positive _ _ _))

/-- The product of ratios along the path x, y, imaginary time through a cube. -/
noncomputable def cube (D : CubeData) (x y t : ℕ) : ℝ :=
  px D.Z x y t ^ 4 * py D.Z (x + 1) y t ^ 2 * pt D.Z (x + 1) (y + 1) t

/-- Three staggered seed inequalities at one scale. -/
def SeedBound (D : CubeData) (ε : ℝ) (x y t : ℕ) : Prop :=
  1 - ε ≤ px D.Z x y t ∧ 1 - ε ≤ py D.Z (x + 1) y t ∧
    1 - ε ≤ pt D.Z (x + 1) (y + 1) t

/-- Cyclically permute the two spatial coordinates and imaginary time. -/
def CubeData.rotate (D : CubeData) : CubeData where
  Z := fun x y t => D.Z t x y
  positive := fun x y t => D.positive t x y
  x_le := fun x y t => D.y_le t x y
  y_le := fun x y t => D.t_le t x y
  t_le := fun x y t => D.x_le t x y
  x_square := fun x y t => D.y_square t x y
  y_square := fun x y t => D.t_square t x y
  t_square := fun x y t => D.x_square t x y

theorem cube_lower (D : CubeData) {ε : ℝ} (_hε : 0 ≤ ε) (hε' : ε ≤ 1)
    {x y t : ℕ} (h : SeedBound D ε x y t) : 1 - 7 * ε ≤ cube D x y t := by
  have hp : (1 - ε) ^ 7 ≤ cube D x y t := by
    rw [show (1 - ε) ^ 7 = (1 - ε) ^ 4 * (1 - ε) ^ 2 * (1 - ε) by ring]
    exact mul_le_mul (mul_le_mul (pow_le_pow_left₀ (sub_nonneg.mpr hε') h.1 4)
      (pow_le_pow_left₀ (sub_nonneg.mpr hε') h.2.1 2) (sq_nonneg _)
      (pow_nonneg (px_pos D _ _ _).le _)) h.2.2 (sub_nonneg.mpr hε')
      (mul_nonneg (pow_nonneg (px_pos D _ _ _).le _) (sq_nonneg _))
  simpa using (one_add_mul_le_pow (a := -ε) (by linarith : -2 ≤ -ε) 7).trans hp

theorem cube_le_x_last (D : CubeData) (x y t : ℕ) :
    cube D x y t ≤ px D.Z x (y + 1) (t + 1) := by
  have heq : cube D x y t =
      py D.Z x y t ^ 4 * pt D.Z x (y + 1) t ^ 2 * px D.Z x (y + 1) (t + 1) := by
    unfold cube px py pt
    field_simp [(D.positive x y t).ne', (D.positive (x + 1) y t).ne',
      (D.positive (x + 1) (y + 1) t).ne', (D.positive x (y + 1) t).ne',
      (D.positive x (y + 1) (t + 1)).ne']
  rw [heq]
  exact mul_le_of_le_one_left (px_pos D _ _ _).le
    ((mul_le_mul_of_nonneg_right (pow_le_one₀ (py_pos D _ _ _).le (D.y_le _ _ _))
      (sq_nonneg _)).trans (by simpa using
        (pow_le_one₀ (n := 2) (pt_pos D _ _ _).le (D.t_le _ _ _))))

theorem squaring_small {ε z : ℝ} (hε : 0 ≤ ε) (hε' : ε ≤ 1 / 1000)
    (hz : 0 ≤ z) (hze : z ≤ 7 * ε) : squaringDefect z ≤ 25 * ε ^ 2 := by
  have hd : 0 ≤ 50 * (1 - z) ^ 2 - 49 := by nlinarith
  rw [squaringDefect, div_le_iff₀ (by nlinarith : 0 < 2 * (1 - z) ^ 2)]
  nlinarith [mul_nonneg (sq_nonneg ε) hd,
    mul_nonneg (sub_nonneg.mpr hze) (by positivity : 0 ≤ 7 * ε + z)]

theorem advance_x (D : CubeData) {ε : ℝ} (hε : 0 ≤ ε) (hε' : ε ≤ 1 / 1000)
    {x y t : ℕ} (h : SeedBound D ε x y t) :
    1 - 25 * ε ^ 2 ≤ px D.Z (x + 1) (y + 1) (t + 1) := by
  have hb := (cube_lower D hε (by linarith) h).trans (cube_le_x_last D x y t)
  have hs := (D.x_square x (y + 1) (t + 1)).trans
    (squaring_small hε hε' (sub_nonneg.mpr (D.x_le _ _ _)) (by linarith))
  linarith

/-- One complete cycle of directional updates squares the common small-defect bound. -/
theorem cube_growth (D : CubeData) {ε : ℝ} (hε : 0 ≤ ε) (hε' : ε ≤ 1 / 1000)
    {x y t : ℕ} (h : SeedBound D ε x y t) :
    SeedBound D (25 * ε ^ 2) (x + 1) (y + 1) (t + 1) := by
  have hr : 25 * ε ^ 2 ≤ ε := by nlinarith
  have hx := advance_x D hε hε' h
  have hy := advance_x D.rotate hε hε' (x := y) (y := t) (t := x + 1)
    ⟨h.2.1, h.2.2, (sub_le_sub_left hr 1).trans hx⟩
  have ht := advance_x D.rotate.rotate hε hε' (x := t) (y := x + 1) (t := y + 1)
    ⟨h.2.2, (sub_le_sub_left hr 1).trans hx, (sub_le_sub_left hr 1).trans hy⟩
  exact ⟨hx, hy, ht⟩

/-- The three staggered defects decrease doubly exponentially with the doubling index. -/
theorem cascade (D : CubeData) {ε : ℝ} (hε : 0 < ε) (hε' : ε ≤ 1 / 1000)
    (h : SeedBound D ε 0 0 0) (j : ℕ) :
    SeedBound D (coefficientDefect 25 ε j) j j j := by
  induction j with
  | zero => simpa [coefficientDefect] using h
  | succ j ih =>
    rw [coefficientDefect_succ (by norm_num : (25 : ℝ) ≠ 0)]
    obtain ⟨hu0, hu⟩ := coefficientDefect_mem (by norm_num : 0 < (25 : ℝ)) hε
      (by linarith : 25 * ε ≤ 1) j
    exact cube_growth D hu0.le (hu.trans hε') ih

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

/-- The spatial moment interpolation input on all even square periods in a dyadic window. -/
def WindowInput (D : CubeData) (E : (L : ℕ) → ι L → ℝ) (n₀ : ℕ) (β₀ : ℝ) : Prop :=
  ∀ j L, n₀ * 2 ^ j ≤ L → L ≤ 2 * (n₀ * 2 ^ j) → Even L →
    cube D j j j ≤ thermalPurity E L (β₀ * 2 ^ j)

/-- Three finite seeds and scalar moment inequalities give a uniform gap on even square tori. -/
theorem gapGT_of_cube_seeds (D : CubeData) (E : (L : ℕ) → ι L → ℝ)
    {n₀ : ℕ} (hn₀ : 0 < n₀) {β₀ ε : ℝ} (hβ₀ : 0 < β₀)
    (hε : 0 < ε) (hε' : ε ≤ 1 / 1000) (hW : WindowInput D E n₀ β₀)
    (hseed : SeedBound D ε 0 0 0) {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) :
    GapGT (E L) (log (1 / (25 * ε)) / β₀) := by
  obtain ⟨j, hjL, hLj⟩ := exists_dyadic_window hn₀ hL
  obtain ⟨hu0, hu⟩ := coefficientDefect_mem (by norm_num : 0 < (25 : ℝ)) hε
    (by linarith : 25 * ε ≤ 1) j
  have hP := (cube_lower D hu0.le (by linarith : coefficientDefect 25 ε j ≤ 1)
    (cascade D hε hε' hseed j)).trans (hW j L hjL hLj hLe)
  refine gapGT_of_thermalPurity E L (by positivity) (by linarith) hP ?_
  have hexp : exp (β₀ * 2 ^ j * (log (1 / (25 * ε)) / β₀)) *
      (25 * coefficientDefect 25 ε j) = 1 := by
    have heq : β₀ * 2 ^ j * (log (1 / (25 * ε)) / β₀) =
        (2 ^ j : ℕ) * log (1 / (25 * ε)) := by
      field_simp
      simp [mul_comm]
    rw [heq, exp_nat_mul, exp_log (by positivity)]
    simp only [coefficientDefect, mul_div_cancel₀ _ (by norm_num : (25 : ℝ) ≠ 0),
      ← mul_pow, one_div_mul_cancel (by positivity : (25 : ℝ) * ε ≠ 0), one_pow]
  nlinarith

end CertifiedGaps.Spacetime
