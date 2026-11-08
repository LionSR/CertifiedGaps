/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.CubeCriterion
import CertifiedGaps.Chain

/-!
# Three-direction bootstrap with unequal defects

A finite chain bounds the staggered defects in the two spatial directions and imaginary
time. Each directional update is checked separately. The finite interpolation windows and
the quadratic tail give a uniform spectral gap. Transfer inequalities and the initial
spectral seeds are explicit inputs, as in the two-purity chain criterion.
-/

open Real

namespace CertifiedGaps.Spacetime

/-- Lower bound on the cube purity from three possibly unequal defects. -/
def cubeProduct (a b c : ℝ) : ℝ := (1 - a) ^ 4 * (1 - b) ^ 2 * (1 - c)

/-- Staggered purity bounds along the path through a spacetime cube. -/
def SeedBounds (D : CubeData) (a b c : ℝ) (x y t : ℕ) : Prop :=
  1 - a ≤ px D.Z x y t ∧ 1 - b ≤ py D.Z (x + 1) y t ∧
    1 - c ≤ pt D.Z (x + 1) (y + 1) t

/-- A complete cycle of three certified directional updates. -/
def CubeStep (a b c a' b' c' : ℝ) : Prop :=
  0 < cubeProduct a b c ∧ squaringDefect (1 - cubeProduct a b c) ≤ a' ∧
  0 < cubeProduct b c a' ∧ squaringDefect (1 - cubeProduct b c a') ≤ b' ∧
  0 < cubeProduct c a' b' ∧ squaringDefect (1 - cubeProduct c a' b') ≤ c'

/-- Multiply the three staggered purity bounds with their cube exponents. -/
theorem cube_lower_unequal (D : CubeData) {a b c : ℝ}
    (ha : a ≤ 1) (hb : b ≤ 1) (hc : c ≤ 1) {x y t : ℕ}
    (h : SeedBounds D a b c x y t) : cubeProduct a b c ≤ cube D x y t := by
  exact mul_le_mul (mul_le_mul (pow_le_pow_left₀ (sub_nonneg.mpr ha) h.1 4)
    (pow_le_pow_left₀ (sub_nonneg.mpr hb) h.2.1 2) (sq_nonneg _)
    (pow_nonneg (px_pos D _ _ _).le _)) h.2.2 (sub_nonneg.mpr hc)
    (mul_nonneg (pow_nonneg (px_pos D _ _ _).le _) (sq_nonneg _))

/-- One directional update, before cyclic permutation. -/
theorem advance_x_unequal (D : CubeData) {a b c a' : ℝ}
    (ha : a ≤ 1) (hb : b ≤ 1) (hc : c ≤ 1)
    (hq : 0 < cubeProduct a b c) (hs : squaringDefect (1 - cubeProduct a b c) ≤ a')
    {x y t : ℕ} (h : SeedBounds D a b c x y t) :
    1 - a' ≤ px D.Z (x + 1) (y + 1) (t + 1) := by
  have hq' := (cube_lower_unequal D ha hb hc h).trans (cube_le_x_last D x y t)
  have hnext := (D.x_square x (y + 1) (t + 1)).trans
    ((squaringDefect_mono (sub_nonneg.mpr (D.x_le _ _ _))
      (sub_le_sub_left hq' 1) (sub_lt_self 1 hq)).trans hs)
  linarith

/-- Three cyclic directional updates advance all staggered bounds by one scale. -/
theorem cube_step (D : CubeData) {a b c a' b' c' : ℝ}
    (ha : a ≤ 1) (hb : b ≤ 1) (hc : c ≤ 1) (ha' : a' ≤ 1) (hb' : b' ≤ 1)
    (hs : CubeStep a b c a' b' c') {x y t : ℕ} (h : SeedBounds D a b c x y t) :
    SeedBounds D a' b' c' (x + 1) (y + 1) (t + 1) := by
  have hx := advance_x_unequal D ha hb hc hs.1 hs.2.1 h
  have hy := advance_x_unequal D.rotate hb hc ha' hs.2.2.1 hs.2.2.2.1
    (x := y) (y := t) (t := x + 1) ⟨h.2.1, h.2.2, hx⟩
  have ht := advance_x_unequal D.rotate.rotate hc ha' hb' hs.2.2.2.2.1 hs.2.2.2.2.2
    (x := t) (y := x + 1) (t := y + 1) ⟨h.2.2, hx, hy⟩
  exact ⟨hx, hy, ht⟩

/-- The purity bounds obtained by replaying a finite chain. -/
theorem cube_purities_ge_chain (D : CubeData) (a b c : ℕ → ℝ) (J : ℕ)
    (hrange : ∀ j ≤ J, a j ≤ 1 ∧ b j ≤ 1 ∧ c j ≤ 1)
    (hstep : ∀ j < J, CubeStep (a j) (b j) (c j) (a (j + 1)) (b (j + 1)) (c (j + 1)))
    (hseed : SeedBounds D (a 0) (b 0) (c 0) 0 0 0) :
    ∀ j ≤ J, SeedBounds D (a j) (b j) (c j) j j j := by
  intro j hj
  induction j with
  | zero => exact hseed
  | succ j ih =>
    exact cube_step D (hrange j (by omega)).1 (hrange j (by omega)).2.1
      (hrange j (by omega)).2.2 (hrange (j + 1) hj).1 (hrange (j + 1) hj).2.1
      (hstep j (by omega)) (ih (by omega))

/-- Begin the same scalar moment data at a later diagonal scale. -/
def CubeData.shift (D : CubeData) (J : ℕ) : CubeData where
  Z := fun x y t => D.Z (x + J) (y + J) (t + J)
  positive := fun x y t => D.positive (x + J) (y + J) (t + J)
  x_le := fun x y t => by
    simpa [px, Nat.add_right_comm] using D.x_le (x + J) (y + J) (t + J)
  y_le := fun x y t => by
    simpa [py, Nat.add_right_comm] using D.y_le (x + J) (y + J) (t + J)
  t_le := fun x y t => by
    simpa [pt, Nat.add_right_comm] using D.t_le (x + J) (y + J) (t + J)
  x_square := fun x y t => by
    simpa [px, Nat.add_right_comm] using D.x_square (x + J) (y + J) (t + J)
  y_square := fun x y t => by
    simpa [py, Nat.add_right_comm] using D.y_square (x + J) (y + J) (t + J)
  t_square := fun x y t => by
    simpa [pt, Nat.add_right_comm] using D.t_square (x + J) (y + J) (t + J)

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

omit [∀ L, Nonempty (ι L)] in
/-- Spatial window interpolation is preserved when the initial scale is shifted. -/
theorem window_shift (D : CubeData) (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} {β₀ : ℝ}
    (hW : WindowInput D E n₀ β₀) (J : ℕ) :
    WindowInput (D.shift J) E (n₀ * 2 ^ J) (β₀ * 2 ^ J) := by
  intro j L hlow hhigh heven
  have h := hW (j + J) L
    (by simpa [pow_add, mul_assoc, mul_comm, mul_left_comm] using hlow)
    (by simpa [pow_add, mul_assoc, mul_comm, mul_left_comm] using hhigh) heven
  simpa [cube, CubeData.shift, px, py, pt, Nat.add_right_comm, pow_add,
    mul_assoc, mul_comm, mul_left_comm] using h

/-- The three-direction criterion with a finite certified chain. All conditions on the chain
and the rate are rational inequalities when the supplied numbers are rational. -/
theorem gapGT_of_cube_chain (D : CubeData) (E : (L : ℕ) → ι L → ℝ)
    {n₀ : ℕ} (hn₀ : 0 < n₀) {β₀ : ℝ} (hβ₀ : 0 < β₀)
    (hW : WindowInput D E n₀ β₀) (a b c : ℕ → ℝ) (J : ℕ) {r ε : ℝ}
    (hr : 0 < r) (hrange : ∀ j ≤ J, a j ≤ 1 ∧ b j ≤ 1 ∧ c j ≤ 1)
    (hstep : ∀ j < J, CubeStep (a j) (b j) (c j) (a (j + 1)) (b (j + 1)) (c (j + 1)))
    (hrate : ∀ j < J, r ^ 2 ^ j * (1 - cubeProduct (a j) (b j) (c j)) <
      cubeProduct (a j) (b j) (c j))
    (hε : 0 < ε) (hε' : ε ≤ 1 / 1000)
    (hJ : a J ≤ ε ∧ b J ≤ ε ∧ c J ≤ ε) (htail : r ^ 2 ^ J * (25 * ε) ≤ 1)
    (hseed : SeedBounds D (a 0) (b 0) (c 0) 0 0 0)
    {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) : GapGT (E L) (log r / β₀) := by
  have hchain := cube_purities_ge_chain D a b c J hrange hstep hseed
  obtain ⟨j, hjL, hLj⟩ := exists_dyadic_window hn₀ hL
  rcases lt_or_ge j J with hjJ | hjJ
  · have hP := (cube_lower_unequal D (hrange j hjJ.le).1
      (hrange j hjJ.le).2.1 (hrange j hjJ.le).2.2 (hchain j hjJ.le)).trans
        (hW j L hjL hLj hLe)
    refine gapGT_of_thermalPurity E L (β := β₀ * 2 ^ j) (by positivity)
      (sub_lt_self 1 (hstep j hjJ).1) (by simpa only [sub_sub_cancel] using hP) ?_
    have heq : β₀ * 2 ^ j * (log r / β₀) = (2 ^ j : ℕ) * log r := by
      field_simp
      simp [mul_comm]
    simpa only [heq, exp_nat_mul, exp_log hr, sub_sub_cancel] using hrate j hjJ
  · have hnJ : n₀ * 2 ^ J ≤ L :=
      (Nat.mul_le_mul_left _ (Nat.pow_le_pow_right two_pos hjJ)).trans hjL
    have hseedJ : SeedBound (D.shift J) ε 0 0 0 := by
      simpa [SeedBound, CubeData.shift, px, py, pt, Nat.add_comm] using
        (show SeedBound D ε J J J from
          ⟨(sub_le_sub_left hJ.1 1).trans (hchain J le_rfl).1,
            (sub_le_sub_left hJ.2.1 1).trans (hchain J le_rfl).2.1,
            (sub_le_sub_left hJ.2.2 1).trans (hchain J le_rfl).2.2⟩)
    have h := gapGT_of_cube_seeds (D.shift J) E (n₀ := n₀ * 2 ^ J)
      (by positivity) (by positivity : (0 : ℝ) < β₀ * 2 ^ J) hε hε'
      (window_shift D E hW J) hseedJ hnJ hLe
    refine h.mono ?_
    have hpow : r ^ 2 ^ J ≤ 1 / (25 * ε) :=
      (le_div_iff₀ (by positivity)).mpr htail
    have hlog : (2 ^ J : ℝ) * log r ≤ log (1 / (25 * ε)) := by
      simpa only [log_pow, Nat.cast_pow, Nat.cast_ofNat] using
        log_le_log (pow_pos hr (2 ^ J)) hpow
    rw [div_le_div_iff₀ hβ₀ (by positivity)]
    nlinarith

end CertifiedGaps.Spacetime
