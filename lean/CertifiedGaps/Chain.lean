/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.PurityCriterion

/-!
# The two-purity criterion with unequal seed defects

`gapGT_of_seeds` needs both seed defects below `1/100`. Certified seeds are usually larger, for
example `(εₛ, εₜ) ≈ (0.064, 0.041)`, and reach that range only after a few doubling steps. This
file runs the exact doubling recurrence of the thermal-MPO note for finitely many steps and then
hands over to `gapGT_of_seeds`.

With `z(a, b) = 1 - (1 - b)(1 - a)²` (`couplingDefect`) and the squaring bound `f`
(`squaringDefect`), one doubling step maps upper bounds `(s, t)` on the defects of
`P_s(n, β), P_τ(2n, β)` to upper bounds `s' ≥ f(z(s, t))`, `t' ≥ f(z(t, s'))` on those of
`P_s(2n, 2β), P_τ(4n, 2β)` (`doubling_step`). At each scale, every even `L ∈ [n, 2n]` has
`P_τ(L, β) ≥ 1 - z(s, t)` (`thermalPurity_ge_of_window`), which bounds the gap
(`gapGT_of_thermalPurity`).

The conclusion is stated as `Δ_L > log r / β₀` for a number `r > 0`. Every remaining condition is
an inequality between rational numbers once the chain `(s j, t j)` and `r` are rational: the
recurrence steps, `r^{2^j} z_j < 1 - z_j` at the scales `j < J`, and `r^{2^J} · 5ε ≤ 1` at the
scale `J` where `gapGT_of_seeds` takes over.

## Main results

* `doubling_step`: one step of the exact recurrence, for any defects with `z < 1`.
* `gapGT_of_chain`: the uniform gap bound `log r / β₀` from a finite certified chain.
-/

open Filter Finset Real

namespace ThermalPurityGap

/-- The coupling defect `z(a, b) = 1 - (1 - b)(1 - a)²`: if `P_s ≥ 1 - a` and `P_τ ≥ 1 - b` then
the coupled purity is at least `1 - z(a, b)`. -/
def couplingDefect (a b : ℝ) : ℝ := 1 - (1 - b) * (1 - a) ^ 2

/-- A gap bound is weakened by lowering the constant. -/
theorem GapGT.mono {κ : Type*} {e : κ → ℝ} {c c' : ℝ} (h : GapGT e c) (hc : c' ≤ c) :
    GapGT e c' := by
  obtain ⟨m₀, hm₀⟩ := h
  exact ⟨m₀, fun m hm ↦ by linarith [hm₀ m hm]⟩

/-- `f(1 - P) ≤ f(z)` from `0 ≤ 1 - P ≤ z < 1`, as an upper bound on a defect. -/
theorem squaringDefect_le_of_le {P z : ℝ} (hP : P ≤ 1) (hPz : 1 - z ≤ P) (hz : z < 1) :
    squaringDefect (1 - P) ≤ squaringDefect z :=
  squaringDefect_mono (by linarith) (by linarith) hz

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

section Steps

variable (E : (L : ℕ) → ι L → ℝ)

/-- Lemma 1 applied to the Gibbs spectrum: `1 - P_τ(L, 2β) ≤ f(1 - P_τ(L, β))`. -/
theorem one_sub_thermalPurity_two_mul_le (L : ℕ) (β : ℝ) :
    1 - thermalPurity E L (2 * β) ≤ squaringDefect (1 - thermalPurity E L β) := by
  have h₁ := hasSum_fintype fun m ↦ exp (-(β * E L m))
  have h₂ := hasSum_fintype fun m ↦ exp (-(β * E L m)) ^ 2
  have h₄ := hasSum_fintype fun m ↦ exp (-(β * E L m)) ^ 4
  rw [← partitionFunction_two_mul] at h₂
  rw [← partitionFunction_four_mul] at h₄
  obtain ⟨-, hmain⟩ := sq_sub_sum_sq_of_hasSum (fun m ↦ (exp_pos _).le) h₁ h₂ h₄
  exact one_sub_div_le_squaringDefect (partitionFunction_pos E L β)
    (partitionFunction_pos E L (2 * β)) hmain

/-- Lemma 1 applied to the spatial spectrum: `1 - P_s(2n, β) ≤ f(1 - P_s(n, β))`. -/
theorem one_sub_spatialPurity_two_mul_le {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ}
    (hn : n₀ ≤ n) (hne : Even n) {β : ℝ} (hβ : 0 < β) :
    1 - spatialPurity E (2 * n) β ≤ squaringDefect (1 - spatialPurity E n β) :=
  one_sub_div_le_squaringDefect (partitionFunction_pos E n β) (partitionFunction_pos E _ β)
    ((hZ β hβ).squaring n hn hne)

/-- **One step of the exact doubling recurrence.** If `P_s(n, β) ≥ 1 - s` and
`P_τ(2n, β) ≥ 1 - t`, then `P_s(2n, 2β) ≥ 1 - s'` and `P_τ(4n, 2β) ≥ 1 - t'` for any
`s' ≥ f(z(s, t))` and `t' ≥ f(z(t, s'))`, provided both arguments of `f` are below `1`. -/
theorem doubling_step {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ} (hn : n₀ ≤ n)
    (hne : Even n) {β : ℝ} (hβ : 0 < β) {s t s' t' : ℝ} (hs1 : s ≤ 1) (ht1 : t ≤ 1)
    (hs'1 : s' ≤ 1) (hz₁ : couplingDefect s t < 1)
    (hs' : squaringDefect (couplingDefect s t) ≤ s') (hz₂ : couplingDefect t s' < 1)
    (ht' : squaringDefect (couplingDefect t s') ≤ t') (hs : 1 - s ≤ spatialPurity E n β)
    (ht : 1 - t ≤ thermalPurity E (2 * n) β) :
    1 - s' ≤ spatialPurity E (2 * n) (2 * β) ∧ 1 - t' ≤ thermalPurity E (4 * n) (2 * β) := by
  have h1s : 0 ≤ 1 - s := by linarith
  have h1t : 0 ≤ 1 - t := by linarith
  have h1s' : 0 ≤ 1 - s' := by linarith
  -- Eq. (9) and `P_τ(n, β) ≤ 1`: `P_s(n, 2β) ≥ (1 - s)² (1 - t)`.
  have hsn : 1 - couplingDefect s t ≤ spatialPurity E n (2 * β) := by
    have hτn := thermalPurity_le_one E n β
    have hτn0 := thermalPurity_pos E n β
    rw [couplingDefect, sub_sub_cancel, spatialPurity_two_mul, le_div_iff₀ (by positivity)]
    have hsq : (1 - s) ^ 2 ≤ spatialPurity E n β ^ 2 := pow_le_pow_left₀ h1s hs 2
    have hτsq : thermalPurity E n β ^ 2 ≤ 1 := pow_le_one₀ hτn0.le hτn
    have hprod := mul_le_mul ht hsq (by positivity) (thermalPurity_pos E (2 * n) β).le
    have : 0 ≤ (1 - t) * (1 - s) ^ 2 := by positivity
    nlinarith [mul_le_mul_of_nonneg_left hτsq this]
  have hsn1 := spatialPurity_le_one E hZ hn hne (by positivity : (0 : ℝ) < 2 * β)
  have hsp : 1 - s' ≤ spatialPurity E (2 * n) (2 * β) := by
    have := one_sub_spatialPurity_two_mul_le E hZ hn hne (by positivity : (0 : ℝ) < 2 * β)
    have := squaringDefect_le_of_le hsn1 hsn hz₁
    linarith
  refine ⟨hsp, ?_⟩
  -- Eq. (10) and `P_s(2n, β) ≤ 1`: `P_τ(4n, β) ≥ (1 - t)² (1 - s')`.
  have hτ4 : 1 - couplingDefect t s' ≤ thermalPurity E (4 * n) β := by
    rw [couplingDefect, sub_sub_cancel, thermalPurity_four_mul,
      le_div_iff₀ (pow_pos (spatialPurity_pos E _ _) 2)]
    have h2n := spatialPurity_le_one E hZ (by omega : n₀ ≤ 2 * n) (by simp) hβ
    have hs2 : spatialPurity E (2 * n) β ^ 2 ≤ 1 :=
      pow_le_one₀ (spatialPurity_pos E _ _).le h2n
    have htsq : (1 - t) ^ 2 ≤ thermalPurity E (2 * n) β ^ 2 := pow_le_pow_left₀ h1t ht 2
    have hprod := mul_le_mul htsq hsp h1s' (by positivity)
    have : 0 ≤ (1 - s') * (1 - t) ^ 2 := by positivity
    nlinarith [mul_le_mul_of_nonneg_left hs2 this]
  have := one_sub_thermalPurity_two_mul_le E (4 * n) β
  have := squaringDefect_le_of_le (thermalPurity_le_one E (4 * n) β) hτ4 hz₂
  linarith

/-- **Interpolation**, eqs. (14) and (15) of the thermal-MPO note, for unequal defects: if
`P_s(n, β) ≥ 1 - s` and `P_τ(2n, β) ≥ 1 - t`, then every even `L` with `n ≤ L ≤ 2n` has
`P_τ(L, β) ≥ 1 - z(s, t)`. -/
theorem thermalPurity_ge_of_window {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ}
    (hn₀n : n₀ ≤ n) (hne : Even n) (hn0 : 0 < n) {β : ℝ} (hβ : 0 < β) {s t : ℝ} (hs1 : s ≤ 1)
    (hs : 1 - s ≤ spatialPurity E n β) (ht : 1 - t ≤ thermalPurity E (2 * n) β) {L : ℕ}
    (hnL : n ≤ L) (hL2 : L ≤ 2 * n) (hLe : Even L) :
    1 - couplingDefect s t ≤ thermalPurity E L β := by
  set A := partitionFunction E L β
  set B := partitionFunction E n β
  set C := partitionFunction E L (2 * β)
  set D := partitionFunction E (2 * n) (2 * β)
  have hA := partitionFunction_pos E L β
  have hB := partitionFunction_pos E n β
  have hC := partitionFunction_pos E L (2 * β)
  have hD := partitionFunction_pos E (2 * n) (2 * β)
  set r : ℝ := (L : ℝ) / (2 * n)
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn0
  have hLpos : (0 : ℝ) < L := by exact_mod_cast lt_of_lt_of_le hn0 hnL
  have hr0 : 0 < r := by positivity
  have hr1 : r ≤ 1 := (div_le_one (by positivity)).2 (by exact_mod_cast hL2)
  have hAB : A ≤ B ^ (2 * r) := by
    have := (hZ β hβ).le_rpow n L hn₀n hnL hne hLe
    convert this using 2; simp only [r]; field_simp
  have hDC : D ≤ C ^ r⁻¹ := by
    have := (hZ (2 * β) (by positivity)).le_rpow L (2 * n) (hn₀n.trans hnL) hL2 hLe (by simp)
    convert this using 2; simp [r]
  have hDr : D ^ r ≤ C := by
    calc D ^ r ≤ (C ^ r⁻¹) ^ r := Real.rpow_le_rpow hD.le hDC hr0.le
      _ = C := Real.rpow_inv_rpow hC.le hr0.ne'
  have hA2 : A ^ 2 ≤ (B ^ 4) ^ r := by
    calc A ^ 2 ≤ (B ^ (2 * r)) ^ 2 := pow_le_pow_left₀ hA.le hAB 2
      _ = (B ^ 4) ^ r := by
        rw [← Real.rpow_mul_natCast hB.le, ← Real.rpow_natCast_mul hB.le]; ring_nf
  set x := D / B ^ 4
  have hx : x = thermalPurity E (2 * n) β * spatialPurity E n β ^ 2 := by
    have := partitionFunction_pos E (2 * n) β
    simp only [x, D, B, thermalPurity, spatialPurity]
    field_simp
  have hx0 : 0 < x := by positivity
  have hx1 : x ≤ 1 := by
    rw [hx]
    have h1 := thermalPurity_le_one E (2 * n) β
    have h2 := spatialPurity_le_one E hZ hn₀n hne hβ
    have h3 := thermalPurity_pos E (2 * n) β
    nlinarith [pow_le_one₀ (n := 2) (spatialPurity_pos E n β).le h2]
  have hxz : 1 - couplingDefect s t ≤ x := by
    rw [hx, couplingDefect, sub_sub_cancel]
    have h1s : 0 ≤ 1 - s := by linarith
    have : (1 - s) ^ 2 ≤ spatialPurity E n β ^ 2 := pow_le_pow_left₀ h1s hs 2
    rcases le_or_gt 0 (1 - t) with h1t | h1t
    · nlinarith [mul_le_mul ht this (by positivity) (thermalPurity_pos E (2 * n) β).le]
    · have : 0 ≤ thermalPurity E (2 * n) β * spatialPurity E n β ^ 2 := by
        have := thermalPurity_pos E (2 * n) β; positivity
      nlinarith [sq_nonneg (1 - s)]
  calc 1 - couplingDefect s t ≤ x := hxz
    _ = x ^ (1 : ℝ) := (Real.rpow_one x).symm
    _ ≤ x ^ r := Real.rpow_le_rpow_of_exponent_ge hx0 hx1 hr1
    _ = D ^ r / (B ^ 4) ^ r := Real.div_rpow hD.le (by positivity) r
    _ ≤ C / A ^ 2 := div_le_div₀ hC.le hDr (by positivity) hA2
    _ = thermalPurity E L β := rfl

/-- **From a thermal purity to a gap**, eq. (4) of the thermal-MPO note: if
`P_τ(L, β) ≥ 1 - z` with `z < 1` and `e^{βc} z < 1 - z`, then `Δ_L > c`. -/
theorem gapGT_of_thermalPurity (L : ℕ) {β : ℝ} (hβ : 0 < β) {z c : ℝ} (hz : z < 1)
    (hP : 1 - z ≤ thermalPurity E L β) (hc : exp (β * c) * z < 1 - z) : GapGT (E L) c := by
  obtain ⟨m₀, -, hm₀⟩ := exists_min_image univ (E L) univ_nonempty
  refine ⟨m₀, fun m hm ↦ ?_⟩
  have h4 := exp_neg_mul_sub_le_of_isMin E L hβ.le (fun m' ↦ hm₀ m' (mem_univ _)) hm
  have hPpos := thermalPurity_pos E L β
  have h1z : 0 < 1 - z := by linarith
  -- `(1 - P)/P ≤ z/(1 - z) < e^{-βc}`.
  have hle : (1 - thermalPurity E L β) / thermalPurity E L β ≤ z / (1 - z) := by
    rw [div_le_div_iff₀ hPpos h1z]; nlinarith
  have hlt : z / (1 - z) < exp (-(β * c)) := by
    rw [div_lt_iff₀ h1z, exp_neg, ← div_eq_inv_mul, lt_div_iff₀ (exp_pos _)]
    linarith [mul_comm z (exp (β * c))]
  have := lt_of_le_of_lt (h4.trans hle) hlt
  rw [exp_lt_exp] at this
  nlinarith

end Steps

/-! ### The chain criterion -/

/-- The purities along a certified chain: `P_s(n₀ 2^j, β₀ 2^j) ≥ 1 - s j` and
`P_τ(2 n₀ 2^j, β₀ 2^j) ≥ 1 - t j` for `j ≤ J`. -/
theorem purities_ge_chain (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀e : Even n₀) {β₀ : ℝ}
    (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀) (s t : ℕ → ℝ) (J : ℕ)
    (hrange : ∀ j ≤ J, s j ≤ 1 ∧ t j ≤ 1)
    (hstep : ∀ j < J, couplingDefect (s j) (t j) < 1 ∧
      squaringDefect (couplingDefect (s j) (t j)) ≤ s (j + 1) ∧
      couplingDefect (t j) (s (j + 1)) < 1 ∧
      squaringDefect (couplingDefect (t j) (s (j + 1))) ≤ t (j + 1))
    (hs : 1 - s 0 ≤ spatialPurity E n₀ β₀) (ht : 1 - t 0 ≤ thermalPurity E (2 * n₀) β₀) :
    ∀ j ≤ J, 1 - s j ≤ spatialPurity E (n₀ * 2 ^ j) (β₀ * 2 ^ j) ∧
      1 - t j ≤ thermalPurity E (2 * (n₀ * 2 ^ j)) (β₀ * 2 ^ j) := by
  intro j hj
  induction j with
  | zero => simp only [pow_zero, mul_one]; exact ⟨hs, ht⟩
  | succ j ih =>
    obtain ⟨ihs, iht⟩ := ih (by omega)
    obtain ⟨hz₁, hs', hz₂, ht'⟩ := hstep j (by omega)
    have h := doubling_step E hZ (n := n₀ * 2 ^ j) (Nat.le_mul_of_pos_right _ (by positivity))
      (hn₀e.mul_right _) (by positivity) (hrange j (by omega)).1 (hrange j (by omega)).2
      (hrange (j + 1) hj).1 hz₁ hs' hz₂ ht' ihs iht
    have e1 : n₀ * 2 ^ (j + 1) = 2 * (n₀ * 2 ^ j) := by ring
    have e2 : 2 * (n₀ * 2 ^ (j + 1)) = 4 * (n₀ * 2 ^ j) := by ring
    have e3 : β₀ * 2 ^ (j + 1) = 2 * (β₀ * 2 ^ j) := by ring
    rw [e1, e3]; rw [e1] at e2; rw [e2]
    exact h

/-- **The two-purity criterion with a certified chain.** Assume the moment inequalities above an
even `n₀ > 0`, seeds `P_s(n₀, β₀) ≥ 1 - s 0` and `P_τ(2n₀, β₀) ≥ 1 - t 0`, and a chain
`(s j, t j)`, `j ≤ J`, of upper bounds for the exact doubling recurrence. If
`r^{2^j} z(s j, t j) < 1 - z(s j, t j)` for `j < J`, and the defects at scale `J` are at most
some `ε ≤ 1/100` with `r^{2^J} · 5ε ≤ 1`, then every even `L ≥ n₀` has a nondegenerate ground
level and gap `Δ_L > log r / β₀`. -/
theorem gapGT_of_chain (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀) (hn₀e : Even n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀) (s t : ℕ → ℝ) (J : ℕ)
    {r ε : ℝ} (hr : 0 < r) (hrange : ∀ j ≤ J, s j ≤ 1 ∧ t j ≤ 1)
    (hstep : ∀ j < J, couplingDefect (s j) (t j) < 1 ∧
      squaringDefect (couplingDefect (s j) (t j)) ≤ s (j + 1) ∧
      couplingDefect (t j) (s (j + 1)) < 1 ∧
      squaringDefect (couplingDefect (t j) (s (j + 1))) ≤ t (j + 1))
    (hrate : ∀ j < J, r ^ 2 ^ j * couplingDefect (s j) (t j) < 1 - couplingDefect (s j) (t j))
    (hε : 0 < ε) (hε₁ : ε ≤ 1 / 100) (hsJ : s J ≤ ε) (htJ : t J ≤ ε)
    (htail : r ^ 2 ^ J * (5 * ε) ≤ 1)
    (hs : 1 - s 0 ≤ spatialPurity E n₀ β₀) (ht : 1 - t 0 ≤ thermalPurity E (2 * n₀) β₀)
    {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) : GapGT (E L) (log r / β₀) := by
  have hchain := purities_ge_chain E hn₀e hβ₀ hZ s t J hrange hstep hs ht
  obtain ⟨j, hjL, hLj⟩ := exists_dyadic_window hn₀ hL
  have h2j : (0 : ℝ) < 2 ^ j := by positivity
  rcases lt_or_ge j J with hjJ | hjJ
  · -- A scale before the hand-over: interpolation and the purity-to-gap bound.
    obtain ⟨hsj, htj⟩ := hchain j hjJ.le
    have hP := thermalPurity_ge_of_window E hZ (n := n₀ * 2 ^ j)
      (Nat.le_mul_of_pos_right _ (by positivity)) (hn₀e.mul_right _) (by positivity)
      (by positivity : (0 : ℝ) < β₀ * 2 ^ j) (hrange j hjJ.le).1 hsj htj hjL hLj hLe
    refine gapGT_of_thermalPurity E L (by positivity) (hstep j hjJ).1 hP ?_
    have : exp (β₀ * 2 ^ j * (log r / β₀)) = r ^ 2 ^ j := by
      rw [show β₀ * 2 ^ j * (log r / β₀) = (2 ^ j : ℕ) * log r by field_simp; push_cast; ring,
        exp_nat_mul, exp_log hr]
    rw [this]
    exact hrate j hjJ
  · -- From scale `J` on: `gapGT_of_seeds` with seed scale `(n₀ 2^J, β₀ 2^J)`.
    obtain ⟨hsJ', htJ'⟩ := hchain J le_rfl
    have hnJ : n₀ * 2 ^ J ≤ L :=
      (Nat.mul_le_mul_left _ (Nat.pow_le_pow_right two_pos hjJ)).trans hjL
    have h := gapGT_of_seeds E (n₀ := n₀ * 2 ^ J) (by positivity) (hn₀e.mul_right _)
      (by positivity : (0 : ℝ) < β₀ * 2 ^ J)
      (fun β hβ ↦ (hZ β hβ).mono (Nat.le_mul_of_pos_right _ (by positivity))) hε hε₁
      (by linarith) (by linarith) hnJ hLe
    refine h.mono ?_
    -- `log r / β₀ ≤ log (1/(5ε)) / (β₀ 2^J)` from `r^{2^J} ≤ 1/(5ε)`.
    have h2J : (0 : ℝ) < 2 ^ J := by positivity
    have hpow : r ^ 2 ^ J ≤ 1 / (5 * ε) := by
      rw [le_div_iff₀ (by positivity)]; exact htail
    have hlog : (2 ^ J : ℕ) * log r ≤ log (1 / (5 * ε)) := by
      rw [← log_pow]; exact log_le_log (by positivity) hpow
    rw [div_le_div_iff₀ hβ₀ (by positivity)]
    push_cast at hlog
    nlinarith

end ThermalPurityGap
