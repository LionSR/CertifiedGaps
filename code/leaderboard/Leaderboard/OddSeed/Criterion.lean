/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.OddSeed.Moment

/-!
# The two-purity criterion at odd seed lengths

The proofs of `CertifiedGaps.PurityCriterion` and `CertifiedGaps.Chain` use the evenness of the
seed length `n₀` only to access the moment inequalities, which hold at even lengths for a merely
self-adjoint transfer operator. When the transfer operator is positive semidefinite
(`HasSpatialMomentBoundsAll`, see `Leaderboard.OddSeed.Moment`), the same inequalities hold at
**every** length, and the whole criterion — the doubling step `doubling_step_all`, the window
interpolation `thermalPurity_ge_of_window_all`, the chain `gapGT_of_chain_all` and its rational
wrapper `gapGT_of_chainCertificate_all` — goes through for an odd seed `n₀` (the final ring
lengths `L` stay even, as the records require).

## Main results

* `doubling_step_all`, `thermalPurity_ge_of_window_all`: one doubling step and the window
  interpolation at any length `n` (no evenness).
* `gapGT_of_seeds_all`, `gapGT_of_chain_all`, `gapGT_of_chainCertificate_all`: the gap bounds from
  seeds at an arbitrary `n₀ > 0`.
-/

open Filter Finset Real

namespace CertifiedGaps

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

section StepsAll

variable (E : (L : ℕ) → ι L → ℝ)

/-- Under the all-lengths moment inequalities, `P_s(n, β) ≤ 1` for every `n ≥ n₀`. -/
theorem spatialPurity_le_one_all {n₀ : ℕ} (hZ : HasSpatialMomentBoundsAll E n₀) {n : ℕ}
    (hn : n₀ ≤ n) {β : ℝ} (hβ : 0 < β) : spatialPurity E n β ≤ 1 := by
  rw [spatialPurity, div_le_one (pow_pos (partitionFunction_pos E _ _) 2)]
  exact (hZ β hβ).sq_le n hn

/-- Lemma 1 applied to the spatial spectrum at any length `n ≥ n₀`. -/
theorem one_sub_spatialPurity_two_mul_le_all {n₀ : ℕ} (hZ : HasSpatialMomentBoundsAll E n₀)
    {n : ℕ} (hn : n₀ ≤ n) {β : ℝ} (hβ : 0 < β) :
    1 - spatialPurity E (2 * n) β ≤ squaringDefect (1 - spatialPurity E n β) :=
  one_sub_div_le_squaringDefect (partitionFunction_pos E n β) (partitionFunction_pos E _ β)
    ((hZ β hβ).squaring n hn)

/-- **One step of the exact doubling recurrence at any length `n`.** The all-lengths variant of
`doubling_step`: no evenness hypothesis on `n`. -/
theorem doubling_step_all {n₀ : ℕ} (hZ : HasSpatialMomentBoundsAll E n₀) {n : ℕ} (hn : n₀ ≤ n)
    {β : ℝ} (hβ : 0 < β) {s t s' t' : ℝ} (hs1 : s ≤ 1) (ht1 : t ≤ 1)
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
  have hsn1 := spatialPurity_le_one_all E hZ hn (by positivity : (0 : ℝ) < 2 * β)
  have hsp : 1 - s' ≤ spatialPurity E (2 * n) (2 * β) := by
    have := one_sub_spatialPurity_two_mul_le_all E hZ hn (by positivity : (0 : ℝ) < 2 * β)
    have := squaringDefect_le_of_le hsn1 hsn hz₁
    linarith
  refine ⟨hsp, ?_⟩
  -- Eq. (10) and `P_s(2n, β) ≤ 1`: `P_τ(4n, β) ≥ (1 - t)² (1 - s')`.
  have hτ4 : 1 - couplingDefect t s' ≤ thermalPurity E (4 * n) β := by
    rw [couplingDefect, sub_sub_cancel, thermalPurity_four_mul,
      le_div_iff₀ (pow_pos (spatialPurity_pos E _ _) 2)]
    have h2n := spatialPurity_le_one_all E hZ (by omega : n₀ ≤ 2 * n) hβ
    have hs2 : spatialPurity E (2 * n) β ^ 2 ≤ 1 :=
      pow_le_one₀ (spatialPurity_pos E _ _).le h2n
    have htsq : (1 - t) ^ 2 ≤ thermalPurity E (2 * n) β ^ 2 := pow_le_pow_left₀ h1t ht 2
    have hprod := mul_le_mul htsq hsp h1s' (by positivity)
    have : 0 ≤ (1 - s') * (1 - t) ^ 2 := by positivity
    nlinarith [mul_le_mul_of_nonneg_left hs2 this]
  have := one_sub_thermalPurity_two_mul_le E (4 * n) β
  have := squaringDefect_le_of_le (thermalPurity_le_one E (4 * n) β) hτ4 hz₂
  linarith

/-- **Interpolation** at any length `n`: the all-lengths variant of
`thermalPurity_ge_of_window`, with no evenness hypothesis on `n` or `L`. -/
theorem thermalPurity_ge_of_window_all {n₀ : ℕ} (hZ : HasSpatialMomentBoundsAll E n₀) {n : ℕ}
    (hn₀n : n₀ ≤ n) (hn0 : 0 < n) {β : ℝ} (hβ : 0 < β) {s t : ℝ} (hs1 : s ≤ 1)
    (hs : 1 - s ≤ spatialPurity E n β) (ht : 1 - t ≤ thermalPurity E (2 * n) β) {L : ℕ}
    (hnL : n ≤ L) (hL2 : L ≤ 2 * n) :
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
    have := (hZ β hβ).le_rpow n L hn₀n hnL
    convert this using 2; simp only [r]; field_simp
  have hDC : D ≤ C ^ r⁻¹ := by
    have := (hZ (2 * β) (by positivity)).le_rpow L (2 * n) (hn₀n.trans hnL) hL2
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
    have h2 := spatialPurity_le_one_all E hZ hn₀n hβ
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

end StepsAll

/-! ### The doubling induction with a coefficient, at any seed length -/

/-- **The doubling step with coefficient `C` at any length `n`.** All-lengths variant of
`doubling_of_coefficient`. -/
theorem doubling_of_coefficient_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ}
    (hZ : HasSpatialMomentBoundsAll E n₀) {n : ℕ} (hn : n₀ ≤ n) {β : ℝ}
    (hβ : 0 < β) {u C : ℝ} (hu0 : 0 ≤ u) (hu1 : u < 1) (hG : doublingCoefficient u ≤ C)
    (hCu : C * u ≤ 1) (hs : 1 - u ≤ spatialPurity E n β)
    (ht : 1 - u ≤ thermalPurity E (2 * n) β) :
    1 - C * u ^ 2 ≤ spatialPurity E (2 * n) (2 * β) ∧
      1 - C * u ^ 2 ≤ thermalPurity E (4 * n) (2 * β) := by
  have hCu2 : C * u ^ 2 ≤ u := by nlinarith
  have hCu20 : 0 ≤ C * u ^ 2 := by
    have : 0 ≤ doublingCoefficient u := by unfold doublingCoefficient; positivity
    nlinarith [sq_nonneg u]
  have hz : couplingDefect u u = 1 - (1 - u) ^ 3 := by unfold couplingDefect; ring
  have hz1 : couplingDefect u u < 1 := by
    rw [hz]; linarith [pow_pos (by linarith : (0 : ℝ) < 1 - u) 3]
  have hf : squaringDefect (couplingDefect u u) ≤ C * u ^ 2 := by
    rw [squaringDefect_couplingDefect_self hu1]
    exact mul_le_mul_of_nonneg_right hG (sq_nonneg u)
  have hsq : 0 ≤ (1 - u) ^ 2 := sq_nonneg _
  have hz2 : couplingDefect u (C * u ^ 2) ≤ couplingDefect u u := by
    unfold couplingDefect; nlinarith [mul_le_mul_of_nonneg_right hCu2 hsq]
  have hz20 : 0 ≤ couplingDefect u (C * u ^ 2) := by
    have hsq1 : (1 - u) ^ 2 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    unfold couplingDefect; nlinarith [mul_nonneg hCu20 hsq]
  exact doubling_step_all E hZ hn hβ hu1.le hu1.le (by linarith) hz1 hf (hz2.trans_lt hz1)
    ((squaringDefect_mono hz20 hz2 hz1).trans hf) hs ht

/-- The induction on scale with coefficient `C` at any seed length: all-lengths variant of
`purities_ge_coefficientDefect`. -/
theorem purities_ge_coefficientDefect_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ}
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBoundsAll E n₀) {u₀ C : ℝ} (hu₀ : 0 < u₀)
    (hu₁ : u₀ < 1) (hG : doublingCoefficient u₀ ≤ C) (hr : C * u₀ ≤ 1)
    (hs : 1 - u₀ ≤ spatialPurity E n₀ β₀) (ht : 1 - u₀ ≤ thermalPurity E (2 * n₀) β₀)
    (j : ℕ) :
    1 - coefficientDefect C u₀ j ≤ spatialPurity E (n₀ * 2 ^ j) (β₀ * 2 ^ j) ∧
      1 - coefficientDefect C u₀ j ≤ thermalPurity E (2 * (n₀ * 2 ^ j)) (β₀ * 2 ^ j) := by
  have hC : 0 < C := lt_of_lt_of_le (doublingCoefficient_pos (by linarith)) hG
  induction j with
  | zero =>
    simp only [coefficientDefect, pow_zero, pow_one, mul_one]
    rw [mul_div_cancel_left₀ _ hC.ne']
    exact ⟨hs, ht⟩
  | succ j ih =>
    obtain ⟨hu0, hu⟩ := coefficientDefect_mem hC hu₀ hr j
    have h := doubling_of_coefficient_all E hZ (n := n₀ * 2 ^ j)
      (Nat.le_mul_of_pos_right _ (by positivity)) (by positivity) hu0.le
      (hu.trans_lt hu₁) ((doublingCoefficient_mono hu hu₁).trans hG)
      ((mul_le_mul_of_nonneg_left hu hC.le).trans hr) ih.1 ih.2
    rw [coefficientDefect_succ hC.ne']
    have e1 : n₀ * 2 ^ (j + 1) = 2 * (n₀ * 2 ^ j) := by ring
    have e2 : 2 * (n₀ * 2 ^ (j + 1)) = 4 * (n₀ * 2 ^ j) := by ring
    have e3 : β₀ * 2 ^ (j + 1) = 2 * (β₀ * 2 ^ j) := by ring
    rw [e1, e3]; rw [e1] at e2; rw [e2]
    exact h

/-- **The two-purity criterion with coefficient `C` at any seed length `n₀ > 0`.** All-lengths
variant of `gapGT_of_coefficient`: no evenness hypothesis on `n₀`. -/
theorem gapGT_of_coefficient_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBoundsAll E n₀) {u₀ C : ℝ} (hu₀ : 0 < u₀)
    (hu₁ : u₀ < 1 / 6) (hG : doublingCoefficient u₀ ≤ C) (hr : C * u₀ < 1)
    (hC3 : 3 < C * (1 - 3 * u₀)) (hs : 1 - u₀ ≤ spatialPurity E n₀ β₀)
    (ht : 1 - u₀ ≤ thermalPurity E (2 * n₀) β₀) {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) :
    GapGT (E L) (log (1 / (C * u₀)) / β₀) := by
  have hC : 0 < C := lt_of_lt_of_le (doublingCoefficient_pos (by linarith)) hG
  obtain ⟨j, hjL, hLj⟩ := exists_dyadic_window hn₀ hL
  obtain ⟨hsj, htj⟩ := purities_ge_coefficientDefect_all E hβ₀ hZ hu₀ (by linarith) hG hr.le
    hs ht j
  set u := coefficientDefect C u₀ j
  obtain ⟨hu0, hu⟩ := coefficientDefect_mem hC hu₀ hr.le j
  have hP := thermalPurity_ge_of_window_all E hZ (n := n₀ * 2 ^ j)
    (Nat.le_mul_of_pos_right _ (by positivity)) (by positivity)
    (by positivity : (0 : ℝ) < β₀ * 2 ^ j) (by linarith) hsj htj hjL hLj
  have hz3 : couplingDefect u u ≤ 3 * u := by
    unfold couplingDefect; nlinarith [sq_nonneg u, mul_nonneg hu0.le (sq_nonneg u)]
  have hz0 : 0 ≤ couplingDefect u u := by
    unfold couplingDefect
    nlinarith [pow_le_one₀ (n := 2) (by linarith : (0 : ℝ) ≤ 1 - u) (by linarith),
      sq_nonneg (1 - u)]
  refine gapGT_of_thermalPurity E L (by positivity) (by linarith) hP ?_
  have hexp : exp (β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀)) * (C * u) = 1 := by
    rw [show β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀) = (2 ^ j : ℕ) * log (1 / (C * u₀)) by
        field_simp; push_cast; ring, exp_nat_mul, exp_log (by positivity)]
    simp only [u, coefficientDefect]
    rw [mul_div_cancel₀ _ hC.ne', ← mul_pow]
    rw [one_div_mul_cancel (by positivity), one_pow]
  have hE := exp_pos (β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀))
  set e := exp (β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀))
  have h1 : e * couplingDefect u u ≤ 3 / C := by
    rw [le_div_iff₀ hC]
    nlinarith [mul_le_mul_of_nonneg_left hz3 hE.le]
  have h2 : 3 / C < 1 - 3 * u := by
    rw [div_lt_iff₀ hC]; nlinarith
  linarith

/-- **Proposition 1 (Two-purity criterion) at any seed length `n₀ > 0`.** All-lengths variant of
`gapGT_of_seeds`: no evenness hypothesis on `n₀`; the ring lengths `L` stay even. -/
theorem gapGT_of_seeds_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBoundsAll E n₀) {ε₀ : ℝ} (hε₀ : 0 < ε₀)
    (hε₁ : ε₀ ≤ 1 / 100) (hs : 1 - ε₀ ≤ spatialPurity E n₀ β₀)
    (ht : 1 - ε₀ ≤ thermalPurity E (2 * n₀) β₀) {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) :
    GapGT (E L) (log (1 / (5 * ε₀)) / β₀) :=
  gapGT_of_coefficient_all E hn₀ hβ₀ hZ hε₀ (by linarith)
    (doublingCoefficient_le_five hε₁) (by linarith) (by linarith) hs ht hL hLe

/-! ### The chain criterion at any seed length -/

/-- The purities along a certified chain at any seed length: all-lengths variant of
`purities_ge_chain`. -/
theorem purities_ge_chain_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ}
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBoundsAll E n₀) (s t : ℕ → ℝ) (J : ℕ)
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
    have h := doubling_step_all E hZ (n := n₀ * 2 ^ j) (Nat.le_mul_of_pos_right _ (by positivity))
      (by positivity) (hrange j (by omega)).1 (hrange j (by omega)).2
      (hrange (j + 1) hj).1 hz₁ hs' hz₂ ht' ihs iht
    have e1 : n₀ * 2 ^ (j + 1) = 2 * (n₀ * 2 ^ j) := by ring
    have e2 : 2 * (n₀ * 2 ^ (j + 1)) = 4 * (n₀ * 2 ^ j) := by ring
    have e3 : β₀ * 2 ^ (j + 1) = 2 * (β₀ * 2 ^ j) := by ring
    rw [e1, e3]; rw [e1] at e2; rw [e2]
    exact h

/-- **The two-purity criterion with a certified chain at any seed length `n₀ > 0`.** All-lengths
variant of `gapGT_of_chain`: no evenness hypothesis on `n₀`; the ring lengths `L` stay even. -/
theorem gapGT_of_chain_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBoundsAll E n₀) (s t : ℕ → ℝ) (J : ℕ)
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
  have hchain := purities_ge_chain_all E hβ₀ hZ s t J hrange hstep hs ht
  obtain ⟨j, hjL, hLj⟩ := exists_dyadic_window hn₀ hL
  have h2j : (0 : ℝ) < 2 ^ j := by positivity
  rcases lt_or_ge j J with hjJ | hjJ
  · -- A scale before the hand-over: interpolation and the purity-to-gap bound.
    obtain ⟨hsj, htj⟩ := hchain j hjJ.le
    have hP := thermalPurity_ge_of_window_all E hZ (n := n₀ * 2 ^ j)
      (Nat.le_mul_of_pos_right _ (by positivity)) (by positivity)
      (by positivity : (0 : ℝ) < β₀ * 2 ^ j) (hrange j hjJ.le).1 hsj htj hjL hLj
    refine gapGT_of_thermalPurity E L (by positivity) (hstep j hjJ).1 hP ?_
    have : exp (β₀ * 2 ^ j * (log r / β₀)) = r ^ 2 ^ j := by
      rw [show β₀ * 2 ^ j * (log r / β₀) = (2 ^ j : ℕ) * log r by field_simp; push_cast; ring,
        exp_nat_mul, exp_log hr]
    rw [this]
    exact hrate j hjJ
  · -- From scale `J` on: `gapGT_of_seeds_all` with seed scale `(n₀ 2^J, β₀ 2^J)`.
    obtain ⟨hsJ', htJ'⟩ := hchain J le_rfl
    have hnJ : n₀ * 2 ^ J ≤ L :=
      (Nat.mul_le_mul_left _ (Nat.pow_le_pow_right two_pos hjJ)).trans hjL
    have h := gapGT_of_seeds_all E (n₀ := n₀ * 2 ^ J) (by positivity)
      (by positivity : (0 : ℝ) < β₀ * 2 ^ J)
      (fun β hβ ↦ (hZ β hβ).mono (Nat.le_mul_of_pos_right _ (by positivity))) hε hε₁
      (by linarith) (by linarith) hnJ hLe
    refine h.mono ?_
    have h2J : (0 : ℝ) < 2 ^ J := by positivity
    have hpow : r ^ 2 ^ J ≤ 1 / (5 * ε) := by
      rw [le_div_iff₀ (by positivity)]; exact htail
    have hlog : (2 ^ J : ℕ) * log r ≤ log (1 / (5 * ε)) := by
      rw [← log_pow]; exact log_le_log (by positivity) hpow
    rw [div_le_div_iff₀ hβ₀ (by positivity)]
    push_cast at hlog
    nlinarith

/-- **The two-purity criterion from a rational chain certificate at any seed length.** All-lengths
variant of `gapGT_of_chainCertificate`: `n₀` may be odd. -/
theorem gapGT_of_chainCertificate_all (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBoundsAll E n₀)
    (c : ChainCertificate) (hc : c.Valid)
    (hs : 1 - (c.sAt 0 : ℝ) ≤ spatialPurity E n₀ β₀)
    (ht : 1 - (c.tAt 0 : ℝ) ≤ thermalPurity E (2 * n₀) β₀)
    {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) : GapGT (E L) (log (c.r : ℝ) / β₀) := by
  obtain ⟨hr, hε, hε₁, hsJ, htJ, htail, hrange, hstep⟩ := hc
  open ChainCertificate in
  refine gapGT_of_chain_all E hn₀ hβ₀ hZ (fun j ↦ c.sAt j) (fun j ↦ c.tAt j) c.J
    (r := c.r) (ε := c.ε) (by exact_mod_cast hr) (fun j hj ↦ ?_) (fun j hj ↦ ?_) (fun j hj ↦ ?_)
    (by exact_mod_cast hε)
    (by rw [show (1 / 100 : ℝ) = ((1 / 100 : ℚ) : ℝ) by norm_num]; exact_mod_cast hε₁)
    (by exact_mod_cast hsJ)
    (by exact_mod_cast htJ) (by exact_mod_cast htail) hs ht hL hLe
  · obtain ⟨h1, h2⟩ := hrange j (by omega)
    exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩
  · obtain ⟨h1, h2, h3, h4, -⟩ := hstep j hj
    simp only [cast_zQ, cast_fQ]
    exact ⟨by exact_mod_cast h1, by exact_mod_cast h2, by exact_mod_cast h3,
      by exact_mod_cast h4⟩
  · obtain ⟨-, -, -, -, h5⟩ := hstep j hj
    simp only [cast_zQ]
    exact_mod_cast h5

end CertifiedGaps
