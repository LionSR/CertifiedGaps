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
(`gapGT_of_thermalPurity`). These three steps are proved in `CertifiedGaps.PurityCriterion`, where
the equal-defect criterion uses them too.

The conclusion is stated as `Δ_L > log r / β₀` for a number `r > 0`. Every remaining condition is
an inequality between rational numbers once the chain `(s j, t j)` and `r` are rational: the
recurrence steps, `r^{2^j} z_j < 1 - z_j` at the scales `j < J`, and `r^{2^J} · 5ε ≤ 1` at the
scale `J` where `gapGT_of_seeds` takes over.

## Main results

* `purities_ge_chain`: the purity bounds along a certified chain.
* `gapGT_of_chain`: the uniform gap bound `log r / β₀` from a finite certified chain.
-/

open Filter Finset Real

namespace CertifiedGaps

/-- A gap bound is weakened by lowering the constant. -/
theorem GapGT.mono {κ : Type*} {e : κ → ℝ} {c c' : ℝ} (h : GapGT e c) (hc : c' ≤ c) :
    GapGT e c' := by
  obtain ⟨m₀, hm₀⟩ := h
  exact ⟨m₀, fun m hm ↦ by linarith [hm₀ m hm]⟩

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

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

end CertifiedGaps
