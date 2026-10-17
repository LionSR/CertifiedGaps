/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing

/-!
# Records from a chain and a bridge seed

A chain certificate from the seed `(n₀, β)` covers the even `L ≥ n₀`. A second seed
`(n₁, β₁)` with `n₁ ≤ n₀ ≤ 2n₁` covers the even `n₁ ≤ L < n₀` in one step: by interpolation
(`thermalPurity_ge_of_window`), `P_τ(L, β₁) ≥ 1 - z(s₁, t₁)`, and a thermal purity bounds the gap
(`gapGT_of_thermalPurity`). For the bound `c = log r / β` this step needs
`e^{β₁ c} z < 1 - z`, that is `r^{β₁/β} z < 1 - z`. With `β₁ / β = a / b` it is the rational
inequality `r^a z^b < (1 - z)^b`, which `decide` checks.
-/

namespace Leaderboard

open CertifiedGaps CertifiedGaps.ChainCertificate Real

/-- `e^{β₁ log r / β} z < 1 - z` from `r^a z^b < (1 - z)^b` when `β₁ b = β a`. -/
theorem exp_mul_lt_of_pow {r β β₁ z : ℝ} {a b : ℕ} (hr : 0 < r) (hβ : 0 < β)
    (hab : β₁ * b = β * a) (hz1 : z < 1) (h : r ^ a * z ^ b < (1 - z) ^ b) :
    exp (β₁ * (log r / β)) * z < 1 - z := by
  have hx : exp (β₁ * (log r / β)) ^ b = r ^ a := by
    rw [← exp_nat_mul, ← exp_log hr, ← exp_nat_mul, exp_log hr]
    congr 1
    field_simp
    linear_combination log r * hab
  refine lt_of_pow_lt_pow_left₀ b (by linarith) ?_
  rwa [mul_pow, hx]

/-- **Records from a chain and a bridge seed.** A valid chain certificate `c` from the even seed
`(n₀, β)` and the bridge seed `(n₁, β₁, s₁, t₁)` with `n₁ ≤ n₀ ≤ 2n₁` and
`r^a z(s₁, t₁)^b < (1 - z(s₁, t₁))^b`, `β₁ b = β a`, give `Δ_L > log r / β` for every even
`L ≥ n₁`. -/
theorem Record.TwoSeeded.of_chain {d : Type} [Fintype d] [DecidableEq d] [Nonempty d]
    (spectrum : (L : ℕ) → (Fin L → d) → ℝ) {n₀ n₁ : ℕ} (hZ : HasSpatialMomentBounds spectrum n₁)
    (c : ChainCertificate) (hc : c.Valid) {r β s t β₁ s₁ t₁ : ℚ} (a b : ℕ)
    (hn₁ : 2 ≤ n₁ := by norm_num) (hn₁e : Even n₁ := by decide) (hn₀e : Even n₀ := by decide)
    (hn : n₁ ≤ n₀ := by norm_num) (hn₂ : n₀ ≤ 2 * n₁ := by norm_num) (hr : c.r = r := by rfl)
    (hs : c.sAt 0 = s := by rfl) (ht : c.tAt 0 = t := by rfl) (hβ : 0 < β := by norm_num)
    (hβ₁ : 0 < β₁ := by norm_num) (hr₁ : 1 < r := by norm_num)
    (hab : β₁ * b = β * a := by norm_num) (hs₁ : s₁ ≤ 1 := by norm_num)
    (hz : zQ s₁ t₁ < 1 := by decide +kernel)
    (hwin : r ^ a * zQ s₁ t₁ ^ b < (1 - zQ s₁ t₁) ^ b := by decide +kernel) :
    Record.TwoSeeded spectrum n₁ r β n₀ s t n₁ β₁ s₁ t₁ := by
  refine ⟨⟨hβ, hr₁⟩, fun hseed ⟨hbs, hbt⟩ L hL hLe ↦ ?_⟩
  by_cases hL₀ : n₀ ≤ L
  · exact (Record.Seeded.of_chain spectrum (fun β hβ ↦ (hZ β hβ).mono hn) c hc (by omega) hn₀e
      hr hs ht hβ hr₁).2 hseed L hL₀ hLe
  · have hP := thermalPurity_ge_of_window spectrum hZ le_rfl hn₁e (by omega)
      (by exact_mod_cast hβ₁) (by exact_mod_cast hs₁) hbs hbt hL (by omega) hLe
    rw [cast_zQ] at hP
    refine gapGT_of_thermalPurity spectrum L (by exact_mod_cast hβ₁) (by exact_mod_cast hz) hP
      (exp_mul_lt_of_pow (by exact_mod_cast (zero_lt_one.trans hr₁)) (by exact_mod_cast hβ)
        (by exact_mod_cast hab) (by exact_mod_cast hz)
        (by exact_mod_cast hwin))

end Leaderboard
