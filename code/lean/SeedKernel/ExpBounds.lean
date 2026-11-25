/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Complex.Exponential

/-!
# Rational upper bounds on the exponential

`expUpQ k y` is a rational upper bound on `e^y` for every rational `y < k` (`exp_le_expUpQ`):

* for `y ≤ 0`, `1 / ∑_{i<k} (-y)^i / i!`, since the Taylor sums of `e^{-y}` lie below it;
* for `0 < y < k`, `(1 - y/k)^{-k}`, since `1 - y/k ≤ e^{-y/k}`.

Both are exact rationals, so an inequality between rationals and `expUpQ` values is decided by
`decide +kernel`. This replaces hand-made decimal bounds on `e` in seed proofs: the seed
conditions of `SeedKernel.Seeds` are decidable propositions in the inverse temperature.
-/

open Real Finset

namespace SeedKernel

/-- The Taylor sum `∑_{i<k} y^i / i!`. -/
def expTaylorQ (k : ℕ) (y : ℚ) : ℚ := ∑ i ∈ range k, y ^ i / (i.factorial : ℚ)

/-- A rational upper bound on `e^y`, valid for `y < k` (`exp_le_expUpQ`). -/
def expUpQ (k : ℕ) (y : ℚ) : ℚ :=
  if y ≤ 0 then 1 / expTaylorQ k (-y) else (1 / (1 - y / k)) ^ k

theorem expTaylorQ_cast (k : ℕ) (y : ℚ) :
    (expTaylorQ k y : ℝ) = ∑ i ∈ range k, (y : ℝ) ^ i / (i.factorial : ℝ) := by
  rw [expTaylorQ]
  exact (map_sum (Rat.castHom ℝ) _ _).trans (by simp)

theorem one_le_expTaylorQ {k : ℕ} (hk : 0 < k) {y : ℚ} (hy : 0 ≤ y) : 1 ≤ expTaylorQ k y := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_lt hk
  rw [expTaylorQ, zero_add, sum_range_succ', pow_zero, Nat.factorial_zero, Nat.cast_one, div_one]
  have : 0 ≤ ∑ i ∈ range k, y ^ (i + 1) / ((i + 1).factorial : ℚ) :=
    sum_nonneg fun i _ ↦ by positivity
  linarith

/-- **`expUpQ` bounds the exponential.** `e^y ≤ expUpQ k y` for `y < k`. -/
theorem exp_le_expUpQ {k : ℕ} (hk : 0 < k) {y : ℚ} (hy : y < k) :
    exp (y : ℝ) ≤ (expUpQ k y : ℝ) := by
  unfold expUpQ
  split_ifs with h0
  · have hT := one_le_expTaylorQ hk (neg_nonneg.2 h0)
    have hle : ((expTaylorQ k (-y) : ℚ) : ℝ) ≤ exp (-(y : ℝ)) := by
      rw [expTaylorQ_cast]
      push_cast
      exact sum_le_exp_of_nonneg (by exact_mod_cast neg_nonneg.2 h0) k
    have hT' : (1 : ℝ) ≤ (expTaylorQ k (-y) : ℝ) := by exact_mod_cast hT
    push_cast
    rw [le_div_iff₀ (by linarith)]
    calc exp (y : ℝ) * (expTaylorQ k (-y) : ℝ) ≤ exp (y : ℝ) * exp (-(y : ℝ)) := by gcongr
      _ = 1 := by rw [← exp_add]; simp
  · push Not at h0
    have hk' : (0 : ℝ) < k := by exact_mod_cast hk
    have hyk : (y : ℝ) / k < 1 := by
      rw [div_lt_one hk']; exact_mod_cast hy
    have h1 : 0 < 1 - (y : ℝ) / k := by linarith
    -- `e^{y/k} ≤ 1 / (1 - y/k)` from `1 - y/k ≤ e^{-y/k}`
    have hstep : exp ((y : ℝ) / k) ≤ 1 / (1 - (y : ℝ) / k) := by
      rw [le_div_iff₀ h1]
      have := add_one_le_exp (-((y : ℝ) / k))
      calc exp ((y : ℝ) / k) * (1 - (y : ℝ) / k) ≤ exp ((y : ℝ) / k) * exp (-((y : ℝ) / k)) := by
            gcongr; linarith
        _ = 1 := by rw [← exp_add]; simp
    have hpow : exp (y : ℝ) = exp ((y : ℝ) / k) ^ k := by
      rw [← exp_nat_mul]; congr 1; field_simp
    push_cast
    rw [hpow]
    exact pow_le_pow_left₀ (exp_pos _).le hstep k

end SeedKernel
