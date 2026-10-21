/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Basic.Real.Basic
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.GCongr

/-!
# Power sums of bounded weights

Bounds on `∑ⱼ tⱼⁿ` for nonnegative weights `tⱼ` from their total mass and a cap on each weight.
In a seed, the weights are `λ^q` for the transfer eigenvalues `λ` of one symmetry sector: the
filters cap each `|λ|`, the table bounds the mass `∑ λ^q`, and the bounds below control the
higher moments `∑ λ^(q r)` (eq. (c-capped) of the certificate note).

## Main results

* `sum_pow_succ_le`: the chord bound `∑ tⱼ^(r+1) ≤ cʳ ∑ tⱼ` for weights in `[0, c]`.
* `sum_pow_le_pow_sum`: superadditivity `∑ tⱼⁿ ≤ (∑ tⱼ)ⁿ`.
* `pow_add_add_pow_le`: moving mass from the smaller to the larger of two weights raises the sum
  of `n`th powers.
* `sum_pow_le_capped`: the capped-mass bound with one full bar, `∑ tⱼⁿ ≤ hⁿ + (m - h)ⁿ` for
  weights in `[0, h]` of mass at most `m ∈ [h, 2h]`.
-/

open Finset

namespace CertifiedGaps

section PowerSums

variable {κ : Type*} (s : Finset κ) (t : κ → ℝ)

/-- The chord bound: weights in `[0, c]` have `∑ tⱼ^(r+1) ≤ cʳ ∑ tⱼ`. -/
theorem sum_pow_succ_le (ht0 : ∀ j ∈ s, 0 ≤ t j) {c : ℝ} (htc : ∀ j ∈ s, t j ≤ c) (r : ℕ) :
    ∑ j ∈ s, t j ^ (r + 1) ≤ c ^ r * ∑ j ∈ s, t j := by
  rw [mul_sum]
  refine sum_le_sum fun j hj ↦ ?_
  rw [pow_succ]
  exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (ht0 j hj) (htc j hj) r) (ht0 j hj)

/-- Superadditivity of `tⁿ`: nonnegative weights have `∑ tⱼⁿ ≤ (∑ tⱼ)ⁿ` for `n ≥ 1`. -/
theorem sum_pow_le_pow_sum (ht0 : ∀ j ∈ s, 0 ≤ t j) {n : ℕ} (hn : 1 ≤ n) :
    ∑ j ∈ s, t j ^ n ≤ (∑ j ∈ s, t j) ^ n := by
  obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le' hn
  have := sum_pow_succ_le s t ht0 (fun j hj ↦ single_le_sum ht0 hj) r
  rwa [← pow_succ] at this

end PowerSums

/-- Increments of `xⁿ` grow with `x`: `(y + d)ⁿ - yⁿ ≤ (x + d)ⁿ - xⁿ` for `0 ≤ y ≤ x`, `d ≥ 0`.
Equivalently, moving mass `d` from the smaller to the larger of two weights raises the sum of
`n`th powers. -/
theorem pow_add_add_pow_le {x y d : ℝ} (hy : 0 ≤ y) (hyx : y ≤ x) (hd : 0 ≤ d) (n : ℕ) :
    x ^ n + (y + d) ^ n ≤ (x + d) ^ n + y ^ n := by
  have key : ∀ z : ℝ, (z + d) ^ n - z ^ n =
      ∑ k ∈ range n, z ^ k * d ^ (n - k) * (n.choose k : ℝ) := fun z ↦ by
    rw [add_pow, sum_range_succ, Nat.sub_self, pow_zero, Nat.choose_self, Nat.cast_one,
      mul_one, mul_one, add_sub_cancel_right]
  have h := sum_le_sum fun k (_ : k ∈ range n) ↦
    mul_le_mul_of_nonneg_right
      (mul_le_mul_of_nonneg_right (pow_le_pow_left₀ hy hyx k) (pow_nonneg hd (n - k)))
      (Nat.cast_nonneg (α := ℝ) (n.choose k))
  rw [← key, ← key] at h
  linarith

/-- **Capped masses** (eq. (c-capped) with one full bar). Weights `tⱼ ∈ [0, h]` of total mass at
most `m`, with `h ≤ m ≤ 2h`, have `∑ tⱼⁿ ≤ hⁿ + (m - h)ⁿ` for `n ≥ 1`. -/
theorem sum_pow_le_capped {κ : Type*} (s : Finset κ) (t : κ → ℝ) {h m : ℝ}
    (ht0 : ∀ j ∈ s, 0 ≤ t j) (hth : ∀ j ∈ s, t j ≤ h) (hsum : ∑ j ∈ s, t j ≤ m)
    (hhm : h ≤ m) (hm2 : m ≤ 2 * h) {n : ℕ} (hn : 1 ≤ n) :
    ∑ j ∈ s, t j ^ n ≤ h ^ n + (m - h) ^ n := by
  classical
  have hmh : 0 ≤ m - h := sub_nonneg.2 hhm
  -- The two-point comparison `aⁿ + (m - a)ⁿ ≤ hⁿ + (m - h)ⁿ` for `m - h ≤ a ≤ h`.
  have two {a : ℝ} (ha : m - h ≤ a) (hah : a ≤ h) : a ^ n + (m - a) ^ n ≤ h ^ n + (m - h) ^ n := by
    have := pow_add_add_pow_le hmh ha (sub_nonneg.2 hah) n
    rwa [show m - h + (h - a) = m - a by ring, show a + (h - a) = h by ring] at this
  rcases s.eq_empty_or_nonempty with hs | hs
  · simp only [hs, sum_empty]
    exact add_nonneg (pow_nonneg (by linarith) n) (pow_nonneg hmh n)
  obtain ⟨j₀, hj₀, hmax⟩ := s.exists_max_image t hs
  rcases le_or_gt (t j₀) (m / 2) with hu | hu
  · -- All weights are at most `m/2`: the chord bound with `c = m/2`.
    obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le' hn
    have hm0 : 0 ≤ m / 2 := by linarith
    calc ∑ j ∈ s, t j ^ (r + 1) ≤ (m / 2) ^ r * ∑ j ∈ s, t j :=
          sum_pow_succ_le s t ht0 (fun j hj ↦ (hmax j hj).trans hu) r
      _ ≤ (m / 2) ^ r * m := mul_le_mul_of_nonneg_left hsum (pow_nonneg hm0 r)
      _ = (m / 2) ^ (r + 1) + (m - m / 2) ^ (r + 1) := by ring
      _ ≤ h ^ (r + 1) + (m - h) ^ (r + 1) := two (by linarith) (by linarith)
  · -- One weight exceeds `m/2`: it stays, the others obey superadditivity.
    have hrest : ∑ j ∈ s.erase j₀, t j ≤ m - t j₀ := by
      have := add_sum_erase s t hj₀; linarith
    have hrest0 : 0 ≤ ∑ j ∈ s.erase j₀, t j :=
      sum_nonneg fun j hj ↦ ht0 j (mem_of_mem_erase hj)
    calc ∑ j ∈ s, t j ^ n = t j₀ ^ n + ∑ j ∈ s.erase j₀, t j ^ n :=
          (add_sum_erase s (fun j ↦ t j ^ n) hj₀).symm
      _ ≤ t j₀ ^ n + (m - t j₀) ^ n := by
          gcongr
          exact (sum_pow_le_pow_sum _ t (fun j hj ↦ ht0 j (mem_of_mem_erase hj)) hn).trans
            (pow_le_pow_left₀ hrest0 hrest n)
      _ ≤ h ^ n + (m - h) ^ n := two (by linarith) (hth j₀ hj₀)

end CertifiedGaps
