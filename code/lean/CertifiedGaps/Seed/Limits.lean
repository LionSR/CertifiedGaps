/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic.FieldSimp

/-!
# Limits and energy shifts of seeds

Two model-independent steps of every seed certificate:

* `seed_limit`: the seed inequalities pass from approximations within `δ` to the traces
  themselves as `δ → 0`;
* `purities_of_shifted`: the spatial and thermal purities are invariant under an energy shift
  `H ↦ H + κ` per bond, which rescales the traces as `Z_n(β) ↦ e^{-βκn} Z_n(β)`
  (`shiftedTrace`). Tables are stated for the shifted traces, the purities for the traces.
-/

open Real

namespace CertifiedGaps

/-! ### The limit `δ → 0` -/

/-- If `a ≤ b + K δ` for all small `δ > 0`, then `a ≤ b`. -/
theorem le_of_forall_small {a b K δ₀ : ℝ} (hK : 0 < K) (hδ₀ : 0 < δ₀)
    (h : ∀ δ, 0 < δ → δ < δ₀ → a ≤ b + K * δ) : a ≤ b := by
  refine le_of_forall_pos_le_add fun ε hε ↦ ?_
  set δ := min (ε / K) (δ₀ / 2)
  have hδ : 0 < δ := lt_min (div_pos hε hK) (by linarith)
  have h1 := h δ hδ (lt_of_le_of_lt (min_le_right _ _) (by linarith))
  have h2 : K * δ ≤ ε := by
    have := min_le_left (ε / K) (δ₀ / 2)
    rwa [le_div_iff₀ hK, mul_comm] at this
  linarith

/-- **The limit `δ → 0` of a seed.** If for every small `δ > 0` there are approximations within
`δ` of `Z₁ > 0` and `Z₂` that satisfy the seed inequalities, in the form
`(1 - p) (Z₁ - δ)² ≤ Z₂ + δ` and `Z₂ ≤ B + δ`, then `(1 - p) Z₁² ≤ Z₂` and `Z₂ ≤ B`. -/
theorem seed_limit {Z₁ Z₂ B p δ₀ : ℝ} (hp : p ≤ 1) (hZ₁ : 0 < Z₁) (hδ₀ : 0 < δ₀)
    (h : ∀ δ, 0 < δ → δ < δ₀ → (1 - p) * (Z₁ - δ) ^ 2 ≤ Z₂ + δ ∧ Z₂ ≤ B + δ) :
    (1 - p) * Z₁ ^ 2 ≤ Z₂ ∧ Z₂ ≤ B := by
  have hp0 : 0 ≤ 1 - p := by linarith
  refine ⟨le_of_forall_small (K := 1 + 2 * (1 - p) * Z₁) (by positivity) hδ₀
      fun δ hδ hδ' ↦ ?_,
    le_of_forall_small (K := 1) one_pos hδ₀ fun δ hδ hδ' ↦ by linarith [(h δ hδ hδ').2]⟩
  have h1 := (h δ hδ hδ').1
  have e : (1 - p) * Z₁ ^ 2 = (1 - p) * (Z₁ - δ) ^ 2 + (1 - p) * (2 * Z₁ * δ) -
      (1 - p) * δ ^ 2 := by ring
  have e2 : (1 + 2 * (1 - p) * Z₁) * δ = δ + (1 - p) * (2 * Z₁ * δ) := by ring
  have t2 : 0 ≤ (1 - p) * δ ^ 2 := mul_nonneg hp0 (sq_nonneg δ)
  linarith

/-! ### Energy shifts -/

/-- The shifted trace `e^{-βκn} z`. -/
noncomputable def shiftedTrace (β κ : ℝ) (n : ℕ) (z : ℝ) : ℝ := exp (-(β * κ * n)) * z

theorem shiftedTrace_eq_pow (β κ : ℝ) (n : ℕ) (z : ℝ) :
    shiftedTrace β κ n z = exp (-(β * κ)) ^ n * z := by
  rw [shiftedTrace, ← exp_nat_mul]; ring_nf

/-- The trial-state bound `e^{βκn} < Z` in shifted form. -/
theorem one_lt_shiftedTrace_of_exp_lt {β κ z : ℝ} {n : ℕ} (h : exp (β * κ * n) < z) :
    1 < shiftedTrace β κ n z := by
  calc (1 : ℝ) = exp (-(β * κ * n)) * exp (β * κ * n) := by rw [← exp_add]; simp
    _ < _ := mul_lt_mul_of_pos_left h (exp_pos _)

/-- **Both purities are invariant under an energy shift.** If the shifted traces
`s^n Z₁`, `s^{2n} Z₂` satisfy the seed inequalities `(1 - p) (sⁿ Z₁)² ≤ s^{2n} Z₂ ≤ B`, the
shifted trace at double temperature `s^{4n} Z₂' > 1`, and `(1 - q) B² ≤ 1`, then
`1 - p ≤ Z₂ / Z₁²` and `1 - q ≤ Z₂' / Z₂²`. -/
theorem purities_of_shifted {s Z₁ Z₂ Z₂' p q B : ℝ} {n : ℕ} (hs : 0 < s) (hZ₁ : 0 < Z₁)
    (hp : (1 - p) * (s ^ n * Z₁) ^ 2 ≤ s ^ (2 * n) * Z₂) (hB : s ^ (2 * n) * Z₂ ≤ B)
    (htrial : 1 < s ^ (2 * n) * Z₂) (htrial₂ : 1 < (s ^ (2 * n)) ^ 2 * Z₂')
    (hq : (1 - q) * B ^ 2 ≤ 1) (hq1 : q ≤ 1) :
    1 - p ≤ Z₂ / Z₁ ^ 2 ∧ 1 - q ≤ Z₂' / Z₂ ^ 2 := by
  have hsn : 0 < s ^ n := pow_pos hs _
  have hs2 : s ^ (2 * n) = (s ^ n) ^ 2 := by ring
  have hZ₂ : 0 < Z₂ := by
    rw [hs2] at htrial
    exact pos_of_mul_pos_right (by linarith : (0 : ℝ) < (s ^ n) ^ 2 * Z₂) (by positivity)
  constructor
  · rw [le_div_iff₀ (by positivity)]
    rw [mul_pow, hs2] at hp
    exact (mul_le_mul_iff_of_pos_left (by positivity : (0 : ℝ) < (s ^ n) ^ 2)).1
      (by linarith : (s ^ n) ^ 2 * ((1 - p) * Z₁ ^ 2) ≤ (s ^ n) ^ 2 * Z₂)
  · have hq0 : 0 ≤ 1 - q := by linarith
    have hsq : (s ^ (2 * n) * Z₂) ^ 2 ≤ B ^ 2 := pow_le_pow_left₀ (by linarith) hB 2
    have hmain : (1 - q) * (s ^ (2 * n) * Z₂) ^ 2 ≤ (s ^ (2 * n)) ^ 2 * Z₂' := by nlinarith
    rw [le_div_iff₀ (by positivity)]
    rw [mul_pow] at hmain
    exact (mul_le_mul_iff_of_pos_left (by positivity : (0 : ℝ) < (s ^ (2 * n)) ^ 2)).1
      (by linarith : (s ^ (2 * n)) ^ 2 * ((1 - q) * Z₂ ^ 2) ≤ (s ^ (2 * n)) ^ 2 * Z₂')

/-- The shifted trace at double temperature: `e^{-2βκm} z = (e^{-βκ})^{2m} z`. -/
theorem shiftedTrace_two_mul (β κ : ℝ) (m : ℕ) (z : ℝ) :
    shiftedTrace (2 * β) κ m z = (exp (-(β * κ)) ^ m) ^ 2 * z := by
  rw [shiftedTrace, ← pow_mul, ← exp_nat_mul]
  push_cast; ring_nf

end CertifiedGaps
