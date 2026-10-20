/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.PolynomialFilter
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.Pow

/-!
# Monotone power-to-polynomial ratios

The boundary initialization of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain*
(Lemma in:seedbounds) bounds a leading contribution `c λᵐ` from a bound on `c P(λ)²` for a cubic
`P`, using only an enclosure `ℓ ≤ λ ≤ U`: if `x ↦ xᵐ / P(x)²` increases on `[ℓ, U]`, then
`c λᵐ = c P(λ)² · λᵐ/P(λ)² ≤ c P(λ)² · Uᵐ/P(U)²`. With `m = 2k`, the ratio `x^k / P(x)` has
derivative `x^{k-1} (k P(x) - x P'(x)) / P(x)²`, so it increases on `[ℓ, U]` as soon as `P > 0`
and `x P'(x) ≤ k P(x)` there (`pow_div_sq_le`).

For a rational polynomial both conditions on a whole ray `x ≥ v` follow from finitely many
signs, in the style of `RayCheck`: the Taylor coefficients at `v` of `P` (with positive constant
term) and of `k P - x P'` are nonnegative (`RatioCheck`, decidable). Since `P > 0`, the condition
for `k₀` implies the one for every `k ≥ k₀`.

## Main definitions

* `derivL p`: the coefficients of the derivative of the polynomial `p`.
* `RatioCheck p k v`: the Taylor-coefficient certificate for `P > 0` and `x P' ≤ k P` on `x ≥ v`.

## Main results

* `hasDerivAt_eval`: `eval (derivL p)` is the derivative of `eval p`.
* `pow_div_sq_le`: `P > 0` and `x P' ≤ k P` on `[ℓ, U]` (`0 ≤ ℓ`) give
  `x^{2k} / P(x)² ≤ U^{2k} / P(U)²` for `x ∈ [ℓ, U]`.
* `pow_div_sq_le_of_ratioCheck`: the same from `RatioCheck p k₀ v` and `k₀ ≤ k`, for
  `v ≤ x ≤ U`.
-/

open Set

namespace CertifiedGaps.PolynomialFilter

/-! ### Derivatives of coefficient lists -/

/-- The derivative of a coefficient list: `(a + x p(x))' = p(x) + x p'(x)`. -/
def derivL : List ℚ → List ℚ
  | [] => []
  | _ :: p => add p (0 :: derivL p)

/-- `eval (derivL p)` is the derivative of `eval p`. -/
theorem hasDerivAt_eval : ∀ (p : List ℚ) (x : ℝ), HasDerivAt (eval p) (eval (derivL p) x) x
  | [], x => by
    have h : eval [] = fun _ ↦ (0 : ℝ) := funext fun y ↦ eval_nil y
    rw [h, derivL, eval_nil]
    exact hasDerivAt_const x 0
  | a :: p, x => by
    have h : eval (a :: p) = fun y ↦ (a : ℝ) + y * eval p y := funext fun y ↦ eval_cons a p y
    rw [h]
    convert ((hasDerivAt_id' x).mul (hasDerivAt_eval p x)).const_add (a : ℝ) using 1
    simp only [derivL, eval_add, eval_cons, Rat.cast_zero, zero_add, one_mul]

/-! ### Monotone ratios -/

/-- **A monotone power-to-polynomial ratio.** If `0 ≤ ℓ`, `P > 0` and `x P'(x) ≤ k P(x)` on
`[ℓ, U]` with `k ≥ 1`, then `x ↦ x^{2k} / P(x)²` increases on `[ℓ, U]`; in particular it is at
most its value at `U`. -/
theorem pow_div_sq_le {p : List ℚ} {k : ℕ} (hk : 1 ≤ k) {ℓ U : ℝ} (hℓ : 0 ≤ ℓ)
    (hpos : ∀ x ∈ Icc ℓ U, 0 < eval p x)
    (hder : ∀ x ∈ Icc ℓ U, x * eval (derivL p) x ≤ k * eval p x) {x : ℝ} (hx : x ∈ Icc ℓ U) :
    x ^ (2 * k) / eval p x ^ 2 ≤ U ^ (2 * k) / eval p U ^ 2 := by
  set f : ℝ → ℝ := fun y ↦ y ^ k / eval p y
  -- The derivative of `f = y^k / P`, from the quotient rule.
  have hf : ∀ y ∈ Icc ℓ U, HasDerivAt f ((k * y ^ (k - 1) * eval p y -
      y ^ k * eval (derivL p) y) / eval p y ^ 2) y := fun y hy ↦
    (hasDerivAt_pow k y).div (hasDerivAt_eval p y) (hpos y hy).ne'
  -- It is `y^{k-1} (k P - y P') / P² ≥ 0`, since `y ≥ ℓ ≥ 0`.
  have hf0 : ∀ y ∈ Icc ℓ U, 0 ≤ (k * y ^ (k - 1) * eval p y - y ^ k * eval (derivL p) y) /
      eval p y ^ 2 := fun y hy ↦ by
    have hy0 : 0 ≤ y := hℓ.trans hy.1
    have e : (k : ℝ) * y ^ (k - 1) * eval p y - y ^ k * eval (derivL p) y =
        y ^ (k - 1) * (k * eval p y - y * eval (derivL p) y) := by
      obtain ⟨j, rfl⟩ : ∃ j, k = j + 1 := ⟨k - 1, by omega⟩
      rw [Nat.add_sub_cancel, pow_succ]; push_cast; ring
    rw [e]
    exact div_nonneg (mul_nonneg (pow_nonneg hy0 _) (sub_nonneg.2 (hder y hy))) (sq_nonneg _)
  have hmono : MonotoneOn f (Icc ℓ U) := by
    refine monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc ℓ U) (f' := fun y ↦
      (k * y ^ (k - 1) * eval p y - y ^ k * eval (derivL p) y) / eval p y ^ 2)
      (fun y hy ↦ (hf y hy).continuousAt.continuousWithinAt) (fun y hy ↦ ?_) (fun y hy ↦ ?_)
    · rw [interior_Icc] at hy
      exact (hf y (Ioo_subset_Icc_self hy)).hasDerivWithinAt
    · rw [interior_Icc] at hy
      exact hf0 y (Ioo_subset_Icc_self hy)
  have hU : U ∈ Icc ℓ U := ⟨hx.1.trans hx.2, le_rfl⟩
  have hle : f x ≤ f U := hmono hx hU hx.2
  have hfx : 0 ≤ f x := div_nonneg (pow_nonneg (hℓ.trans hx.1) k) (hpos x hx).le
  have hsq := pow_le_pow_left₀ hfx hle 2
  simp only [f, div_pow, ← pow_mul, mul_comm k 2] at hsq
  exact hsq

/-! ### Certificates on a ray -/

/-- The ratio certificate at `v` for the exponent `k`: the Taylor coefficients at `v` of `P` are
nonnegative with a positive constant term (`RayCheck p v 0`), and those of `k P - x P'` are
nonnegative. It gives `P > 0` and `x P' ≤ k P` on `x ≥ v`. -/
def RatioCheck (p : List ℚ) (k : ℕ) (v : ℚ) : Prop :=
  RayCheck p v 0 ∧ ∀ c ∈ shift v (add (smul k p) (smul (-1) (0 :: derivL p))), 0 ≤ c

instance (p : List ℚ) (k : ℕ) (v : ℚ) : Decidable (RatioCheck p k v) := by
  unfold RatioCheck; infer_instance

/-- The positivity half of the ratio certificate: `P > 0` on `x ≥ v`. -/
theorem RatioCheck.eval_pos {p : List ℚ} {k : ℕ} {v : ℚ} (h : RatioCheck p k v) {x : ℝ}
    (hx : (v : ℝ) ≤ x) : 0 < eval p x := by
  simpa using lt_eval_of_rayCheck h.1 hx

/-- The derivative half of the ratio certificate: `x P' ≤ k P` on `x ≥ v`. -/
theorem RatioCheck.mul_deriv_le {p : List ℚ} {k : ℕ} {v : ℚ} (h : RatioCheck p k v) {x : ℝ}
    (hx : (v : ℝ) ≤ x) : x * eval (derivL p) x ≤ k * eval p x := by
  have h0 := eval_nonneg_of_coeffs _ h.2 (x - v) (sub_nonneg.2 hx)
  rw [eval_shift, add_sub_cancel, eval_add, eval_smul, eval_smul, eval_cons] at h0
  push_cast at h0
  linarith

/-- **The monotone ratio from a certificate.** If `RatioCheck p k₀ v` with `0 ≤ v`, then for every
`k ≥ max k₀ 1` and `v ≤ x ≤ U`, `x^{2k} / P(x)² ≤ U^{2k} / P(U)²`. -/
theorem pow_div_sq_le_of_ratioCheck {p : List ℚ} {k₀ : ℕ} {v : ℚ} (h : RatioCheck p k₀ v)
    (hv : 0 ≤ v) {k : ℕ} (hk : 1 ≤ k) (hk₀ : k₀ ≤ k) {x U : ℝ} (hx : (v : ℝ) ≤ x) (hxU : x ≤ U) :
    x ^ (2 * k) / eval p x ^ 2 ≤ U ^ (2 * k) / eval p U ^ 2 := by
  refine pow_div_sq_le hk (by exact_mod_cast hv) (fun y hy ↦ h.eval_pos hy.1)
    (fun y hy ↦ ?_) ⟨hx, hxU⟩
  have hk' : (k₀ : ℝ) ≤ k := by exact_mod_cast hk₀
  nlinarith [h.mul_deriv_le hy.1, h.eval_pos hy.1]

end CertifiedGaps.PolynomialFilter
