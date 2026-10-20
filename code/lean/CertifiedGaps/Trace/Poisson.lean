/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.SpecificLimits.Normed
import CertifiedGaps.Spectral.FunctionalCalculus
import CertifiedGaps.Spectral.OperatorBounds

/-!
# Conditioned Poisson weights and exponential tails

The thermal tables of the spin-1 Haldane-gap paper (`thermal72.tex`, `thermal120.tex` in
`references/periodic-gap-source`) evaluate `F = e^{m(R - I)}` for a Hermitian contraction `R`
through the polynomial `∑_{j ≤ J} pⱼ Rʲ` with rational Poisson weights of mean `m`, conditioned
to `j ≤ J`. This file proves the truncation estimate of that step (the coefficient argument of
Fox and Glynn, Proposition 1) and the Taylor enclosure of `e^x` used to turn the integer traces
into rational bounds.

Both operators are functions of `R`, so the matrix estimate reduces to a scalar one on each
eigenvalue through the functional calculus `hR.cfc f = U diag(f(λᵢ)) Uᴴ` of
`Spectral/FunctionalCalculus.lean` (`nsq_hermFun_mulVec_le` bounds `‖f(R) v‖` by
`max |f(λᵢ)| ‖v‖`).

## Definitions

* `expTail y K = e^y - ∑_{j<K} y^j/j!`, the tail of the exponential series.
* `poissonMass m J = ∑_{ℓ ≤ J} m^ℓ/ℓ!` and `poissonWeight m J j = (m^j/j!) / poissonMass m J`
  (rational and computable).
* `poissonError m J = 2 e^{-m} ∑_{j > J} m^j/j!`.

## Main results

* `sum_poissonWeight`, `poissonWeight_nonneg`: the weights form a probability vector;
  `poissonWeight_eq`: the paper's form `vⱼ / ∑ v_ℓ` with `vⱼ = m^j J!/j!`.
* `abs_exp_sub_poissonWeight_le`: `|e^{m(x-1)} - ∑_{j ≤ J} pⱼ xʲ| ≤ η(m, J)` on `[-1, 1]`.
* `poissonError_le_exp_div`, `poissonError_le_three_pow_div`: `η(m, J) ≤ e^m/2^J ≤ 3^m/2^J`.
* `nsq_exp_sub_poisson_mulVec_le` (and `_of_qform`): the matrix estimate
  `‖(e^{m(R - 1)} - ∑_{j ≤ J} pⱼ Rʲ) v‖² ≤ η(m, J)² ‖v‖²` for Hermitian `R` with spectrum in
  `[-1, 1]`.
* `nsq_exp_sub_poisson_98_280`, `nsq_exp_sub_poisson_84_256`: the instances with
  `η < 2^{-84}` and `η < 2^{-116}`.
* `exp_le_sum_add_geometric`, `exp_mem_Icc_taylor`, `exp_mem_Icc_taylor_100`: the enclosure
  `s ≤ e^x ≤ s + t` with `s = ∑_{j ≤ N} x^j/j!` and `t = x^N/N!` for `0 ≤ 2x ≤ N + 1`.
-/

open Finset

open scoped Nat

namespace CertifiedGaps.Trace

/-! ### Tails of the exponential series -/

/-- The tail `e^y - ∑_{j<K} y^j/j!` of the exponential series after `K` terms. -/
noncomputable def expTail (y : ℝ) (K : ℕ) : ℝ := Real.exp y - ∑ j ∈ range K, y ^ j / j !

/-- The tail as a series. -/
theorem expTail_eq_tsum (y : ℝ) (K : ℕ) :
    expTail y K = ∑' i, y ^ (i + K) / (i + K)! := by
  have h : Real.exp y = ∑' n, y ^ n / n ! := by
    rw [Real.exp_eq_exp_ℝ, NormedSpace.exp_eq_tsum_div]
  rw [expTail, h, ← (Real.summable_pow_div_factorial y).sum_add_tsum_nat_add K]
  ring

/-- The shifted exponential series is summable. -/
theorem summable_tail (y : ℝ) (K : ℕ) : Summable fun i ↦ y ^ (i + K) / (i + K)! :=
  (summable_nat_add_iff K).2 (Real.summable_pow_div_factorial y)

/-- For `y ≥ 0` every omitted term is nonnegative. -/
theorem expTail_nonneg {y : ℝ} (hy : 0 ≤ y) (K : ℕ) : 0 ≤ expTail y K := by
  rw [expTail_eq_tsum]
  exact tsum_nonneg fun i ↦ by positivity

/-- The tail at `y` is bounded by the tail at any `m ≥ |y|`. -/
theorem abs_expTail_le {y m : ℝ} (hy : |y| ≤ m) (K : ℕ) : |expTail y K| ≤ expTail m K := by
  rw [expTail_eq_tsum, expTail_eq_tsum]
  have habs : ∀ i, ‖y ^ (i + K) / ((i + K)! : ℝ)‖ = |y| ^ (i + K) / (i + K)! := by
    intro i
    rw [Real.norm_eq_abs, abs_div, abs_pow, Nat.abs_cast]
  calc |∑' i, y ^ (i + K) / ((i + K)! : ℝ)|
      ≤ ∑' i, ‖y ^ (i + K) / ((i + K)! : ℝ)‖ := by
        rw [← Real.norm_eq_abs]
        refine norm_tsum_le_tsum_norm ?_
        simpa only [habs] using summable_tail |y| K
    _ = ∑' i, |y| ^ (i + K) / ((i + K)! : ℝ) := tsum_congr habs
    _ ≤ ∑' i, m ^ (i + K) / ((i + K)! : ℝ) := by
        refine (summable_tail |y| K).tsum_mono (summable_tail m K) fun i ↦ ?_
        dsimp only
        gcongr

/-- Doubling the argument: `2^K (e^m - ∑_{j<K} m^j/j!) ≤ e^{2m}` for `m ≥ 0`, because every
omitted index `k ≥ K` has `2^K ≤ 2^k`. -/
theorem expTail_mul_two_pow_le {m : ℝ} (hm : 0 ≤ m) (K : ℕ) :
    expTail m K * 2 ^ K ≤ Real.exp (2 * m) := by
  rw [expTail_eq_tsum, ← tsum_mul_right]
  calc ∑' i, m ^ (i + K) / ((i + K)! : ℝ) * 2 ^ K
      ≤ ∑' i, (2 * m) ^ (i + K) / ((i + K)! : ℝ) := by
        refine ((summable_tail m K).mul_right _).tsum_mono (summable_tail (2 * m) K)
          fun i ↦ ?_
        dsimp only
        rw [mul_pow, div_mul_eq_mul_div, mul_comm (m ^ (i + K))]
        gcongr
        · norm_num
        · omega
    _ = expTail (2 * m) K := (expTail_eq_tsum _ _).symm
    _ ≤ Real.exp (2 * m) := by
        rw [expTail, sub_le_self_iff]
        exact sum_nonneg fun j _ ↦ by positivity

/-! ### Conditioned Poisson weights -/

/-- The normalization `∑_{ℓ ≤ J} m^ℓ/ℓ!` of the Poisson weights of mean `m` conditioned to
indices at most `J`. -/
def poissonMass (m : ℚ) (J : ℕ) : ℚ := ∑ ℓ ∈ range (J + 1), m ^ ℓ / ℓ !

/-- The Poisson weight `p_j = (m^j/j!) / ∑_{ℓ ≤ J} m^ℓ/ℓ!` of mean `m`, conditioned to indices
at most `J`. It is used only for `j ≤ J`; the formula is not truncated beyond `J`. -/
def poissonWeight (m : ℚ) (J j : ℕ) : ℚ := m ^ j / j ! / poissonMass m J

/-- The `ℓ = 0` term of the normalization is `1` and the others are nonnegative. -/
theorem one_le_poissonMass {m : ℚ} (hm : 0 ≤ m) (J : ℕ) : 1 ≤ poissonMass m J := by
  rw [poissonMass, sum_range_succ', pow_zero, Nat.factorial_zero, Nat.cast_one, div_one]
  exact le_add_of_nonneg_left (sum_nonneg fun j _ ↦ by positivity)

/-- The normalization is positive, so the weights are well defined. -/
theorem poissonMass_pos {m : ℚ} (hm : 0 ≤ m) (J : ℕ) : 0 < poissonMass m J :=
  one_pos.trans_le (one_le_poissonMass hm J)

/-- The weights are nonnegative for `m ≥ 0`. -/
theorem poissonWeight_nonneg {m : ℚ} (hm : 0 ≤ m) (J j : ℕ) : 0 ≤ poissonWeight m J j :=
  div_nonneg (by positivity) (poissonMass_pos hm J).le

/-- The conditioned weights `p_0, …, p_J` sum to one. -/
theorem sum_poissonWeight {m : ℚ} (hm : 0 ≤ m) (J : ℕ) :
    ∑ j ∈ range (J + 1), poissonWeight m J j = 1 := by
  simp only [poissonWeight, ← sum_div]
  exact div_self (poissonMass_pos hm J).ne'

/-- The paper's form `p_j = v_j / ∑_{ℓ ≤ J} v_ℓ` with the integer numerators
`v_j = m^j J!/j!`. -/
theorem poissonWeight_eq (m : ℚ) (J j : ℕ) :
    poissonWeight m J j =
      m ^ j * J ! / j ! / ∑ ℓ ∈ range (J + 1), m ^ ℓ * J ! / ℓ ! := by
  have h : ∑ ℓ ∈ range (J + 1), m ^ ℓ * J ! / ℓ ! = poissonMass m J * J ! := by
    rw [poissonMass, sum_mul]
    exact sum_congr rfl fun ℓ _ ↦ by ring
  have hJ : (J ! : ℚ) ≠ 0 := by positivity
  rw [h, poissonWeight, mul_div_right_comm, mul_div_mul_right _ _ hJ]

/-- The real cast of a Poisson weight. -/
theorem poissonWeight_cast (m : ℚ) (J j : ℕ) :
    (poissonWeight m J j : ℝ) =
      (m : ℝ) ^ j / j ! / ∑ ℓ ∈ range (J + 1), (m : ℝ) ^ ℓ / ℓ ! := by
  simp [poissonWeight, poissonMass]

/-! ### The scalar Fox–Glynn estimate -/

/-- The error bound `η(m, J) = 2 e^{-m} ∑_{j > J} m^j/j!`: twice the probability that a Poisson
variable of mean `m` exceeds `J`. -/
noncomputable def poissonError (m : ℝ) (J : ℕ) : ℝ := 2 * Real.exp (-m) * expTail m (J + 1)

/-- The error bound is nonnegative for `m ≥ 0`. -/
theorem poissonError_nonneg {m : ℝ} (hm : 0 ≤ m) (J : ℕ) : 0 ≤ poissonError m J :=
  mul_nonneg (by positivity) (expTail_nonneg hm _)

/-- Fox–Glynn for a point of `[-1, 1]`: with `S = ∑_{j ≤ J} m^j/j!`, the conditioned Poisson
average `∑_{j ≤ J} (m^j/j!/S) x^j` approximates `e^{m(x-1)}` within `η(m, J)`.

Write `e^{mx} = A + T` with `A = ∑_{j ≤ J} (mx)^j/j!` and the tail `T`, and `τ = e^m - S`.
Then `e^{m(x-1)} - A/S = e^{-m} T - e^{-m} τ (A/S)`: the truncated tail plus the
renormalization, each at most `e^{-m} τ` because `|T| ≤ τ` and `|A| ≤ S`. -/
theorem abs_exp_sub_poisson_le {m x : ℝ} (hm : 0 ≤ m) (hx : |x| ≤ 1) (J : ℕ) :
    |Real.exp (m * (x - 1)) -
        ∑ j ∈ range (J + 1), m ^ j / j ! / (∑ ℓ ∈ range (J + 1), m ^ ℓ / ℓ !) * x ^ j|
      ≤ poissonError m J := by
  set S := ∑ ℓ ∈ range (J + 1), m ^ ℓ / (ℓ ! : ℝ) with hSdef
  set A := ∑ j ∈ range (J + 1), (m * x) ^ j / (j ! : ℝ) with hAdef
  set T := expTail (m * x) (J + 1) with hTdef
  set τ := expTail m (J + 1) with hτdef
  have hS1 : 1 ≤ S := by
    rw [hSdef, sum_range_succ', pow_zero, Nat.factorial_zero, Nat.cast_one, div_one]
    exact le_add_of_nonneg_left (sum_nonneg fun j _ ↦ by positivity)
  have hS0 : S ≠ 0 := (one_pos.trans_le hS1).ne'
  have hmx : |m * x| ≤ m := by
    rw [abs_mul, abs_of_nonneg hm]
    exact mul_le_of_le_one_right hm hx
  -- `|A| ≤ S` termwise, since `|mx| ≤ m`.
  have hA : |A / S| ≤ 1 := by
    rw [abs_div, abs_of_pos (one_pos.trans_le hS1), div_le_one (one_pos.trans_le hS1)]
    refine (abs_sum_le_sum_abs _ _).trans (sum_le_sum fun j _ ↦ ?_)
    rw [abs_div, abs_pow, Nat.abs_cast]
    gcongr
  have hT : |T| ≤ τ := abs_expTail_le hmx _
  have hτ : 0 ≤ τ := expTail_nonneg hm _
  have hexp : Real.exp (m * (x - 1)) = Real.exp (-m) * (A + T) := by
    have hAT : A + T = Real.exp (m * x) := by simp only [hTdef, expTail, hAdef]; ring
    rw [hAT, ← Real.exp_add]
    ring_nf
  have hpoly : ∑ j ∈ range (J + 1), m ^ j / j ! / S * x ^ j = A / S := by
    rw [hAdef, sum_div]
    exact sum_congr rfl fun j _ ↦ by rw [mul_pow]; ring
  have hE : Real.exp (-m) * Real.exp m = 1 := by rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  have hSτ : S = Real.exp m - τ := by simp only [hτdef, expTail, hSdef]; ring
  have key : Real.exp (-m) * (A + T) - A / S =
      Real.exp (-m) * T - Real.exp (-m) * τ * (A / S) := by
    have hτ' : τ = Real.exp m - S := by linarith
    rw [hτ']
    field_simp
    linear_combination A * hE
  rw [hexp, hpoly, key, poissonError]
  calc |Real.exp (-m) * T - Real.exp (-m) * τ * (A / S)|
      ≤ |Real.exp (-m) * T| + |Real.exp (-m) * τ * (A / S)| := abs_sub _ _
    _ ≤ Real.exp (-m) * τ + Real.exp (-m) * τ * 1 := by
        rw [abs_mul, abs_mul, abs_mul, abs_of_pos (Real.exp_pos _), abs_of_nonneg hτ]
        gcongr
    _ = 2 * Real.exp (-m) * τ := by ring

/-- The scalar estimate for rational `m` and the rational weights `poissonWeight`. -/
theorem abs_exp_sub_poissonWeight_le {m : ℚ} (hm : 0 ≤ m) {x : ℝ} (hx : |x| ≤ 1) (J : ℕ) :
    |Real.exp (m * (x - 1)) - ∑ j ∈ range (J + 1), (poissonWeight m J j : ℝ) * x ^ j|
      ≤ poissonError m J := by
  simp only [poissonWeight_cast]
  exact abs_exp_sub_poisson_le (by exact_mod_cast hm) hx J

/-- The clean form `η(m, J) ≤ e^m / 2^J`: bound the tail by `2^{-(J+1)} e^{2m}`
(`expTail_mul_two_pow_le`). -/
theorem poissonError_le_exp_div {m : ℝ} (hm : 0 ≤ m) (J : ℕ) :
    poissonError m J ≤ Real.exp m / 2 ^ J := by
  have h := expTail_mul_two_pow_le hm (J + 1)
  have h2 : Real.exp (2 * m) = Real.exp m * Real.exp m := by rw [← Real.exp_add]; ring_nf
  have hE : Real.exp (-m) * Real.exp m = 1 := by rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  rw [poissonError, le_div_iff₀ (by positivity)]
  calc 2 * Real.exp (-m) * expTail m (J + 1) * 2 ^ J
      = Real.exp (-m) * (expTail m (J + 1) * 2 ^ (J + 1)) := by ring
    _ ≤ Real.exp (-m) * Real.exp (2 * m) := by gcongr
    _ = Real.exp m := by rw [h2, ← mul_assoc, hE, one_mul]

/-- For natural `m`, `η(m, J) ≤ 3^m / 2^J`, from `e < 3`. -/
theorem poissonError_le_three_pow_div (m J : ℕ) :
    poissonError m J ≤ 3 ^ m / 2 ^ J := by
  refine (poissonError_le_exp_div (Nat.cast_nonneg m) J).trans ?_
  gcongr
  rw [← Real.exp_one_pow]
  gcongr
  exact (Real.exp_one_lt_d9.trans (by norm_num)).le

/-! ### The numerical instances -/

/-- The instance `(m, J) = (98, 280)` of `thermal72.tex`: `η < 2^{-84}`. -/
theorem poissonError_98_280 : poissonError 98 280 < 1 / 2 ^ 84 := by
  have h := poissonError_le_three_pow_div 98 280
  simp only [Nat.cast_ofNat] at h
  refine h.trans_lt ?_
  rw [div_lt_div_iff₀ (by positivity) (by positivity), one_mul]
  -- `3^98 < 4^98 = 2^196`; `norm_num` does not expand the power `2^280` itself.
  calc (3 : ℝ) ^ 98 * 2 ^ 84 < (2 ^ 2) ^ 98 * 2 ^ 84 := by gcongr; norm_num
    _ = 2 ^ 280 := by rw [← pow_mul, ← pow_add]

/-- The instance `(m, J) = (84, 256)` of `thermal120.tex`: `η < 2^{-116}`. The paper bounds
`η ≤ 2 e^{84}/2^{257} < 2 · 27^{28}/2^{257}`; `e < 3` gives the same `3^{84}/2^{256}`. -/
theorem poissonError_84_256 : poissonError 84 256 < 1 / 2 ^ 116 := by
  have h := poissonError_le_three_pow_div 84 256
  simp only [Nat.cast_ofNat] at h
  refine h.trans_lt ?_
  rw [div_lt_div_iff₀ (by positivity) (by positivity), one_mul]
  -- `3^84 < 2^140`, that is `84 log₂ 3 ≈ 133.1 < 140`.
  calc (3 : ℝ) ^ 84 * 2 ^ 116 < 2 ^ 140 * 2 ^ 116 := by gcongr; norm_num
    _ = 2 ^ 256 := by rw [← pow_add]

/-! ### Taylor enclosure of the exponential -/

/-- Lower Taylor bound: `∑_{j ≤ N} x^j/j! ≤ e^x` for `x ≥ 0`. -/
theorem sum_le_exp {x : ℝ} (hx : 0 ≤ x) (N : ℕ) :
    ∑ j ∈ range (N + 1), x ^ j / j ! ≤ Real.exp x :=
  Real.sum_le_exp_of_nonneg hx _

/-- Upper Taylor bound with a geometric tail: for `0 ≤ x < N + 1`,
`e^x ≤ ∑_{j ≤ N} x^j/j! + (x^N/N!) · x/(N + 1 - x)`. The omitted term of index `N + k`
(`k ≥ 1`) is at most `(x^N/N!) (x/(N+1))^k`, since `N! (N+1)^k ≤ (N+k)!`. -/
theorem exp_le_sum_add_geometric {x : ℝ} (hx : 0 ≤ x) {N : ℕ} (hxN : x < N + 1) :
    Real.exp x ≤ ∑ j ∈ range (N + 1), x ^ j / j ! + x ^ N / N ! * (x / (N + 1 - x)) := by
  set r := x / (N + 1) with hr
  have hN : (0 : ℝ) < N + 1 := by positivity
  have hr0 : 0 ≤ r := by positivity
  have hr1 : r < 1 := by rw [hr, div_lt_one hN]; exact hxN
  have hterm : ∀ i : ℕ,
      x ^ (i + (N + 1)) / ((i + (N + 1))! : ℝ) ≤ x ^ N / N ! * r ^ i * r := by
    intro i
    have hfac : (N ! : ℝ) * (N + 1) ^ (i + 1) ≤ ((i + (N + 1))! : ℝ) := by
      rw [show i + (N + 1) = N + (i + 1) by ring]
      exact_mod_cast Nat.factorial_mul_pow_le_factorial
    have heq : x ^ N / N ! * r ^ i * r = x ^ (i + (N + 1)) / (N ! * (N + 1) ^ (i + 1)) := by
      rw [hr, div_pow]
      field_simp
      ring
    rw [heq]
    exact div_le_div_of_nonneg_left (by positivity) (by positivity) hfac
  have hgeom : Summable fun i : ℕ ↦ x ^ N / N ! * r ^ i * r :=
    ((summable_geometric_of_lt_one hr0 hr1).mul_left _).mul_right _
  have htail : expTail x (N + 1) ≤ x ^ N / N ! * (x / (N + 1 - x)) := by
    rw [expTail_eq_tsum]
    refine ((summable_tail x (N + 1)).tsum_mono hgeom hterm).trans_eq ?_
    rw [tsum_mul_right, tsum_mul_left, tsum_geometric_of_lt_one hr0 hr1, hr]
    have : (N : ℝ) + 1 - x ≠ 0 := by linarith
    field_simp
  have : Real.exp x = ∑ j ∈ range (N + 1), x ^ j / j ! + expTail x (N + 1) := by
    rw [expTail]; ring
  linarith

/-- The clean enclosure `s ≤ e^x ≤ s + t` with `s = ∑_{j ≤ N} x^j/j!` and `t = x^N/N!`,
valid when `2x ≤ N + 1` (then the geometric tail factor `x/(N + 1 - x)` is at most `1`). -/
theorem exp_mem_Icc_taylor {x : ℝ} (hx : 0 ≤ x) {N : ℕ} (hxN : 2 * x ≤ N + 1) :
    ∑ j ∈ range (N + 1), x ^ j / j ! ≤ Real.exp x ∧
      Real.exp x ≤ ∑ j ∈ range (N + 1), x ^ j / j ! + x ^ N / N ! := by
  refine ⟨sum_le_exp hx N, ?_⟩
  have hlt : x < N + 1 := by
    rcases hx.lt_or_eq with h | h
    · linarith
    · rw [← h]; positivity
  refine (exp_le_sum_add_geometric hx hlt).trans ?_
  have hq : x / (N + 1 - x) ≤ 1 := by
    rw [div_le_one (by linarith)]; linarith
  have ht : 0 ≤ x ^ N / (N ! : ℝ) := by positivity
  nlinarith

/-- The enclosure of `thermal72.tex` (`N = 100`, `0 ≤ x < 15`): `s ≤ e^x ≤ s + t`. -/
theorem exp_mem_Icc_taylor_100 {x : ℝ} (hx : 0 ≤ x) (hx15 : x < 15) :
    ∑ j ∈ range 101, x ^ j / j ! ≤ Real.exp x ∧
      Real.exp x ≤ ∑ j ∈ range 101, x ^ j / j ! + x ^ 100 / (100 ! : ℝ) :=
  exp_mem_Icc_taylor hx (N := 100) (by push_cast; linarith)

/-! ### The matrix Fox–Glynn estimate -/

section PoissonMatrix

open Matrix NormedSpace CertifiedGaps.Spectral

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Fox–Glynn for a Hermitian contraction.** If the eigenvalues of the Hermitian matrix `R`
lie in `[-1, 1]` and `m ≥ 0`, then for every vector `v`
`‖(e^{m(R - 1)} - ∑_{j ≤ J} p_j R^j) v‖ ≤ η(m, J) ‖v‖`, stated with squared norms. Both
operators are functions of `R`, so the scalar estimate `abs_exp_sub_poissonWeight_le` on each
eigenvalue bounds their difference. -/
theorem nsq_exp_sub_poisson_mulVec_le {R : Matrix n n ℂ} (hR : R.IsHermitian)
    (hspec : ∀ i, |hR.eigenvalues i| ≤ 1) {m : ℚ} (hm : 0 ≤ m) (J : ℕ) (v : n → ℂ) :
    nsq ((exp ((m : ℂ) • (R - 1)) -
        ∑ j ∈ range (J + 1), (poissonWeight m J j : ℂ) • R ^ j) *ᵥ v)
      ≤ poissonError m J ^ 2 * nsq v := by
  have hexp : exp ((m : ℂ) • (R - 1)) = hR.cfc (fun x ↦ Real.exp (m * (x - 1))) := by
    have hR1 : R - 1 = hR.cfc (fun x ↦ x - 1) := by
      conv_lhs => rw [← hermFun_id hR, ← hermFun_one hR]
      rw [hermFun_sub]
      rfl
    rw [hR1, show (m : ℂ) = ((m : ℝ) : ℂ) by push_cast; rfl, exp_smul_hermFun]
  have hpoly : ∑ j ∈ range (J + 1), (poissonWeight m J j : ℂ) • R ^ j =
      hR.cfc (fun x ↦ ∑ j ∈ range (J + 1), (poissonWeight m J j : ℝ) * x ^ j) := by
    rw [← hermFun_sum]
    refine sum_congr rfl fun j _ ↦ ?_
    rw [hermFun_pow hR j, ← hermFun_smul]
    congr 1
  rw [hexp, hpoly, hermFun_sub]
  exact nsq_hermFun_mulVec_le hR _ (fun i ↦ abs_exp_sub_poissonWeight_le hm (hspec i) J) v

/-- The same estimate under the contraction hypothesis `|⟨x, R x⟩| ≤ ‖x‖²`. -/
theorem nsq_exp_sub_poisson_mulVec_le_of_qform {R : Matrix n n ℂ} (hR : R.IsHermitian)
    (hcon : ∀ x, |qform R x| ≤ nsq x) {m : ℚ} (hm : 0 ≤ m) (J : ℕ) (v : n → ℂ) :
    nsq ((exp ((m : ℂ) • (R - 1)) -
        ∑ j ∈ range (J + 1), (poissonWeight m J j : ℂ) • R ^ j) *ᵥ v)
      ≤ poissonError m J ^ 2 * nsq v :=
  nsq_exp_sub_poisson_mulVec_le hR (abs_eigenvalues_le_one_of_qform hR hcon) hm J v

/-- The instance `(m, J) = (98, 280)`: `‖(e^{98(R - 1)} - ∑_{j ≤ 280} p_j R^j) v‖ ≤
2^{-84} ‖v‖`. -/
theorem nsq_exp_sub_poisson_98_280 {R : Matrix n n ℂ} (hR : R.IsHermitian)
    (hspec : ∀ i, |hR.eigenvalues i| ≤ 1) (v : n → ℂ) :
    nsq ((exp ((98 : ℂ) • (R - 1)) -
        ∑ j ∈ range 281, (poissonWeight 98 280 j : ℂ) • R ^ j) *ᵥ v)
      ≤ (1 / 2 ^ 84) ^ 2 * nsq v := by
  have h := nsq_exp_sub_poisson_mulVec_le hR hspec (m := 98) (by norm_num) 280 v
  simp only [Rat.cast_ofNat] at h
  refine h.trans (mul_le_mul_of_nonneg_right ?_ (nsq_nonneg v))
  exact pow_le_pow_left₀ (poissonError_nonneg (by norm_num) _) poissonError_98_280.le 2

/-- The instance `(m, J) = (84, 256)`: `‖(e^{84(R - 1)} - ∑_{j ≤ 256} p_j R^j) v‖ ≤
2^{-116} ‖v‖`. -/
theorem nsq_exp_sub_poisson_84_256 {R : Matrix n n ℂ} (hR : R.IsHermitian)
    (hspec : ∀ i, |hR.eigenvalues i| ≤ 1) (v : n → ℂ) :
    nsq ((exp ((84 : ℂ) • (R - 1)) -
        ∑ j ∈ range 257, (poissonWeight 84 256 j : ℂ) • R ^ j) *ᵥ v)
      ≤ (1 / 2 ^ 116) ^ 2 * nsq v := by
  have h := nsq_exp_sub_poisson_mulVec_le hR hspec (m := 84) (by norm_num) 256 v
  simp only [Rat.cast_ofNat] at h
  refine h.trans (mul_le_mul_of_nonneg_right ?_ (nsq_nonneg v))
  exact pow_le_pow_left₀ (poissonError_nonneg (by norm_num) _) poissonError_84_256.le 2

end PoissonMatrix

end CertifiedGaps.Trace
