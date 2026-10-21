/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Basic.Real.Basic
import Mathlib.Data.Rat.Cast.Order
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Polynomial filters: eigenvalue caps from finitely many moments

Let `x i` (`i ∈ S`) be real numbers, such as the eigenvalues of a self-adjoint transfer
operator in one symmetry sector, whose power sums `M k = ∑ i ∈ S, x i ^ k` are known to within
`e` of rational centers `m k`. A polynomial `F = ∑ₖ fₖ Xᵏ` that is nonnegative on all of `ℝ`
caps every single `x i`: each term of `∑ i ∈ S, F (x i) = ∑ₖ fₖ M k` is nonnegative, so
`F (x i) ≤ ∑ₖ fₖ M k ≤ B := ∑ₖ (fₖ m k + |fₖ| e)`. If moreover `F > B` on both rays `|x| ≥ v`,
then `|x i| < v`. This is the filter step of the spin-one seed (Appendix B of the certificate
note, eq. (c-filter); the term is from OpenAI's *The periodic spin-one Haldane gap*).

Polynomials are lists of rational coefficients, constant term first, so that every finite
condition is decidable and a concrete filter is checked by `decide`:

* `CapCheck p v B`: the Taylor coefficients of `p` at `v` and of `p (-·)` at `v` are
  nonnegative apart from the constant term, which exceeds `B`. This gives `F > B` on both rays.
* `momentBound p m e`: the rational upper bound `B` on `∑ₖ fₖ M k` from centers `m k ± e`.
* `pairingBound p lo hi`: the same bound from optional lower and upper bounds on the moments.
  A term `fₖ M k` needs only `hi k` if `fₖ > 0` and only `lo k` if `fₖ < 0`, so one-sided
  tables suffice; the bound is `none` if a needed side is missing.

## Main results

* `CertifiedGaps.PolynomialFilter.abs_lt_of_filter`: the cap `|x i| < v` from a filter.
* `CertifiedGaps.PolynomialFilter.abs_lt_of_filterPoly`: the same for the filters
  `F = X^(2a) q²`, which are nonnegative by construction.
-/

open Finset

namespace CertifiedGaps.PolynomialFilter

/-! ### Rational polynomials as coefficient lists -/

/-- The value at a real point of the rational polynomial with coefficients `p` (constant term
first), by Horner's rule. -/
def eval : List ℚ → ℝ → ℝ
  | [], _ => 0
  | a :: p, x => a + x * eval p x

@[simp]
theorem eval_nil (x : ℝ) : eval [] x = 0 := rfl

@[simp]
theorem eval_cons (a : ℚ) (p : List ℚ) (x : ℝ) : eval (a :: p) x = a + x * eval p x := rfl

/-- The sum of two coefficient lists. -/
def add : List ℚ → List ℚ → List ℚ
  | [], q => q
  | p, [] => p
  | a :: p, b :: q => (a + b) :: add p q

theorem eval_add : ∀ (p q : List ℚ) (x : ℝ), eval (add p q) x = eval p x + eval q x
  | [], q, x => by simp [add]
  | a :: p, [], x => by simp [add]
  | a :: p, b :: q, x => by simp only [add, eval_cons, eval_add p q x]; push_cast; ring

/-- A rational multiple of a coefficient list. -/
def smul (c : ℚ) : List ℚ → List ℚ
  | [] => []
  | a :: p => c * a :: smul c p

theorem eval_smul (c : ℚ) : ∀ (p : List ℚ) (x : ℝ), eval (smul c p) x = c * eval p x
  | [], x => by simp [smul]
  | a :: p, x => by simp only [smul, eval_cons, eval_smul c p x]; push_cast; ring

/-- The product of two coefficient lists. -/
def mul : List ℚ → List ℚ → List ℚ
  | [], _ => []
  | a :: p, q => add (smul a q) (0 :: mul p q)

theorem eval_mul : ∀ (p q : List ℚ) (x : ℝ), eval (mul p q) x = eval p x * eval q x
  | [], q, x => by simp [mul]
  | a :: p, q, x => by
    simp only [mul, eval_add, eval_smul, eval_cons, eval_mul p q x]; push_cast; ring

/-- The Taylor shift: the coefficients of `y ↦ p (v + y)`. -/
def shift (v : ℚ) : List ℚ → List ℚ
  | [] => []
  | a :: p => add [a] (mul [v, 1] (shift v p))

theorem eval_shift (v : ℚ) : ∀ (p : List ℚ) (y : ℝ), eval (shift v p) y = eval p (v + y)
  | [], y => by simp [shift]
  | a :: p, y => by
    simp only [shift, eval_add, eval_mul, eval_cons, eval_nil, eval_shift v p y]; push_cast; ring

/-- The reflection: the coefficients of `x ↦ p (-x)`. -/
def reflect : List ℚ → List ℚ
  | [] => []
  | a :: p => a :: smul (-1) (reflect p)

theorem eval_reflect : ∀ (p : List ℚ) (x : ℝ), eval (reflect p) x = eval p (-x)
  | [], x => by simp [reflect]
  | a :: p, x => by
    simp only [reflect, eval_cons, eval_smul, eval_reflect p x]; push_cast; ring

/-- `X ^ n * r`: the list `r` after `n` zero coefficients. -/
theorem eval_replicate_append (r : List ℚ) (x : ℝ) :
    ∀ n : ℕ, eval (List.replicate n 0 ++ r) x = x ^ n * eval r x
  | 0 => by simp
  | n + 1 => by
    rw [List.replicate_succ, List.cons_append, eval_cons, eval_replicate_append r x n]
    push_cast; ring

/-! ### Both rays from Taylor coefficients -/

/-- A polynomial with nonnegative coefficients is nonnegative on `[0, ∞)`. -/
theorem eval_nonneg_of_coeffs :
    ∀ (p : List ℚ), (∀ c ∈ p, 0 ≤ c) → ∀ y : ℝ, 0 ≤ y → 0 ≤ eval p y
  | [], _, y, _ => by simp
  | a :: p, h, y, hy => by
    have ha : (0 : ℝ) ≤ a := by exact_mod_cast h a (by simp)
    have hp := eval_nonneg_of_coeffs p (fun c hc ↦ h c (by simp [hc])) y hy
    rw [eval_cons]; positivity

/-- A polynomial whose coefficients are nonnegative apart from the constant term is at least
that constant term on `[0, ∞)`. -/
theorem headD_le_eval :
    ∀ (g : List ℚ), (∀ c ∈ g.tail, 0 ≤ c) → ∀ y : ℝ, 0 ≤ y → (g.headD 0 : ℝ) ≤ eval g y
  | [], _, y, _ => by simp
  | a :: t, h, y, hy => by
    have ht := eval_nonneg_of_coeffs t h y hy
    simp only [List.headD_cons, eval_cons, le_add_iff_nonneg_right]
    positivity

/-- The ray condition at `v`: the Taylor coefficients of `p` at `v` are nonnegative apart from
the constant term `p v`, which exceeds `B`. It gives `p x > B` for all `x ≥ v`. -/
def RayCheck (p : List ℚ) (v B : ℚ) : Prop :=
  B < (shift v p).headD 0 ∧ ∀ c ∈ (shift v p).tail, 0 ≤ c

instance (p : List ℚ) (v B : ℚ) : Decidable (RayCheck p v B) := by
  unfold RayCheck; infer_instance

/-- The cap condition at `v`: the ray condition for `p` and for `p (-·)`. It gives `p x > B`
for all `|x| ≥ v`. -/
def CapCheck (p : List ℚ) (v B : ℚ) : Prop :=
  RayCheck p v B ∧ RayCheck (reflect p) v B

instance (p : List ℚ) (v B : ℚ) : Decidable (CapCheck p v B) := by
  unfold CapCheck; infer_instance

theorem lt_eval_of_rayCheck {p : List ℚ} {v B : ℚ} (h : RayCheck p v B) {x : ℝ}
    (hx : (v : ℝ) ≤ x) : (B : ℝ) < eval p x := by
  have hle := headD_le_eval (shift v p) h.2 (x - v) (sub_nonneg.2 hx)
  rw [eval_shift, add_sub_cancel] at hle
  exact lt_of_lt_of_le (by exact_mod_cast h.1) hle

theorem lt_eval_of_capCheck {p : List ℚ} {v B : ℚ} (h : CapCheck p v B) {x : ℝ}
    (hx : (v : ℝ) ≤ |x|) : (B : ℝ) < eval p x := by
  rcases le_or_gt 0 x with h0 | h0
  · exact lt_eval_of_rayCheck h.1 (by rwa [abs_of_nonneg h0] at hx)
  · have := lt_eval_of_rayCheck h.2 (x := -x) (by rwa [abs_of_neg h0] at hx)
    rwa [eval_reflect, neg_neg] at this

/-! ### Moments -/

/-- The pairing `∑ₖ pₖ M k` of a coefficient list with a moment sequence. -/
def pairing : List ℚ → (ℕ → ℝ) → ℝ
  | [], _ => 0
  | a :: p, M => a * M 0 + pairing p (fun k ↦ M (k + 1))

/-- A weighted sum of values of `p` is the pairing of `p` with the weighted power sums. -/
theorem sum_mul_eval {ι : Type*} (S : Finset ι) (x : ι → ℝ) :
    ∀ (p : List ℚ) (w : ι → ℝ),
      ∑ i ∈ S, w i * eval p (x i) = pairing p (fun k ↦ ∑ i ∈ S, w i * x i ^ k)
  | [], w => by simp [pairing]
  | a :: p, w => by
    have ih := sum_mul_eval S x p (fun i ↦ w i * x i)
    have h1 : (fun k ↦ ∑ i ∈ S, w i * x i ^ (k + 1)) = fun k ↦ ∑ i ∈ S, w i * x i * x i ^ k := by
      funext k; refine sum_congr rfl fun i _ ↦ ?_; ring
    simp only [pairing, pow_zero, mul_one]
    rw [h1, ← ih, mul_sum, ← sum_add_distrib]
    refine sum_congr rfl fun i _ ↦ ?_
    rw [eval_cons]; ring

/-- The sum of values of `p` is the pairing of `p` with the power sums. -/
theorem sum_eval {ι : Type*} (S : Finset ι) (x : ι → ℝ) (p : List ℚ) :
    ∑ i ∈ S, eval p (x i) = pairing p (fun k ↦ ∑ i ∈ S, x i ^ k) := by
  simpa using sum_mul_eval S x p 1

/-- The rational upper bound `∑ₖ (pₖ m k + |pₖ| e)` on the pairing of `p` with any moment
sequence within `e` of `m`. -/
def momentBound : List ℚ → (ℕ → ℚ) → ℚ → ℚ
  | [], _, _ => 0
  | a :: p, m, e => a * m 0 + |a| * e + momentBound p (fun k ↦ m (k + 1)) e

/-- The moment orders that a filter uses: the indices of its nonzero coefficients. -/
def support (p : List ℚ) : List ℕ :=
  (List.range p.length).filter fun k ↦ p.getD k 0 ≠ 0

theorem mem_support {p : List ℚ} {k : ℕ} (h : p.getD k 0 ≠ 0) : k ∈ support p := by
  have hk : k < p.length := by
    by_contra hk
    exact h (by rw [List.getD_eq_getElem?_getD, List.getElem?_eq_none (not_lt.1 hk)]; rfl)
  simpa [support, hk, List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hk] using h

/-- Only the moments of the orders in the support enter the pairing. -/
theorem pairing_le_momentBound :
    ∀ (p : List ℚ) (M : ℕ → ℝ) (m : ℕ → ℚ) (e : ℚ),
      (∀ k, p.getD k 0 ≠ 0 → |M k - m k| ≤ e) → pairing p M ≤ momentBound p m e
  | [], _, _, _, _ => by simp [pairing, momentBound]
  | a :: p, M, m, e, h => by
    have ih := pairing_le_momentBound p (fun k ↦ M (k + 1)) (fun k ↦ m (k + 1)) e
      (fun k hk ↦ h (k + 1) (by simpa using hk))
    have h0 : (a : ℝ) * M 0 ≤ a * m 0 + |(a : ℝ)| * e := by
      rcases eq_or_ne a 0 with ha | ha
      · simp [ha]
      · have hM := h 0 (by simpa using ha)
        have := (le_abs_self ((a : ℝ) * (M 0 - m 0))).trans_eq (abs_mul _ _)
        nlinarith [abs_nonneg (a : ℝ)]
    simp only [pairing, momentBound]
    push_cast
    linarith

/-! ### Moments with one-sided bounds -/

/-- The bound on one term `a M`: `a hi` if `a > 0`, `a lo` if `a < 0`, and `0` if `a = 0`. -/
def termBound (a : ℚ) (lo hi : Option ℚ) : Option ℚ :=
  if 0 < a then hi.map (a * ·) else if a < 0 then lo.map (a * ·) else some 0

/-- The upper bound on the pairing `∑ₖ pₖ M k` from optional bounds `lo k ≤ M k ≤ hi k`. -/
def pairingBound : List ℚ → (ℕ → Option ℚ) → (ℕ → Option ℚ) → Option ℚ
  | [], _, _ => some 0
  | a :: p, lo, hi => (termBound a (lo 0) (hi 0)).bind fun t ↦
      (pairingBound p (fun k ↦ lo (k + 1)) (fun k ↦ hi (k + 1))).map (t + ·)

/-- The moments `M` obey the optional bounds `lo`, `hi` wherever they are given. -/
def MomentBounds (M : ℕ → ℝ) (lo hi : ℕ → Option ℚ) : Prop :=
  ∀ k b, (lo k = some b → (b : ℝ) ≤ M k) ∧ (hi k = some b → M k ≤ b)

theorem termBound_spec {a : ℚ} {lo hi : Option ℚ} {t : ℚ} (h : termBound a lo hi = some t)
    {M : ℝ} (hlo : ∀ b, lo = some b → (b : ℝ) ≤ M) (hhi : ∀ b, hi = some b → M ≤ b) :
    (a : ℝ) * M ≤ t := by
  unfold termBound at h
  split_ifs at h with ha ha'
  · obtain ⟨b, hb, rfl⟩ := Option.map_eq_some_iff.1 h
    push_cast
    exact mul_le_mul_of_nonneg_left (hhi b hb) (by exact_mod_cast ha.le)
  · obtain ⟨b, hb, rfl⟩ := Option.map_eq_some_iff.1 h
    push_cast
    exact mul_le_mul_of_nonpos_left (hlo b hb) (by exact_mod_cast ha'.le)
  · cases h
    have : a = 0 := le_antisymm (not_lt.1 ha) (not_lt.1 ha')
    simp [this]

/-- The pairing is at most `pairingBound` when the moments obey the bounds. -/
theorem pairing_le_pairingBound :
    ∀ (p : List ℚ) (lo hi : ℕ → Option ℚ) (M : ℕ → ℝ) {B : ℚ},
      pairingBound p lo hi = some B → MomentBounds M lo hi → pairing p M ≤ B
  | [], _, _, _, B, h, _ => by
    simp only [pairingBound, Option.some.injEq] at h; simp [pairing, ← h]
  | a :: p, lo, hi, M, B, h, hM => by
    obtain ⟨t, ht, h'⟩ := Option.bind_eq_some_iff.1 h
    obtain ⟨r, hr, rfl⟩ := Option.map_eq_some_iff.1 h'
    have h1 := termBound_spec ht (fun b hb ↦ (hM 0 b).1 hb) (fun b hb ↦ (hM 0 b).2 hb)
    have h2 := pairing_le_pairingBound p _ _ (fun k ↦ M (k + 1)) hr fun k b ↦ hM (k + 1) b
    simp only [pairing]
    push_cast
    linarith

/-- **Polynomial filter with one-sided moment bounds.** If `p ≥ 0` on `ℝ`, the power sums obey
optional bounds `lo`, `hi`, `pairingBound p lo hi = some B` and `p > B` on `|x| ≥ v`, then every
`|x i| < v`. -/
theorem abs_lt_of_filter_bounds {ι : Type*} (S : Finset ι) (x : ι → ℝ) {p : List ℚ}
    {lo hi : ℕ → Option ℚ} {B v : ℚ} (hp : ∀ y, 0 ≤ eval p y)
    (hM : MomentBounds (fun k ↦ ∑ i ∈ S, x i ^ k) lo hi) (hB : pairingBound p lo hi = some B)
    (hcap : CapCheck p v B) {i : ι} (hiS : i ∈ S) : |x i| < v := by
  by_contra h
  have h1 := lt_eval_of_capCheck hcap (not_lt.1 h)
  have h2 : eval p (x i) ≤ ∑ j ∈ S, eval p (x j) := single_le_sum (fun j _ ↦ hp (x j)) hiS
  have h3 := pairing_le_pairingBound p lo hi _ hB hM
  rw [← sum_eval] at h3
  linarith

/-- The cap condition for an optional bound: the bound exists and passes `CapCheck`. -/
def CapCheckOpt (p : List ℚ) (v : ℚ) : Option ℚ → Prop
  | some B => CapCheck p v B
  | none => False

instance (p : List ℚ) (v : ℚ) : (o : Option ℚ) → Decidable (CapCheckOpt p v o)
  | some B => inferInstanceAs (Decidable (CapCheck p v B))
  | none => isFalse id

theorem CapCheckOpt.exists {p : List ℚ} {v : ℚ} {o : Option ℚ} (h : CapCheckOpt p v o) :
    ∃ B, o = some B ∧ CapCheck p v B := by
  cases o with
  | none => exact h.elim
  | some B => exact ⟨B, rfl, h⟩

/-! ### The cap -/

/-- **Polynomial filter.** If `p ≥ 0` on `ℝ`, the power sums of the orders that `p` uses are
within `e` of `m`, and `p > momentBound p m e` on `|x| ≥ v`, then every `|x i| < v`. -/
theorem abs_lt_of_filter {ι : Type*} (S : Finset ι) (x : ι → ℝ) {p : List ℚ} {m : ℕ → ℚ}
    {e v : ℚ} (hp : ∀ y, 0 ≤ eval p y)
    (hm : ∀ k ∈ support p, |∑ i ∈ S, x i ^ k - m k| ≤ e)
    (hcap : CapCheck p v (momentBound p m e)) {i : ι} (hi : i ∈ S) : |x i| < v := by
  by_contra h
  have h1 := lt_eval_of_capCheck hcap (not_lt.1 h)
  have h2 : eval p (x i) ≤ ∑ j ∈ S, eval p (x j) := single_le_sum (fun j _ ↦ hp (x j)) hi
  have h3 := pairing_le_momentBound p _ m e fun k hk ↦ hm k (mem_support hk)
  rw [← sum_eval] at h3
  linarith

/-- The filter `X^(2a) q²`, nonnegative on `ℝ` by construction. -/
def filterPoly (a : ℕ) (q : List ℚ) : List ℚ :=
  List.replicate (2 * a) 0 ++ mul q q

theorem eval_filterPoly_nonneg (a : ℕ) (q : List ℚ) (y : ℝ) : 0 ≤ eval (filterPoly a q) y := by
  rw [filterPoly, eval_replicate_append, eval_mul, pow_mul']
  exact mul_nonneg (by positivity) (mul_self_nonneg _)

/-- **Polynomial filter** for `X^(2a) q²`. -/
theorem abs_lt_of_filterPoly {ι : Type*} (S : Finset ι) (x : ι → ℝ) {a : ℕ} {q : List ℚ}
    {m : ℕ → ℚ} {e v : ℚ}
    (hm : ∀ k ∈ support (filterPoly a q), |∑ i ∈ S, x i ^ k - m k| ≤ e)
    (hcap : CapCheck (filterPoly a q) v (momentBound (filterPoly a q) m e)) {i : ι}
    (hi : i ∈ S) : |x i| < v :=
  abs_lt_of_filter S x (eval_filterPoly_nonneg a q) hm hcap hi

end CertifiedGaps.PolynomialFilter
