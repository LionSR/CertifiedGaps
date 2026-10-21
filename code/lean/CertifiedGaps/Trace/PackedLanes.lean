/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

/-!
# Packed lanes

A vector `(f 0, …, f (d - 1))` of naturals is stored as the single natural

  `pack w d f = ∑_{i < d} f i · 2^(w i)`,

lane `i` occupying the bits `w i, …, w i + w - 1`. The Lean kernel evaluates natural-number
arithmetic and bitwise operations with GMP, so one operation on a packed number acts on all
lanes in a single kernel step. This file proves the lane-wise meaning of the operations used by
the packed Horner recursion (`Trace/PackedHorner.lean`):

* sums, scalar multiples and lane-wise dominated differences of packed numbers are packed
  (`pack_add`, `pack_mul`, `pack_sub`); these hold without any bound on the lanes;
* when every lane is below `2^w` (`Canonical`), the bits of `pack w d f` are those of the lanes
  (`testBit_pack`), hence lane extraction (`pack_lane`), lane-wise `&&&` (`pack_and`) and the
  lane-wise right shift `(P >>> S) &&& (2^(w-S) - 1, …)` (`pack_shiftRight_and`);
* `pack_const`, `pack_single`: the constant vectors `c · (1, …, 1)` (computed through the repunit
  `(2^(w d) - 1)/(2^w - 1)`) and the unit vectors `t · 2^(w k)`;
* `lanes_mem_of_and_eq`: one `&&&` test checks `2^b ≤ f i < 2^(b + 1)` in every lane.
-/

namespace CertifiedGaps.Trace.Packed

open Finset

/-- The packed number `∑_{i < d} f i · 2^(w i)` of the lanes `f 0, …, f (d - 1)`. -/
def pack (w d : ℕ) (f : ℕ → ℕ) : ℕ := ∑ i ∈ range d, f i * 2 ^ (w * i)

/-- Every lane is below `2^w`, so that lanes do not overlap. -/
def Canonical (w d : ℕ) (f : ℕ → ℕ) : Prop := ∀ i < d, f i < 2 ^ w

/-- The repunit `∑_{i < d} 2^(w i) = (2^(w d) - 1)/(2^w - 1)`, in the form the kernel computes. -/
def repunit (w d : ℕ) : ℕ := (2 ^ (w * d) - 1) / (2 ^ w - 1)

variable {w d : ℕ}

@[simp] theorem pack_zero (f : ℕ → ℕ) : pack w 0 f = 0 := by simp [pack]

/-- Peeling off the lowest lane. -/
theorem pack_succ (f : ℕ → ℕ) :
    pack w (d + 1) f = f 0 + 2 ^ w * pack w d (fun i ↦ f (i + 1)) := by
  rw [pack, sum_range_succ', pack, mul_sum, add_comm]
  congr 1
  · simp
  · refine sum_congr rfl fun i _ ↦ ?_
    rw [mul_add, mul_one, pow_add]; ring

theorem pack_congr {f g : ℕ → ℕ} (h : ∀ i < d, f i = g i) : pack w d f = pack w d g :=
  sum_congr rfl fun i hi ↦ by rw [h i (mem_range.1 hi)]

theorem pack_add (f g : ℕ → ℕ) : pack w d (fun i ↦ f i + g i) = pack w d f + pack w d g := by
  simp only [pack, add_mul, sum_add_distrib]

theorem pack_mul (c : ℕ) (f : ℕ → ℕ) : pack w d (fun i ↦ c * f i) = c * pack w d f := by
  simp only [pack, mul_sum, mul_assoc]

/-- A lane-wise dominated difference is the difference of the packed numbers. -/
theorem pack_sub {f g : ℕ → ℕ} (h : ∀ i < d, g i ≤ f i) :
    pack w d (fun i ↦ f i - g i) = pack w d f - pack w d g := by
  have e : pack w d f = pack w d (fun i ↦ f i - g i) + pack w d g := by
    rw [← pack_add]; exact pack_congr fun i hi ↦ (Nat.sub_add_cancel (h i hi)).symm
  omega

theorem pack_sum {ι : Type*} (s : Finset ι) (f : ι → ℕ → ℕ) :
    pack w d (fun i ↦ ∑ a ∈ s, f a i) = ∑ a ∈ s, pack w d (f a) := by
  simp only [pack, sum_mul]
  exact sum_comm

/-- The unit vector `t e_k` is the shift `t · 2^(w k)`. -/
theorem pack_single {k : ℕ} (hk : k < d) (t : ℕ) :
    pack w d (fun i ↦ if i = k then t else 0) = t * 2 ^ (w * k) := by
  rw [pack, sum_eq_single k (fun i _ hi ↦ by simp [hi]) (fun h ↦ (h (mem_range.2 hk)).elim)]
  simp

/-- The constant vector `c (1, …, 1)` is `c` times the repunit. -/
theorem pack_const (hw : 0 < w) (c : ℕ) : pack w d (fun _ ↦ c) = c * repunit w d := by
  have h1 : pack w d (fun _ ↦ c) = c * pack w d (fun _ ↦ 1) := by
    rw [← pack_mul]; simp
  have h2 : 2 ≤ 2 ^ w := by
    calc 2 = 2 ^ 1 := by norm_num
      _ ≤ 2 ^ w := Nat.pow_le_pow_right (by norm_num) hw
  have h3 : pack w d (fun _ ↦ 1) = repunit w d := by
    rw [pack, repunit, pow_mul, ← Nat.geomSum_eq h2 d]
    exact sum_congr rfl fun i _ ↦ by rw [one_mul, pow_mul]
  rw [h1, h3]

/-! ### Canonical lanes -/

theorem Canonical.tail {f : ℕ → ℕ} (hf : Canonical w (d + 1) f) :
    Canonical w d (fun i ↦ f (i + 1)) := fun i hi ↦ hf (i + 1) (by omega)

/-- **The bits of a canonical packed number** are those of its lanes: bit `n` is bit `n % w` of
lane `n / w` (and `0` beyond the last lane). -/
theorem testBit_pack (hw : 0 < w) {f : ℕ → ℕ} (hf : Canonical w d f) (n : ℕ) :
    (pack w d f).testBit n = (decide (n / w < d) && (f (n / w)).testBit (n % w)) := by
  induction d generalizing f n with
  | zero => simp
  | succ d ih =>
    rw [pack_succ, add_comm, Nat.testBit_two_pow_mul_add _ (hf 0 (by omega))]
    by_cases hn : n < w
    · rw [ite_eq_left hn, Nat.div_eq_of_lt hn, Nat.mod_eq_of_lt hn]
      simp
    · rw [ite_eq_right hn, ih hf.tail]
      obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le (not_lt.1 hn)
      rw [Nat.add_sub_cancel_left, Nat.add_div_left _ hw, Nat.add_mod_left]
      simp only [Nat.add_comm (m / w) 1]
      congr 1
      simp only [decide_eq_decide]; omega

/-- A canonical packed number has at most `w d` bits. -/
theorem pack_lt (hw : 0 < w) {f : ℕ → ℕ} (hf : Canonical w d f) : pack w d f < 2 ^ (w * d) := by
  refine Nat.lt_pow_two_of_testBit _ fun n hn ↦ ?_
  rw [testBit_pack hw hf]
  have : ¬ n / w < d := by
    rw [not_lt, Nat.le_div_iff_mul_le hw]; rw [mul_comm]; exact hn
  simp [this]

/-- **Lane extraction.** -/
theorem pack_lane (hw : 0 < w) {f : ℕ → ℕ} (hf : Canonical w d f) {i : ℕ} (hi : i < d) :
    (pack w d f >>> (w * i)) % 2 ^ w = f i := by
  refine Nat.eq_of_testBit_eq fun n ↦ ?_
  rw [Nat.testBit_mod_two_pow, Nat.testBit_shiftRight, testBit_pack hw hf]
  by_cases hn : n < w
  · have h1 : (w * i + n) / w = i := by
      rw [Nat.add_comm, Nat.add_mul_div_left _ _ hw, Nat.div_eq_of_lt hn, zero_add]
    have h2 : (w * i + n) % w = n := by
      rw [Nat.add_comm, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hn]
    simp [hn, h1, h2, hi]
  · have := Nat.testBit_lt_two_pow (lt_of_lt_of_le (hf i hi)
      (Nat.pow_le_pow_right (by norm_num) (not_lt.1 hn)))
    simp [hn, this]

/-- Two canonical packed numbers are equal only if their lanes are. -/
theorem lane_eq_of_pack_eq (hw : 0 < w) {f g : ℕ → ℕ} (hf : Canonical w d f)
    (hg : Canonical w d g) (h : pack w d f = pack w d g) {i : ℕ} (hi : i < d) : f i = g i := by
  rw [← pack_lane hw hf hi, h, pack_lane hw hg hi]

/-- **Lane-wise `&&&`.** -/
theorem pack_and (hw : 0 < w) {f g : ℕ → ℕ} (hf : Canonical w d f) (hg : Canonical w d g) :
    pack w d f &&& pack w d g = pack w d (fun i ↦ f i &&& g i) := by
  have hfg : Canonical w d (fun i ↦ f i &&& g i) := fun i hi ↦ Nat.and_lt_two_pow _ (hg i hi)
  refine Nat.eq_of_testBit_eq fun n ↦ ?_
  rw [Nat.testBit_and, testBit_pack hw hf, testBit_pack hw hg, testBit_pack hw hfg,
    Nat.testBit_and]
  cases decide (n / w < d) <;> simp

/-- **Lane-wise right shift.** Shifting a canonical packed number right by `S ≤ w` bits and
masking each lane to its low `w - S` bits shifts every lane right by `S`. -/
theorem pack_shiftRight_and (hw : 0 < w) {f : ℕ → ℕ} (hf : Canonical w d f) {S : ℕ}
    (hS : S ≤ w) :
    (pack w d f >>> S) &&& pack w d (fun _ ↦ 2 ^ (w - S) - 1) = pack w d (fun i ↦ f i >>> S) := by
  have hm : Canonical w d (fun _ ↦ 2 ^ (w - S) - 1) := fun _ _ ↦
    lt_of_lt_of_le (Nat.sub_lt (Nat.two_pow_pos _) one_pos)
      (Nat.pow_le_pow_right (by norm_num) (Nat.sub_le _ _))
  have hfs : Canonical w d (fun i ↦ f i >>> S) := fun i hi ↦
    lt_of_le_of_lt (by simp only [Nat.shiftRight_eq_div_pow]; exact Nat.div_le_self _ _)
      (hf i hi)
  refine Nat.eq_of_testBit_eq fun n ↦ ?_
  rw [Nat.testBit_and, Nat.testBit_shiftRight, testBit_pack hw hf, testBit_pack hw hm,
    testBit_pack hw hfs, Nat.testBit_two_pow_sub_one, Nat.testBit_shiftRight]
  set q := n / w
  set ρ := n % w
  have hρ : ρ < w := Nat.mod_lt _ hw
  have hn : n = ρ + w * q := (Nat.mod_add_div n w).symm
  by_cases hs : ρ < w - S
  · have h1 : (S + n) / w = q := by
      rw [hn, ← Nat.add_assoc, Nat.add_mul_div_left _ _ hw, Nat.div_eq_of_lt (by omega),
        zero_add]
    have h2 : (S + n) % w = S + ρ := by
      rw [hn, ← Nat.add_assoc, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt (by omega)]
    rw [h1, h2]
    by_cases hq : q < d <;> simp [hs, hq]
  · rw [decide_eq_false hs]
    by_cases hq : q < d
    · have := Nat.testBit_lt_two_pow (lt_of_lt_of_le (hf q hq)
        (Nat.pow_le_pow_right (by norm_num) (by omega : w ≤ S + ρ)))
      simp [this]
    · simp [hq]

/-- **One test for all lanes.** If every lane of `f` is below `2^w`, `b < w`, and keeping the
bits `b, …, w - 1` of every lane leaves exactly `2^b` in every lane, then
`2^b ≤ f i < 2^(b + 1)` for every lane. -/
theorem lanes_mem_of_and_eq (hw : 0 < w) {b : ℕ} (hb : b < w) {f : ℕ → ℕ}
    (hf : Canonical w d f)
    (h : pack w d f &&& pack w d (fun _ ↦ 2 ^ w - 2 ^ b) = pack w d (fun _ ↦ 2 ^ b)) :
    ∀ i < d, 2 ^ b ≤ f i ∧ f i < 2 ^ (b + 1) := by
  have hm : Canonical w d (fun _ ↦ 2 ^ w - 2 ^ b) := fun _ _ ↦
    Nat.sub_lt (Nat.two_pow_pos _) (Nat.two_pow_pos _)
  have hc : Canonical w d (fun _ ↦ 2 ^ b) := fun _ _ ↦ Nat.pow_lt_pow_right (by norm_num) hb
  rw [pack_and hw hf hm] at h
  have hfm : Canonical w d (fun i ↦ f i &&& (2 ^ w - 2 ^ b)) := fun i hi ↦
    Nat.and_lt_two_pow _ (hm i hi)
  intro i hi
  have hl := lane_eq_of_pack_eq hw hfm hc h hi
  -- The bits of the mask `2^w - 2^b = (2^(w-b) - 1) 2^b` are `b, …, w - 1`.
  have hmask : ∀ j, (2 ^ w - 2 ^ b).testBit j = decide (b ≤ j ∧ j < w) := by
    intro j
    have e : 2 ^ w - 2 ^ b = (2 ^ (w - b) - 1) * 2 ^ b := by
      rw [Nat.sub_mul, one_mul, ← pow_add, Nat.sub_add_cancel hb.le]
    rw [e, Nat.testBit_mul_two_pow, Nat.testBit_two_pow_sub_one]
    by_cases hj : b ≤ j
    · simp only [hj, decide_true, Bool.true_and, decide_eq_decide, true_and]
      omega
    · simp [hj]
  have hbit : ∀ j, (f i &&& (2 ^ w - 2 ^ b)).testBit j = (decide (j = b)) := by
    intro j
    rw [hl, Nat.testBit_two_pow]
    simp only [decide_eq_decide]
    exact eq_comm
  constructor
  · refine Nat.ge_two_pow_of_testBit ?_
    have := hbit b
    rw [Nat.testBit_and, hmask] at this
    simpa [hb] using this
  · refine Nat.lt_pow_two_of_testBit _ fun j hj ↦ ?_
    by_cases hjw : j < w
    · have := hbit j
      rw [Nat.testBit_and, hmask] at this
      have hjb : j ≠ b := by omega
      simpa [hjw, show b ≤ j by omega, hjb] using this
    · exact Nat.testBit_lt_two_pow (lt_of_lt_of_le (hf i hi)
        (Nat.pow_le_pow_right (by norm_num) (not_lt.1 hjw)))

end CertifiedGaps.Trace.Packed
