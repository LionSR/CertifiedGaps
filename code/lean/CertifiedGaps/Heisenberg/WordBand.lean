/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.SzBridge
import CertifiedGaps.Trace.PackedHorner

/-!
# Spin-one words as numbers, and exchange bonds as bands

The integer Horner recursions of both spin-one papers act on the `S^z` basis of `n` spins one,
the words `σ ∈ {0, 1, 2}ⁿ` (letter `σⱼ` is `m = σⱼ - 1`). Numbering the words in base `3` with
site `0` least significant, `wordNum σ = ∑ⱼ σⱼ 3ʲ`, turns the spin-flip hop
`σⱼ ↦ σⱼ - 1`, `σₖ ↦ σₖ + 1` of an exchange bond into a band of constant offset `3ᵏ - 3ʲ`. This
file collects the word numbering and the bands shared by the boundary-field recurrence
(`Heisenberg/BoundaryData.lean`, `Heisenberg/BoundaryThermal.lean`) and the periodic recurrence
(`Heisenberg/ThermalPacked.lean`).

## Main definitions

* `digit k j`, `letter k j`: the letter of site `j` of the word number `k`; `wordNum n` the
  numbering.
* `hopUp j k σ`: the word after the hop `σⱼ ↦ σⱼ - 1`, `σₖ ↦ σₖ + 1`, if it is allowed.
* `hopBand n j k c`: the band of that hop on the word numbers, with entry magnitude `c`.

## Main results

* `digit_wordNum`, `wordNum_update₂`: the digits and the place values of the numbering.
* `bondOp_apply_of_flip`, `bondOp_exchangeZ_apply`, `bondOp_twistedExchangeZ_apply`: the
  entries of an embedded exchange bond in terms of the two hops.
* `hopBand_term`: the band `hopBand n j k c` has the entry `c` exactly at the hop `j → k`.
-/

open Matrix

namespace CertifiedGaps.SpinOneHeisenberg

open Trace.Packed

/-! ### Word numbers -/

/-- Digit `j` of the word number `k`: the letter `wⱼ + 1 ∈ {0, 1, 2}`. -/
def digit (k j : ℕ) : ℕ := k / 3 ^ j % 3

/-- The letter `wⱼ = digit k j - 1 ∈ {-1, 0, 1}`. -/
def letter (k j : ℕ) : ℤ := (digit k j : ℤ) - 1

/-- The word number `e σ = ∑ᵢ σᵢ 3ⁱ` (site `0` least significant). -/
abbrev wordNum (n : ℕ) : (Fin n → Fin 3) ≃ Fin (3 ^ n) := finFunctionFinEquiv

theorem wordNum_val {n : ℕ} (σ : Fin n → Fin 3) :
    ((wordNum n σ : ℕ) : ℤ) = ∑ i, ((σ i : ℕ) : ℤ) * 3 ^ (i : ℕ) := by
  rw [wordNum, finFunctionFinEquiv_apply]; push_cast; rfl

/-- The digits of a word number are the letters. -/
theorem digit_wordNum {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) :
    digit (wordNum n σ) j = σ j := by
  have := finFunctionFinEquiv_symm_apply_val (m := 3) (wordNum n σ) j
  rw [wordNum, Equiv.symm_apply_apply] at this
  exact this.symm

theorem letter_wordNum {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) :
    letter (wordNum n σ) j = mz (σ j) := by
  rw [letter, digit_wordNum, mz]

/-- Changing two letters changes a weighted letter sum by the weighted differences. -/
theorem sum_update₂ {n : ℕ} (w : Fin n → ℤ) (σ : Fin n → Fin 3) {j k : Fin n} (hjk : j ≠ k)
    (a b : Fin 3) :
    ∑ i, ((Function.update (Function.update σ j a) k b i : ℕ) : ℤ) * w i =
      ∑ i, ((σ i : ℕ) : ℤ) * w i + ((a : ℤ) - σ j) * w j + ((b : ℤ) - σ k) * w k := by
  have h : ∀ i, ((Function.update (Function.update σ j a) k b i : ℕ) : ℤ) * w i =
      ((σ i : ℕ) : ℤ) * w i + (if i = j then ((a : ℤ) - σ j) * w j else 0) +
        (if i = k then ((b : ℤ) - σ k) * w k else 0) := by
    intro i
    by_cases hk : i = k
    · subst hk; simp [Ne.symm hjk]; ring
    · by_cases hj : i = j
      · subst hj; simp [hk]; ring
      · simp [hk, hj]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

/-- Changing one letter changes a weighted letter sum by the weighted difference. -/
theorem sum_update {n : ℕ} (w : Fin n → ℤ) (σ : Fin n → Fin 3) (j : Fin n) (a : Fin 3) :
    ∑ i, ((Function.update σ j a i : ℕ) : ℤ) * w i =
      ∑ i, ((σ i : ℕ) : ℤ) * w i + ((a : ℤ) - σ j) * w j := by
  have h : ∀ i, ((Function.update σ j a i : ℕ) : ℤ) * w i =
      ((σ i : ℕ) : ℤ) * w i + (if i = j then ((a : ℤ) - σ j) * w j else 0) := by
    intro i
    by_cases hj : i = j
    · subst hj; simp; ring
    · simp [hj]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

/-- Changing two letters changes the word number by the differences of the place values. -/
theorem wordNum_update₂ {n : ℕ} (σ : Fin n → Fin 3) {j k : Fin n} (hjk : j ≠ k) (a b : Fin 3) :
    ((wordNum n (Function.update (Function.update σ j a) k b) : ℕ) : ℤ) =
      wordNum n σ + ((a : ℤ) - σ j) * 3 ^ (j : ℕ) + ((b : ℤ) - σ k) * 3 ^ (k : ℕ) := by
  rw [wordNum_val, wordNum_val, sum_update₂ (fun i ↦ 3 ^ (i : ℕ)) σ hjk]

/-- Changing one letter changes the word number by the difference of the place value. -/
theorem wordNum_update {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) (a : Fin 3) :
    ((wordNum n (Function.update σ j a) : ℕ) : ℤ) =
      wordNum n σ + ((a : ℤ) - σ j) * 3 ^ (j : ℕ) := by
  rw [wordNum_val, wordNum_val, sum_update (fun i ↦ 3 ^ (i : ℕ)) σ j]

/-- A word with a nonzero letter at site `j` has number at least `3ʲ`. -/
theorem le_wordNum {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) (hj : 1 ≤ (σ j : ℕ)) :
    3 ^ (j : ℕ) ≤ (wordNum n σ : ℕ) := by
  have h := wordNum_val σ
  have h1 : ((σ j : ℕ) : ℤ) * 3 ^ (j : ℕ) ≤ ∑ i, ((σ i : ℕ) : ℤ) * 3 ^ (i : ℕ) :=
    Finset.single_le_sum (f := fun i ↦ ((σ i : ℕ) : ℤ) * 3 ^ (i : ℕ))
      (fun i _ ↦ by positivity) (Finset.mem_univ j)
  have h2 : (3 : ℤ) ^ (j : ℕ) ≤ ((σ j : ℕ) : ℤ) * 3 ^ (j : ℕ) :=
    le_mul_of_one_le_left (by positivity) (by exact_mod_cast hj)
  have : ((3 ^ (j : ℕ) : ℕ) : ℤ) ≤ ((wordNum n σ : ℕ) : ℤ) := by push_cast; linarith
  exact_mod_cast this

theorem wordNum_eq_iff {n : ℕ} {σ τ : Fin n → Fin 3} :
    ((wordNum n σ : ℕ) : ℤ) = wordNum n τ ↔ σ = τ := by
  constructor
  · intro h; exact (wordNum n).injective (Fin.ext (by exact_mod_cast h))
  · rintro rfl; rfl

/-! ### Sums over word numbers -/

theorem list_sum_range {M : Type*} [AddCommMonoid M] (f : ℕ → M) (m : ℕ) :
    ((List.range m).map f).sum = ∑ i ∈ Finset.range m, f i := by
  induction m with
  | zero => simp
  | succ m ih => rw [List.sum_range_succ, ih, Finset.sum_range_succ]

theorem getD_map_range {α : Type*} (f : ℕ → α) (d k : ℕ) (a : α) (hk : k < d) :
    ((List.range d).map f).getD k a = f k := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_range hk]; rfl

theorem sum_map_flatMap {α β M : Type*} [AddCommMonoid M] (L : List α) (f : α → List β)
    (g : β → M) : ((L.flatMap f).map g).sum = (L.map fun a ↦ ((f a).map g).sum).sum := by
  induction L with
  | nil => simp
  | cons a L ih => simp [List.flatMap_cons, List.map_append, List.sum_append, ih]

/-- A sum over words is a sum over word numbers. -/
theorem sum_words {n : ℕ} {M : Type*} [AddCommMonoid M] (f : ℕ → M) :
    ∑ σ : Fin n → Fin 3, f (wordNum n σ) = ∑ k ∈ Finset.range (3 ^ n), f k := by
  rw [(wordNum n).sum_comp (fun k ↦ f k), Fin.sum_univ_eq_sum_range]

/-! ### Hops between words -/

/-- The bond hop `σⱼ ↦ σⱼ - 1`, `σₖ ↦ σₖ + 1`, allowed when `σⱼ ≥ 1` and `σₖ ≤ 1`. -/
def hopUp {n : ℕ} (j k : Fin n) (σ : Fin n → Fin 3) : Option (Fin n → Fin 3) :=
  if h : 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1 then
    some (Function.update (Function.update σ j ⟨σ j - 1, by have := (σ j).isLt; omega⟩) k
      ⟨σ k + 1, by omega⟩)
  else none

theorem hopUp_eq_some {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (σ τ : Fin n → Fin 3) :
    hopUp j k σ = some τ ↔ (∀ i, i ≠ j → i ≠ k → σ i = τ i) ∧
      (σ j : ℕ) = τ j + 1 ∧ (σ k : ℕ) + 1 = τ k := by
  unfold hopUp
  constructor
  · intro h
    split_ifs at h with hc
    cases h
    refine ⟨fun i hj hk ↦ by simp [Function.update_of_ne hk, Function.update_of_ne hj], ?_, ?_⟩
    · rw [Function.update_of_ne hjk, Function.update_self]; dsimp only; omega
    · simp
  · rintro ⟨hag, h1, h2⟩
    have hτk := (τ k).isLt
    rw [dite_eq_left ⟨by omega, by omega⟩]
    congr 1
    funext i
    by_cases hk : i = k
    · subst hk; rw [Function.update_self]; exact Fin.ext (by dsimp only; omega)
    · by_cases hj : i = j
      · subst hj; rw [Function.update_of_ne hk, Function.update_self]
        exact Fin.ext (by dsimp only; omega)
      · simp [Function.update_of_ne hk, Function.update_of_ne hj, hag i hj hk]

theorem hopUp_ne_self {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (σ : Fin n → Fin 3) :
    hopUp j k σ ≠ some σ := fun h ↦ by have := ((hopUp_eq_some hjk σ σ).1 h).2.1; omega

/-- **A bond of spin-flip type in the `S^z` basis, entry by entry**: if the two-site matrix `B`
is `e mz mz` on the diagonal, `c₁` on the hop that lowers the first site and `c₂` on the hop that
raises it, the bond embedded on sites `j, k` has the diagonal `e mz mz`, the entry `c₁` for the
hop `j → k` and `c₂` for the hop `k → j`. (The `C`-twisted seam has `c₁ ≠ c₂`.) -/
theorem bondOp_apply_of_flip₂ {n : ℕ} {j k : Fin n} (hjk : j ≠ k)
    {B : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℤ} {e c₁ c₂ : ℤ}
    (hB : ∀ p q, B p q = (if p = q then e * (mz p.1 * mz p.2) else 0) +
      c₁ * (if p.1.val = q.1.val + 1 ∧ p.2.val + 1 = q.2.val then 1 else 0) +
      c₂ * (if p.1.val + 1 = q.1.val ∧ p.2.val = q.2.val + 1 then 1 else 0))
    (σ τ : Fin n → Fin 3) :
    bondOp j k B σ τ = (if σ = τ then e * (mz (σ j) * mz (σ k)) else 0) +
      c₁ * (if hopUp j k σ = some τ then 1 else 0) +
        c₂ * (if hopUp k j σ = some τ then 1 else 0) := by
  simp only [bondOp, of_apply, hB]
  by_cases hag : ∀ i, i ≠ j → i ≠ k → σ i = τ i
  · rw [ite_eq_left hag]
    have e1 : (σ j, σ k) = (τ j, τ k) ↔ σ = τ := by
      constructor
      · intro h
        simp only [Prod.mk.injEq] at h
        funext i
        by_cases hj : i = j
        · subst hj; exact h.1
        · by_cases hk : i = k
          · subst hk; exact h.2
          · exact hag i hj hk
      · rintro rfl; rfl
    have e2 : ((σ j : ℕ) = τ j + 1 ∧ (σ k : ℕ) + 1 = τ k) ↔ hopUp j k σ = some τ := by
      rw [hopUp_eq_some hjk]; exact ⟨fun h ↦ ⟨hag, h⟩, fun h ↦ h.2⟩
    have e3 : ((σ j : ℕ) + 1 = τ j ∧ (σ k : ℕ) = τ k + 1) ↔ hopUp k j σ = some τ := by
      rw [hopUp_eq_some (Ne.symm hjk)]
      exact ⟨fun h ↦ ⟨fun i hk hj ↦ hag i hj hk, h.2, h.1⟩, fun h ↦ ⟨h.2.2, h.2.1⟩⟩
    simp only [e1, e2, e3]
  · rw [ite_eq_right hag]
    have h1 : σ ≠ τ := fun h ↦ hag fun i _ _ ↦ by rw [h]
    have h2 : hopUp j k σ ≠ some τ := fun h ↦ hag ((hopUp_eq_some hjk σ τ).1 h).1
    have h3 : hopUp k j σ ≠ some τ := fun h ↦
      hag fun i hj hk ↦ ((hopUp_eq_some (Ne.symm hjk) σ τ).1 h).1 i hk hj
    simp [h1, h2, h3]

/-- **An exchange-type bond in the `S^z` basis, entry by entry**: if the two-site matrix `B` is
`mz mz` on the diagonal and `c` on the two spin-flip hops, the bond embedded on sites `j, k` has
the diagonal `mz mz` and the entry `c` for each of the two hops. -/
theorem bondOp_apply_of_flip {n : ℕ} {j k : Fin n} (hjk : j ≠ k)
    {B : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℤ} {c : ℤ}
    (hB : ∀ p q, B p q = (if p = q then mz p.1 * mz p.2 else 0) +
      c * (if p.1.val = q.1.val + 1 ∧ p.2.val + 1 = q.2.val then 1 else 0) +
      c * (if p.1.val + 1 = q.1.val ∧ p.2.val = q.2.val + 1 then 1 else 0))
    (σ τ : Fin n → Fin 3) :
    bondOp j k B σ τ = (if σ = τ then mz (σ j) * mz (σ k) else 0) +
      c * (if hopUp j k σ = some τ then 1 else 0) +
        c * (if hopUp k j σ = some τ then 1 else 0) := by
  simpa only [one_mul] using bondOp_apply_of_flip₂ hjk (e := 1) (c₁ := c) (c₂ := c)
    (fun p q ↦ by rw [hB p q, one_mul]) σ τ

/-- **The integer bond in the `S^z` basis, entry by entry**: the diagonal `mz mz` and one unit
entry for each of the two bond hops. -/
theorem bondOp_exchangeZ_apply {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (σ τ : Fin n → Fin 3) :
    bondOp j k exchangeZ σ τ = (if σ = τ then mz (σ j) * mz (σ k) else 0) +
      (if hopUp j k σ = some τ then 1 else 0) + (if hopUp k j σ = some τ then 1 else 0) := by
  rw [bondOp_apply_of_flip hjk (c := 1) (fun p q ↦ by simp only [exchangeZ, one_mul]) σ τ,
    one_mul, one_mul]

/-- The twisted bond in the `S^z` basis: the diagonal `mz mz` and the entry `-1` for each of the
two bond hops. -/
theorem bondOp_twistedExchangeZ_apply {n : ℕ} {j k : Fin n} (hjk : j ≠ k)
    (σ τ : Fin n → Fin 3) :
    bondOp j k twistedExchangeZ σ τ = (if σ = τ then mz (σ j) * mz (σ k) else 0) -
      (if hopUp j k σ = some τ then 1 else 0) - (if hopUp k j σ = some τ then 1 else 0) := by
  rw [bondOp_apply_of_flip hjk (c := -1) (fun p q ↦ by simp only [twistedExchangeZ]; ring) σ τ]
  ring

/-! ### Hops as bands of the word numbers -/

/-- The band of the hop `σⱼ ↦ σⱼ - 1`, `σₖ ↦ σₖ + 1` (`j ≠ k`) on the numbers of the words of
length `n`, with entry magnitude `c`: the hop changes the word number by `3ᵏ - 3ʲ`. -/
def hopBand (n j k c : ℕ) : Hop where
  up := decide (j < k)
  off := if j < k then 3 ^ k - 3 ^ j else 3 ^ j - 3 ^ k
  mask := (List.range (3 ^ n)).map fun x ↦ if 1 ≤ digit x j ∧ digit x k ≤ 1 then c else 0

/-- **The band of a hop, entry by entry**: `hopBand n j k c` has the entry `c` from the number of
`σ` to the number of `τ` exactly when the hop `j → k` takes `σ` to `τ`. -/
theorem hopBand_term {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (c : ℕ) (σ τ : Fin n → Fin 3) :
    (if (wordNum n τ : ℕ) = (hopBand n j k c).nb (wordNum n σ) then
      (hopBand n j k c).mask.getD (wordNum n σ) 0 else 0) =
      if hopUp j k σ = some τ then c else 0 := by
  have hx := (wordNum n σ).isLt
  have hm : (hopBand n j k c).mask.getD (wordNum n σ) 0 =
      if 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1 then c else 0 := by
    rw [hopBand, getD_map_range _ _ _ _ hx, digit_wordNum, digit_wordNum]
  rw [hm]
  by_cases hc : 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1
  · have hh : hopUp j k σ = some _ := dite_eq_left hc
    rw [hh, ite_eq_left hc]
    set τ' := Function.update (Function.update σ j ⟨σ j - 1, by have := (σ j).isLt; omega⟩) k
      ⟨σ k + 1, by omega⟩
    have he : ((wordNum n τ' : ℕ) : ℤ) = wordNum n σ + 3 ^ (k : ℕ) - 3 ^ (j : ℕ) := by
      rw [wordNum_update₂ σ hjk]
      push_cast [Nat.cast_sub hc.1]
      ring
    have hnb : (((hopBand n j k c).nb (wordNum n σ) : ℕ) : ℤ) =
        wordNum n σ + 3 ^ (k : ℕ) - 3 ^ (j : ℕ) := by
      simp only [Hop.nb, hopBand]
      by_cases hlt : (j : ℕ) < k
      · have hp : 3 ^ (j : ℕ) ≤ 3 ^ (k : ℕ) := Nat.pow_le_pow_right (by norm_num) hlt.le
        simp only [hlt, decide_true, ite_true]
        push_cast [hp]
        ring
      · have hp : 3 ^ (k : ℕ) ≤ 3 ^ (j : ℕ) :=
          Nat.pow_le_pow_right (by norm_num) (Nat.le_of_not_lt hlt)
        have hle : 3 ^ (j : ℕ) - 3 ^ (k : ℕ) ≤ (wordNum n σ : ℕ) :=
          le_trans (Nat.sub_le _ _) (le_wordNum σ j hc.1)
        simp only [hlt, decide_false, ite_false, Bool.false_eq_true]
        push_cast [hle, hp]
        ring
    have : (wordNum n τ : ℕ) = (hopBand n j k c).nb (wordNum n σ) ↔ some τ' = some τ := by
      rw [Option.some_inj, ← wordNum_eq_iff (σ := τ') (τ := τ), he, ← Nat.cast_inj (R := ℤ),
        hnb]
      constructor <;> intro h <;> linarith
    simp only [this]
  · have hh : hopUp j k σ = none := dite_eq_right hc
    rw [hh, ite_eq_right hc]
    simp

/-- **A hop band stays in range**: a nonzero entry in row `x < 3ⁿ` couples to a word number
`< 3ⁿ`, and a band below the diagonal does not reach below `0` (`Band.InRange`). -/
theorem hopBand_inRange {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (c x : ℕ) (hx : x < 3 ^ n)
    (hm : (hopBand n j k c).mask.getD x 0 ≠ 0) :
    ((hopBand n j k c).up = true → x + (hopBand n j k c).off < 3 ^ n) ∧
      ((hopBand n j k c).up = false → (hopBand n j k c).off ≤ x) := by
  set σ := (wordNum n).symm ⟨x, hx⟩
  have hσ : (wordNum n σ : ℕ) = x := by simp [σ]
  have hmask : (hopBand n j k c).mask.getD x 0 =
      if 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1 then c else 0 := by
    rw [hopBand, getD_map_range _ _ _ _ hx, ← hσ, digit_wordNum, digit_wordNum]
  rw [hmask] at hm
  have hc : 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1 := by by_contra h; exact hm (ite_eq_right h)
  obtain ⟨τ, hτ⟩ : ∃ τ, hopUp j k σ = some τ := ⟨_, dite_eq_left hc⟩
  have ht := hopBand_term hjk c σ τ
  rw [hτ, ite_eq_left rfl, hσ] at ht
  have hnb : (wordNum n τ : ℕ) = (hopBand n j k c).nb x := by
    by_contra h; rw [ite_eq_right h] at ht; exact hm (by rw [ite_eq_left hc]; exact ht.symm)
  have hlt := (wordNum n τ).isLt
  refine ⟨fun hu ↦ ?_, fun hu ↦ ?_⟩
  · simp only [Hop.nb, hu, ite_true] at hnb; omega
  · have h3 := le_wordNum σ j hc.1
    rw [hσ] at h3
    simp only [hopBand, decide_eq_false_iff_not] at hu ⊢
    simp only [hu, ite_false]
    exact le_trans (Nat.sub_le _ _) h3

end CertifiedGaps.SpinOneHeisenberg
