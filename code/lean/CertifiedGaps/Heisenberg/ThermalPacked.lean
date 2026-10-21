/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalTable
import CertifiedGaps.Heisenberg.WordBand

/-!
# A fast evaluator of the periodic integer totals

`HornerTotals` (`Heisenberg/ThermalTable.lean`) says that `hornerTotal n g`, the integer total
`W_{n,g}` of `thermal72.tex`, equals the printed `paperTotal n g`. `hornerTotal` is written with
the generic `Trace.hornerFloor`, whose vectors are functions, so it cannot be evaluated in
practice. This file defines `packedHornerTotal n g`, the same total computed by the packed
big-integer recursion `Trace.Packed.runCols`, and proves that the two agree
(`hornerTotal_eq_of_packed`).

## The evaluator

* **The matrix.** Words are numbered `∑ⱼ σⱼ 3ʲ` (`wordNum`); `D = (32 - 3n) I - 2 H_n(g)` is the
  banded matrix `ringBand n g`: the diagonal `ringDiag`, the entries `-2` of the hops of the
  bonds `j, j + 1` (`hopBand`), and the entries `∓2` of the seam hops (`seamHops`; `+2` for the
  twisted ring `g = 1`, whose spin-flip entries on the seam bond are negated). `ringBand_entry`
  and `thermalMatrix_apply` compute both sides entry by entry.
* **The coefficients.** `c_j = ⌊2⁵⁵ p_j⌋ = ⌊2⁵⁵ v_j / ∑ v⌋` with `v_j = 98ʲ 280!/j!`
  (`thermalCoeff_eq`); they vanish for `j ≥ 192`, so the recursion starts at `j = 191`
  (`specCol_append_zeros`), as the paper notes.
* **The columns.** `repNum n x` is the number of the representative `thermalRep` of the word
  numbered `x`, computed on numbers (the shift is `x ↦ ⌊x/3⌋ + (x mod 3) 3ⁿ⁻¹`, the reversal is
  the code `wordCode`, the negation is `3ⁿ - 1 - x`; `repNum_wordNum`). `countArr` counts the
  words of every representative, and the representatives are the numbers with a nonzero count
  (`repCols`). One packed run computes the representative columns only (lane `r` is column
  `repCols[r]`), and `numerator` weighs each by its count: `∑_ρ ν_ρ ∑_w X_ρ(w)²`, the paper's
  eq. (thermal72-total). This step, `packedTotal`, takes the band, the representative map and
  the coefficients as arguments, so that `Heisenberg/ThermalPacked60.lean` reuses it for the
  `b = 21/2` table.
* **The parameters** (`thermalParams`): lanes of `64` bits, bias `2⁵⁷`, scale `32`, and the
  exact division `⌊u / 32⌋ = u >>> 5` (magic multiplier `1`). The range tests of the packed run
  check `|X| < 2⁵⁶`; the paper's bound is `|X| < 2Q = 2⁵⁶`.

## Main results

* `ringBand_entry`, `thermalMatrix_apply`, `thermalMatrix_eq_submatrix`: `D` is the banded
  matrix relabelled by the word numbers.
* `thermalColumn_eq`: the columns `X_ρ` are the column recursions of the banded matrix.
* `repNum_wordNum`: the search for the representative on numbers.
* `hornerTotal_eq_sum_rep`: `W = ∑_x ∑_w X_{rep x}(w)²`, the sum over representatives with
  multiplicities written as a sum over all words.
* `packedTotal_spec`: a packed total is the sum over all rows of the squared norms of their
  representative columns.
* `hornerTotal_eq_of_packed : packedHornerTotal n g = some W → hornerTotal n g = W`.

## Comparison with the paper

* The paper runs the recursion in signed 64-bit integers; here the lanes are 64 bits wide with
  the bias `2⁵⁷`, and every step checks the range instead of relying on the paper's a priori
  bound, so a run either returns the exact integers or `none`.
* The paper omits the steps `j > 191` because their coefficients vanish; so does
  `packedHornerTotal`, and `coeffNat_tail` checks the vanishing.
-/

open Finset Matrix
open scoped Nat

namespace CertifiedGaps.SpinOneHeisenberg

open Trace Trace.Packed

namespace ThermalPacked

/-! ### The banded matrix -/

/-- The diagonal `32 - 3n - 2 ∑ⱼ wⱼ wⱼ₊₁` of `D` at the word numbered `x` (the sum includes the
seam bond `n - 1, 0`). -/
def ringDiag (n x : ℕ) : ℤ :=
  32 - 3 * n - 2 * (((List.range (n - 1)).map fun j ↦ letter x j * letter x (j + 1)).sum +
    letter x (n - 1) * letter x 0)

/-- The two hops of the seam bond `n - 1, 0`, with entry magnitude `2`. -/
def seamHops (n : ℕ) : List Hop := [hopBand n (n - 1) 0 2, hopBand n 0 (n - 1) 2]

/-- **The integer matrix `D = (32 - 3n) I - 2 H_n(g)` as a banded matrix** on the word numbers:
the seam hops have the entry `-2` on the untwisted ring and `+2` on the twisted ring. -/
def ringBand (n : ℕ) (g : Fin 2) : Band where
  pos := (List.range (3 ^ n)).map fun x ↦ (ringDiag n x).toNat
  neg := (List.range (3 ^ n)).map fun x ↦ (-ringDiag n x).toNat
  hops := (List.range (n - 1)).flatMap (fun j ↦ [hopBand n j (j + 1) 2, hopBand n (j + 1) j 2]) ++
    (if g = 0 then seamHops n else [])
  phops := if g = 0 then [] else seamHops n

/-- The sign of the spin-flip entries of the seam bond: `+1` untwisted, `-1` twisted. -/
def seamSign : Fin 2 → ℤ := ![1, -1]

/-- The diagonal of `D` in terms of the letters of the word. -/
theorem ringDiag_wordNum (m : ℕ) (σ : Fin (m + 2) → Fin 3) :
    ringDiag (m + 2) (wordNum (m + 2) σ) = 32 - 3 * ((m + 2 : ℕ) : ℤ) -
      2 * (∑ i : Fin (m + 1), mz (σ i.castSucc) * mz (σ i.succ) +
        mz (σ (Fin.last (m + 1))) * mz (σ 0)) := by
  have hm : m + 2 - 1 = m + 1 := by omega
  rw [ringDiag, hm, list_sum_range, ← Fin.sum_univ_eq_sum_range]
  have h0 : letter (wordNum (m + 2) σ) 0 = mz (σ 0) := letter_wordNum σ 0
  have hl : letter (wordNum (m + 2) σ) (m + 1) = mz (σ (Fin.last (m + 1))) :=
    letter_wordNum σ (Fin.last (m + 1))
  have hj : ∀ j : Fin (m + 1), letter (wordNum (m + 2) σ) j *
      letter (wordNum (m + 2) σ) (j + 1) = mz (σ j.castSucc) * mz (σ j.succ) := fun j ↦ by
    rw [show (j : ℕ) = (j.castSucc : Fin (m + 2)) from rfl, letter_wordNum,
      show (j.castSucc : ℕ) + 1 = (j.succ : Fin (m + 2)) by simp, letter_wordNum]
  simp only [h0, hl, hj]

/-- **The entries of the banded matrix in terms of hops.** -/
theorem ringBand_entry (g : Fin 2) (m : ℕ) (σ τ : Fin (m + 2) → Fin 3) :
    (ringBand (m + 2) g).matrix (3 ^ (m + 2)) (wordNum (m + 2) σ) (wordNum (m + 2) τ) =
      (if σ = τ then ringDiag (m + 2) (wordNum (m + 2) σ) else 0) -
      ∑ i : Fin (m + 1), ((if hopUp i.castSucc i.succ σ = some τ then (2 : ℤ) else 0) +
        (if hopUp i.succ i.castSucc σ = some τ then (2 : ℤ) else 0)) -
      seamSign g * ((if hopUp (Fin.last (m + 1)) 0 σ = some τ then (2 : ℤ) else 0) +
        (if hopUp 0 (Fin.last (m + 1)) σ = some τ then (2 : ℤ) else 0)) := by
  have hk := (wordNum (m + 2) σ).isLt
  have hm : m + 2 - 1 = m + 1 := by omega
  have hbond : ((List.range (m + 1)).map fun a ↦ (([hopBand (m + 2) a (a + 1) 2,
      hopBand (m + 2) (a + 1) a 2]).map fun h ↦
        if ((wordNum (m + 2) τ : Fin _) : ℕ) = h.nb (wordNum (m + 2) σ) then
          (h.mask.getD (wordNum (m + 2) σ) 0 : ℤ) else 0).sum).sum =
      ∑ i : Fin (m + 1), ((if hopUp i.castSucc i.succ σ = some τ then (2 : ℤ) else 0) +
        (if hopUp i.succ i.castSucc σ = some τ then (2 : ℤ) else 0)) := by
    rw [list_sum_range, ← Fin.sum_univ_eq_sum_range]
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
    have h1 := congrArg (fun x : ℕ ↦ (x : ℤ))
      (hopBand_term (j := i.castSucc) (k := i.succ) (by simp [Fin.ext_iff]) 2 σ τ)
    have h2 := congrArg (fun x : ℕ ↦ (x : ℤ))
      (hopBand_term (j := i.succ) (k := i.castSucc) (by simp [Fin.ext_iff]) 2 σ τ)
    simp only [Fin.val_castSucc, Fin.val_succ, Nat.cast_ite, Nat.cast_zero,
      Nat.cast_ofNat] at h1 h2 ⊢
    rw [h1, h2]
  have hseam : ((seamHops (m + 2)).map fun h ↦
        if ((wordNum (m + 2) τ : Fin _) : ℕ) = h.nb (wordNum (m + 2) σ) then
          (h.mask.getD (wordNum (m + 2) σ) 0 : ℤ) else 0).sum =
      (if hopUp (Fin.last (m + 1)) 0 σ = some τ then (2 : ℤ) else 0) +
        (if hopUp 0 (Fin.last (m + 1)) σ = some τ then (2 : ℤ) else 0) := by
    simp only [seamHops, hm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil,
      add_zero]
    have h1 := congrArg (fun x : ℕ ↦ (x : ℤ))
      (hopBand_term (n := m + 2) (j := Fin.last (m + 1)) (k := 0) (by simp [Fin.ext_iff]) 2 σ τ)
    have h2 := congrArg (fun x : ℕ ↦ (x : ℤ))
      (hopBand_term (n := m + 2) (j := 0) (k := Fin.last (m + 1)) (by simp [Fin.ext_iff]) 2 σ τ)
    push_cast at h1 h2
    exact congrArg₂ (· + ·) h1 h2
  have hdiag : (if wordNum (m + 2) σ = wordNum (m + 2) τ then
      (((List.range (3 ^ (m + 2))).map fun x ↦ (ringDiag (m + 2) x).toNat).getD
        (wordNum (m + 2) σ) 0 : ℤ) -
        (((List.range (3 ^ (m + 2))).map fun x ↦ (-ringDiag (m + 2) x).toNat).getD
          (wordNum (m + 2) σ) 0 : ℤ) else 0) =
      if σ = τ then ringDiag (m + 2) (wordNum (m + 2) σ) else 0 := by
    rw [getD_map_range _ _ _ _ hk, getD_map_range _ _ _ _ hk, Int.toNat_sub_toNat_neg]
    simp only [(wordNum (m + 2)).injective.eq_iff]
  simp only [Band.matrix, ringBand, List.map_append, List.sum_append, sum_map_flatMap, hm]
  rw [hdiag, hbond]
  fin_cases g
  · simp only [Fin.zero_eta, ↓reduceIte, List.map_nil, List.sum_nil, add_zero, hseam,
      seamSign, Matrix.cons_val_zero, one_mul]
    ring
  · simp only [Fin.mk_one, Fin.isValue, one_ne_zero, ↓reduceIte, List.map_nil, List.sum_nil,
      add_zero, hseam, seamSign, Matrix.cons_val_one, Matrix.cons_val_fin_one]
    ring

/-- The cyclic successor of `i < m + 1` in a ring of `m + 2` sites. -/
theorem finRotate_castSucc (m : ℕ) (i : Fin (m + 1)) :
    finRotate (m + 2) i.castSucc = i.succ :=
  Fin.ext (by rw [coe_finRotate_of_ne_last (Fin.castSucc_lt_last i).ne]; simp)

/-- **The entries of `D = (32 - 3n) I - 2 H_n(g)` in terms of hops.** -/
theorem thermalMatrix_apply (g : Fin 2) (m : ℕ) (σ τ : Fin (m + 2) → Fin 3) :
    thermalMatrix (m + 2) g σ τ =
      (if σ = τ then ringDiag (m + 2) (wordNum (m + 2) σ) else 0) -
      ∑ i : Fin (m + 1), ((if hopUp i.castSucc i.succ σ = some τ then (2 : ℤ) else 0) +
        (if hopUp i.succ i.castSucc σ = some τ then (2 : ℤ) else 0)) -
      seamSign g * ((if hopUp (Fin.last (m + 1)) 0 σ = some τ then (2 : ℤ) else 0) +
        (if hopUp 0 (Fin.last (m + 1)) σ = some τ then (2 : ℤ) else 0)) := by
  have hne : ∀ i : Fin (m + 1), (i.castSucc : Fin (m + 2)) ≠ i.succ := fun i ↦ by
    simp [Fin.ext_iff]
  have hl0 : (Fin.last (m + 1) : Fin (m + 2)) ≠ 0 := by simp [Fin.ext_iff]
  have hseam : bondOp (Fin.last (m + 1)) 0 (seamBond g) σ τ =
      (if σ = τ then mz (σ (Fin.last (m + 1))) * mz (σ 0) else 0) +
        seamSign g * ((if hopUp (Fin.last (m + 1)) 0 σ = some τ then 1 else 0) +
          (if hopUp 0 (Fin.last (m + 1)) σ = some τ then 1 else 0)) := by
    fin_cases g
    · simp only [Fin.zero_eta, seamBond, Matrix.cons_val_zero, seamSign, one_mul]
      rw [bondOp_exchangeZ_apply hl0]
      ring
    · simp only [Fin.mk_one, seamBond, Matrix.cons_val_one, Matrix.cons_val_fin_one, seamSign]
      rw [bondOp_twistedExchangeZ_apply hl0]
      ring
  have hring : ringZ (m + 2) g σ τ =
      ∑ i : Fin (m + 1), bondOp i.castSucc i.succ exchangeZ σ τ +
        bondOp (Fin.last (m + 1)) 0 (seamBond g) σ τ := by
    simp only [ringZ, twistedRingHamiltonian, Matrix.sum_apply]
    rw [Fin.sum_univ_castSucc]
    congr 1
    · refine Finset.sum_congr rfl fun i _ ↦ ?_
      have hi : ¬((i.castSucc : Fin (m + 2)) : ℕ) + 1 = m + 2 := by simp; omega
      rw [finRotate_castSucc]
      simp only [hi, ite_false]
    · rw [finRotate_last]
      simp only [Fin.val_last, ite_true]
  have ite2 : ∀ (c : Prop) [Decidable c], (if c then (2 : ℤ) else 0) = 2 * if c then 1 else 0 :=
    fun c _ ↦ by split_ifs <;> simp
  have hb : ∀ i : Fin (m + 1), bondOp i.castSucc i.succ exchangeZ σ τ =
      (if σ = τ then mz (σ i.castSucc) * mz (σ i.succ) else 0) +
        (if hopUp i.castSucc i.succ σ = some τ then 1 else 0) +
          (if hopUp i.succ i.castSucc σ = some τ then 1 else 0) := fun i ↦
    bondOp_exchangeZ_apply (hne i) σ τ
  rw [thermalMatrix, Matrix.sub_apply, Matrix.smul_apply, Matrix.smul_apply, hring, hseam]
  simp only [hb, Matrix.one_apply, smul_eq_mul, nsmul_eq_mul, ite2, Finset.sum_add_distrib,
    ← Finset.mul_sum]
  by_cases hστ : σ = τ
  · subst hστ
    have hu : ∀ i : Fin (m + 1), hopUp i.castSucc i.succ σ ≠ some σ := fun i ↦
      hopUp_ne_self (hne i) σ
    have hd : ∀ i : Fin (m + 1), hopUp i.succ i.castSucc σ ≠ some σ := fun i ↦
      hopUp_ne_self (hne i).symm σ
    simp only [ite_true, hu, hd, hopUp_ne_self hl0, hopUp_ne_self hl0.symm, ite_false,
      Finset.sum_const_zero, mul_zero, add_zero, ringDiag_wordNum]
    push_cast
    ring
  · simp only [hστ, ite_false, Finset.sum_const_zero]
    push_cast
    ring

/-- **`D` is the banded matrix relabelled by the word numbers.** -/
theorem thermalMatrix_eq_submatrix (g : Fin 2) (m : ℕ) :
    thermalMatrix (m + 2) g =
      ((ringBand (m + 2) g).matrix (3 ^ (m + 2))).submatrix (wordNum (m + 2)) (wordNum (m + 2)) :=
  by
  ext σ τ
  rw [submatrix_apply, ringBand_entry, thermalMatrix_apply]

/-! ### The coefficients -/

/-- The numerators `v_j = 98ʲ 280!/j!` of the conditioned Poisson weights. -/
def poissonNum (j : ℕ) : ℕ := 98 ^ j * (280 ! / j !)

/-- `∑_{j ≤ 280} v_j`. -/
def poissonDen : ℕ := ((List.range 281).map poissonNum).sum

/-- The coefficient `c_j = ⌊2⁵⁵ v_j / ∑ v⌋` in natural-number arithmetic. -/
def coeffNat (j : ℕ) : ℕ := 2 ^ 55 * poissonNum j / poissonDen

/-- The coefficients `c_0, …, c_191`; the later ones vanish. -/
def thermalCoeffs : List ℕ := (List.range 192).map coeffNat

/-- `p_j = v_j / ∑ v`: the factor `280!` cancels. -/
theorem poissonWeight_eq {j : ℕ} (hj : j ≤ 280) :
    poissonWeight 98 280 j = (poissonNum j : ℚ) / poissonDen := by
  have hcast : ∀ l ≤ 280, (poissonNum l : ℚ) = 98 ^ l / l ! * 280 ! := fun l hl ↦ by
    rw [poissonNum, Nat.cast_mul, Nat.cast_div (Nat.factorial_dvd_factorial hl)
      (by positivity)]
    push_cast
    ring
  have hden : (poissonDen : ℚ) = poissonMass 98 280 * 280 ! := by
    rw [poissonDen, list_sum_range, Nat.cast_sum, poissonMass, Finset.sum_mul]
    exact Finset.sum_congr rfl fun l hl ↦ hcast l (by simp at hl; omega)
  have h280 : (280 ! : ℚ) ≠ 0 := by positivity
  rw [poissonWeight, hcast j hj, hden, mul_div_mul_right _ _ h280]

/-- `c_j = ⌊2⁵⁵ p_j⌋` is the natural-number quotient. -/
theorem thermalCoeff_eq {j : ℕ} (hj : j ≤ 280) : thermalCoeff j = coeffNat j := by
  rw [thermalCoeff, poissonWeight_eq hj, thermalScale, ← mul_div_assoc,
    show (2 : ℚ) ^ 55 * poissonNum j = ((2 ^ 55 * poissonNum j : ℕ) : ℚ) by push_cast; rfl,
    Rat.floor_natCast_div_natCast, coeffNat]
  push_cast
  rfl

set_option maxRecDepth 10000 in
/-- The coefficients `c_j` with `192 ≤ j ≤ 280` vanish. -/
theorem coeffNat_tail : ∀ j ∈ List.range' 192 89, coeffNat j = 0 := by decide +kernel

/-- The coefficient list of `thermalColumn` is `thermalCoeffs` followed by `89` zeros. -/
theorem coeff_list :
    (List.range (280 + 1)).map thermalCoeff =
      thermalCoeffs.map (↑) ++ List.replicate 89 (0 : ℤ) := by
  refine List.ext_getElem (by simp [thermalCoeffs]) fun j h1 h2 ↦ ?_
  have hj : j < 281 := by simpa using h1
  rw [List.getElem_map, List.getElem_range, thermalCoeff_eq (by omega), List.getElem_append]
  split_ifs with h192
  · simp [thermalCoeffs]
  · have h192' : 192 ≤ j := by simpa [thermalCoeffs] using h192
    rw [List.getElem_replicate, coeffNat_tail j (List.mem_range'_1.2 ⟨h192', by omega⟩)]
    rfl

/-- Zero coefficients at the end of the list (the first steps of the recursion) leave the
column at zero. -/
theorem specCol_append_zeros (D : Band) (s r : ℕ) (l : List ℤ) (z : ℕ) :
    specCol D s r (l ++ List.replicate z 0) = specCol D s r l := by
  have h0 : (List.replicate z (0 : ℤ)).foldr (specStep D s r) 0 = 0 := by
    induction z with
    | zero => rfl
    | succ z ih =>
      rw [List.replicate_succ, List.foldr_cons, ih]
      funext k
      simp [specStep, Band.apply]
  rw [specCol, List.foldr_append, h0]
  rfl

/-- The Horner columns of `D`, for any coefficients, are the column recursions of the banded
matrix. -/
theorem hornerFloor_thermalMatrix {g : Fin 2} {m : ℕ}
    (hD : (ringBand (m + 2) g).InRange (3 ^ (m + 2))) (c : ℕ → ℤ) (J : ℕ)
    (ρ w : Fin (m + 2) → Fin 3) :
    hornerFloor (thermalMatrix (m + 2) g) 32 c ρ J w = specCol (ringBand (m + 2) g) 32
      (wordNum (m + 2) ρ) ((List.range (J + 1)).map c) (wordNum (m + 2) w) := by
  rw [thermalMatrix_eq_submatrix, hornerFloor_submatrix, hornerFloor_eq_specCol hD]

/-- **The columns of `D` are the column recursions of the banded matrix.** -/
theorem thermalColumn_eq {g : Fin 2} {m : ℕ} (hD : (ringBand (m + 2) g).InRange (3 ^ (m + 2)))
    (ρ w : Fin (m + 2) → Fin 3) :
    thermalColumn (m + 2) g ρ w = specCol (ringBand (m + 2) g) 32 (wordNum (m + 2) ρ)
      (thermalCoeffs.map (↑)) (wordNum (m + 2) w) := by
  rw [thermalColumn, hornerFloor_thermalMatrix hD, coeff_list, specCol_append_zeros]

/-! ### Representatives on word numbers

`thermalRep` composes permutations of `Fin n → Fin 3`, which is slow to evaluate. The same
search runs here on word numbers, with the shift, the reversal and the negation as arithmetic
(`repNum_wordNum`). -/

/-- The cyclic shift of the word numbered `x < 3ⁿ`: digit `0` moves to site `n - 1`. -/
def rotNum (n x : ℕ) : ℕ := x / 3 + x % 3 * 3 ^ (n - 1)

/-- The code `wordCode` of the word numbered `x`, its digits with site `0` most significant;
also the number of the reversed word. -/
def codeNum (n x : ℕ) : ℕ := ((List.range n).map fun i ↦ digit x i * 3 ^ (n - 1 - i)).sum

/-- The number of the negated word. -/
def negNum (n x : ℕ) : ℕ := 3 ^ n - 1 - x

/-- The weight `∑ⱼ wⱼ` of the word numbered `x`. -/
def weightNum (n x : ℕ) : ℤ := ((List.range n).map fun i ↦ letter x i).sum

/-- The numbers of the cyclic shifts of the word numbered `y` and of their reversals. -/
def dihedralNums (n y : ℕ) : List ℕ :=
  (List.range n).flatMap fun k ↦ [(rotNum n)^[k] y, codeNum n ((rotNum n)^[k] y)]

/-- The numbers of the words of `orbitWords`. -/
def orbitNums (n x : ℕ) : List ℕ :=
  (if 0 ≤ weightNum n x then dihedralNums n x else []) ++
    (if weightNum n x ≤ 0 then dihedralNums n (negNum n x) else [])

/-- `lexMinStep` on word numbers. -/
def lexMinNum (n : ℕ) (a : ℕ × ℕ) (y : ℕ) : ℕ × ℕ :=
  if codeNum n y < a.1 then (codeNum n y, y) else a

/-- **The representative on word numbers**: the number of `thermalRep` of the word numbered `x`
(`repNum_wordNum`). -/
def repNum (n x : ℕ) : ℕ :=
  let x₀ := if 0 ≤ weightNum n x then x else negNum n x
  ((orbitNums n x).foldl (lexMinNum n) (codeNum n x₀, x₀)).2

theorem wordNum_coe {n : ℕ} (σ : Fin n → Fin 3) :
    (wordNum n σ : ℕ) = ∑ i, (σ i : ℕ) * 3 ^ (i : ℕ) := by
  rw [wordNum, finFunctionFinEquiv_apply]

/-- The cyclic successor of `i < m` in a ring of `m + 1` sites. -/
theorem finRotate_castSucc' (m : ℕ) (i : Fin m) : finRotate (m + 1) i.castSucc = i.succ :=
  Fin.ext (by rw [coe_finRotate_of_ne_last (Fin.castSucc_lt_last i).ne]; simp)

/-- The shift is a division by `3` and a carry of digit `0` to the top. -/
theorem wordNum_rotWord {n : ℕ} (σ : Fin n → Fin 3) :
    (wordNum n (rotWord n σ) : ℕ) = rotNum n (wordNum n σ) := by
  cases n with
  | zero => simp [rotNum]
  | succ m =>
    have hT : ∑ i : Fin m, (σ i.succ : ℕ) * 3 ^ ((i.succ : Fin (m + 1)) : ℕ) =
        3 * ∑ i : Fin m, (σ i.succ : ℕ) * 3 ^ (i : ℕ) := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ ↦ by rw [Fin.val_succ, pow_succ]; ring
    have h0 := (σ 0).isLt
    rw [wordNum_coe, wordNum_coe, rotNum, Fin.sum_univ_castSucc, Fin.sum_univ_succ, hT]
    simp only [rotWord, relabel_apply, finRotate_castSucc', finRotate_last, Fin.val_castSucc,
      Fin.val_last, Fin.val_zero, pow_zero, mul_one, Nat.add_sub_cancel]
    rw [Nat.add_mul_div_left _ _ (by norm_num), Nat.div_eq_of_lt h0, zero_add,
      Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt h0]

/-- The code of a word is the code of its number. -/
theorem codeNum_wordNum {n : ℕ} (σ : Fin n → Fin 3) : codeNum n (wordNum n σ) = wordCode σ := by
  rw [codeNum, list_sum_range,
    ← Fin.sum_univ_eq_sum_range (fun i ↦ digit (wordNum n σ) i * 3 ^ (n - 1 - i)), wordCode]
  simp only [digit_wordNum]

/-- The number of the reversed word is the code of the word. -/
theorem wordNum_revWord {n : ℕ} (σ : Fin n → Fin 3) :
    (wordNum n (revWord n σ) : ℕ) = wordCode σ := by
  rw [wordNum_coe, wordCode]
  refine Fintype.sum_equiv Fin.revPerm _ _ fun i ↦ ?_
  simp only [revWord, relabel_apply, Fin.revPerm_apply, Fin.val_rev]
  congr 2
  omega

/-- The number of the negated word. -/
theorem wordNum_negWord {n : ℕ} (σ : Fin n → Fin 3) :
    (wordNum n (negWord n σ) : ℕ) = negNum n (wordNum n σ) := by
  have hg : ∀ m : ℕ, ∑ i ∈ range m, 2 * 3 ^ i + 1 = 3 ^ m := fun m ↦ by
    induction m with
    | zero => simp
    | succ m ih => rw [Finset.sum_range_succ, pow_succ]; omega
  have hs : (wordNum n (negWord n σ) : ℕ) + wordNum n σ + 1 = 3 ^ n := by
    rw [wordNum_coe, wordNum_coe, ← Finset.sum_add_distrib, ← hg n,
      ← Fin.sum_univ_eq_sum_range (fun i ↦ 2 * 3 ^ i)]
    congr 1
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    simp only [negWord, onsite_apply, Fin.revPerm_apply, Fin.val_rev]
    have := (σ i).isLt
    rw [← add_mul]
    congr 1
    omega
  rw [negNum]
  omega

/-- The weight of a word is the weight of its number. -/
theorem weightNum_wordNum {n : ℕ} (σ : Fin n → Fin 3) :
    weightNum n (wordNum n σ) = wordWeight σ := by
  rw [weightNum, list_sum_range,
    ← Fin.sum_univ_eq_sum_range (fun i ↦ letter (wordNum n σ) i), wordWeight]
  simp only [letter_wordNum]

/-- The shifts by `k` sites, on numbers. -/
theorem wordNum_rotWord_pow {n : ℕ} (σ : Fin n → Fin 3) (k : ℕ) :
    (wordNum n ((rotWord n ^ k) σ) : ℕ) = (rotNum n)^[k] (wordNum n σ) := by
  induction k with
  | zero => rfl
  | succ k ih =>
    rw [pow_succ', Equiv.Perm.mul_apply, wordNum_rotWord, ih, Function.iterate_succ_apply']

theorem map_dihedralPerms {n : ℕ} (σ : Fin n → Fin 3) :
    ((dihedralPerms n).map (· σ)).map (fun v ↦ (wordNum n v : ℕ)) =
      dihedralNums n (wordNum n σ) := by
  simp only [dihedralPerms, dihedralNums, List.map_flatMap, List.map_cons,
    List.map_nil, Equiv.Perm.mul_apply, wordNum_revWord, wordNum_rotWord_pow,
    ← codeNum_wordNum]

/-- The numbers of `orbitWords`. -/
theorem map_orbitWords {n : ℕ} (σ : Fin n → Fin 3) :
    (orbitWords σ).map (fun v ↦ (wordNum n v : ℕ)) = orbitNums n (wordNum n σ) := by
  have hneg : negNum n (wordNum n σ) = wordNum n (negWord n σ) := (wordNum_negWord σ).symm
  rw [orbitWords, orbitNums, weightNum_wordNum, List.map_append, hneg]
  split_ifs <;>
    simp only [List.map_nil, List.nil_append, List.append_nil, map_dihedralPerms]

/-- `lexMinStep` and `lexMinNum` agree through the numbering. -/
theorem foldl_lexMinNum {n : ℕ} (l : List (Fin n → Fin 3)) (c : ℕ) (w : Fin n → Fin 3) :
    ((l.map fun v ↦ (wordNum n v : ℕ)).foldl (lexMinNum n) (c, wordNum n w)).2 =
      wordNum n ((l.foldl lexMinStep (c, w)).2) := by
  induction l generalizing c w with
  | nil => rfl
  | cons b l ih =>
    simp only [List.map_cons, List.foldl_cons, lexMinNum, lexMinStep, codeNum_wordNum]
    split_ifs
    · exact ih _ _
    · exact ih _ _

/-- **The representative on numbers is the number of the representative.** -/
theorem repNum_wordNum {n : ℕ} (σ : Fin n → Fin 3) :
    repNum n (wordNum n σ) = wordNum n (thermalRep σ) := by
  have hx₀ : (if 0 ≤ weightNum n (wordNum n σ) then (wordNum n σ : ℕ)
      else negNum n (wordNum n σ)) =
      wordNum n (if 0 ≤ wordWeight σ then σ else negWord n σ) := by
    rw [weightNum_wordNum]
    split_ifs
    · rfl
    · rw [wordNum_negWord]
  simp only [repNum, thermalRep]
  rw [hx₀, ← map_orbitWords, codeNum_wordNum, foldl_lexMinNum]

/-- The representative of a word number is a word number. -/
theorem repNum_lt {n x : ℕ} (hx : x < 3 ^ n) : repNum n x < 3 ^ n := by
  have h := repNum_wordNum ((wordNum n).symm ⟨x, hx⟩)
  simp only [Equiv.apply_symm_apply] at h
  rw [h]
  exact Fin.isLt _

/-! ### Multiplicities -/

/-- `cnt[j]` = the number of `x < d` with `f x = j`, accumulated in one pass. -/
def countArr (d : ℕ) (f : ℕ → ℕ) : Array ℕ :=
  (List.range d).foldl (fun a x ↦ a.modify (f x) (· + 1)) (Array.replicate d 0)

/-- The numbers below `d` with a nonzero count: the representatives. -/
def repCols (cnt : Array ℕ) (d : ℕ) : List ℕ := (List.range d).filter fun j ↦ cnt.getD j 0 != 0

/-- The units of the columns `cols`: `2^(w r)` in row `cols[r]`. -/
def colUnits (w : ℕ) (cols : List ℕ) (d : ℕ) : List ℕ :=
  (List.range d).map fun k ↦ if k ∈ cols then 2 ^ (w * cols.idxOf k) else 0

theorem countArr_spec (d : ℕ) (f : ℕ → ℕ) : ∀ m : ℕ,
    ((List.range m).foldl (fun a x ↦ a.modify (f x) (· + 1)) (Array.replicate d 0)).size = d ∧
      ∀ j < d, ((List.range m).foldl (fun a x ↦ a.modify (f x) (· + 1))
        (Array.replicate d 0)).getD j 0 = #{x ∈ range m | f x = j}
  | 0 => ⟨by simp, fun j hj ↦ by simp [hj]⟩
  | m + 1 => by
    obtain ⟨hs, hc⟩ := countArr_spec d f m
    rw [List.range_succ, List.foldl_append, List.foldl_cons, List.foldl_nil]
    refine ⟨by rw [Array.size_modify, hs], fun j hj ↦ ?_⟩
    rw [Finset.range_add_one, Finset.filter_insert]
    have hj' : j < ((List.range m).foldl (fun a x ↦ a.modify (f x) (· + 1))
        (Array.replicate d 0)).size := by rw [hs]; exact hj
    rw [Array.getD_eq_getD_getElem?, Array.getElem?_modify]
    by_cases hfm : f m = j
    · simp only [hfm, ite_true]
      rw [Finset.card_insert_of_notMem (by simp), ← hc j hj, Array.getD_eq_getD_getElem?]
      simp [Array.getElem?_eq_getElem hj']
    · simp only [hfm, ite_false]
      rw [← hc j hj, Array.getD_eq_getD_getElem?]

/-- The count array counts. -/
theorem countArr_getD {d : ℕ} (f : ℕ → ℕ) {j : ℕ} (hj : j < d) :
    (countArr d f).getD j 0 = #{x ∈ range d | f x = j} :=
  (countArr_spec d f d).2 j hj

/-- A sum over `x < d` of a function of `f x < d`, grouped by the value of `f`. -/
theorem sum_comp_eq_sum_count {d : ℕ} {f : ℕ → ℕ} (hf : ∀ x < d, f x < d) (F : ℕ → ℤ) :
    ∑ x ∈ range d, F (f x) = ∑ j ∈ range d, ((countArr d f).getD j 0 : ℤ) * F j := by
  rw [← Finset.sum_fiberwise_of_maps_to (g := f) (t := range d)
    (fun x hx ↦ mem_range.2 (hf x (mem_range.1 hx)))]
  refine Finset.sum_congr rfl fun j hj ↦ ?_
  rw [countArr_getD f (mem_range.1 hj), Finset.sum_congr rfl
    (g := fun _ ↦ F j) (fun x hx ↦ by rw [(Finset.mem_filter.1 hx).2]),
    Finset.sum_const, nsmul_eq_mul]

/-- A list sum over the entries of a list, as a sum over positions. -/
theorem list_sum_map_eq_sum_getD (l : List ℕ) (F : ℕ → ℤ) :
    (l.map F).sum = ∑ r ∈ range l.length, F (l.getD r 0) := by
  induction l with
  | nil => simp
  | cons a l ih =>
    rw [List.map_cons, List.sum_cons, ih, List.length_cons, Finset.sum_range_succ']
    simp [add_comm]

/-- Dropping the zero terms of a sum over `j < d`: the sum over the representatives. -/
theorem sum_range_eq_sum_repCols (cnt : Array ℕ) (d : ℕ) (G : ℕ → ℤ) :
    ∑ j ∈ range d, (cnt.getD j 0 : ℤ) * G j =
      ∑ r ∈ range (repCols cnt d).length, (((repCols cnt d).map (cnt.getD · 0)).getD r 0 : ℤ) *
        G ((repCols cnt d).getD r 0) := by
  have hfilter : ∀ l : List ℕ, (((l.filter fun j ↦ cnt.getD j 0 != 0).map
      fun j ↦ (cnt.getD j 0 : ℤ) * G j).sum) = (l.map fun j ↦ (cnt.getD j 0 : ℤ) * G j).sum := by
    intro l
    induction l with
    | nil => rfl
    | cons a l ih =>
      by_cases ha : cnt.getD a 0 = 0
      · rw [List.filter_cons_of_neg (by simpa using ha), ih, List.map_cons, List.sum_cons, ha]
        simp
      · rw [List.filter_cons_of_pos (by simpa using ha), List.map_cons, List.sum_cons, ih,
          List.map_cons, List.sum_cons]
  rw [← list_sum_range, ← hfilter, list_sum_map_eq_sum_getD]
  refine Finset.sum_congr rfl fun r hr ↦ ?_
  have hr' : r < (repCols cnt d).length := mem_range.1 hr
  rw [repCols] at hr' ⊢
  simp only [List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_eq_getElem hr',
    Option.map_some, Option.getD_some]

/-- **The units of a list of distinct columns.** -/
theorem colUnits_ok (P : Params) {cols : List ℕ} (hnd : cols.Nodup) (d : ℕ) :
    UnitsOk P d cols.length (fun r ↦ cols.getD r 0) (colUnits P.w cols d) := by
  refine ⟨by simp [colUnits], fun k hk ↦ ?_⟩
  rw [colUnits, getD_map_range _ _ _ _ hk]
  by_cases hm : k ∈ cols
  · have hi := List.idxOf_lt_length_of_mem hm
    rw [ite_eq_left hm, pack_congr (g := fun r ↦ if r = cols.idxOf k then 1 else 0),
      pack_single hi,
      one_mul]
    intro r hr
    simp only [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hr, Option.getD_some]
    congr 1
    apply propext
    constructor
    · rintro rfl; exact (hnd.idxOf_getElem r hr).symm
    · rintro rfl; exact List.getElem_idxOf hi
  · rw [ite_eq_right hm]
    symm
    rw [pack_congr (g := fun _ ↦ 0 * 0), pack_mul, zero_mul]
    intro r hr
    simp only [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hr, Option.getD_some]
    rw [ite_eq_right (fun h : cols[r] = k ↦ hm (h ▸ List.getElem_mem hr))]

/-- **The total as a sum over all words**: the representatives weighted by their multiplicities
are the representatives of all words, `W = ∑_x ∑_w X_{rep x}(w)²`. -/
theorem hornerTotal_eq_sum_rep (n : ℕ) (g : Fin 2) :
    hornerTotal n g = ∑ x : Fin n → Fin 3, ∑ w, thermalColumn n g (thermalRep x) w ^ 2 := by
  have h := Finset.sum_comp (s := univ) (fun ρ ↦ ∑ w, thermalColumn n g ρ w ^ 2)
    (thermalRep (n := n))
  simp only [nsmul_eq_mul] at h
  rw [h, hornerTotal]

/-! ### The packed total -/

/-- **A packed total** `∑_ρ ν_ρ ∑_k X_ρ(k)²` of a banded matrix `D` on `d` rows: `f` maps a row
number to the number of its representative, `ν_ρ` counts the rows with representative `ρ`, and
`X_ρ` is the column recursion of `D` with the coefficients `cs`. One packed run computes the
representative columns (the numbers with a nonzero count, lane `r` for column `cols[r]`), and
`numerator` weighs each by its count. The result is `none` if the band check or a range test
fails. -/
def packedTotal (P : Params) (D : Band) (f : ℕ → ℕ) (d : ℕ) (cs : List ℕ) : Option ℕ :=
  let cst := D.cst P d
  if D.check P cst d then
    let cnt := countArr d f
    let cols := repCols cnt d
    (runCols P D cst d (colUnits P.w cols d) cols.length cs).map fun ys ↦
      numerator P (List.replicate d 1) (cols.map (cnt.getD · 0)) ys
  else none

/-- **The packed total is the sum over all rows of the squared norm of their representative
column**, `W = ∑_{x < d} ∑_k X_{f x}(k)²`; the band check also gives `D.InRange d`. -/
theorem packedTotal_spec {P : Params} {D : Band} {f : ℕ → ℕ} {d : ℕ} {cs : List ℕ} {W : ℕ}
    (hP : P.Ok) (hcs : ∀ t ∈ cs, t ≤ P.B) (hf : ∀ x < d, f x < d)
    (h : packedTotal P D f d cs = some W) :
    D.InRange d ∧
      (W : ℤ) = ∑ x ∈ range d, ∑ k ∈ range d, specCol D P.s (f x) (cs.map (↑)) k ^ 2 := by
  simp only [packedTotal] at h
  split_ifs at h with hchk
  obtain ⟨ys, hrun, rfl⟩ := Option.map_eq_some_iff.1 h
  obtain ⟨hD, hR⟩ := Band.check_spec hchk
  refine ⟨hD, ?_⟩
  set cnt := countArr d f
  set cols := repCols cnt d
  have hnd : cols.Nodup := List.nodup_range.filter _
  have hN := runCols_numerator hP hD hR (colUnits_ok P hnd _) hcs hrun
    (u := List.replicate d 1) (v := cols.map (cnt.getD · 0)) (by simp) (by simp)
  -- `G j = ∑_k X_j(k)²`, the squared norm of the column `j` of the banded matrix.
  let G : ℕ → ℤ := fun j ↦ ∑ k ∈ range d, specCol D P.s j (cs.map (↑)) k ^ 2
  -- The numerator: `∑_r ν_r G(cols r)`.
  rw [Finset.sum_comm] at hN
  have hN' : ∀ r ∈ range cols.length,
      ∑ k ∈ range d, (((List.replicate d 1).getD k 0 : ℕ) : ℤ) *
        ((cols.map (cnt.getD · 0)).getD r 0 : ℤ) *
        specCol D P.s (cols.getD r 0) (cs.map (↑)) k ^ 2 =
      ((cols.map (cnt.getD · 0)).getD r 0 : ℤ) * G (cols.getD r 0) := fun r _ ↦ by
    simp only [G, Finset.mul_sum]
    refine Finset.sum_congr rfl fun k hk ↦ ?_
    rw [List.getD_eq_getElem?_getD, List.getElem?_replicate_of_lt (mem_range.1 hk),
      Option.getD_some]
    simp only [Nat.cast_one, one_mul]
  rw [Finset.sum_congr rfl hN', ← sum_range_eq_sum_repCols, ← sum_comp_eq_sum_count hf G] at hN
  exact hN

/-- The packing parameters: lanes of `64` bits, bias `2⁵⁷`, scale `32`, and the exact division
`⌊u / 32⌋ = ⌊u · 1 / 2⁵⌋`. -/
def thermalParams : Params := ⟨64, 57, 32, 5, 1⟩

/-- **The packed total** `W_{n,g} = ∑_ρ ν_ρ ∑_w X_ρ(w)²`: one packed run on the representative
columns of the banded matrix, or `none` if the band check or a range test fails. -/
def packedHornerTotal (n : ℕ) (g : Fin 2) : Option ℕ :=
  packedTotal thermalParams (ringBand n g) (repNum n) (3 ^ n) thermalCoeffs

theorem thermalParams_ok : thermalParams.Ok := by decide +kernel

theorem thermalCoeffs_le : ∀ t ∈ thermalCoeffs, t ≤ thermalParams.B := by decide +kernel

/-- **The packed total is the integer total.** If the packed run succeeds with the value `W`,
then `hornerTotal n g = W`. -/
theorem hornerTotal_eq_of_packed {n : ℕ} (hn : 2 ≤ n) {g : Fin 2} {W : ℕ}
    (h : packedHornerTotal n g = some W) : hornerTotal n g = W := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  obtain ⟨hD, hW⟩ := packedTotal_spec thermalParams_ok thermalCoeffs_le
    (fun _ hx ↦ repNum_lt hx) h
  rw [hornerTotal_eq_sum_rep, hW, ← sum_words (fun x ↦ ∑ k ∈ range (3 ^ (m + 2)),
    specCol (ringBand (m + 2) g) thermalParams.s (repNum (m + 2) x) (thermalCoeffs.map (↑)) k ^ 2)]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [repNum_wordNum, ← sum_words (fun k ↦ specCol (ringBand (m + 2) g) thermalParams.s
    (wordNum (m + 2) (thermalRep x)) (thermalCoeffs.map (↑)) k ^ 2)]
  exact Finset.sum_congr rfl fun w _ ↦ by rw [thermalColumn_eq hD]; rfl

end ThermalPacked

end CertifiedGaps.SpinOneHeisenberg
