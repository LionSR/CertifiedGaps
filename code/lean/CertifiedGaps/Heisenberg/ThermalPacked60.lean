/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalPacked
import CertifiedGaps.Heisenberg.ThermalTable60

/-!
# A fast evaluator of the `b = 21/2` totals

`HornerTotals60` (`Heisenberg/ThermalTable60.lean`) says that `hornerTotal60 n g`, the integer
total `t_{n,g}` of `thermal120.tex`, equals the printed `paperTotal60 n g` for `n = 6, 8, 10` and
the twists `g = 1, P, C`. `hornerTotal60` runs the recursion over `ℤ[ω]` with vectors stored as
functions, so it cannot be evaluated in practice. This file computes the same totals with the
packed evaluator of `Heisenberg/ThermalPacked.lean` and proves that they agree
(`hornerTotal60_eq_of_packedTotal60`, `hornerTotals60_of_packed`).

## The reduction

* **The `b = 49/4` band.** For even `n`, `2B = (32 - 3n) I - 2H_n(g)` and `⌊2y / 32⌋ = ⌊y / 16⌋`
  (`hornerFloor_smul`), so the recursion may run on `2B` with divisor `32`, the scale of the
  `b = 49/4` evaluator.
* **The twists `1` and `P`.** The ring `H_n(g)` over `ℤ[ω]` is the integer ring `ringZ n g`
  (`ringE_castSucc`), so the recursion keeps the `ω` coordinate zero
  (`hornerFloorE_map_intCast`), and `2B` is the matrix `D = ringBand n g` of the `b = 49/4`
  table (`thermalMatrix_eq_two_smul`, `thermalColumn60_eq`); only the coefficients change.
* **The twist `C`.** `2B` has the entries `-2ω` and `2 + 2ω` on the two seam hops. Writing
  `z = u + vω`, the recursion is the integer recursion of the real form `[[A, -I], [I, A - I]]`
  of `2B = A + Iω` on `(u, v)` (`hornerFloorE_realify`); the floor acts on each coordinate.
  The four blocks are ring bands (`seamBand`, `thermalMatrix60_two_entry`), placed as one band
  `ringBandC n` on `2 · 3ⁿ` rows (`realBand`, `realify_thermalMatrix60_two`). The total sums
  `|u + vω|² = u² - uv + v²` over the rows `x` and `3ⁿ + x` of each column (`packedTotalE`).
* **The coefficients.** `q_j = ⌊2⁵⁶ p_j⌋ = ⌊2⁵⁶ v_j / ∑ v⌋` with `v_j = 84ʲ 256!/j!`
  (`thermalCoeff60_eq`); they vanish for `j ≥ 173` (`thermalCoeff60_support`), so the recursion
  starts at `j = 172`, as in the paper.
* **The representatives.** `repNum60 n x` is the number of `thermalRep60` of the word numbered
  `x` (`repNum60_wordNum`): the least code among the shifts and reversals of the word, or of its
  negation if its weight is negative.
* **The parameters** (`thermalParams60`): lanes of `65` bits, bias `2⁵⁸`, scale `32`. The range
  tests check `|u|, |v| < 2⁵⁷`.

## Main results

* `thermalColumn60_eq`, `thermalColumn60_two_eq`: the columns `z_ρ` are column recursions of
  `ringBand` (twists `1`, `P`) and of `ringBandC` (twist `C`).
* `hornerTotal60_eq_of_packedTotal60`:
  `2 ≤ n → Even n → packedTotal60 n g = some W → hornerTotal60 n g = W`.
* `hornerTotals60_of_packed`: `HornerTotals60` follows from the nine packed totals.

## Comparison with the paper

* The paper evaluates `B z / 16` in signed 64-bit integers, with `u` and `v` as two integer
  vectors. Here the band is `2B` with the divisor `32`, which gives the same floors, so that the
  banded matrix of the `b = 49/4` table and the proof of its entries serve both tables.
-/

open Finset Matrix
open scoped Nat

namespace CertifiedGaps.SpinOneHeisenberg

open Trace Trace.Packed Trace.EisensteinInt

namespace ThermalPacked

/-! ### The matrix of the twists `1` and `P` -/

/-- For the twists `1` and `P` the ring over `ℤ[ω]` is the integer ring. -/
theorem ringE_castSucc (n : ℕ) (g : Fin 2) : ringE n g.castSucc = (ringZ n g).map (↑) := by
  have h := twistedRingHamiltonian_map exchangeZ (seamBond g) n (Int.castAddHom EisensteinInt)
  rw [Int.coe_castAddHom] at h
  rw [ringZ, h, ringE]
  congr 1
  refine Matrix.ext fun ⟨a, b⟩ ⟨c, e⟩ ↦ ?_
  rw [Matrix.map_apply]
  fin_cases g
  · exact (seamE_zero_one a b c e).1
  · exact (seamE_zero_one a b c e).2

/-- The matrix `B` of the twists `1` and `P` is an integer matrix. -/
theorem thermalMatrix60_castSucc (n : ℕ) (g : Fin 2) :
    thermalMatrix60 n g.castSucc =
      (thermalDiag60 n • (1 : Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℤ) - ringZ n g).map (↑) := by
  rw [thermalMatrix60, ringE_castSucc]
  refine Matrix.ext fun σ τ ↦ ?_
  simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply, Matrix.map_apply,
    smul_eq_mul]
  split_ifs <;> push_cast <;> ring

/-- **For even `n`, `2B` is the matrix `D` of the `b = 49/4` table.** -/
theorem thermalMatrix_eq_two_smul {n : ℕ} (he : Even n) (g : Fin 2) :
    thermalMatrix n g =
      ((2 : ℕ) : ℤ) • (thermalDiag60 n • (1 : Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℤ) -
        ringZ n g) := by
  obtain ⟨k, rfl⟩ := he
  have h : 32 - 3 * ((k + k : ℕ) : ℤ) = 2 * thermalDiag60 (k + k) := by
    rw [thermalDiag60]; push_cast; omega
  ext σ τ
  simp only [thermalMatrix, Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul,
    h]
  simp only [nsmul_eq_mul]
  push_cast
  ring

/-- `|u + 0ω|² = u²`. -/
theorem normSq_intCast (a : ℤ) : normSq (a : EisensteinInt) = a ^ 2 := by
  simp [normSq]

/-! ### The coefficients -/

/-- The numerators `v_j = 84ʲ 256!/j!` of the conditioned Poisson weights. -/
def poissonNum60 (j : ℕ) : ℕ := 84 ^ j * (256 ! / j !)

/-- `∑_{j ≤ 256} v_j`. -/
def poissonDen60 : ℕ := ((List.range 257).map poissonNum60).sum

/-- The coefficient `q_j = ⌊2⁵⁶ v_j / ∑ v⌋` in natural-number arithmetic. -/
def coeffNat60 (j : ℕ) : ℕ := 2 ^ 56 * poissonNum60 j / poissonDen60

/-- The coefficients `q_0, …, q_172`; the later ones vanish. -/
def thermalCoeffs60 : List ℕ := (List.range 173).map coeffNat60

/-- `p_j = v_j / ∑ v`: the factor `256!` cancels. -/
theorem poissonWeight_eq60 {j : ℕ} (hj : j ≤ 256) :
    poissonWeight 84 256 j = (poissonNum60 j : ℚ) / poissonDen60 := by
  have hcast : ∀ l ≤ 256, (poissonNum60 l : ℚ) = 84 ^ l / l ! * 256 ! := fun l hl ↦ by
    rw [poissonNum60, Nat.cast_mul, Nat.cast_div (Nat.factorial_dvd_factorial hl)
      (by positivity)]
    push_cast
    ring
  have hden : (poissonDen60 : ℚ) = poissonMass 84 256 * 256 ! := by
    rw [poissonDen60, list_sum_range, Nat.cast_sum, poissonMass, Finset.sum_mul]
    exact Finset.sum_congr rfl fun l hl ↦ hcast l (by simp at hl; omega)
  have h256 : (256 ! : ℚ) ≠ 0 := by positivity
  rw [poissonWeight, hcast j hj, hden, mul_div_mul_right _ _ h256]

/-- `q_j = ⌊2⁵⁶ p_j⌋` is the natural-number quotient. -/
theorem thermalCoeff60_eq {j : ℕ} (hj : j ≤ 256) : thermalCoeff60 j = coeffNat60 j := by
  rw [thermalCoeff60, poissonWeight_eq60 hj, thermalScale60, ← mul_div_assoc,
    show (2 : ℚ) ^ 56 * poissonNum60 j = ((2 ^ 56 * poissonNum60 j : ℕ) : ℚ) by push_cast; rfl,
    Rat.floor_natCast_div_natCast, coeffNat60]
  push_cast
  rfl

/-- The coefficient list of `thermalColumn60` is `thermalCoeffs60` followed by `84` zeros. -/
theorem coeff_list60 :
    (List.range (256 + 1)).map thermalCoeff60 =
      thermalCoeffs60.map (↑) ++ List.replicate 84 (0 : ℤ) := by
  refine List.ext_getElem (by simp [thermalCoeffs60]) fun j h1 h2 ↦ ?_
  have hj : j < 257 := by simpa using h1
  rw [List.getElem_map, List.getElem_range, List.getElem_append]
  split_ifs with h173
  · simp [thermalCoeffs60, thermalCoeff60_eq (by omega : j ≤ 256)]
  · have h173' : 173 ≤ j := by simpa [thermalCoeffs60] using h173
    rw [List.getElem_replicate,
      thermalCoeff60_support.2 j (List.mem_range'_1.2 ⟨h173', by omega⟩)]

/-- **The columns of the twists `1` and `P` are column recursions of the banded matrix**, with
the band `D = 2B` of the `b = 49/4` table and the divisor `32`. -/
theorem thermalColumn60_eq {g : Fin 2} {m : ℕ} (he : Even (m + 2))
    (hD : (ringBand (m + 2) g).InRange (3 ^ (m + 2))) (ρ w : Fin (m + 2) → Fin 3) :
    thermalColumn60 (m + 2) g.castSucc ρ w = ((specCol (ringBand (m + 2) g) 32
      (wordNum (m + 2) ρ) (thermalCoeffs60.map (↑)) (wordNum (m + 2) w) : ℤ) :
        EisensteinInt) := by
  rw [thermalColumn60, thermalMatrix60_castSucc, hornerFloorE_map_intCast,
    ← hornerFloor_smul _ (k := 2) (by norm_num), ← thermalMatrix_eq_two_smul he,
    show (2 * 16 : ℕ) = 32 from rfl]
  beta_reduce
  rw [hornerFloor_thermalMatrix hD, coeff_list60, specCol_append_zeros]

/-! ### Representatives on word numbers -/

/-- **The representative of `thermalRep60` on word numbers** (`repNum60_wordNum`): no negation
within a sector, unlike `repNum`. -/
def repNum60 (n x : ℕ) : ℕ :=
  let x₀ := if 0 ≤ weightNum n x then x else negNum n x
  ((dihedralNums n x₀).foldl (lexMinNum n) (codeNum n x₀, x₀)).2

/-- The representative on numbers is the number of the representative. -/
theorem repNum60_wordNum {n : ℕ} (σ : Fin n → Fin 3) :
    repNum60 n (wordNum n σ) = wordNum n (thermalRep60 σ) := by
  have hx₀ : (if 0 ≤ weightNum n (wordNum n σ) then (wordNum n σ : ℕ)
      else negNum n (wordNum n σ)) =
      wordNum n (if 0 ≤ wordWeight σ then σ else negWord n σ) := by
    rw [weightNum_wordNum]
    split_ifs
    · rfl
    · rw [wordNum_negWord]
  simp only [repNum60, thermalRep60]
  rw [hx₀, ← map_dihedralPerms, codeNum_wordNum, foldl_lexMinNum]

/-- The representative of a word number is a word number. -/
theorem repNum60_lt {n x : ℕ} (hx : x < 3 ^ n) : repNum60 n x < 3 ^ n := by
  have h := repNum60_wordNum ((wordNum n).symm ⟨x, hx⟩)
  simp only [Equiv.apply_symm_apply] at h
  rw [h]
  exact Fin.isLt _

/-- **The total as a sum over all words**, `t = ∑_x ∑_w |z_{rep x}(w)|²`. -/
theorem hornerTotal60_eq_sum_rep (n : ℕ) (g : Fin 3) :
    hornerTotal60 n g =
      ∑ x : Fin n → Fin 3, ∑ w, normSq (thermalColumn60 n g (thermalRep60 x) w) := by
  have h := Finset.sum_comp (s := univ) (fun ρ ↦ ∑ w, normSq (thermalColumn60 n g ρ w))
    (thermalRep60 (n := n))
  simp only [nsmul_eq_mul] at h
  rw [h, hornerTotal60]

/-! ### The packed total -/

/-- The packing parameters: lanes of `65` bits, bias `2⁵⁸`, scale `32`, and the exact division
`⌊u / 32⌋ = ⌊u · 1 / 2⁵⌋`. The columns reach `Q = 2⁵⁶`, one bit more than at `b = 49/4`. -/
def thermalParams60 : Params := ⟨65, 58, 32, 5, 1⟩

/-- **The packed total** `t_{n,g}` for the twists `1` and `P`: one packed run on the
representative columns of the banded matrix, or `none` if the band check or a range test
fails. -/
def packedHornerTotal60 (n : ℕ) (g : Fin 2) : Option ℕ :=
  packedTotal thermalParams60 (ringBand n g) (repNum60 n) (3 ^ n) thermalCoeffs60

theorem thermalParams60_ok : thermalParams60.Ok := by decide +kernel

theorem thermalCoeffs60_le : ∀ t ∈ thermalCoeffs60, t ≤ thermalParams60.B := by decide +kernel

/-- **The packed total is the integer total.** If the packed run for the twist `1` (`g = 0`) or
`P` (`g = 1`) succeeds with the value `W`, then `hornerTotal60 n g = W`. -/
theorem hornerTotal60_eq_of_packed {n : ℕ} (hn : 2 ≤ n) (he : Even n) {g : Fin 2} {W : ℕ}
    (h : packedHornerTotal60 n g = some W) : hornerTotal60 n g.castSucc = W := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  obtain ⟨hD, hW⟩ := packedTotal_spec thermalParams60_ok thermalCoeffs60_le
    (fun _ hx ↦ repNum60_lt hx) h
  rw [hornerTotal60_eq_sum_rep, hW, ← sum_words (fun x ↦ ∑ k ∈ range (3 ^ (m + 2)),
    specCol (ringBand (m + 2) g) thermalParams60.s (repNum60 (m + 2) x)
      (thermalCoeffs60.map (↑)) k ^ 2)]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [repNum60_wordNum, ← sum_words (fun k ↦ specCol (ringBand (m + 2) g) thermalParams60.s
    (wordNum (m + 2) (thermalRep60 x)) (thermalCoeffs60.map (↑)) k ^ 2)]
  exact Finset.sum_congr rfl fun w _ ↦ by rw [thermalColumn60_eq he hD, normSq_intCast]; rfl

/-! ### The twist `C`: the real form of `2B`

For the twist `C`, `2B = (32 - 3n) I - 2H_n(C)` has the entry `-2ω` on the seam hop
`n - 1 → 0` and `-2ω̄ = 2 + 2ω` on `0 → n - 1`; elsewhere it is the integer matrix of the open
chain. Its real form `[[A, -I], [I, A - I]]` (`realify`, with `2B = A + Iω`) has four blocks of
the form `seamBand`. -/

/-- The block `A` of the real form of `2B`: the integer part, with `+2` on the hop `0 → n - 1`. -/
def bandUU (n : ℕ) : Band := seamBand n true true 0 2

/-- The block `-I`: `+2` on the hop `n - 1 → 0`, `-2` on `0 → n - 1`. -/
def bandUV (n : ℕ) : Band := seamBand n false false 2 (-2)

/-- The block `I`, the `ω` part of `2B`: `-2` on the hop `n - 1 → 0`, `+2` on `0 → n - 1`. -/
def bandVU (n : ℕ) : Band := seamBand n false false (-2) 2

/-- The block `A - I`: `+2` on the hop `n - 1 → 0`. -/
def bandVV (n : ℕ) : Band := seamBand n true true 2 0

/-- **The real form of `2B` for the twist `C` as a band** on `2 · 3ⁿ` rows: row `x` holds the
coordinate `u` and row `3ⁿ + x` the coordinate `v` of the word numbered `x`. -/
def ringBandC (n : ℕ) : Band := realBand (3 ^ n) (bandUU n) (bandUV n) (bandVU n) (bandVV n)

/-- The bands of a ring band are hop bands. -/
theorem mem_seamBand_hops {n : ℕ} (hn : 2 ≤ n) {e f : Bool} {a b : ℤ} {h : Hop}
    (hh : h ∈ (seamBand n e f a b).hops ++ (seamBand n e f a b).phops) :
    ∃ (j k : Fin n) (c : ℕ), j ≠ k ∧ h = hopBand n j k c := by
  have hseam : ∀ c, (∃ (j k : Fin n) (c' : ℕ), j ≠ k ∧ hopBand n (n - 1) 0 c = hopBand n j k c') ∧
      ∃ (j k : Fin n) (c' : ℕ), j ≠ k ∧ hopBand n 0 (n - 1) c = hopBand n j k c' := fun c ↦
    ⟨⟨⟨n - 1, by omega⟩, ⟨0, by omega⟩, c, by simp [Fin.ext_iff]; omega, rfl⟩,
      ⟨⟨0, by omega⟩, ⟨n - 1, by omega⟩, c, by simp [Fin.ext_iff]; omega, rfl⟩⟩
  simp only [seamBand, negHop, posHop, List.mem_append] at hh
  rcases hh with ((hb | hA) | hB) | (hA | hB)
  · split_ifs at hb
    · obtain ⟨j, hj, hh⟩ := List.mem_flatMap.1 hb
      simp only [List.mem_range, List.mem_cons, List.not_mem_nil, or_false] at hj hh
      rcases hh with rfl | rfl
      · exact ⟨⟨j, by omega⟩, ⟨j + 1, by omega⟩, 2, by simp [Fin.ext_iff], rfl⟩
      · exact ⟨⟨j + 1, by omega⟩, ⟨j, by omega⟩, 2, by simp [Fin.ext_iff], rfl⟩
    · simp at hb
  all_goals split_ifs at * <;> simp only [List.mem_cons, List.not_mem_nil, or_false] at *
  · subst hA; exact (hseam _).1
  · subst hB; exact (hseam _).2
  · subst hA; exact (hseam _).1
  · subst hB; exact (hseam _).2

/-- **A ring band stays in range.** -/
theorem seamBand_inRange {n : ℕ} (hn : 2 ≤ n) (e f : Bool) (a b : ℤ) :
    (seamBand n e f a b).InRange (3 ^ n) := by
  intro h hh x hx hm
  obtain ⟨j, k, c, hjk, rfl⟩ := mem_seamBand_hops hn hh
  exact hopBand_inRange hjk c x hx hm

/-- A ring band without diagonal has zero diagonal. -/
theorem seamBand_diag_zero (n : ℕ) (f : Bool) (a b : ℤ) (k : ℕ) :
    (seamBand n false f a b).pos.getD k 0 = (seamBand n false f a b).neg.getD k 0 := by
  simp [seamBand]

/-- The numbering of the coordinates: `u` of word `σ` is row `wordNum σ`, `v` is row
`3ⁿ + wordNum σ`. -/
def coordNum (n : ℕ) : (Fin n → Fin 3) ⊕ (Fin n → Fin 3) ≃ Fin (3 ^ n + 3 ^ n) :=
  (Equiv.sumCongr (wordNum n) (wordNum n)).trans finSumFinEquiv

theorem coordNum_inl {n : ℕ} (σ : Fin n → Fin 3) :
    coordNum n (Sum.inl σ) = ⟨blk (3 ^ n) false (wordNum n σ), blk_lt false (wordNum n σ).isLt⟩ :=
  Fin.ext (by simp [coordNum, blk])

theorem coordNum_inr {n : ℕ} (σ : Fin n → Fin 3) :
    coordNum n (Sum.inr σ) = ⟨blk (3 ^ n) true (wordNum n σ), blk_lt true (wordNum n σ).isLt⟩ :=
  Fin.ext (by simp [coordNum, blk, add_comm])

@[simp] theorem uHom_apply (z : EisensteinInt) : uHom z = z.u := rfl

@[simp] theorem vHom_apply (z : EisensteinInt) : vHom z = z.v := rfl

/-- The coordinate `u` of the `C`-twisted seam bond: `mz mz` on the diagonal, `0` on the hop
lowering the first site, `-1` on the hop raising it. -/
theorem seamE_two_u (p q : Fin 3 × Fin 3) :
    (seamE 2).map uHom p q = (if p = q then 1 * (mz p.1 * mz p.2) else 0) +
      0 * (if p.1.val = q.1.val + 1 ∧ p.2.val + 1 = q.2.val then 1 else 0) +
      (-1) * (if p.1.val + 1 = q.1.val ∧ p.2.val = q.2.val + 1 then 1 else 0) := by
  obtain ⟨a, b⟩ := p
  obtain ⟨c, e⟩ := q
  revert a b c e
  decide +kernel

/-- The coordinate `v` of the `C`-twisted seam bond: `1` on the hop lowering the first site,
`-1` on the hop raising it. -/
theorem seamE_two_v (p q : Fin 3 × Fin 3) :
    (seamE 2).map vHom p q = (if p = q then 0 * (mz p.1 * mz p.2) else 0) +
      1 * (if p.1.val = q.1.val + 1 ∧ p.2.val + 1 = q.2.val then 1 else 0) +
      (-1) * (if p.1.val + 1 = q.1.val ∧ p.2.val = q.2.val + 1 then 1 else 0) := by
  obtain ⟨a, b⟩ := p
  obtain ⟨c, e⟩ := q
  revert a b c e
  decide +kernel

/-- `2 · (16 - 3n/2) = 32 - 3n` for even `n`. -/
theorem two_mul_thermalDiag60 {n : ℕ} (he : Even n) : 2 * thermalDiag60 n = 32 - 3 * (n : ℤ) := by
  obtain ⟨k, rfl⟩ := he
  rw [thermalDiag60]
  push_cast
  omega

/-- **The entries of `2B` for the twist `C`**, coordinate by coordinate. -/
theorem thermalMatrix60_two_entry {m : ℕ} (he : Even (m + 2)) (σ τ : Fin (m + 2) → Fin 3) :
    2 * (thermalMatrix60 (m + 2) 2 σ τ).u =
        (bandUU (m + 2)).matrix (3 ^ (m + 2)) (wordNum (m + 2) σ) (wordNum (m + 2) τ) ∧
      2 * (thermalMatrix60 (m + 2) 2 σ τ).v =
        (bandVU (m + 2)).matrix (3 ^ (m + 2)) (wordNum (m + 2) σ) (wordNum (m + 2) τ) := by
  have hne : ∀ i : Fin (m + 1), (i.castSucc : Fin (m + 2)) ≠ i.succ := fun i ↦ by
    simp [Fin.ext_iff]
  have hl0 : (Fin.last (m + 1) : Fin (m + 2)) ≠ 0 := by simp [Fin.ext_iff]
  -- The two coordinates of the ring.
  have hmap := fun (f : EisensteinInt →+ ℤ) ↦
    congrFun (congrFun (twistedRingHamiltonian_map exchangeE (seamE 2) (m + 2) f) σ) τ
  have hxu : exchangeE.map uHom = exchangeZ := by ext p q; simp [exchangeE]
  have hxv : exchangeE.map vHom = 0 := by ext p q; simp [exchangeE]
  have hu : (ringE (m + 2) 2 σ τ).u =
      twistedRingHamiltonian exchangeZ ((seamE 2).map uHom) (m + 2) σ τ := by
    rw [← hxu]; exact hmap uHom
  have hv : (ringE (m + 2) 2 σ τ).v =
      twistedRingHamiltonian 0 ((seamE 2).map vHom) (m + 2) σ τ := by
    rw [← hxv]; exact hmap vHom
  rw [twistedRing_apply, bondOp_apply_of_flip₂ hl0 seamE_two_u] at hu
  rw [twistedRing_apply, bondOp_apply_of_flip₂ hl0 seamE_two_v] at hv
  simp only [bondOp_exchangeZ_apply (hne _)] at hu
  have h0 : ∀ i : Fin (m + 1), bondOp i.castSucc i.succ
      (0 : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℤ) σ τ = 0 := fun i ↦ by
    simp [bondOp]
  simp only [h0, Finset.sum_const_zero, zero_add] at hv
  have hB : thermalMatrix60 (m + 2) 2 σ τ =
      (if σ = τ then (thermalDiag60 (m + 2) : EisensteinInt) else 0) - ringE (m + 2) 2 σ τ := by
    simp only [thermalMatrix60, Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply,
      smul_eq_mul]
    split_ifs <;> simp
  have ite2 : ∀ (c : Prop) [Decidable c], (if c then (2 : ℤ) else 0) = 2 * if c then 1 else 0 :=
    fun c _ ↦ by split_ifs <;> simp
  have h2 := two_mul_thermalDiag60 he
  rw [hB, bandUU, bandVU, seamBand_entry, seamBand_entry]
  constructor
  · rw [sub_u, hu]
    by_cases hστ : σ = τ
    · subst hστ
      simp only [ite_true, hopUp_ne_self (hne _), hopUp_ne_self (hne _).symm, hopUp_ne_self hl0,
        hopUp_ne_self hl0.symm, ite_false, add_zero, mul_zero, Finset.sum_const_zero, sub_zero,
        intCast_u, ringDiag_wordNum]
      push_cast at h2 ⊢
      linarith
    · simp only [hστ, ite_false, ite_true, zero_u, zero_add, Finset.sum_add_distrib, ite2,
        ← Finset.mul_sum]
      ring
  · rw [sub_v, hv]
    by_cases hστ : σ = τ
    · subst hστ
      simp [hopUp_ne_self hl0, hopUp_ne_self hl0.symm]
    · simp only [hστ, ite_false, zero_v, Bool.false_eq_true, sub_zero]
      ring

/-- **The real form of `2B` for the twist `C` is the band `ringBandC`**, relabelled by the
coordinate numbering. -/
theorem realify_thermalMatrix60_two {m : ℕ} (he : Even (m + 2)) :
    ((2 : ℕ) : ℤ) • realify (thermalMatrix60 (m + 2) 2) =
      ((ringBandC (m + 2)).matrix (3 ^ (m + 2) + 3 ^ (m + 2))).submatrix (coordNum (m + 2))
        (coordNum (m + 2)) := by
  have hr := seamBand_inRange (n := m + 2) (by omega)
  have hentry := realBand_entry (hr true true 0 2) (hr false false 2 (-2)) (hr false false (-2) 2)
    (hr true true 2 0) (seamBand_diag_zero _ false 2 (-2)) (seamBand_diag_zero _ false (-2) 2)
  ext (σ | σ) (τ | τ) <;>
    simp only [Matrix.smul_apply, realify, Matrix.fromBlocks_apply₁₁, Matrix.fromBlocks_apply₁₂,
      Matrix.fromBlocks_apply₂₁, Matrix.fromBlocks_apply₂₂, Matrix.submatrix_apply, coordNum_inl,
      coordNum_inr, Matrix.map_apply, Matrix.neg_apply, Matrix.sub_apply, uHom_apply, vHom_apply,
      smul_eq_mul, ringBandC, bandUU, bandUV, bandVU, bandVV] <;>
    rw [hentry] <;>
    simp only [pickBlock, ite_true, ite_false, Bool.false_eq_true] <;>
    obtain ⟨hu, hv⟩ := thermalMatrix60_two_entry he σ τ
  · push_cast; exact hu
  · rw [seamBand_entry]
    rw [bandVU, seamBand_entry] at hv
    push_cast
    simp only [Bool.false_eq_true, ite_false, ite_self, sub_zero] at hv ⊢
    linarith
  · push_cast; exact hv
  · rw [seamBand_entry]
    rw [bandUU, seamBand_entry] at hu
    rw [bandVU, seamBand_entry] at hv
    push_cast
    simp only [Bool.false_eq_true, ite_false, ite_true, ite_self, sub_zero] at hu hv ⊢
    linarith

/-- **The columns of the twist `C` are column recursions of `ringBandC`**: the coordinates
`u, v` of `z_ρ(w)` are rows `wordNum w` and `3ⁿ + wordNum w`. -/
theorem thermalColumn60_two_eq {m : ℕ} (he : Even (m + 2))
    (hD : (ringBandC (m + 2)).InRange (3 ^ (m + 2) + 3 ^ (m + 2))) (ρ w : Fin (m + 2) → Fin 3) :
    thermalColumn60 (m + 2) 2 ρ w =
      ⟨specCol (ringBandC (m + 2)) 32 (wordNum (m + 2) ρ) (thermalCoeffs60.map (↑))
          (wordNum (m + 2) w),
        specCol (ringBandC (m + 2)) 32 (wordNum (m + 2) ρ) (thermalCoeffs60.map (↑))
          (3 ^ (m + 2) + wordNum (m + 2) w)⟩ := by
  rw [thermalColumn60, hornerFloorE_realify,
    ← hornerFloor_smul _ (k := 2) (by norm_num), show (2 * 16 : ℕ) = 32 from rfl,
    realify_thermalMatrix60_two he, hornerFloor_submatrix, hornerFloor_submatrix,
    hornerFloor_eq_specCol hD, hornerFloor_eq_specCol hD, coeff_list60, specCol_append_zeros,
    coordNum_inl, coordNum_inr, coordNum_inl]
  rfl

/-- The packed total for the twist `C`: the norms of the representative columns of `ringBandC`. -/
def packedHornerTotal60C (n : ℕ) : Option ℤ :=
  packedTotalE thermalParams60 (ringBandC n) (repNum60 n) (3 ^ n) thermalCoeffs60

/-- **The packed total for the twist `C` is the integer total.** -/
theorem hornerTotal60_two_eq_of_packed {n : ℕ} (hn : 2 ≤ n) (he : Even n) {W : ℤ}
    (h : packedHornerTotal60C n = some W) : hornerTotal60 n 2 = W := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  obtain ⟨hD, hW⟩ := packedTotalE_spec thermalParams60_ok thermalCoeffs60_le
    (fun _ hx ↦ repNum60_lt hx) h
  rw [hornerTotal60_eq_sum_rep, hW, ← sum_words (fun x ↦ ∑ k ∈ range (3 ^ (m + 2)),
    normSqZ (specCol (ringBandC (m + 2)) thermalParams60.s (repNum60 (m + 2) x)
      (thermalCoeffs60.map (↑)) k) (specCol (ringBandC (m + 2)) thermalParams60.s
        (repNum60 (m + 2) x) (thermalCoeffs60.map (↑)) (3 ^ (m + 2) + k)))]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [repNum60_wordNum, ← sum_words (fun k ↦ normSqZ (specCol (ringBandC (m + 2))
    thermalParams60.s (wordNum (m + 2) (thermalRep60 x)) (thermalCoeffs60.map (↑)) k)
    (specCol (ringBandC (m + 2)) thermalParams60.s (wordNum (m + 2) (thermalRep60 x))
      (thermalCoeffs60.map (↑)) (3 ^ (m + 2) + k)))]
  exact Finset.sum_congr rfl fun w _ ↦ by rw [thermalColumn60_two_eq he hD]; rfl

/-! ### All three twists -/

/-- **The packed total for each twist** `1`, `P`, `C`. -/
def packedTotal60 (n : ℕ) : Fin 3 → Option ℤ :=
  ![(packedHornerTotal60 n 0).map (↑), (packedHornerTotal60 n 1).map (↑), packedHornerTotal60C n]

/-- **The packed total is the integer total**, for every twist. -/
theorem hornerTotal60_eq_of_packedTotal60 {n : ℕ} (hn : 2 ≤ n) (he : Even n) {g : Fin 3}
    {W : ℤ} (h : packedTotal60 n g = some W) : hornerTotal60 n g = W := by
  fin_cases g
  · obtain ⟨W', h', rfl⟩ := Option.map_eq_some_iff.1 h
    exact hornerTotal60_eq_of_packed (g := 0) hn he h'
  · obtain ⟨W', h', rfl⟩ := Option.map_eq_some_iff.1 h
    exact hornerTotal60_eq_of_packed (g := 1) hn he h'
  · exact hornerTotal60_two_eq_of_packed hn he h

/-- **`HornerTotals60` from the packed evaluator**: if the packed totals are the printed ones at
the nine points, the hypothesis holds. -/
theorem hornerTotals60_of_packed
    (h : ∀ n ∈ tableLengths60, ∀ g, packedTotal60 n g = some (paperTotal60 n g)) :
    HornerTotals60 := by
  intro n hn g
  have hn' : 2 ≤ n ∧ Even n := by
    simp only [tableLengths60, List.mem_cons, List.not_mem_nil, or_false] at hn
    rcases hn with rfl | rfl | rfl <;> decide
  exact hornerTotal60_eq_of_packedTotal60 hn'.1 hn'.2 (h n hn g)

end ThermalPacked

end CertifiedGaps.SpinOneHeisenberg
