/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.ThermalPacked
import CertifiedGaps.Heisenberg.ThermalTable60

/-!
# A fast evaluator of the `b = 21/2` totals for the twists `1` and `P`

`HornerTotals60` (`Heisenberg/ThermalTable60.lean`) says that `hornerTotal60 n g`, the integer
total `t_{n,g}` of `thermal120.tex`, equals the printed `paperTotal60 n g` for `n = 6, 8, 10` and
the twists `g = 1, P, C`. `hornerTotal60` runs the recursion over `ℤ[ω]` with vectors stored as
functions, so it cannot be evaluated in practice. For the twists `1` and `P` the matrix has
integer entries, and this file computes the same totals with the packed evaluator of
`Heisenberg/ThermalPacked.lean` (`hornerTotal60_eq_of_packed`).

## The reduction

* **Integer arithmetic.** For `g = 1, P` the ring `H_n(g)` over `ℤ[ω]` is the integer ring
  `ringZ n g` (`ringE_castSucc`), so `B = (16 - 3n/2) I - H_n(g)` is an integer matrix. On an
  integer matrix the recursion over `ℤ[ω]` keeps the `ω` coordinate zero
  (`hornerFloorE_map_intCast`), and `|u + 0ω|² = u²`.
* **The `b = 49/4` band.** For even `n`, `2B = (32 - 3n) I - 2H_n(g)` is the matrix `D` of the
  `b = 49/4` table (`thermalMatrix_eq_two_smul`), and `⌊2y / 32⌋ = ⌊y / 16⌋`
  (`hornerFloor_smul`). So the columns are the column recursions of the banded matrix
  `ringBand n g` with divisor `32` (`thermalColumn60_eq`); only the coefficients change.
* **The coefficients.** `q_j = ⌊2⁵⁶ p_j⌋ = ⌊2⁵⁶ v_j / ∑ v⌋` with `v_j = 84ʲ 256!/j!`
  (`thermalCoeff60_eq`); they vanish for `j ≥ 173` (`thermalCoeff60_support`), so the recursion
  starts at `j = 172`, as in the paper.
* **The representatives.** `repNum60 n x` is the number of `thermalRep60` of the word numbered
  `x` (`repNum60_wordNum`): the least code among the shifts and reversals of the word, or of its
  negation if its weight is negative.
* **The parameters** (`thermalParams60`): lanes of `65` bits, bias `2⁵⁸`, scale `32`. The range
  tests check `|z| < 2⁵⁷ = 2Q`.

## Main results

* `thermalColumn60_eq`: the columns `z_ρ` for `g = 1, P` are column recursions of `ringBand`.
* `hornerTotal60_eq_of_packed`:
  `2 ≤ n → Even n → packedHornerTotal60 n g = some W → hornerTotal60 n g.castSucc = W`.

## Comparison with the paper

* The paper evaluates `B z / 16` in signed 64-bit integers. Here the band is `2B` with the
  divisor `32`, which gives the same floors, so that one banded matrix and one proof of its
  entries serve both tables.
* The `C` twist has entries `±ω` on the seam bond and is not covered here; its totals remain in
  `HornerTotals60`.
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

end ThermalPacked

end CertifiedGaps.SpinOneHeisenberg
