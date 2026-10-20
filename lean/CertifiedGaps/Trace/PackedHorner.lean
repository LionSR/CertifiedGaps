/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Trace.Horner
import CertifiedGaps.Trace.PackedLanes

/-!
# Integer Horner recursions of banded matrices, all columns at once

The integer Horner recursion `X ← ⌊D X / s⌋ + t I` of `Trace/Horner.lean` (`hornerFloor`, one
column per call) is evaluated here for all `d` columns simultaneously, for a banded integer
matrix `D` on `{0, …, d - 1}`:

  `(D x)_k = (p_k - q_k) x_k - ∑_h m_h(k) x_{k ± o_h}`,

with diagonal `p - q` (`p, q ≥ 0`) and off-diagonal bands `h` (a `Hop`) of offset `o_h` and
nonpositive entries `-m_h(k)`. Row `k` of the `d × d` matrix `X` is stored as one packed number
(`Trace/PackedLanes.lean`) with lane `r` holding `X_{k r} + B`, `B = 2^b`. Then

* the product `D X` needs one scalar multiple per band and row: the band `h` of `D X` in row `k`
  is `m_h(k)` times row `k ± o_h` (`negAcc`, built from shifted lists of rows);
* the floor division by `s` is one multiplication by a magic constant `M`, one shift and one
  mask per row (`magic_div`, `pack_shiftRight_and`);
* the identity term adds `t` to lane `k` of row `k`;
* one `&&&` test per row checks that every new entry lies in `[-2^(b-1), 2^(b-1))`
  (`lanes_mem_of_and_eq`); the bias `B` and the lane width are chosen so that no lane overflows
  as long as these tests pass.

The recursion may run on a window of `nc` consecutive columns `col₀, …, col₀ + nc - 1` (the
lanes of a row are then these columns); the whole matrix is the window `col₀ = 0`, `nc = d`.
It returns `none` as soon as a test fails. `run_spec` proves that if it returns `some ys`, then
lane `r` of row `k` is `hornerFloor D s t (col₀ + r) K k + B` for every `k < d`, `r < nc`
(through the column recursion `specCol`), so that `decide +kernel` on a packed run certifies all
columns of `hornerFloor` with a number of kernel steps proportional to `d K` instead of `d² K`.
The kernel keeps every intermediate packed row until the end of a declaration, so the memory of
a check grows like the number of kernel steps (about `10` GB for `d = 243`, `K = 192`).

## Main definitions

* `Hop`, `Band`: the banded matrix, and `Band.matrix` its `Fin d × Fin d` integer matrix.
* `Params`: lane width `w`, bias exponent `b`, scale `s`, magic shift `S` and multiplier `M`.
* `run`: the packed recursion; `numerator`: `∑_k u_k ∑_r v_r (lane_r - B)²` over the rows.
* `Band.check`: the conditions on the band data as one Boolean, linear in the data.

## Main results

* `Band.matrix_mulVec`: the matrix acts by the band formula; `hornerFloor_eq_specCol`.
* `run_spec`, `run_lane`: the packed recursion computes `hornerFloor` of `Band.matrix`.
* `run_numerator`: the weighted sums of squares of the computed columns.
* `hornerFloor_submatrix`: relabelling the index set.
-/

namespace CertifiedGaps.Trace.Packed

open Finset Matrix

/-! ### Banded matrices -/

/-- One band of a banded matrix: row `k` couples to row `k + off` (`up`) or `k - off`, with the
nonpositive entry `-(mask k)`. -/
structure Hop where
  /-- Whether the band lies above the diagonal. -/
  up : Bool
  /-- The distance of the band from the diagonal. -/
  off : ℕ
  /-- The magnitudes of the band entries, row by row. -/
  mask : List ℕ

/-- The column coupled to row `k` by the band. -/
def Hop.nb (h : Hop) (k : ℕ) : ℕ := if h.up then k + h.off else k - h.off

/-- A banded integer matrix: diagonal `pos - neg` and bands `hops`. -/
structure Band where
  /-- The positive parts of the diagonal entries. -/
  pos : List ℕ
  /-- The negative parts of the diagonal entries. -/
  neg : List ℕ
  /-- The off-diagonal bands. -/
  hops : List Hop

/-- `(D x)_k = (p_k - q_k) x_k - ∑_h m_h(k) x_{nb_h(k)}`. -/
def Band.apply (D : Band) (x : ℕ → ℤ) (k : ℕ) : ℤ :=
  ((D.pos.getD k 0 : ℤ) - D.neg.getD k 0) * x k -
    (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * x (h.nb k)).sum

/-- The banded matrix on `Fin d`. -/
def Band.matrix (D : Band) (d : ℕ) : Matrix (Fin d) (Fin d) ℤ := fun k l ↦
  (if k = l then (D.pos.getD k 0 : ℤ) - D.neg.getD k 0 else 0) -
    (D.hops.map fun h ↦ if (l : ℕ) = h.nb k then (h.mask.getD k 0 : ℤ) else 0).sum

/-- Every band stays inside `{0, …, d - 1}`: a nonzero entry in row `k < d` couples to a
column `< d`. -/
def Band.InRange (D : Band) (d : ℕ) : Prop :=
  ∀ h ∈ D.hops, ∀ k < d, h.mask.getD k 0 ≠ 0 →
    (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)

/-- A vector on `Fin d`, extended by zero to `ℕ`. -/
def extend {d : ℕ} (X : Fin d → ℤ) (n : ℕ) : ℤ := if h : n < d then X ⟨n, h⟩ else 0

theorem extend_fin {d : ℕ} (X : Fin d → ℤ) (k : Fin d) : extend X k = X k := by
  simp [extend, k.isLt]

/-- The neighbour of a row in range is in range. -/
theorem Band.InRange.nb_lt {D : Band} {d : ℕ} (hD : D.InRange d) {h : Hop} (hh : h ∈ D.hops)
    {k : ℕ} (hk : k < d) (hm : h.mask.getD k 0 ≠ 0) : h.nb k < d := by
  have := hD h hh k hk hm
  unfold Hop.nb
  cases hu : h.up
  · simp only [Bool.false_eq_true, ite_false]; omega
  · simp only [ite_true]; exact this.1 hu

/-- Summing against the indicator of one index `n` picks `extend X n`. -/
theorem sum_ite_val_eq {d : ℕ} (X : Fin d → ℤ) (n : ℕ) (c : ℤ) :
    ∑ l : Fin d, (if (l : ℕ) = n then c else 0) * X l = c * extend X n := by
  by_cases hn : n < d
  · rw [Finset.sum_eq_single ⟨n, hn⟩ (fun l _ hl ↦ by
      have : (l : ℕ) ≠ n := fun e ↦ hl (Fin.ext e)
      simp [this]) (by simp)]
    simp [extend, hn]
  · have : ∀ l : Fin d, (l : ℕ) ≠ n := fun l e ↦ hn (e ▸ l.isLt)
    simp [this, extend, hn]

/-- Exchanging a finite sum with the sum over a list. -/
theorem sum_list_sum_mul {d : ℕ} (L : List Hop) (g : Hop → Fin d → ℤ) (X : Fin d → ℤ) :
    ∑ l : Fin d, (L.map fun h ↦ g h l).sum * X l = (L.map fun h ↦ ∑ l, g h l * X l).sum := by
  induction L with
  | nil => simp
  | cons h L ih =>
    simp only [List.map_cons, List.sum_cons, add_mul, Finset.sum_add_distrib, ih]

/-- **The banded matrix acts by the band formula.** -/
theorem Band.matrix_mulVec (D : Band) {d : ℕ} (X : Fin d → ℤ) (k : Fin d) :
    (D.matrix d *ᵥ X) k = D.apply (extend X) k := by
  simp only [mulVec, dotProduct, Band.matrix, sub_mul, Finset.sum_sub_distrib, ite_mul,
    zero_mul, Finset.sum_ite_eq, Finset.mem_univ, ite_true, Band.apply, extend_fin]
  rw [sum_list_sum_mul]
  simp only [sum_ite_val_eq]

/-- Changing a vector away from `{0, …, d - 1}` does not change the band formula in a row
`k < d`, when the bands stay in range. -/
theorem Band.apply_congr {D : Band} {d : ℕ} (hD : D.InRange d) {x y : ℕ → ℤ}
    (hxy : ∀ n < d, x n = y n) {k : ℕ} (hk : k < d) : D.apply x k = D.apply y k := by
  unfold Band.apply
  rw [hxy k hk]
  congr 1
  refine congrArg List.sum (List.map_congr_left fun h hh ↦ ?_)
  by_cases hm : h.mask.getD k 0 = 0
  · rw [hm, Nat.cast_zero, zero_mul, zero_mul]
  · rw [hxy _ (hD.nb_lt hh hk hm)]

/-! ### The column recursion -/

/-- One step of column `r` of the integer Horner recursion: `x ↦ ⌊D x / s⌋ + t e_r`. -/
def specStep (D : Band) (s r : ℕ) (t : ℤ) (x : ℕ → ℤ) : ℕ → ℤ :=
  fun k ↦ D.apply x k / s + if k = r then t else 0

/-- Column `r` of the integer Horner recursion with coefficients `cs` (last one first). -/
def specCol (D : Band) (s r : ℕ) (cs : List ℤ) : ℕ → ℤ := cs.foldr (specStep D s r) 0

/-- **The column recursion is `hornerList` of the banded matrix.** -/
theorem hornerList_eq_specCol {D : Band} {d : ℕ} (hD : D.InRange d) (s : ℕ) (r : Fin d)
    (cs : List ℤ) (k : Fin d) : hornerList (D.matrix d) s r cs k = specCol D s r cs k := by
  induction cs generalizing k with
  | nil => rfl
  | cons t cs ih =>
    simp only [hornerList_cons, hornerStep, specCol, List.foldr_cons, specStep]
    rw [Band.matrix_mulVec]
    congr 2
    · refine Band.apply_congr hD (fun n hn ↦ ?_) k.isLt
      simp only [extend, hn, ↓reduceDIte, ih]; rfl
    · simp [Fin.ext_iff]

/-- `hornerFloor` of the banded matrix in terms of the column recursion. -/
theorem hornerFloor_eq_specCol {D : Band} {d : ℕ} (hD : D.InRange d) (s : ℕ) (c : ℕ → ℤ)
    (r : Fin d) (K : ℕ) (k : Fin d) :
    hornerFloor (D.matrix d) s c r K k = specCol D s r ((List.range (K + 1)).map c) k :=
  hornerList_eq_specCol hD s r _ k

/-- `hornerList` commutes with relabelling the index set. -/
theorem hornerList_submatrix {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] (M : Matrix κ κ ℤ) (e : ι ≃ κ) (s : ℕ) (r : ι) (l : List ℤ) (w : ι) :
    hornerList (M.submatrix e e) s r l w = hornerList M s (e r) l (e w) := by
  induction l generalizing w with
  | nil => rfl
  | cons c l ih =>
    simp only [hornerList_cons, hornerStep]
    have hX : hornerList (M.submatrix e e) s r l = hornerList M s (e r) l ∘ e := funext ih
    rw [hX, submatrix_mulVec_equiv]
    simp [Function.comp_def, e.injective.eq_iff]

/-- `hornerFloor` commutes with relabelling the index set. -/
theorem hornerFloor_submatrix {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] (M : Matrix κ κ ℤ) (e : ι ≃ κ) (s : ℕ) (c : ℕ → ℤ) (r : ι) (K : ℕ)
    (w : ι) : hornerFloor (M.submatrix e e) s c r K w = hornerFloor M s c (e r) K (e w) :=
  hornerList_submatrix M e s r _ w

/-! ### List helpers -/

/-- `a + m y`, entrywise. -/
def axpy : List ℕ → List ℕ → List ℕ → List ℕ
  | a :: as, m :: ms, y :: ys => (a + m * y) :: axpy as ms ys
  | _, _, _ => []

/-- `m y`, entrywise. -/
def mulList : List ℕ → List ℕ → List ℕ
  | m :: ms, y :: ys => (m * y) :: mulList ms ys
  | _, _ => []

/-- The rows coupled by a band: entry `k` is row `nb k` (zero outside the list). -/
def Hop.shift (h : Hop) (ys : List ℕ) : List ℕ :=
  if h.up then ys.drop h.off ++ List.replicate h.off 0 else List.replicate h.off 0 ++ ys

theorem length_axpy : ∀ a m y : List ℕ,
    (axpy a m y).length = min (min a.length m.length) y.length
  | _ :: _, _ :: _, _ :: _ => by simp [axpy, length_axpy]
  | [], _, _ => by simp [axpy]
  | _ :: _, [], _ => by simp [axpy]
  | _ :: _, _ :: _, [] => by simp [axpy]

theorem getD_axpy : ∀ (a m y : List ℕ) (k : ℕ), k < a.length → k < m.length → k < y.length →
    (axpy a m y).getD k 0 = a.getD k 0 + m.getD k 0 * y.getD k 0
  | _ :: _, _ :: _, _ :: _, 0, _, _, _ => by simp [axpy]
  | _ :: as, _ :: ms, _ :: ys, k + 1, h1, h2, h3 => by
    simp only [axpy, List.getD_cons_succ]
    exact getD_axpy as ms ys k (by simpa using h1) (by simpa using h2) (by simpa using h3)

theorem length_mulList : ∀ m y : List ℕ, (mulList m y).length = min m.length y.length
  | _ :: _, _ :: _ => by simp [mulList, length_mulList]
  | [], _ => by simp [mulList]
  | _ :: _, [] => by simp [mulList]

theorem getD_mulList : ∀ (m y : List ℕ) (k : ℕ), k < m.length → k < y.length →
    (mulList m y).getD k 0 = m.getD k 0 * y.getD k 0
  | _ :: _, _ :: _, 0, _, _ => by simp [mulList]
  | _ :: ms, _ :: ys, k + 1, h1, h2 => by
    simp only [mulList, List.getD_cons_succ]
    exact getD_mulList ms ys k (by simpa using h1) (by simpa using h2)

theorem length_shift (h : Hop) (ys : List ℕ) : ys.length ≤ (h.shift ys).length := by
  unfold Hop.shift
  split_ifs
  · simp; omega
  · simp

/-- The shifted list holds the neighbouring rows wherever the band entry is nonzero. -/
theorem mask_mul_getD_shift {d : ℕ} {h : Hop} (ys : List ℕ) (hys : ys.length = d) {k : ℕ}
    (hr : h.mask.getD k 0 ≠ 0 → (h.up = true → k + h.off < d) ∧
      (h.up = false → h.off ≤ k)) :
    h.mask.getD k 0 * (h.shift ys).getD k 0 = h.mask.getD k 0 * ys.getD (h.nb k) 0 := by
  by_cases hm : h.mask.getD k 0 = 0
  · rw [hm, zero_mul, zero_mul]
  obtain ⟨hup, hdn⟩ := hr hm
  congr 1
  unfold Hop.shift Hop.nb
  cases hu : h.up
  · have := hdn hu
    simp only [Bool.false_eq_true, ite_false, List.getD_eq_getElem?_getD]
    rw [List.getElem?_append_right (by simp; omega)]
    simp only [List.length_replicate]
  · have := hup hu
    simp only [ite_true, List.getD_eq_getElem?_getD]
    rw [List.getElem?_append_left (by simp; omega), List.getElem?_drop, Nat.add_comm]

/-- The accumulated band products `q_k y_k + ∑_h m_h(k) y_{nb_h(k)}` of all rows. -/
def negAcc (D : Band) (ys : List ℕ) : List ℕ :=
  D.hops.foldl (fun acc h ↦ axpy acc h.mask (h.shift ys)) (mulList D.neg ys)

theorem foldl_axpy_spec (ys : List ℕ) : ∀ (L : List Hop) (acc : List ℕ),
    (∀ h ∈ L, acc.length ≤ h.mask.length) → acc.length ≤ ys.length →
    (L.foldl (fun acc h ↦ axpy acc h.mask (h.shift ys)) acc).length = acc.length ∧
      ∀ k < acc.length, (L.foldl (fun acc h ↦ axpy acc h.mask (h.shift ys)) acc).getD k 0 =
        acc.getD k 0 + (L.map fun h ↦ h.mask.getD k 0 * (h.shift ys).getD k 0).sum
  | [], acc, _, _ => by simp
  | h :: L, acc, hm, hy => by
    have hsh := length_shift h ys
    have hlen : (axpy acc h.mask (h.shift ys)).length = acc.length := by
      rw [length_axpy]
      have := hm h (by simp)
      omega
    have ih := foldl_axpy_spec ys L (axpy acc h.mask (h.shift ys))
      (fun h' hh' ↦ hlen ▸ hm h' (by simp [hh'])) (hlen ▸ hy)
    rw [hlen] at ih
    refine ⟨ih.1, fun k hk ↦ ?_⟩
    rw [List.foldl_cons, ih.2 k hk, getD_axpy _ _ _ k hk
      (lt_of_lt_of_le hk (hm h (by simp))) (by omega), List.map_cons, List.sum_cons]
    ring

/-! ### Lanes of the band products -/

/-- The packed recursion's parameters: lane width `w`, bias `B = 2^b`, scale `s`, and the magic
shift `S` and multiplier `M` for the floor division by `s`. -/
structure Params where
  /-- The lane width. -/
  w : ℕ
  /-- The bias exponent: lanes hold `x + 2^b`. -/
  b : ℕ
  /-- The scale `s` of the recursion `X ← ⌊D X / s⌋ + t I`. -/
  s : ℕ
  /-- The shift of the magic division. -/
  S : ℕ
  /-- The multiplier of the magic division. -/
  M : ℕ

namespace Params

variable (P : Params)

/-- The bias `B = 2^b`. -/
def B : ℕ := 2 ^ P.b

/-- Entries stay in `[-2^(b-1), 2^(b-1))`. -/
def half : ℕ := 2 ^ (P.b - 1)

/-- The arithmetic conditions under which no lane overflows and the magic division is exact
for lanes up to `2 s B`. -/
def Ok : Prop :=
  0 < P.s ∧ 1 ≤ P.b ∧ P.b < P.w ∧ P.S ≤ P.w ∧ 2 ^ P.S ≤ P.M * P.s ∧
    2 * P.s * P.B * (P.M * P.s - 2 ^ P.S) < 2 ^ P.S ∧ 2 * P.s * P.B * P.M < 2 ^ P.w ∧
    4 * P.B < 2 ^ P.w

instance : Decidable P.Ok := by unfold Ok; infer_instance

theorem B_eq (h : 1 ≤ P.b) : P.B = 2 * P.half := by
  rw [B, half, ← pow_succ']; congr 1; omega

end Params

/-- **Floor division by a magic multiplication**: `⌊u M / 2^S⌋ = ⌊u / s⌋` for `u ≤ U` when
`M s ≥ 2^S` and `U (M s - 2^S) < 2^S`. -/
theorem magic_div {s S M U u : ℕ} (hs : 0 < s) (h1 : 2 ^ S ≤ M * s)
    (h2 : U * (M * s - 2 ^ S) < 2 ^ S) (hu : u ≤ U) : u * M / 2 ^ S = u / s := by
  refine Nat.div_eq_of_lt_le ?_ ?_
  · calc u / s * 2 ^ S ≤ u / s * (M * s) := Nat.mul_le_mul_left _ h1
      _ = u / s * s * M := by ring
      _ ≤ u * M := Nat.mul_le_mul_right _ (Nat.div_mul_le_self u s)
  · -- `u M s = u 2^S + u e < (u + 1) 2^S ≤ (⌊u/s⌋ + 1) s 2^S`.
    set e := M * s - 2 ^ S
    have hMs : M * s = 2 ^ S + e := by omega
    have hue : u * e < 2 ^ S := lt_of_le_of_lt (Nat.mul_le_mul_right _ hu) h2
    have hq : u + 1 ≤ (u / s + 1) * s := by
      have := Nat.lt_div_mul_add (a := u) hs
      rw [add_mul, one_mul]; omega
    have key : u * M * s < (u / s + 1) * 2 ^ S * s := by
      calc u * M * s = u * (M * s) := by ring
        _ = u * 2 ^ S + u * e := by rw [hMs, mul_add]
        _ < u * 2 ^ S + 2 ^ S := by omega
        _ = (u + 1) * 2 ^ S := by ring
        _ ≤ (u / s + 1) * s * 2 ^ S := Nat.mul_le_mul_right _ hq
        _ = (u / s + 1) * 2 ^ S * s := by ring
    exact Nat.lt_of_mul_lt_mul_right key

/-- The packed sum over a list. -/
theorem pack_list_sum {w d : ℕ} {α : Type*} (L : List α) (f : α → ℕ → ℕ) :
    pack w d (fun r ↦ (L.map fun a ↦ f a r).sum) = (L.map fun a ↦ pack w d (f a)).sum := by
  induction L with
  | nil => simp [pack]
  | cons a L ih => simp only [List.map_cons, List.sum_cons, pack_add, ih]

/-! ### One step of the packed recursion -/

/-- The new row `k`: `⌊(p y + c - a) / s⌋` lane by lane (by the magic division), plus the
identity term `tk`. Here `p y + c - a` holds `(D X)_k + s B` in every lane. -/
def rowVal (P : Params) (LM tk p c y a : ℕ) : ℕ :=
  ((((p * y + c - a) * P.M) >>> P.S) &&& LM) + tk

/-- The identity term `t e_k` of row `k` in the window of columns `col₀, …, col₀ + nc - 1`:
`t` in lane `k - col₀` if column `k` is in the window, else zero. -/
def idTerm (P : Params) (t col₀ nc k : ℕ) : ℕ :=
  if col₀ ≤ k ∧ k < col₀ + nc then t <<< (P.w * (k - col₀)) else 0

/-- The new rows, row `k₀` first; `none` as soon as one row fails the range test
`(v + hR) &&& HM = T`. -/
def finish (P : Params) (LM HM T hR t col₀ nc : ℕ) :
    ℕ → List ℕ → List ℕ → List ℕ → List ℕ → Option (List ℕ)
  | k, p :: ps, c :: cs, y :: ys, a :: as =>
    if (rowVal P LM (idTerm P t col₀ nc k) p c y a + hR) &&& HM = T then
      (finish P LM HM T hR t col₀ nc (k + 1) ps cs ys as).map
        (rowVal P LM (idTerm P t col₀ nc k) p c y a :: ·)
    else none
  | _, _, _, _, _ => some []

theorem finish_spec (P : Params) (LM HM T hR t col₀ nc : ℕ) :
    ∀ (k₀ : ℕ) (ps cs ys as out : List ℕ),
    finish P LM HM T hR t col₀ nc k₀ ps cs ys as = some out →
    out.length = min (min (min ps.length cs.length) ys.length) as.length ∧
      ∀ i < out.length, out.getD i 0 = rowVal P LM (idTerm P t col₀ nc (k₀ + i)) (ps.getD i 0)
          (cs.getD i 0) (ys.getD i 0) (as.getD i 0) ∧
        (rowVal P LM (idTerm P t col₀ nc (k₀ + i)) (ps.getD i 0) (cs.getD i 0) (ys.getD i 0)
          (as.getD i 0) + hR) &&& HM = T
  | k₀, p :: ps, c :: cs, y :: ys, a :: as, out, h => by
    simp only [finish] at h
    split_ifs at h with htest
    obtain ⟨out', h', rfl⟩ := Option.map_eq_some_iff.1 h
    obtain ⟨hl, hv⟩ := finish_spec P LM HM T hR t col₀ nc (k₀ + 1) ps cs ys as out' h'
    refine ⟨by simp [hl], fun i hi ↦ ?_⟩
    cases i with
    | zero => simpa using htest
    | succ i =>
      have := hv i (by simpa using hi)
      simpa [Nat.add_assoc, Nat.add_comm 1 i] using this
  | _, [], _, _, _, out, h => by simp only [finish, Option.some.injEq] at h; subst h; simp
  | _, _ :: _, [], _, _, out, h => by simp only [finish, Option.some.injEq] at h; subst h; simp
  | _, _ :: _, _ :: _, [], _, out, h => by
    simp only [finish, Option.some.injEq] at h; subst h; simp
  | _, _ :: _, _ :: _, _ :: _, [], out, h => by
    simp only [finish, Option.some.injEq] at h; subst h; simp

/-- The row conditions, checked once per matrix: list lengths, absolute row sums at most `2 s`,
and the constants `c_k = B (s - p_k + q_k + ∑_h m_h(k))` that make every lane of `p y + c - a`
equal to `(D X)_k + s B`. -/
def Band.RowsOk (D : Band) (P : Params) (cst : List ℕ) (d : ℕ) : Prop :=
  D.pos.length = d ∧ D.neg.length = d ∧ cst.length = d ∧ (∀ h ∈ D.hops, h.mask.length = d) ∧
    ∀ k < d, D.pos.getD k 0 + D.neg.getD k 0 + (D.hops.map fun h ↦ h.mask.getD k 0).sum ≤
        2 * P.s ∧
      cst.getD k 0 + P.B * D.pos.getD k 0 =
        P.B * (P.s + D.neg.getD k 0 + (D.hops.map fun h ↦ h.mask.getD k 0).sum)

/-- The invariant of the packed recursion on the window of columns `col₀, …, col₀ + nc - 1`:
row `k` packs the lanes `x (col₀ + r) k + B`, and every entry lies in `[-2^(b-1), 2^(b-1))`. -/
def Inv (P : Params) (d col₀ nc : ℕ) (ys : List ℕ) (x : ℕ → ℕ → ℤ) : Prop :=
  ys.length = d ∧ ∀ k < d, ys.getD k 0 = pack P.w nc (fun r ↦ (x (col₀ + r) k + P.B).toNat) ∧
    ∀ r < nc, -(P.half : ℤ) ≤ x (col₀ + r) k ∧ x (col₀ + r) k < P.half

/-- One step of the packed recursion, all columns of the window at once. -/
def step (P : Params) (D : Band) (cstP : List ℕ) (LM HM T hR t col₀ nc : ℕ) (ys : List ℕ) :
    Option (List ℕ) :=
  finish P LM HM T hR t col₀ nc 0 D.pos cstP ys (negAcc D ys)

section StepSpec

variable {P : Params} {D : Band} {d col₀ nc : ℕ} {cst : List ℕ} {x : ℕ → ℕ → ℤ}
  {ys : List ℕ}

theorem length_negAcc (hR : D.RowsOk P cst d) (hy : ys.length = d) :
    (negAcc D ys).length = d := by
  obtain ⟨-, hn, -, hm, -⟩ := hR
  have hml : (mulList D.neg ys).length = d := by rw [length_mulList, hn, hy, min_self]
  have hf := foldl_axpy_spec ys D.hops (mulList D.neg ys)
    (fun h hh ↦ by rw [hml, hm h hh]) (by rw [hml, hy])
  rw [hml] at hf
  exact hf.1

/-- The band products of row `k`, lane by lane. -/
theorem getD_negAcc (hD : D.InRange d) (hR : D.RowsOk P cst d) (hinv : Inv P d col₀ nc ys x)
    {k : ℕ} (hk : k < d) :
    (negAcc D ys).getD k 0 = pack P.w nc (fun r ↦
      D.neg.getD k 0 * (x (col₀ + r) k + P.B).toNat +
        (D.hops.map fun h ↦ h.mask.getD k 0 * (x (col₀ + r) (h.nb k) + P.B).toNat).sum) := by
  obtain ⟨hp, hn, -, hm, -⟩ := hR
  obtain ⟨hy, hrow⟩ := hinv
  have hml : (mulList D.neg ys).length = d := by rw [length_mulList, hn, hy, min_self]
  have hf := foldl_axpy_spec ys D.hops (mulList D.neg ys)
    (fun h hh ↦ by rw [hml, hm h hh]) (by rw [hml, hy])
  rw [hml] at hf
  unfold negAcc
  rw [hf.2 k hk, getD_mulList _ _ _ (by omega) (by omega), (hrow k hk).1, ← pack_mul,
    pack_add, pack_list_sum]
  congr 1
  refine congrArg List.sum (List.map_congr_left fun h hh ↦ ?_)
  rw [mask_mul_getD_shift ys hy (hD h hh k hk), pack_mul]
  by_cases hm0 : h.mask.getD k 0 = 0
  · rw [hm0, zero_mul, zero_mul]
  · rw [(hrow _ (hD.nb_lt hh hk hm0)).1]

/-- `∑ m (y + B) = ∑ m y + B ∑ m` over a list. -/
theorem list_sum_mul_add {α : Type*} (L : List α) (m : α → ℤ) (y : α → ℤ) (B : ℤ) :
    (L.map fun a ↦ m a * (y a + B)).sum = (L.map fun a ↦ m a * y a).sum + B * (L.map m).sum := by
  induction L with
  | nil => simp
  | cons a L ih => simp only [List.map_cons, List.sum_cons, ih]; ring

/-- **One lane of one row.** If every entry of the column `y` lies in `[-2^(b-1), 2^(b-1))`,
the lane of `p y + c` dominates that of the band products, and their difference is
`(D y)_k + s B ∈ [0, 2 s B]`. -/
theorem lane_identity (hP : P.Ok) (hD : D.InRange d) (hR : D.RowsOk P cst d) {y : ℕ → ℤ}
    (hbd : ∀ n < d, -(P.half : ℤ) ≤ y n ∧ y n < P.half) {k : ℕ} (hk : k < d) :
    let pos := D.pos.getD k 0 * (y k + P.B).toNat + cst.getD k 0
    let neg := D.neg.getD k 0 * (y k + P.B).toNat +
      (D.hops.map fun h ↦ h.mask.getD k 0 * (y (h.nb k) + P.B).toNat).sum
    neg ≤ pos ∧ ((pos - neg : ℕ) : ℤ) = D.apply y k + P.s * P.B ∧
      pos - neg ≤ 2 * P.s * P.B := by
  intro pos neg
  obtain ⟨hs, hb, -⟩ := hP
  obtain ⟨-, -, -, -, hrow⟩ := hR
  obtain ⟨hsum, hcst⟩ := hrow k hk
  have hB := P.B_eq hb
  have hcast : ∀ n < d, (((y n + P.B).toNat : ℕ) : ℤ) = y n + P.B := fun n hn ↦ by
    have := (hbd n hn).1
    rw [Int.toNat_of_nonneg]; push_cast [hB]; omega
  set p : ℕ := D.pos.getD k 0
  set q : ℕ := D.neg.getD k 0
  set sm : ℕ := (D.hops.map fun h ↦ h.mask.getD k 0).sum
  -- The band sum, cast to `ℤ`, with the biases separated.
  have hnb : (((D.hops.map fun h ↦ h.mask.getD k 0 * (y (h.nb k) + P.B).toNat).sum : ℕ) : ℤ)
      = (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * y (h.nb k)).sum + P.B * sm := by
    rw [Nat.cast_list_sum, List.map_map]
    have : (D.hops.map ((↑) ∘ fun h ↦ h.mask.getD k 0 * (y (h.nb k) + P.B).toNat)) =
        D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * (y (h.nb k) + P.B) := by
      refine List.map_congr_left fun h hh ↦ ?_
      simp only [Function.comp_apply, Nat.cast_mul]
      by_cases hm0 : h.mask.getD k 0 = 0
      · rw [hm0]; simp
      · rw [hcast _ (hD.nb_lt hh hk hm0)]
    rw [this, list_sum_mul_add]
    simp [sm, Nat.cast_list_sum, List.map_map, Function.comp_def]
  -- Bounds on the band sum.
  have hlo : -((sm : ℤ) * P.half) ≤ (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * y (h.nb k)).sum
      := by
    have : (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * (-(P.half : ℤ))).sum ≤
        (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * y (h.nb k)).sum := by
      refine List.sum_le_sum fun h hh ↦ ?_
      by_cases hm0 : h.mask.getD k 0 = 0
      · rw [hm0]; simp
      · exact mul_le_mul_of_nonneg_left (hbd _ (hD.nb_lt hh hk hm0)).1 (Nat.cast_nonneg _)
    refine le_trans (le_of_eq ?_) this
    rw [List.sum_map_mul_right]
    simp only [sm, Nat.cast_list_sum, List.map_map, Function.comp_def]
    ring
  have hhi : (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * y (h.nb k)).sum ≤ (sm : ℤ) * P.half
      := by
    have : (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * y (h.nb k)).sum ≤
        (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * P.half).sum := by
      refine List.sum_le_sum fun h hh ↦ ?_
      by_cases hm0 : h.mask.getD k 0 = 0
      · rw [hm0]; simp
      · exact mul_le_mul_of_nonneg_left (hbd _ (hD.nb_lt hh hk hm0)).2.le (Nat.cast_nonneg _)
    refine this.trans (le_of_eq ?_)
    rw [List.sum_map_mul_right]
    simp only [sm, Nat.cast_list_sum, List.map_map, Function.comp_def]
  have hxk := hbd k hk
  have hck := hcast k hk
  have hcst' : (cst.getD k 0 : ℤ) = P.B * (P.s + q + sm) - P.B * p := by
    have := congrArg (fun n : ℕ ↦ (n : ℤ)) hcst
    push_cast at this; linarith
  have hsum' : (p : ℤ) + q + sm ≤ 2 * P.s := by exact_mod_cast hsum
  have happly : D.apply y k = ((p : ℤ) - q) * y k -
      (D.hops.map fun h ↦ (h.mask.getD k 0 : ℤ) * y (h.nb k)).sum := rfl
  have hval : (pos : ℤ) - neg = D.apply y k + P.s * P.B := by
    simp only [pos, neg, Nat.cast_add, Nat.cast_mul, hck, hnb, happly, hcst']
    ring
  have hp0 : (0 : ℤ) ≤ p := Nat.cast_nonneg _
  have hq0 : (0 : ℤ) ≤ q := Nat.cast_nonneg _
  have hB' : (P.B : ℤ) = 2 * P.half := by exact_mod_cast hB
  have hs0 : (0 : ℤ) ≤ P.half := Nat.cast_nonneg _
  have a1 : (p : ℤ) * (-(P.half : ℤ)) ≤ p * y k := mul_le_mul_of_nonneg_left hxk.1 hp0
  have a2 : (q : ℤ) * y k ≤ q * P.half := mul_le_mul_of_nonneg_left hxk.2.le hq0
  have a3 : (p : ℤ) * y k ≤ p * P.half := mul_le_mul_of_nonneg_left hxk.2.le hp0
  have a4 : (q : ℤ) * (-(P.half : ℤ)) ≤ q * y k := mul_le_mul_of_nonneg_left hxk.1 hq0
  have h3 : ((p : ℤ) + q + sm) * P.half ≤ 2 * P.s * P.half :=
    mul_le_mul_of_nonneg_right hsum' hs0
  have hle : (neg : ℤ) ≤ pos := by
    have : (0 : ℤ) ≤ D.apply y k + P.s * P.B := by
      rw [happly, hB']; linarith
    linarith
  have hle' : neg ≤ pos := by exact_mod_cast hle
  refine ⟨hle', ?_, ?_⟩
  · rw [Nat.cast_sub hle', hval]
  · have : ((pos - neg : ℕ) : ℤ) ≤ 2 * P.s * P.B := by
      rw [Nat.cast_sub hle', hval, happly, hB']; linarith
    exact_mod_cast this

/-- The identity term of row `k` is the packed vector with `t` in the lane of column `k`. -/
theorem idTerm_eq (P : Params) (t col₀ nc k : ℕ) :
    idTerm P t col₀ nc k = pack P.w nc (fun r ↦ if col₀ + r = k then t else 0) := by
  unfold idTerm
  split_ifs with hk
  · rw [Nat.shiftLeft_eq, ← pack_single (show k - col₀ < nc by omega)]
    exact pack_congr fun r _ ↦ by
      congr 1; apply propext; omega
  · symm
    rw [pack_congr (g := fun r ↦ 0 * 0) fun r hr ↦ by
      rw [ite_eq_right (by omega)]]
    rw [pack_mul, zero_mul]

/-- **One step of the packed recursion.** If the rows of `ys` pack the columns `x r` of the
window and the step succeeds, the new rows pack the columns `specStep D s r t (x r)`, with all
entries in range. -/
theorem step_spec (hP : P.Ok) (hD : D.InRange d) (hR : D.RowsOk P cst d) {t : ℕ}
    (ht : t ≤ P.B) (hinv : Inv P d col₀ nc ys x) {ys' : List ℕ}
    (h : step P D (cst.map (· * repunit P.w nc)) ((2 ^ (P.w - P.S) - 1) * repunit P.w nc)
      ((2 ^ P.w - 2 ^ P.b) * repunit P.w nc) (2 ^ P.b * repunit P.w nc)
      (P.half * repunit P.w nc) t col₀ nc ys = some ys') :
    Inv P d col₀ nc ys' (fun r k ↦ specStep D P.s r t (x r) k) := by
  obtain ⟨hs, hb, hbw, hSw, hM1, hM2, hM3, hM4⟩ := hP
  have hP : P.Ok := ⟨hs, hb, hbw, hSw, hM1, hM2, hM3, hM4⟩
  have hw : 0 < P.w := by omega
  have hB := P.B_eq hb
  set rep := repunit P.w nc
  have hpl := hR.1
  have hcl := hR.2.2.1
  obtain ⟨hlen, hrows⟩ := finish_spec P _ _ _ _ t col₀ nc 0 D.pos (cst.map (· * rep)) ys
    (negAcc D ys) ys' h
  rw [hpl, List.length_map, hcl, hinv.1, length_negAcc hR hinv.1] at hlen
  simp only [min_self] at hlen
  refine ⟨hlen, fun k hk ↦ ?_⟩
  obtain ⟨hv, htest⟩ := hrows k (hlen ▸ hk)
  rw [zero_add] at hv htest
  have hc : (cst.map (· * rep)).getD k 0 = cst.getD k 0 * rep := by
    simp only [List.getD_eq_getElem?_getD, List.getElem?_map]
    cases cst[k]? <;> simp
  rw [hc] at hv htest
  -- The lanes of `p y + c - a`.
  set ℓ : ℕ → ℕ → ℕ := fun r n ↦ (x (col₀ + r) n + P.B).toNat
  set pos : ℕ → ℕ := fun r ↦ D.pos.getD k 0 * ℓ r k + cst.getD k 0
  set neg : ℕ → ℕ := fun r ↦ D.neg.getD k 0 * ℓ r k +
    (D.hops.map fun h ↦ h.mask.getD k 0 * ℓ r (h.nb k)).sum
  have hid : ∀ r < nc, neg r ≤ pos r ∧ ((pos r - neg r : ℕ) : ℤ) =
      D.apply (x (col₀ + r)) k + P.s * P.B ∧ pos r - neg r ≤ 2 * P.s * P.B := fun r hr ↦
    lane_identity hP hD hR (fun n hn ↦ (hinv.2 n hn).2 r hr) hk
  have hZ : D.pos.getD k 0 * ys.getD k 0 + cst.getD k 0 * rep - (negAcc D ys).getD k 0 =
      pack P.w nc (fun r ↦ pos r - neg r) := by
    rw [(hinv.2 k hk).1, getD_negAcc hD hR hinv hk, ← pack_mul, ← pack_const hw,
      ← pack_add, pack_sub fun r hr ↦ (hid r hr).1]
  -- The magic division, lane by lane.
  have hcanZ : Canonical P.w nc (fun r ↦ (pos r - neg r) * P.M) := fun r hr ↦
    lt_of_le_of_lt (Nat.mul_le_mul_right _ (hid r hr).2.2) hM3
  have hdiv : ∀ r < nc, ((pos r - neg r) * P.M) >>> P.S = (pos r - neg r) / P.s :=
    fun r hr ↦ by
      rw [Nat.shiftRight_eq_div_pow]
      exact magic_div hs hM1 hM2 (hid r hr).2.2
  set L : ℕ → ℕ := fun r ↦ (pos r - neg r) / P.s + if col₀ + r = k then t else 0
  have hrow : rowVal P ((2 ^ (P.w - P.S) - 1) * rep) (idTerm P t col₀ nc k) (D.pos.getD k 0)
      (cst.getD k 0 * rep) (ys.getD k 0) ((negAcc D ys).getD k 0) = pack P.w nc L := by
    rw [rowVal, hZ, mul_comm (pack _ _ _), ← pack_mul, ← pack_const hw,
      pack_congr (g := fun r ↦ (pos r - neg r) * P.M) (fun r _ ↦ mul_comm _ _),
      pack_shiftRight_and hw hcanZ hSw, pack_congr hdiv, idTerm_eq, ← pack_add]
  -- The lanes are the new entries plus the bias.
  have hL : ∀ r < nc, (L r : ℤ) = specStep D P.s (col₀ + r) t (x (col₀ + r)) k + P.B := by
    intro r hr
    have hsB : ((pos r - neg r) / P.s : ℕ) =
        (D.apply (x (col₀ + r)) k + P.s * P.B) / (P.s : ℤ) := by
      rw [Int.natCast_div, (hid r hr).2.1]
    simp only [L, specStep, Nat.cast_add, hsB]
    rw [show (D.apply (x (col₀ + r)) k + P.s * P.B) = D.apply (x (col₀ + r)) k + P.B * P.s by
      ring, Int.add_mul_ediv_right _ _ (by exact_mod_cast hs.ne')]
    have e1 : (if col₀ + r = k then (t : ℤ) else 0) = (if k = col₀ + r then (t : ℤ) else 0) := by
      simp only [eq_comm]
    push_cast
    rw [e1]
    ring
  have hLnat : ∀ r < nc, L r = (specStep D P.s (col₀ + r) t (x (col₀ + r)) k + P.B).toNat :=
    fun r hr ↦ by rw [← hL r hr, Int.toNat_natCast]
  refine ⟨by rw [hv, hrow]; exact pack_congr hLnat, fun r hr ↦ ?_⟩
  -- The range test.
  have hcan : Canonical P.w nc (fun r ↦ L r + P.half) := fun r' hr' ↦ by
    have h1 : (pos r' - neg r') / P.s ≤ 2 * P.B := by
      rw [Nat.div_le_iff_le_mul_add_pred hs]
      have := (hid r' hr').2.2
      nlinarith
    have : L r' ≤ 3 * P.B := by simp only [L]; split_ifs <;> omega
    change L r' + P.half < 2 ^ P.w
    omega
  rw [hrow, ← pack_const hw, ← pack_add, ← pack_const hw, ← pack_const hw] at htest
  have hmem := lanes_mem_of_and_eq hw hbw hcan htest r hr
  have hL' := hL r hr
  have h2b : ((2 ^ P.b : ℕ) : ℤ) = P.B := rfl
  have h2b1 : ((2 ^ (P.b + 1) : ℕ) : ℤ) = 2 * P.B := by
    rw [pow_succ]; push_cast; rw [mul_comm]; rfl
  have m1 : ((2 ^ P.b : ℕ) : ℤ) ≤ L r + P.half := by exact_mod_cast hmem.1
  have m2 : (L r : ℤ) + P.half < ((2 ^ (P.b + 1) : ℕ) : ℤ) := by exact_mod_cast hmem.2
  have hB' : (P.B : ℤ) = 2 * P.half := by exact_mod_cast hB
  rw [h2b] at m1; rw [h2b1] at m2
  constructor <;> linarith

end StepSpec

/-! ### The packed recursion -/

/-- **The packed Horner recursion.** The columns `col₀, …, col₀ + nc - 1` of
`X ← ⌊D X / s⌋ + t I` (coefficients `cs`, last one first, from `X = 0`), row `k` packed with
lanes `X_{k, col₀ + r} + B`; `none` if a range test fails. -/
def run (P : Params) (D : Band) (cst : List ℕ) (d col₀ nc : ℕ) (cs : List ℕ) :
    Option (List ℕ) :=
  cs.foldr (fun t oy ↦ oy.bind (step P D (cst.map (· * repunit P.w nc))
    ((2 ^ (P.w - P.S) - 1) * repunit P.w nc) ((2 ^ P.w - 2 ^ P.b) * repunit P.w nc)
    (2 ^ P.b * repunit P.w nc) (P.half * repunit P.w nc) t col₀ nc))
    (some (List.replicate d (P.B * repunit P.w nc)))

/-- **Correctness of the packed recursion.** If it succeeds, row `k` packs the lanes
`specCol D s (col₀ + r) cs k + B` of the columns of the window. -/
theorem run_spec {P : Params} {D : Band} {d col₀ nc : ℕ} {cst : List ℕ} (hP : P.Ok)
    (hD : D.InRange d) (hR : D.RowsOk P cst d) : ∀ (cs : List ℕ), (∀ t ∈ cs, t ≤ P.B) →
    ∀ {ys : List ℕ}, run P D cst d col₀ nc cs = some ys →
      Inv P d col₀ nc ys (fun r k ↦ specCol D P.s r (cs.map (↑)) k)
  | [], _, ys, h => by
    simp only [run, List.foldr_nil, Option.some.injEq] at h
    subst h
    have hw : 0 < P.w := by have := hP.2.2.1; omega
    refine ⟨by simp, fun k hk ↦ ⟨?_, fun r _ ↦ ?_⟩⟩
    · rw [List.getD_eq_getElem?_getD, List.getElem?_replicate_of_lt hk, Option.getD_some,
        ← pack_const hw]
      exact pack_congr fun r _ ↦ by simp [specCol]
    · have : (0 : ℤ) < P.half := by simp [Params.half]
      simp only [specCol, List.map_nil, List.foldr_nil, Pi.zero_apply]
      constructor <;> omega
  | t :: cs, hcs, ys, h => by
    simp only [run, List.foldr_cons] at h
    obtain ⟨ys₀, h₀, h₁⟩ := Option.bind_eq_some_iff.1 h
    have ih := run_spec hP hD hR cs (fun t' ht' ↦ hcs t' (by simp [ht'])) h₀
    exact step_spec hP hD hR (hcs t (by simp)) ih h₁

/-- Lane `r` of row `k` of a successful packed run, as an integer: entry `(k, col₀ + r)` of the
recursion. -/
theorem run_lane {P : Params} {D : Band} {d col₀ nc : ℕ} {cst : List ℕ} (hP : P.Ok)
    (hD : D.InRange d) (hR : D.RowsOk P cst d) {cs : List ℕ} (hcs : ∀ t ∈ cs, t ≤ P.B)
    {ys : List ℕ} (h : run P D cst d col₀ nc cs = some ys) {k r : ℕ} (hk : k < d)
    (hr : r < nc) :
    (((ys.getD k 0 >>> (P.w * r)) % 2 ^ P.w : ℕ) : ℤ) - P.B =
      specCol D P.s (col₀ + r) (cs.map (↑)) k := by
  have hw : 0 < P.w := by have := hP.2.2.1; omega
  have hB := P.B_eq hP.2.1
  have hM4 := hP.2.2.2.2.2.2.2
  have hinv := run_spec hP hD hR cs hcs h
  obtain ⟨hrow, hbd⟩ := hinv.2 k hk
  set c : ℕ → ℤ := fun r' ↦ specCol D P.s (col₀ + r') (cs.map (↑)) k with hc
  have hbd' : ∀ r' < nc, -(P.half : ℤ) ≤ c r' ∧ c r' < P.half := hbd
  have hrow' : ys.getD k 0 = pack P.w nc (fun r' ↦ (c r' + P.B).toNat) := hrow
  have hB' : (P.B : ℤ) = 2 * P.half := by exact_mod_cast hB
  have hnn : ∀ r' < nc, 0 ≤ c r' + P.B := fun r' hr' ↦ by
    have := (hbd' r' hr').1; omega
  have hcan : Canonical P.w nc (fun r' ↦ (c r' + P.B).toNat) := fun r' hr' ↦ by
    have h1 := (hbd' r' hr').2
    have h3 : (((c r' + P.B).toNat : ℕ) : ℤ) = c r' + P.B := Int.toNat_of_nonneg (hnn r' hr')
    have h2 : (c r' + P.B).toNat < 4 * P.B := by omega
    change (c r' + P.B).toNat < 2 ^ P.w
    omega
  rw [hrow', pack_lane hw hcan hr, Int.toNat_of_nonneg (hnn r hr)]
  ring

/-! ### Weighted sums of squares of the computed columns -/

/-- `∑_r v_r (lane_r(y) - B)²` over the lanes of one packed row, peeling one lane at a time. -/
def laneSq (P : Params) : ℕ → List ℕ → ℕ
  | y, vr :: vs => vr * ((((y % 2 ^ P.w : ℕ) : ℤ) - P.B).natAbs ^ 2) + laneSq P (y >>> P.w) vs
  | _, [] => 0

/-- `∑_k u_k ∑_r v_r (X_{k r})²` over the rows of a packed run. -/
def numerator (P : Params) (u v ys : List ℕ) : ℕ :=
  (List.zipWith (fun uk y ↦ uk * laneSq P y v) u ys).sum

theorem laneSq_eq (P : Params) : ∀ (v : List ℕ) (y : ℕ),
    laneSq P y v = ∑ r ∈ range v.length,
      v.getD r 0 * ((((y >>> (P.w * r)) % 2 ^ P.w : ℕ) : ℤ) - P.B).natAbs ^ 2
  | [], y => by simp [laneSq]
  | vr :: vs, y => by
    rw [laneSq, laneSq_eq P vs, List.length_cons, sum_range_succ']
    simp only [List.getD_cons_succ, List.getD_cons_zero, mul_zero, Nat.shiftRight_zero]
    rw [add_comm]
    congr 1
    refine sum_congr rfl fun r _ ↦ ?_
    rw [← Nat.shiftRight_add, mul_add, mul_one, Nat.add_comm (P.w * r)]

theorem sum_zipWith_eq {f : ℕ → ℕ → ℕ} : ∀ (u ys : List ℕ) (d : ℕ), u.length = d →
    ys.length = d → (List.zipWith f u ys).sum = ∑ k ∈ range d, f (u.getD k 0) (ys.getD k 0)
  | [], [], 0, _, _ => by simp
  | a :: u, y :: ys, d + 1, hu, hy => by
    simp only [List.length_cons, Nat.add_right_cancel_iff] at hu hy
    rw [List.zipWith_cons_cons, List.sum_cons, sum_zipWith_eq u ys d hu hy, sum_range_succ',
      add_comm]
    simp
  | [], _ :: _, _, hu, hy => by simp at hu hy; omega
  | _ :: _, [], _, hu, hy => by simp at hu hy; omega

/-- **The weighted sum of squares of the columns of the window**, from a successful packed
run: `∑_k u_k ∑_{r < nc} v_r X_{k, col₀ + r}²`. -/
theorem run_numerator {P : Params} {D : Band} {d col₀ nc : ℕ} {cst : List ℕ} (hP : P.Ok)
    (hD : D.InRange d) (hR : D.RowsOk P cst d) {cs : List ℕ} (hcs : ∀ t ∈ cs, t ≤ P.B)
    {ys : List ℕ} (h : run P D cst d col₀ nc cs = some ys) {u v : List ℕ} (hu : u.length = d)
    (hv : v.length = nc) :
    (numerator P u v ys : ℤ) = ∑ k ∈ range d, ∑ r ∈ range nc,
      (u.getD k 0 : ℤ) * v.getD r 0 * specCol D P.s (col₀ + r) (cs.map (↑)) k ^ 2 := by
  have hys := (run_spec hP hD hR cs hcs h).1
  rw [numerator, sum_zipWith_eq u ys d hu hys]
  push_cast
  refine sum_congr rfl fun k hk ↦ ?_
  rw [laneSq_eq, hv, Nat.cast_sum, Finset.mul_sum]
  refine sum_congr rfl fun r hr ↦ ?_
  rw [run_lane hP hD hR hcs h (mem_range.1 hk) (mem_range.1 hr)]
  push_cast
  rw [sq_abs]
  ring

/-! ### A linear-time check of the band conditions -/

/-- Entrywise sum of two lists (truncated to the shorter one). -/
def addList : List ℕ → List ℕ → List ℕ
  | a :: as, b :: bs => (a + b) :: addList as bs
  | _, _ => []

/-- The row sums `∑_h m_h(k)` of the band entries. -/
def maskSums (d : ℕ) (hops : List Hop) : List ℕ :=
  hops.foldl (fun acc h ↦ addList acc h.mask) (List.replicate d 0)

/-- The row conditions of `Band.RowsOk`, row by row along the lists. -/
def rowCheck (P : Params) : List ℕ → List ℕ → List ℕ → List ℕ → Bool
  | p :: ps, q :: qs, c :: cs, m :: ms =>
    (decide (p + q + m ≤ 2 * P.s) && c + P.B * p == P.B * (P.s + q + m)) &&
      rowCheck P ps qs cs ms
  | [], [], [], [] => true
  | _, _, _, _ => false

/-- The range condition of one band, along its list of entries from row `k`. -/
def hopCheck (up : Bool) (off d : ℕ) : ℕ → List ℕ → Bool
  | k, m :: ms => (m == 0 || if up then decide (k + off < d) else decide (off ≤ k)) &&
      hopCheck up off d (k + 1) ms
  | _, [] => true

/-- **The band conditions as one Boolean**, evaluated in time linear in the number of entries:
list lengths, `Band.InRange` and `Band.RowsOk`. -/
def Band.check (D : Band) (P : Params) (cst : List ℕ) (d : ℕ) : Bool :=
  D.pos.length == d && D.neg.length == d && cst.length == d &&
    D.hops.all (fun h ↦ h.mask.length == d && hopCheck h.up h.off d 0 h.mask) &&
    rowCheck P D.pos D.neg cst (maskSums d D.hops)

theorem length_addList : ∀ a b : List ℕ, (addList a b).length = min a.length b.length
  | _ :: _, _ :: _ => by simp [addList, length_addList]
  | [], _ => by simp [addList]
  | _ :: _, [] => by simp [addList]

theorem getD_addList : ∀ (a b : List ℕ) (k : ℕ), k < a.length → k < b.length →
    (addList a b).getD k 0 = a.getD k 0 + b.getD k 0
  | _ :: _, _ :: _, 0, _, _ => by simp [addList]
  | _ :: as, _ :: bs, k + 1, h1, h2 => by
    simp only [addList, List.getD_cons_succ]
    exact getD_addList as bs k (by simpa using h1) (by simpa using h2)

theorem maskSums_spec (d : ℕ) : ∀ (L : List Hop) (acc : List ℕ), acc.length = d →
    (∀ h ∈ L, h.mask.length = d) →
    (L.foldl (fun acc h ↦ addList acc h.mask) acc).length = d ∧
      ∀ k < d, (L.foldl (fun acc h ↦ addList acc h.mask) acc).getD k 0 =
        acc.getD k 0 + (L.map fun h ↦ h.mask.getD k 0).sum
  | [], acc, hacc, _ => by simp [hacc]
  | h :: L, acc, hacc, hm => by
    have hl : (addList acc h.mask).length = d := by
      rw [length_addList, hacc, hm h (by simp), min_self]
    have ih := maskSums_spec d L (addList acc h.mask) hl (fun h' hh' ↦ hm h' (by simp [hh']))
    refine ⟨ih.1, fun k hk ↦ ?_⟩
    rw [List.foldl_cons, ih.2 k hk, getD_addList _ _ _ (by omega) (by rw [hm h (by simp)]; omega),
      List.map_cons, List.sum_cons]
    ring

theorem rowCheck_spec (P : Params) : ∀ (ps qs cs ms : List ℕ), rowCheck P ps qs cs ms = true →
    qs.length = ps.length ∧ cs.length = ps.length ∧ ms.length = ps.length ∧
      ∀ k < ps.length, ps.getD k 0 + qs.getD k 0 + ms.getD k 0 ≤ 2 * P.s ∧
        cs.getD k 0 + P.B * ps.getD k 0 = P.B * (P.s + qs.getD k 0 + ms.getD k 0)
  | p :: ps, q :: qs, c :: cs, m :: ms, h => by
    simp only [rowCheck, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    obtain ⟨l1, l2, l3, hk⟩ := rowCheck_spec P ps qs cs ms h3
    refine ⟨by simp [l1], by simp [l2], by simp [l3], fun k hk' ↦ ?_⟩
    cases k with
    | zero => exact ⟨h1, h2⟩
    | succ k => simpa using hk k (by simpa using hk')
  | [], [], [], [], _ => by simp
  | [], [], [], _ :: _, h => by simp [rowCheck] at h
  | [], [], _ :: _, _, h => by simp [rowCheck] at h
  | [], _ :: _, _, _, h => by simp [rowCheck] at h
  | _ :: _, [], _, _, h => by simp [rowCheck] at h
  | _ :: _, _ :: _, [], _, h => by simp [rowCheck] at h
  | _ :: _, _ :: _, _ :: _, [], h => by simp [rowCheck] at h

theorem hopCheck_spec (up : Bool) (off d : ℕ) : ∀ (k₀ : ℕ) (ms : List ℕ),
    hopCheck up off d k₀ ms = true → ∀ i < ms.length, ms.getD i 0 ≠ 0 →
      (up = true → k₀ + i + off < d) ∧ (up = false → off ≤ k₀ + i)
  | k₀, m :: ms, h, i, hi, hm => by
    simp only [hopCheck, Bool.and_eq_true, Bool.or_eq_true, beq_iff_eq] at h
    obtain ⟨h1, h2⟩ := h
    cases i with
    | zero =>
      simp only [List.getD_cons_zero, ne_eq] at hm
      rcases h1 with h1 | h1
      · exact (hm h1).elim
      · cases up <;> simp_all
    | succ i =>
      have := hopCheck_spec up off d (k₀ + 1) ms h2 i (by simpa using hi) (by simpa using hm)
      simpa [Nat.add_assoc, Nat.add_comm 1 i] using this
  | _, [], _, i, hi, _ => by simp at hi

/-- The constants `c_k = B (s + q_k + m_k - p_k)` of `Band.RowsOk`, along the lists. -/
def cstOf (P : Params) : List ℕ → List ℕ → List ℕ → List ℕ
  | p :: ps, q :: qs, m :: ms => P.B * (P.s + q + m - p) :: cstOf P ps qs ms
  | _, _, _ => []

/-- The constants of `Band.RowsOk` for a band. -/
def Band.cst (D : Band) (P : Params) (d : ℕ) : List ℕ := cstOf P D.pos D.neg (maskSums d D.hops)

/-- **The Boolean check gives the band conditions.** -/
theorem Band.check_spec {D : Band} {P : Params} {cst : List ℕ} {d : ℕ}
    (h : D.check P cst d = true) : D.InRange d ∧ D.RowsOk P cst d := by
  simp only [Band.check, Bool.and_eq_true, beq_iff_eq, List.all_eq_true] at h
  obtain ⟨⟨⟨⟨hp, hn⟩, hc⟩, hh⟩, hr⟩ := h
  have hm : ∀ h ∈ D.hops, h.mask.length = d := fun h hh' ↦ (hh h hh').1
  obtain ⟨hml, hms⟩ := maskSums_spec d D.hops (List.replicate d 0) (by simp) hm
  obtain ⟨-, -, -, hrow⟩ := rowCheck_spec P _ _ _ _ hr
  rw [hp] at hrow
  refine ⟨fun h hh' k hk hmk ↦ ?_, hp, hn, hc, hm, fun k hk ↦ ?_⟩
  · have := hopCheck_spec h.up h.off d 0 h.mask (hh h hh').2 k (by rw [hm h hh']; exact hk) hmk
    simpa using this
  · have := hrow k hk
    have hk' : (maskSums d D.hops).getD k 0 = (D.hops.map fun h ↦ h.mask.getD k 0).sum := by
      rw [maskSums, hms k hk]; simp [List.getD_eq_getElem?_getD, hk]
    rwa [hk'] at this

end CertifiedGaps.Trace.Packed
