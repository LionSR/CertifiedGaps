/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Trace.PackedHorner

/-!
# Block bands: four banded matrices as one

A matrix over the Eisenstein integers `ℤ[ω]` acts on `u + vω` as a `2 × 2` block matrix of
integer matrices on `(u, v)`. The packed recursion of `Trace/PackedHorner.lean` runs on one
banded integer matrix, so this file places four banded matrices `B_uu, B_uv, B_vu, B_vv` on `d`
rows as the blocks of one banded matrix `realBand d B_uu B_uv B_vu B_vv` on `2d` rows: row `k` of
block `a` is row `blk d a k` (`k` for the first block, `d + k` for the second).

## Main definitions

* `blk d a k`: the row `k` of block `a`.
* `Hop.lift d a b h`: the band `h` moved to block `(a, b)`. A band within a block keeps its
  offset; a band from the first block to the second gains `d`, and conversely.
* `realBand d B_uu B_uv B_vu B_vv`: the block band.

## Main results

* `lift_term`: the entry of a moved band, given that the band stays in range in its block.
* `realBand_entry`: the entry of the block band at rows `blk d a k`, `blk d b l` is the entry of
  block `(a, b)` at `k, l`, provided the four blocks stay in range and the off-diagonal blocks
  have zero diagonal.
-/

open Finset

namespace CertifiedGaps.Trace.Packed

/-- Row `k` of block `a` of a matrix on `2d` rows: `k` for the first block, `d + k` for the
second. -/
def blk (d : ℕ) (a : Bool) (k : ℕ) : ℕ := if a then d + k else k

theorem blk_lt {d k : ℕ} (a : Bool) (hk : k < d) : blk d a k < d + d := by
  unfold blk; split_ifs <;> omega

/-- The entries `m_0, …, m_{d-1}` of a mask placed in block `a` of `2d` rows. -/
def liftMask (d : ℕ) (a : Bool) (m : List ℕ) : List ℕ :=
  if a then List.replicate d 0 ++ (List.range d).map (m.getD · 0)
  else (List.range d).map (m.getD · 0) ++ List.replicate d 0

/-- **The band `h` of a `d`-row matrix moved to block `(a, b)`** of a `2d`-row matrix: rows of
block `a` couple to columns of block `b`. -/
def Hop.lift (d : ℕ) (a b : Bool) (h : Hop) : Hop where
  up := if a = b then h.up else !a
  off := if a = b then h.off else if h.up = b then h.off + d else d - h.off
  mask := liftMask d a h.mask

theorem getD_append_left' {l m : List ℕ} {k : ℕ} (h : k < l.length) :
    (l ++ m).getD k 0 = l.getD k 0 := by
  simp only [List.getD_eq_getElem?_getD, List.getElem?_append_left h]

theorem getD_append_right' {l m : List ℕ} {k : ℕ} (h : l.length ≤ k) :
    (l ++ m).getD k 0 = m.getD (k - l.length) 0 := by
  simp only [List.getD_eq_getElem?_getD, List.getElem?_append_right h]

theorem getD_replicate_zero (d k : ℕ) : (List.replicate d 0).getD k 0 = 0 := by
  simp only [List.getD_eq_getElem?_getD, List.getElem?_replicate]
  split_ifs <;> rfl

theorem getD_map_range_getD (d j : ℕ) (m : List ℕ) :
    ((List.range d).map (m.getD · 0)).getD j 0 = if j < d then m.getD j 0 else 0 := by
  split_ifs with hj
  · simp [List.getD_eq_getElem?_getD, hj]
  · simp [List.getD_eq_getElem?_getD, hj]

theorem liftMask_getD {d k : ℕ} (hk : k < d) (a a' : Bool) (m : List ℕ) :
    (liftMask d a m).getD (blk d a' k) 0 = if a' = a then m.getD k 0 else 0 := by
  cases a <;> cases a' <;> simp only [liftMask, blk, Bool.false_eq_true, ite_true, ite_false]
  · rw [getD_append_left' (by simpa using hk), getD_map_range_getD]; simp [hk]
  · rw [getD_append_right' (by simp), getD_replicate_zero]; rfl
  · rw [getD_append_left' (by simpa using hk), getD_replicate_zero]
  · rw [getD_append_right' (by simp), List.length_replicate, Nat.add_sub_cancel_left,
      getD_map_range_getD]
    simp [hk]

/-- **The entry of a moved band**, if the band stays in range in its block. -/
theorem lift_term {d : ℕ} (h : Hop) (a b a' b' : Bool) {k l : ℕ} (hk : k < d) (hl : l < d)
    (hin : h.mask.getD k 0 ≠ 0 → (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) :
    (if blk d b' l = (h.lift d a b).nb (blk d a' k) then
        ((h.lift d a b).mask.getD (blk d a' k) 0 : ℤ) else 0) =
      if a' = a ∧ b' = b ∧ l = h.nb k then (h.mask.getD k 0 : ℤ) else 0 := by
  have hm : (h.lift d a b).mask.getD (blk d a' k) 0 = if a' = a then h.mask.getD k 0 else 0 :=
    liftMask_getD hk a a' h.mask
  rw [hm]
  by_cases ha : a' = a
  · subst ha
    by_cases h0 : h.mask.getD k 0 = 0
    · simp only [h0, Nat.cast_zero, ite_self]
    obtain ⟨hup, hdn⟩ := hin h0
    have key : blk d b' l = (h.lift d a' b).nb (blk d a' k) ↔ b' = b ∧ l = h.nb k := by
      cases a' <;> cases b <;> cases b' <;> cases hu : h.up <;>
        simp_all [blk, Hop.lift, Hop.nb] <;> omega
    simp only [true_and, ite_true, key]
  · simp [ha]

/-- The terms of a list of moved bands. -/
theorem sum_lift_term {d : ℕ} (L : List Hop) (a b a' b' : Bool) {k l : ℕ} (hk : k < d)
    (hl : l < d) (hin : ∀ h ∈ L, h.mask.getD k 0 ≠ 0 →
      (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) :
    ((L.map (Hop.lift d a b)).map fun h ↦ if blk d b' l = h.nb (blk d a' k) then
        (h.mask.getD (blk d a' k) 0 : ℤ) else 0).sum =
      if a' = a ∧ b' = b then
        (L.map fun h ↦ if l = h.nb k then (h.mask.getD k 0 : ℤ) else 0).sum else 0 := by
  induction L with
  | nil => simp
  | cons h L ih =>
    simp only [List.map_cons, List.sum_cons]
    rw [lift_term h a b a' b' hk hl (hin h (List.mem_cons_self ..)),
      ih fun h' hh ↦ hin h' (List.mem_cons_of_mem _ hh)]
    by_cases hab : a' = a ∧ b' = b
    · simp only [hab, true_and, ite_true]
    · simp only [hab, ite_false, add_zero]
      rw [ite_eq_right]
      intro hc
      exact (hab ⟨hc.1, hc.2.1⟩).elim

/-- The block `(a, b)` among four bands. -/
def pickBlock (Buu Buv Bvu Bvv : Band) (a b : Bool) : Band :=
  if a then (if b then Bvv else Bvu) else (if b then Buv else Buu)

/-- **The block band** of four `d`-row bands: `B_uu` on the first `d` rows and columns, `B_vv`
on the last, `B_uv` and `B_vu` coupling them. The diagonal is taken from `B_uu` and `B_vv`. -/
def realBand (d : ℕ) (Buu Buv Bvu Bvv : Band) : Band where
  pos := (List.range d).map (Buu.pos.getD · 0) ++ (List.range d).map (Bvv.pos.getD · 0)
  neg := (List.range d).map (Buu.neg.getD · 0) ++ (List.range d).map (Bvv.neg.getD · 0)
  hops := Buu.hops.map (Hop.lift d false false) ++ Buv.hops.map (Hop.lift d false true) ++
    Bvu.hops.map (Hop.lift d true false) ++ Bvv.hops.map (Hop.lift d true true)
  phops := Buu.phops.map (Hop.lift d false false) ++ Buv.phops.map (Hop.lift d false true) ++
    Bvu.phops.map (Hop.lift d true false) ++ Bvv.phops.map (Hop.lift d true true)

/-- **The entries of the block band**: at rows `blk d a k`, `blk d b l` the entry of block
`(a, b)` at `k, l`, provided the blocks stay in range and `B_uv`, `B_vu` have zero diagonal. -/
theorem realBand_entry {d : ℕ} {Buu Buv Bvu Bvv : Band} (huu : Buu.InRange d)
    (huv : Buv.InRange d) (hvu : Bvu.InRange d) (hvv : Bvv.InRange d)
    (huv0 : ∀ k, Buv.pos.getD k 0 = Buv.neg.getD k 0)
    (hvu0 : ∀ k, Bvu.pos.getD k 0 = Bvu.neg.getD k 0) (a b : Bool) (k l : Fin d) :
    (realBand d Buu Buv Bvu Bvv).matrix (d + d) ⟨blk d a k, blk_lt a k.isLt⟩
        ⟨blk d b l, blk_lt b l.isLt⟩ = (pickBlock Buu Buv Bvu Bvv a b).matrix d k l := by
  have hin : ∀ (B : Band), B.InRange d → (∀ h ∈ B.hops, h.mask.getD k 0 ≠ 0 →
      (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) ∧
      (∀ h ∈ B.phops, h.mask.getD k 0 ≠ 0 →
        (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) := fun B hB ↦
    ⟨fun h hh ↦ hB h (List.mem_append_left _ hh) k k.isLt,
      fun h hh ↦ hB h (List.mem_append_right _ hh) k k.isLt⟩
  have hR : ∀ (L₁ L₂ L₃ L₄ : List Hop), (∀ h ∈ L₁, h.mask.getD k 0 ≠ 0 →
      (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) →
      (∀ h ∈ L₂, h.mask.getD k 0 ≠ 0 →
      (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) →
      (∀ h ∈ L₃, h.mask.getD k 0 ≠ 0 →
      (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) →
      (∀ h ∈ L₄, h.mask.getD k 0 ≠ 0 →
      (h.up = true → k + h.off < d) ∧ (h.up = false → h.off ≤ k)) →
      ((L₁.map (Hop.lift d false false) ++ L₂.map (Hop.lift d false true) ++
        L₃.map (Hop.lift d true false) ++ L₄.map (Hop.lift d true true)).map fun h ↦
          if ((⟨blk d b l, blk_lt b l.isLt⟩ : Fin (d + d)) : ℕ) =
            h.nb (⟨blk d a k, blk_lt a k.isLt⟩ : Fin (d + d)) then
          (h.mask.getD (⟨blk d a k, blk_lt a k.isLt⟩ : Fin (d + d)) 0 : ℤ) else 0).sum =
        ((if a then (if b then L₄ else L₃) else (if b then L₂ else L₁)).map fun h ↦
          if (l : ℕ) = h.nb k then (h.mask.getD k 0 : ℤ) else 0).sum := by
    intro L₁ L₂ L₃ L₄ h₁ h₂ h₃ h₄
    simp only [List.map_append, List.sum_append]
    rw [sum_lift_term L₁ _ _ _ _ k.isLt l.isLt h₁, sum_lift_term L₂ _ _ _ _ k.isLt l.isLt h₂,
      sum_lift_term L₃ _ _ _ _ k.isLt l.isLt h₃, sum_lift_term L₄ _ _ _ _ k.isLt l.isLt h₄]
    cases a <;> cases b <;> simp
  have hpos : ∀ (p q : List ℕ), ((List.range d).map (p.getD · 0) ++
      (List.range d).map (q.getD · 0)).getD (blk d a k) 0 = if a then q.getD k 0 else p.getD k 0 :=
    fun p q ↦ by
      have h1 : ∀ m : List ℕ, ((List.range d).map (m.getD · 0)).getD k 0 = m.getD k 0 :=
        fun m ↦ by simp [List.getD_eq_getElem?_getD, k.isLt]
      cases a
      · simp only [blk, Bool.false_eq_true, ite_false]
        rw [getD_append_left' (by simp), h1]
      · simp only [blk, ite_true]
        rw [getD_append_right' (by simp), List.length_map, List.length_range,
          Nat.add_sub_cancel_left, h1]
  have heq : ((⟨blk d a k, blk_lt a k.isLt⟩ : Fin (d + d)) = ⟨blk d b l, blk_lt b l.isLt⟩) ↔
      a = b ∧ k = l := by
    have := k.isLt; have := l.isLt
    rw [Fin.mk.injEq, Fin.ext_iff]
    cases a <;> cases b <;> simp [blk] <;> omega
  simp only [Band.matrix, realBand]
  rw [hR _ _ _ _ (hin _ huu).1 (hin _ huv).1 (hin _ hvu).1 (hin _ hvv).1,
    hR _ _ _ _ (hin _ huu).2 (hin _ huv).2 (hin _ hvu).2 (hin _ hvv).2, hpos, hpos]
  simp only [heq]
  have huv0' : Buv.pos[(k : ℕ)]?.getD 0 = Buv.neg[(k : ℕ)]?.getD 0 := by
    simpa [List.getD_eq_getElem?_getD] using huv0 k
  have hvu0' : Bvu.pos[(k : ℕ)]?.getD 0 = Bvu.neg[(k : ℕ)]?.getD 0 := by
    simpa [List.getD_eq_getElem?_getD] using hvu0 k
  cases a <;> cases b <;> simp [pickBlock, huv0', hvu0']

/-! ### The norm of `u + vω` from the lanes -/

/-- `a² - ab + b²`, the norm `|a + bω|²`. -/
def normSqZ (a b : ℤ) : ℤ := a ^ 2 - a * b + b ^ 2

/-- `∑_r v_r N(a_r, b_r)` with `N = normSqZ`, over the lanes `a_r` of `y` and `b_r` of `y'`. -/
def laneNormSq (P : Params) : ℕ → ℕ → List ℕ → ℤ
  | y, y', vr :: vs => vr * normSqZ (((y % 2 ^ P.w : ℕ) : ℤ) - P.B)
      (((y' % 2 ^ P.w : ℕ) : ℤ) - P.B) + laneNormSq P (y >>> P.w) (y' >>> P.w) vs
  | _, _, [] => 0

/-- **The weighted norms of the computed columns** of a packed run on `2d` rows, whose row `k`
holds the coordinates `u` and row `d + k` the coordinates `v`:
`∑_{k < d} ∑_r v_r N(X_{k, col r}, X_{d + k, col r})`. -/
def numeratorE (P : Params) (d : ℕ) (v ys : List ℕ) : ℤ :=
  (List.zipWith (fun y y' ↦ laneNormSq P y y' v) (ys.take d) (ys.drop d)).sum

theorem laneNormSq_eq (P : Params) : ∀ (v : List ℕ) (y y' : ℕ),
    laneNormSq P y y' v = ∑ r ∈ range v.length, v.getD r 0 *
      normSqZ ((((y >>> (P.w * r)) % 2 ^ P.w : ℕ) : ℤ) - P.B)
        ((((y' >>> (P.w * r)) % 2 ^ P.w : ℕ) : ℤ) - P.B)
  | [], y, y' => by simp [laneNormSq]
  | vr :: vs, y, y' => by
    rw [laneNormSq, laneNormSq_eq P vs, List.length_cons, sum_range_succ']
    simp only [List.getD_cons_succ, List.getD_cons_zero, mul_zero, Nat.shiftRight_zero]
    rw [add_comm]
    congr 1
    refine sum_congr rfl fun r _ ↦ ?_
    rw [← Nat.shiftRight_add, ← Nat.shiftRight_add, mul_add, mul_one, Nat.add_comm (P.w * r)]

theorem sum_zipWith_eqZ {f : ℕ → ℕ → ℤ} : ∀ (u ys : List ℕ) (d : ℕ), u.length = d →
    ys.length = d → (List.zipWith f u ys).sum = ∑ k ∈ range d, f (u.getD k 0) (ys.getD k 0)
  | [], [], 0, _, _ => by simp
  | a :: u, y :: ys, d + 1, hu, hy => by
    simp only [List.length_cons, Nat.add_right_cancel_iff] at hu hy
    rw [List.zipWith_cons_cons, List.sum_cons, sum_zipWith_eqZ u ys d hu hy, sum_range_succ',
      add_comm]
    simp
  | [], _ :: _, _, hu, hy => by simp at hu hy; omega
  | _ :: _, [], _, hu, hy => by simp at hu hy; omega

/-- The two halves of a list of `2d` rows, paired row by row. -/
theorem sum_zipWith_take_drop {f : ℕ → ℕ → ℤ} {ys : List ℕ} {d : ℕ} (h : ys.length = d + d) :
    (List.zipWith f (ys.take d) (ys.drop d)).sum =
      ∑ k ∈ range d, f (ys.getD k 0) (ys.getD (d + k) 0) := by
  rw [sum_zipWith_eqZ _ _ d (by simp [h]) (by simp [h])]
  refine sum_congr rfl fun k hk ↦ ?_
  have hk' : k < d := mem_range.1 hk
  simp [List.getD_eq_getElem?_getD, hk', List.getElem?_drop]

/-- **The weighted norms of the columns**, from a successful packed run on `2d` rows. -/
theorem runCols_numeratorE {P : Params} {D : Band} {d nc : ℕ} {col : ℕ → ℕ} {cst E : List ℕ}
    (hP : P.Ok) (hD : D.InRange (d + d)) (hR : D.RowsOk P cst (d + d))
    (hE : UnitsOk P (d + d) nc col E) {cs : List ℕ} (hcs : ∀ t ∈ cs, t ≤ P.B) {ys : List ℕ}
    (h : runCols P D cst (d + d) E nc cs = some ys) {v : List ℕ} (hv : v.length = nc) :
    numeratorE P d v ys = ∑ k ∈ range d, ∑ r ∈ range nc, (v.getD r 0 : ℤ) *
      normSqZ (specCol D P.s (col r) (cs.map (↑)) k)
        (specCol D P.s (col r) (cs.map (↑)) (d + k)) := by
  have hys := (runCols_spec hP hD hR hE cs hcs h).1
  rw [numeratorE, sum_zipWith_take_drop hys]
  refine sum_congr rfl fun k hk ↦ ?_
  have hk' : k < d := mem_range.1 hk
  rw [laneNormSq_eq, hv]
  refine sum_congr rfl fun r hr ↦ ?_
  rw [runCols_lane hP hD hR hE hcs h (by omega) (mem_range.1 hr),
    runCols_lane hP hD hR hE hcs h (by omega) (mem_range.1 hr)]

end CertifiedGaps.Trace.Packed
