/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.OpenTrialCompute

/-!
# The block iteration computes matrix products

The kernel computations of `OpenTrialCompute` are lists; this file proves what they compute, for
every input that passes the checks.

* `toVec_vecmat`: when its lane condition holds, the packed product `vecmat` is the row times the
  matrix, `toVec x ᵥ* t`.
* `Params.run_spec`: when the iteration passes all its lane conditions, its outputs are
  `(ℓ tᵐ ν, zₘ₋₁ ν, ℓ tᵐ⁻¹ ν′)` for `m = 34, 78, 142, 240`, with
  `z_k = 500000 ∑_{j < k} ℓ tʲ o t^(k-1-j) - 300000 ℓ′ tᵏ` (`Params.zVec`).
-/

namespace SeedKernel.OpenTrial

open Matrix SeedKernel.Trial

/-- The vector of a list, `xᵢ` (zero past its end). -/
def toVec (n : ℕ) (x : List ℤ) : Fin n → ℤ := fun i ↦ x.getD i 0

/-! ### Lists of length `n` -/

/-- The vector of an entrywise sum. -/
theorem toVec_vadd {n : ℕ} {x y : List ℤ} (hx : x.length = n) (hy : y.length = n) :
    toVec n (vadd x y) = toVec n x + toVec n y := by
  ext i
  exact getD_vadd x y hx hy i i.isLt

/-- An entrywise sum of lists of length `n` has length `n`. -/
theorem length_vadd_eq {n : ℕ} {x y : List ℤ} (hx : x.length = n) (hy : y.length = n) :
    (vadd x y).length = n := by
  rw [length_vadd, hx, hy, Nat.min_self]

/-- The vector of a scaled list. -/
theorem toVec_map_mul (n : ℕ) (k : ℤ) (x : List ℤ) :
    toVec n (x.map (k * ·)) = k • toVec n x := by
  ext i
  simp only [toVec, Pi.smul_apply, smul_eq_mul]
  rcases Nat.lt_or_ge (i : ℕ) x.length with h | h
  · rw [List.getD_eq_getElem _ _ (by simpa using h), List.getD_eq_getElem _ _ h,
      List.getElem_map]
  · rw [List.getD_eq_default _ _ (by simpa using h), List.getD_eq_default _ _ h, mul_zero]

/-- The list dot product is the dot product of the vectors. -/
theorem ldot_eq_dotProduct {n : ℕ} {x y : List ℤ} (hx : x.length = n) (hy : y.length = n) :
    ldot x y = toVec n x ⬝ᵥ toVec n y :=
  ldot_eq_sum n x y hx hy

/-! ### The packed vector-matrix product -/

/-- `maxAbs` bounds every entry. -/
theorem natAbs_le_maxAbs : ∀ (x : List ℤ), ∀ v ∈ x, v.natAbs ≤ maxAbs x
  | [], _, h => absurd h List.not_mem_nil
  | y :: ys, v, h => by
    rcases List.mem_cons.1 h with rfl | h
    · exact le_max_left _ _
    · exact (natAbs_le_maxAbs ys v h).trans (le_max_right _ _)

/-- `colsOf` packs `n` columns. -/
theorem length_colsOf (w : ℕ) (Y : List (List ℕ)) (n : ℕ) : (colsOf w Y n).length = n := by
  simp [colsOf]

/-- `vecmat` returns one entry per column. -/
theorem length_vecmat (w off n : ℕ) (Y : List (List ℕ)) (x : List ℤ) :
    (vecmat w off n (colsOf w Y n) x).length = n := by
  simp [vecmat, length_colsOf]

/-- **The packed product is exact**: entry `j` of `vecmat` is the dot product of the row with
column `j` of the offset matrix, when the lane condition holds. -/
theorem vecmat_getD {n off w : ℕ} {Y : List (List ℕ)} {x : List ℤ} (hx : x.length = n)
    (hn : 0 < n) (hY : ∀ r ∈ Y, r.length = n ∧ ∀ y ∈ r, y ≤ 2 * off) (hYn : Y.length = n)
    (hok : vecmatOk w off n x = true) (j : ℕ) (hj : j < n) :
    (vecmat w off n (colsOf w Y n) x).getD j 0 = ldot x (lcol (shiftL off Y) j) := by
  set β := maxAbs x
  set xs := x.map fun v ↦ (v + (β : ℤ)).toNat with hxs
  have hb : ∀ v ∈ x, v.natAbs ≤ β := natAbs_le_maxAbs x
  have hxs_len : xs.length = n := by simp [xs, hx]
  have hxs_le : ∀ y ∈ xs, y ≤ 2 * β := by
    intro y hy
    simp only [xs, List.mem_map] at hy
    obtain ⟨v, hv, rfl⟩ := hy
    have := hb v hv
    omega
  have hxs_sh : xs.map (sh β) = x := by
    simp only [xs, List.map_map]
    conv_rhs => rw [← List.map_id x]
    apply List.map_congr_left
    intro v hv
    have := hb v hv
    simp only [Function.comp, sh, id]
    omega
  have hne : xs ≠ [] := by
    intro h; rw [h] at hxs_len; simp at hxs_len; omega
  have hw : xs.length * xs.length * (2 * β) * (2 * off) < 2 ^ w := by
    rw [hxs_len]; simpa [vecmatOk] using hok
  have hcol_len : (colN Y j).length = xs.length := by rw [length_colN, hYn, hxs_len]
  have hcol_le := colN_le Y j (2 * off) fun r hr ↦ (hY r hr).2
  have he := entryB_eq w β off xs (colN Y j) hcol_len.symm hne hxs_le hcol_le hw
  rw [hxs_len, hxs_sh] at he
  rw [lcol_shiftL off Y j fun r hr ↦ by rw [(hY r hr).1]; exact hj]
  have hj' : j < (colsOf w Y n).length := by rw [length_colsOf]; exact hj
  simp only [vecmat]
  rw [List.getD_eq_getElem _ _ (by simpa using hj'), List.getElem_map]
  simp only [colsOf, List.getElem_map, List.getElem_range]
  rw [← hxs, he]
  push_cast
  ring

/-- **The packed product is the vector-matrix product** `x ↦ x t`. -/
theorem toVec_vecmat {n w : ℕ} {t : OM} (ht : t.Good n) {x : List ℤ} (hx : x.length = n)
    (hok : vecmatOk w t.1 n x = true) :
    toVec n (vecmat w t.1 n (colsOf w t.2 n) x) = toVec n x ᵥ* t.mat n := by
  have hn : 0 < n := by
    have := ht.2.1; rw [← List.length_pos_iff, ht.1] at this; exact this
  ext j
  rw [toVec, vecmat_getD hx hn (fun r hr ↦ ht.2.2 r hr) ht.1 hok j j.isLt,
    ldot_eq_sum n _ _ hx (by rw [length_lcol, shiftL, List.length_map, ht.1])]
  simp only [vecMul, dotProduct, toVec, OM.mat, toMatB, toMat_apply, getD_lcol]

/-! ### The block iteration -/

namespace Params

variable (P : Params)

/-- The block `t` as a matrix. -/
def T : Matrix (Fin P.n) (Fin P.n) ℤ := P.t.mat P.n
/-- The block `o` as a matrix. -/
def O : Matrix (Fin P.n) (Fin P.n) ℤ := P.o.mat P.n
/-- `ℓ`. -/
def lv : Fin P.n → ℤ := toVec P.n P.l
/-- `ℓ′`. -/
def lv' : Fin P.n → ℤ := toVec P.n P.l'
/-- `ν`. -/
def vv : Fin P.n → ℤ := toVec P.n P.v
/-- `ν′`. -/
def vv' : Fin P.n → ℤ := toVec P.n P.v'

/-- The accumulated row `z_{k+1} = 500000 ∑_{j < k} ℓ tʲ o t^(k-1-j) - 300000 ℓ′ tᵏ`. -/
def zVec (k : ℕ) : Fin P.n → ℤ :=
  (500000 : ℤ) • ∑ j ∈ Finset.range k, P.lv ᵥ* (P.T ^ j * P.O * P.T ^ (k - 1 - j)) -
    (300000 : ℤ) • (P.lv' ᵥ* P.T ^ k)

/-- The expected outputs `(ℓ t^(k+1) ν, z_{k+1} ν, ℓ tᵏ ν′)` of the state at `r = k + 1`. -/
def expected (k : ℕ) : ℤ × ℤ × ℤ :=
  ((P.lv ᵥ* P.T ^ (k + 1)) ⬝ᵥ P.vv, P.zVec k ⬝ᵥ P.vv, (P.lv ᵥ* P.T ^ k) ⬝ᵥ P.vv')

/-- The well-formedness of the inputs. -/
structure Good : Prop where
  t : P.t.Good P.n
  o : P.o.Good P.n
  l : P.l.length = P.n
  l' : P.l'.length = P.n
  v : P.v.length = P.n
  v' : P.v'.length = P.n

/-- The state at `r = k + 1`: `(ℓ tᵏ, ℓ t^(k+1), z_{k+1})` as lists of length `n`. -/
def Inv (k : ℕ) (s : List ℤ × List ℤ × List ℤ) : Prop :=
  s.1.length = P.n ∧ s.2.1.length = P.n ∧ s.2.2.length = P.n ∧
    toVec P.n s.1 = P.lv ᵥ* P.T ^ k ∧ toVec P.n s.2.1 = P.lv ᵥ* P.T ^ (k + 1) ∧
      toVec P.n s.2.2 = P.zVec k

variable {P}

/-- The product by `t` of the iteration is `x ↦ x t`. -/
theorem toVec_mt (hP : P.Good) {x : List ℤ} (hx : x.length = P.n)
    (hok : vecmatOk P.wt P.t.1 P.n x = true) : toVec P.n (P.mt x) = toVec P.n x ᵥ* P.T :=
  toVec_vecmat hP.t hx hok

/-- The product by `o` of the iteration is `x ↦ x o`. -/
theorem toVec_mo (hP : P.Good) {x : List ℤ} (hx : x.length = P.n)
    (hok : vecmatOk P.wo P.o.1 P.n x = true) : toVec P.n (P.mo x) = toVec P.n x ᵥ* P.O :=
  toVec_vecmat hP.o hx hok

/-- The product by `t` has length `n`. -/
theorem length_mt (x : List ℤ) : (P.mt x).length = P.n := length_vecmat _ _ _ _ _

/-- The product by `o` has length `n`. -/
theorem length_mo (x : List ℤ) : (P.mo x).length = P.n := length_vecmat _ _ _ _ _

/-- The recursion of the accumulated row: `z_{k+2} = z_{k+1} t + 500000 ℓ tᵏ o`. -/
theorem zVec_succ (k : ℕ) :
    P.zVec (k + 1) = P.zVec k ᵥ* P.T + (500000 : ℤ) • ((P.lv ᵥ* P.T ^ k) ᵥ* P.O) := by
  simp only [zVec, sub_vecMul, smul_vecMul, sum_vecMul, vecMul_vecMul, Finset.sum_range_succ,
    Nat.add_sub_cancel, Nat.sub_self, pow_zero, Matrix.mul_one, smul_add]
  have h : ∀ j ∈ Finset.range k, P.lv ᵥ* (P.T ^ j * P.O * P.T ^ (k - 1 - j) * P.T) =
      P.lv ᵥ* (P.T ^ j * P.O * P.T ^ (k - j)) := by
    intro j hj
    rw [Finset.mem_range] at hj
    rw [Matrix.mul_assoc, ← pow_succ, show k - 1 - j + 1 = k - j by omega]
  rw [Finset.sum_congr rfl h, ← pow_succ]
  abel

/-- One step preserves the invariant when its lane conditions hold. -/
theorem inv_step (hP : P.Good) {k : ℕ} {s : List ℤ × List ℤ × List ℤ} (hs : P.Inv k s)
    (hok : P.stepOk s = true) : P.Inv (k + 1) (P.step s) := by
  obtain ⟨h1, h2, h3, e1, e2, e3⟩ := hs
  simp only [stepOk, Bool.and_eq_true] at hok
  obtain ⟨⟨ok1, ok2⟩, ok3⟩ := hok
  refine ⟨h2, length_mt _, length_vadd_eq (length_mt _) (by simp [length_mo]), e2, ?_, ?_⟩
  · simp only [step]
    rw [toVec_mt hP h2 ok1, e2, vecMul_vecMul, ← pow_succ]
  · simp only [step]
    rw [toVec_vadd (length_mt _) (by simp [length_mo]), toVec_map_mul, toVec_mt hP h3 ok2,
      toVec_mo hP h1 ok3, e3, e1, zVec_succ]

/-- The outputs of a state satisfying the invariant. -/
theorem out_eq (hP : P.Good) {k : ℕ} {s : List ℤ × List ℤ × List ℤ} (hs : P.Inv k s) :
    P.out s = P.expected k := by
  obtain ⟨h1, h2, h3, e1, e2, e3⟩ := hs
  simp only [out, expected, ldot_eq_dotProduct h2 hP.v, ldot_eq_dotProduct h3 hP.v,
    ldot_eq_dotProduct h1 hP.v', e1, e2, e3, vv, vv']

/-- A passed final lane check implies the lane check passed before the loop. -/
theorem loop_ok : ∀ (f r : ℕ) (s : List ℤ × List ℤ × List ℤ) (acc : List (ℤ × ℤ × ℤ))
    (ok : Bool), (P.loop f r s acc ok).2 = true → ok = true
  | 0, _, _, _, _, h => h
  | f + 1, r, s, acc, ok, h => by
    have := loop_ok f (r + 1) (P.step s) _ _ h
    simp only [Bool.and_eq_true] at this
    exact this.1

/-- **The loop computes the expected outputs** at the `lengths` it passes. -/
theorem loop_spec (hP : P.Good) : ∀ (f k : ℕ) (s : List ℤ × List ℤ × List ℤ)
    (acc : List (ℤ × ℤ × ℤ)) (ok : Bool), P.Inv k s → (P.loop f (k + 1) s acc ok).2 = true →
      (P.loop f (k + 1) s acc ok).1 = acc ++ ((List.range f).filter
        fun i ↦ lengths.contains (k + 1 + i)).map fun i ↦ P.expected (k + i)
  | 0, _, _, _, _, _, _ => by simp [loop]
  | f + 1, k, s, acc, ok, hs, h => by
    have hok := loop_ok f (k + 1 + 1) (P.step s) _ _ h
    simp only [Bool.and_eq_true] at hok
    have ih := loop_spec hP f (k + 1) (P.step s)
      (if lengths.contains (k + 1) then acc ++ [P.out s] else acc) (ok && P.stepOk s)
      (inv_step hP hs hok.2) h
    simp only [loop]
    rw [ih, List.range_succ_eq_map, List.filter_cons, out_eq hP hs]
    have e1 : ∀ i, k + 1 + 1 + i = k + 1 + (i + 1) := fun i ↦ by ring
    have e2 : ∀ i, k + 1 + i = k + (i + 1) := fun i ↦ by ring
    simp only [List.filter_map, Function.comp_def, e1, e2]
    split_ifs <;> simp

/-- The invariant of the initial state, when its lane condition holds. -/
theorem inv_init (hP : P.Good) (hok : vecmatOk P.wt P.t.1 P.n P.l = true) : P.Inv 0 P.init := by
  refine ⟨hP.l, length_mt _, by simp [init, hP.l'], ?_, ?_, ?_⟩
  · change toVec P.n P.l = P.lv ᵥ* P.T ^ 0
    rw [pow_zero, vecMul_one]
    rfl
  · simp only [init]; rw [toVec_mt hP hP.l hok, zero_add, pow_one]; rfl
  · simp only [init, zVec, Finset.range_zero, Finset.sum_empty, smul_zero, zero_sub, pow_zero,
      vecMul_one]
    rw [toVec_map_mul, neg_smul]
    rfl

/-- **The block iteration**: when every lane condition holds, the outputs are
`(ℓ tᵐ ν, zₘ ν, ℓ tᵐ⁻¹ ν′)` at `m = 34, 78, 142, 240`. -/
theorem run_spec (hP : P.Good) (h : P.run.2 = true) :
    P.run.1 = [P.expected 33, P.expected 77, P.expected 141, P.expected 239] := by
  have hok := loop_ok 240 1 P.init [] _ h
  have := loop_spec hP 240 0 P.init [] _ (inv_init hP hok) h
  rw [run, this]
  have hr : ((List.range 240).filter fun i ↦ lengths.contains (0 + 1 + i)) = [33, 77, 141, 239] :=
    by decide
  rw [hr]
  simp

end Params

end SeedKernel.OpenTrial
