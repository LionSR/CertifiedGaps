/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.List.GetD
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.GCongr

/-!
# Integer list matrices for kernel computation

The trial-state contractions are computed by the kernel (`decide +kernel`) on integer matrices
stored as lists of rows. This file connects that representation to Mathlib's `Matrix`:

* `toMat n m X`: the `n × m` matrix whose entry `(i, j)` is `X[i][j]`, zero outside `X`.
* `lmul m X Y`: the schoolbook product of two list matrices; `toMat_lmul` identifies it with the
  matrix product.
* `lmulB a b X Y m`: the product of matrices stored as naturals with offsets `a`, `b` (entries
  `x - a`), computed with one big multiplication per entry (Kronecker substitution): row and
  reversed column are packed in radix `2^w` and the dot product is read off the middle digit.
  `lmulB_spec` proves it equal to `lmul` for every input within its offsets, since the radix is
  chosen from the offsets at run time; `toMatB_lmulB` is the matrix form.
* `radix_eq_zero`: digits of absolute value below the radix are determined by the number.
* `decodeN`, `decodeSparse`: the generated data, stored as single `Nat` literals.

The kernel evaluates `Nat` arithmetic on literals with GMP but pays a large fixed cost per
reduction step (and more for `Int`), so one big multiplication per entry is far cheaper than `n`
small ones.
-/

namespace SeedKernel.Trial

open Matrix

/-! ### Dot products and list matrices -/

/-- The dot product of two integer lists, truncated to the shorter one. -/
def ldot : List ℤ → List ℤ → ℤ
  | x :: xs, y :: ys => x * y + ldot xs ys
  | _, _ => 0

theorem ldot_eq_sum : ∀ (k : ℕ) (xs ys : List ℤ), xs.length = k → ys.length = k →
    ldot xs ys = ∑ l : Fin k, xs.getD l 0 * ys.getD l 0
  | 0, [], [], _, _ => by simp [ldot]
  | k + 1, x :: xs, y :: ys, hx, hy => by
    rw [Fin.sum_univ_succ, ldot, ldot_eq_sum k xs ys (by simpa using hx) (by simpa using hy)]
    simp

/-- Entry `(i, j)` of a list of rows, zero outside. -/
def lget (X : List (List ℤ)) (i j : ℕ) : ℤ := (X.getD i []).getD j 0

/-- The `n × m` integer matrix stored as the list of rows `X`. -/
def toMat (n m : ℕ) (X : List (List ℤ)) : Matrix (Fin n) (Fin m) ℤ :=
  Matrix.of fun i j ↦ lget X i j

@[simp] theorem toMat_apply (n m : ℕ) (X : List (List ℤ)) (i : Fin n) (j : Fin m) :
    toMat n m X i j = lget X i j := rfl

/-- Column `j` of a list matrix. -/
def lcol (Y : List (List ℤ)) (j : ℕ) : List ℤ := Y.map fun r ↦ r.getD j 0

theorem length_lcol (Y : List (List ℤ)) (j : ℕ) : (lcol Y j).length = Y.length := by
  simp [lcol]

theorem getD_lcol (Y : List (List ℤ)) (j l : ℕ) : (lcol Y j).getD l 0 = lget Y l j := by
  rcases Nat.lt_or_ge l Y.length with h | h
  · simp [lcol, lget, h]
  · simp [lcol, lget, h]

theorem getD_getD_of_ge (X : List (List ℤ)) (i j : ℕ) (h : X.length ≤ i) :
    (X.getD i []).getD j 0 = 0 := by
  simp [h]

/-- The schoolbook product of list matrices with `m` result columns. -/
def lmul (m : ℕ) (X Y : List (List ℤ)) : List (List ℤ) :=
  X.map fun r ↦ (List.range m).map fun j ↦ ldot r (lcol Y j)

theorem length_lmul (m : ℕ) (X Y : List (List ℤ)) : (lmul m X Y).length = X.length := by
  simp [lmul]

theorem length_of_mem_lmul (m : ℕ) (X Y : List (List ℤ)) :
    ∀ r ∈ lmul m X Y, r.length = m := by
  simp [lmul]

theorem lget_lmul {k : ℕ} (m : ℕ) (X Y : List (List ℤ)) (hX : ∀ r ∈ X, r.length = k)
    (hY : Y.length = k) (i j : ℕ) (hj : j < m) :
    lget (lmul m X Y) i j = ∑ l : Fin k, lget X i l * lget Y l j := by
  rcases Nat.lt_or_ge i X.length with hi | hi
  · have hrow : (X.getD i []).length = k := by
      rw [List.getD_eq_getElem _ _ hi]; exact hX _ (List.getElem_mem hi)
    have h1 : lget (lmul m X Y) i j = ldot (X.getD i []) (lcol Y j) := by
      simp [lget, lmul, hi, hj]
    rw [h1, ldot_eq_sum k _ _ hrow (by rw [length_lcol, hY])]
    simp only [getD_lcol, lget]
  · have h1 : lget (lmul m X Y) i j = 0 :=
      getD_getD_of_ge _ _ _ (by rw [length_lmul]; exact hi)
    rw [h1]
    simp [lget, hi]

/-- The list product is the matrix product. -/
theorem toMat_lmul {n k m : ℕ} (X Y : List (List ℤ)) (hX : ∀ r ∈ X, r.length = k)
    (hY : Y.length = k) : toMat n m (lmul m X Y) = toMat n k X * toMat k m Y := by
  ext i j
  rw [toMat_apply, lget_lmul m X Y hX hY i j j.isLt, mul_apply]
  rfl

/-! ### Packed products of offset matrices

An integer matrix with entries in `[-β, β]` is stored as the natural numbers `x + β` (offset
`β`). The kernel multiplies two such matrices with one big `Nat` product per entry: the row and
the reversed column are packed in radix `2^w`, the middle digit of their product is the dot
product of the stored entries, and a correction by row and column sums gives the stored entry
of the product, with offset `k a b`. -/

/-- The radix-`2^w` value `∑ᵢ xs[i] 2^(w i)`. -/
def packN (w : ℕ) : List ℕ → ℕ
  | [] => 0
  | x :: xs => x + (packN w xs <<< w)

theorem packN_cons (w x : ℕ) (xs : List ℕ) : packN w (x :: xs) = x + 2 ^ w * packN w xs := by
  rw [packN, Nat.shiftLeft_eq, Nat.mul_comm]

theorem packN_append_singleton (w y : ℕ) : ∀ l : List ℕ,
    packN w (l ++ [y]) = packN w l + 2 ^ (w * l.length) * y
  | [] => by simp [packN]
  | x :: xs => by
    rw [List.cons_append, packN_cons, packN_cons, packN_append_singleton w y xs,
      List.length_cons, Nat.mul_succ, pow_add]
    ring

/-- Dot product of lists of naturals. -/
def ldotN : List ℕ → List ℕ → ℕ
  | x :: xs, y :: ys => x * y + ldotN xs ys
  | _, _ => 0

/-- Sum of a list of naturals. -/
def sumN (l : List ℕ) : ℕ := l.foldr (· + ·) 0

theorem two_pow_mul_packN_le (w B : ℕ) : ∀ l : List ℕ, (∀ x ∈ l, x ≤ B) →
    2 ^ w * packN w l ≤ B * l.length * 2 ^ (w * l.length)
  | [], _ => by simp [packN]
  | x :: xs, h => by
    have hx : x ≤ B := h x (by simp)
    have ih := two_pow_mul_packN_le w B xs (fun y hy ↦ h y (by simp [hy]))
    have e : 2 ^ (w * (xs.length + 1)) = 2 ^ w * 2 ^ (w * xs.length) := by
      rw [Nat.mul_succ, pow_add, Nat.mul_comm]
    have h1 : 1 ≤ 2 ^ (w * xs.length) := Nat.one_le_two_pow
    rw [packN_cons, List.length_cons, e]
    calc 2 ^ w * (x + 2 ^ w * packN w xs) = 2 ^ w * x + 2 ^ w * (2 ^ w * packN w xs) := by ring
      _ ≤ 2 ^ w * B * 2 ^ (w * xs.length) + 2 ^ w * (B * xs.length * 2 ^ (w * xs.length)) := by
          apply Nat.add_le_add
          · calc 2 ^ w * x ≤ 2 ^ w * B := Nat.mul_le_mul_left _ hx
              _ ≤ 2 ^ w * B * 2 ^ (w * xs.length) := Nat.le_mul_of_pos_right _ h1
          · exact Nat.mul_le_mul_left _ ih
      _ = B * (xs.length + 1) * (2 ^ w * 2 ^ (w * xs.length)) := by ring

theorem ldotN_le (A B : ℕ) : ∀ xs ys : List ℕ, (∀ x ∈ xs, x ≤ A) → (∀ y ∈ ys, y ≤ B) →
    ldotN xs ys ≤ xs.length * A * B
  | [], _, _, _ => by simp [ldotN]
  | _ :: _, [], _, _ => by simp [ldotN]
  | x :: xs, y :: ys, hx, hy => by
    have ih := ldotN_le A B xs ys (fun a ha ↦ hx a (by simp [ha])) (fun a ha ↦ hy a (by simp [ha]))
    have h1 : x * y ≤ A * B := Nat.mul_le_mul (hx x (by simp)) (hy y (by simp))
    rw [ldotN, List.length_cons]
    nlinarith

/-- The product of a packing and a reversed packing splits into a low part, the dot product at
digit `k - 1`, and a high part. -/
theorem packN_mul_reverse (w A B : ℕ) : ∀ (xs ys : List ℕ), xs.length = ys.length → xs ≠ [] →
    (∀ x ∈ xs, x ≤ A) → (∀ y ∈ ys, y ≤ B) →
    ∃ lo hi, packN w xs * packN w ys.reverse =
        lo + ldotN xs ys * 2 ^ (w * (xs.length - 1)) + 2 ^ (w * xs.length) * hi ∧
      2 ^ w * lo ≤ xs.length * xs.length * A * B * 2 ^ (w * (xs.length - 1))
  | [], _, _, h, _, _ => absurd rfl h
  | [x], [y], _, _, _, _ => ⟨0, 0, by simp [packN, ldotN]⟩
  | [_], [], h, _, _, _ => by simp at h
  | [_], _ :: _ :: _, h, _, _, _ => by simp at h
  | _ :: _ :: _, [], h, _, _, _ => by simp at h
  | x :: x' :: xs, y :: ys, hlen, _, hx, hy => by
    have hlen' : (x' :: xs).length = ys.length := by simpa using hlen
    obtain ⟨lo, hi, heq, hb⟩ := packN_mul_reverse w A B (x' :: xs) ys hlen' (by simp)
      (fun a ha ↦ hx a (by simp_all)) (fun a ha ↦ hy a (by simp [ha]))
    set n := xs.length with hn
    have hys : ys.length = n + 1 := by simpa using hlen'.symm
    have hrev : packN w (y :: ys).reverse = packN w ys.reverse + 2 ^ (w * (n + 1)) * y := by
      rw [List.reverse_cons, packN_append_singleton, List.length_reverse, hys]
    have hxB : x ≤ A := hx x (by simp)
    have hyB : y ≤ B := hy y (by simp)
    have hY : 2 ^ w * packN w ys.reverse ≤ B * (n + 1) * 2 ^ (w * (n + 1)) := by
      have := two_pow_mul_packN_le w B ys.reverse (fun a ha ↦ hy a (by simp_all))
      simpa [hys] using this
    simp only [List.length_cons] at heq hb ⊢
    simp only [Nat.add_sub_cancel] at heq hb ⊢
    refine ⟨x * packN w ys.reverse + 2 ^ w * lo, hi + packN w (x' :: xs) * y, ?_, ?_⟩
    · rw [packN_cons w x (x' :: xs), hrev, ldotN]
      have e1 : 2 ^ (w * (n + 1)) = 2 ^ w * 2 ^ (w * n) := by
        rw [Nat.mul_succ, pow_add, Nat.mul_comm]
      have e2 : 2 ^ (w * (n + 1 + 1)) = 2 ^ w * 2 ^ (w * (n + 1)) := by
        rw [show w * (n + 1 + 1) = w + w * (n + 1) by ring, pow_add]
      rw [e2]
      calc (x + 2 ^ w * packN w (x' :: xs)) * (packN w ys.reverse + 2 ^ (w * (n + 1)) * y)
          = x * packN w ys.reverse + x * y * 2 ^ (w * (n + 1)) +
              2 ^ w * (packN w (x' :: xs) * packN w ys.reverse) +
              2 ^ w * 2 ^ (w * (n + 1)) * (packN w (x' :: xs) * y) := by ring
        _ = _ := by rw [heq, e1]; ring
    · have e1 : 2 ^ (w * (n + 1)) = 2 ^ w * 2 ^ (w * n) := by
        rw [Nat.mul_succ, pow_add, Nat.mul_comm]
      calc 2 ^ w * (x * packN w ys.reverse + 2 ^ w * lo)
          = x * (2 ^ w * packN w ys.reverse) + 2 ^ w * (2 ^ w * lo) := by ring
        _ ≤ A * (B * (n + 1) * 2 ^ (w * (n + 1))) +
              2 ^ w * ((n + 1) * (n + 1) * A * B * 2 ^ (w * n)) := by gcongr
        _ = (A * B * (n + 1) + (n + 1) * (n + 1) * A * B) * 2 ^ (w * (n + 1)) := by
            rw [e1]; ring
        _ ≤ (n + 1 + 1) * (n + 1 + 1) * A * B * 2 ^ (w * (n + 1)) := by
            have h2 : (n + 1) + (n + 1) * (n + 1) ≤ (n + 1 + 1) * (n + 1 + 1) := by nlinarith
            have h3 := Nat.mul_le_mul_right (A * B * 2 ^ (w * (n + 1))) h2
            calc (A * B * (n + 1) + (n + 1) * (n + 1) * A * B) * 2 ^ (w * (n + 1))
                = ((n + 1) + (n + 1) * (n + 1)) * (A * B * 2 ^ (w * (n + 1))) := by ring
              _ ≤ (n + 1 + 1) * (n + 1 + 1) * (A * B * 2 ^ (w * (n + 1))) := h3
              _ = _ := by ring

/-- **Packed dot product.** When every convolution coefficient fits in a digit, digit `k - 1` of
the product of a packing and a reversed packing is the dot product. -/
theorem packN_dot (w A B : ℕ) (xs ys : List ℕ) (hlen : xs.length = ys.length) (hne : xs ≠ [])
    (hx : ∀ x ∈ xs, x ≤ A) (hy : ∀ y ∈ ys, y ≤ B)
    (hw : xs.length * xs.length * A * B < 2 ^ w) :
    ((packN w xs * packN w ys.reverse) >>> (w * (xs.length - 1))) % 2 ^ w = ldotN xs ys := by
  obtain ⟨lo, hi, heq, hb⟩ := packN_mul_reverse w A B xs ys hlen hne hx hy
  have hk : 1 ≤ xs.length := List.length_pos_iff.mpr hne
  set D := 2 ^ (w * (xs.length - 1)) with hD
  have hDpos : 0 < D := Nat.two_pow_pos _
  have hlo : lo < D := by
    by_contra h
    push Not at h
    have : 2 ^ w * D ≤ 2 ^ w * lo := Nat.mul_le_mul_left _ h
    have h2 : xs.length * xs.length * A * B * D < 2 ^ w * D := Nat.mul_lt_mul_of_pos_right hw hDpos
    omega
  have hd : ldotN xs ys < 2 ^ w := by
    have := ldotN_le A B xs ys hx hy
    have h3 : xs.length * A * B ≤ xs.length * xs.length * A * B := by
      have := Nat.mul_le_mul_right (A * B) (Nat.le_mul_of_pos_left xs.length hk)
      simpa [Nat.mul_assoc] using this
    omega
  have hpow : 2 ^ (w * xs.length) = D * 2 ^ w := by
    rw [hD, ← pow_add]
    congr 1
    obtain ⟨k, hk'⟩ : ∃ k, xs.length = k + 1 := ⟨xs.length - 1, by omega⟩
    rw [hk']
    simp [Nat.mul_succ]
  rw [heq, Nat.shiftRight_eq_div_pow, ← hD, hpow]
  have : lo + ldotN xs ys * D + D * 2 ^ w * hi = lo + D * (ldotN xs ys + 2 ^ w * hi) := by ring
  rw [this, Nat.add_mul_div_left _ _ hDpos, Nat.div_eq_of_lt hlo, Nat.zero_add,
    Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hd]

/-- The stored entry of a product of offset matrices: digit `k - 1` of the packed product of a
row and a reversed column, corrected by their sums to the offset `k a b`. -/
def entryB (w k a b px sx py sy : ℕ) : ℕ :=
  ((px * py) >>> (w * (k - 1))) % 2 ^ w + 2 * k * a * b - (b * sx + a * sy)

/-- Column `j` of a list of rows of naturals. -/
def colN (Y : List (List ℕ)) (j : ℕ) : List ℕ := Y.map fun r ↦ r.getD j 0

/-- Packed product of offset matrices: `X` stored with offset `a`, `Y` with offset `b`, and `m`
result columns. The result is stored with offset `Y.length * a * b`; the radix is chosen from
the offsets, so the product is exact for all inputs within their offsets (`lmulB_spec`). -/
def lmulB (a b : ℕ) (X Y : List (List ℕ)) (m : ℕ) : List (List ℕ) :=
  let k := Y.length
  let w := Nat.log2 (4 * k * k * a * b) + 1
  let cols := (List.range m).map fun j ↦ (packN w (colN Y j).reverse, sumN (colN Y j))
  X.map fun r ↦ cols.map fun c ↦ entryB w k a b (packN w r) (sumN r) c.1 c.2

/-- The integer `x - β` stored as `x` with offset `β`. -/
def sh (β x : ℕ) : ℤ := (x : ℤ) - β

/-- The integer matrix stored as `X` with offset `β`. -/
def shiftL (β : ℕ) (X : List (List ℕ)) : List (List ℤ) :=
  X.map (·.map (sh β))

theorem ldot_shift (a b : ℕ) : ∀ r c : List ℕ, r.length = c.length →
    ldot (r.map (sh a)) (c.map (sh b)) =
      ldotN r c - b * sumN r - a * sumN c + r.length * a * b
  | [], [], _ => by simp [ldot, ldotN, sumN]
  | x :: r, y :: c, h => by
    have ih := ldot_shift a b r c (by simpa using h)
    have hr : sumN (x :: r) = x + sumN r := rfl
    have hc : sumN (y :: c) = y + sumN c := rfl
    simp only [List.map_cons, ldot, ldotN, List.length_cons, hr, hc]
    rw [ih]
    unfold sh
    push_cast
    ring

theorem abs_ldot_shift_le (a b : ℕ) : ∀ r c : List ℕ, (∀ x ∈ r, x ≤ 2 * a) →
    (∀ y ∈ c, y ≤ 2 * b) →
    |ldot (r.map (sh a)) (c.map (sh b))| ≤ r.length * a * b
  | [], _, _, _ => by simp [ldot]
  | _ :: _, [], _, _ => by
    simp only [List.map_cons, List.map_nil, ldot, abs_zero]
    positivity
  | x :: r, y :: c, hr, hc => by
    have ih := abs_ldot_shift_le a b r c (fun z hz ↦ hr z (by simp [hz]))
      (fun z hz ↦ hc z (by simp [hz]))
    have hx : |(x : ℤ) - a| ≤ a := by
      have := hr x (by simp); rw [abs_le]; constructor <;> omega
    have hy : |(y : ℤ) - b| ≤ b := by
      have := hc y (by simp); rw [abs_le]; constructor <;> omega
    simp only [List.map_cons, ldot, List.length_cons]
    calc |((x : ℤ) - a) * ((y : ℤ) - b) + ldot (r.map (sh a))
          (c.map (sh b))|
        ≤ |((x : ℤ) - a) * ((y : ℤ) - b)| + |ldot (r.map (sh a))
          (c.map (sh b))| := abs_add_le _ _
      _ ≤ a * b + r.length * a * b := by
          rw [abs_mul]; gcongr
      _ = ((r.length + 1 : ℕ) : ℤ) * a * b := by push_cast; ring

/-- One stored entry of the packed product is the offset dot product. -/
theorem entryB_eq (w a b : ℕ) (r c : List ℕ) (hlen : r.length = c.length) (hne : r ≠ [])
    (hr : ∀ x ∈ r, x ≤ 2 * a) (hc : ∀ y ∈ c, y ≤ 2 * b)
    (hw : r.length * r.length * (2 * a) * (2 * b) < 2 ^ w) :
    (entryB w r.length a b (packN w r) (sumN r) (packN w c.reverse) (sumN c) : ℤ) =
      ldot (r.map (sh a)) (c.map (sh b)) + r.length * a * b := by
  have hmid := packN_dot w (2 * a) (2 * b) r c hlen hne hr hc hw
  have hs := ldot_shift a b r c hlen
  have hb := abs_ldot_shift_le a b r c hr hc
  rw [abs_le] at hb
  have hle : b * sumN r + a * sumN c ≤ ldotN r c + 2 * r.length * a * b := by
    have : ((b * sumN r + a * sumN c : ℕ) : ℤ) ≤ (ldotN r c + 2 * r.length * a * b : ℕ) := by
      push_cast; linarith
    exact_mod_cast this
  unfold entryB
  rw [hmid, Nat.cast_sub hle]
  push_cast
  linarith

theorem entryB_le (w a b : ℕ) (r c : List ℕ) (hlen : r.length = c.length) (hne : r ≠ [])
    (hr : ∀ x ∈ r, x ≤ 2 * a) (hc : ∀ y ∈ c, y ≤ 2 * b)
    (hw : r.length * r.length * (2 * a) * (2 * b) < 2 ^ w) :
    entryB w r.length a b (packN w r) (sumN r) (packN w c.reverse) (sumN c) ≤
      2 * (r.length * a * b) := by
  have h := entryB_eq w a b r c hlen hne hr hc hw
  have hb := abs_ldot_shift_le a b r c hr hc
  rw [abs_le] at hb
  have : ((entryB w r.length a b (packN w r) (sumN r) (packN w c.reverse) (sumN c) : ℕ) : ℤ) ≤
      ((2 * (r.length * a * b) : ℕ) : ℤ) := by push_cast; linarith
  exact_mod_cast this

theorem length_colN (Y : List (List ℕ)) (j : ℕ) : (colN Y j).length = Y.length := by
  simp [colN]

theorem colN_le (Y : List (List ℕ)) (j β : ℕ) (hY : ∀ r ∈ Y, ∀ y ∈ r, y ≤ β) :
    ∀ y ∈ colN Y j, y ≤ β := by
  intro y hy
  simp only [colN, List.mem_map] at hy
  obtain ⟨r, hr, rfl⟩ := hy
  rcases Nat.lt_or_ge j r.length with h | h
  · rw [List.getD_eq_getElem _ _ h]; exact hY r hr _ (List.getElem_mem h)
  · simp [h]

theorem lcol_shiftL (β : ℕ) (Y : List (List ℕ)) (j : ℕ) (hY : ∀ r ∈ Y, j < r.length) :
    lcol (shiftL β Y) j = (colN Y j).map (sh β) := by
  simp only [lcol, shiftL, colN, List.map_map]
  apply List.map_congr_left
  intro r hr
  simp [List.getElem?_eq_getElem (hY r hr)]

theorem radix_bound (k a b : ℕ) :
    k * k * (2 * a) * (2 * b) < 2 ^ (Nat.log2 (4 * k * k * a * b) + 1) := by
  have := @Nat.lt_log2_self (4 * k * k * a * b)
  calc k * k * (2 * a) * (2 * b) = 4 * k * k * a * b := by ring
    _ < _ := this

/-- **The packed product is exact.** -/
theorem lmulB_spec (a b m : ℕ) (X Y : List (List ℕ)) (hk : Y ≠ [])
    (hX : ∀ r ∈ X, r.length = Y.length ∧ ∀ x ∈ r, x ≤ 2 * a)
    (hY : ∀ r ∈ Y, r.length = m ∧ ∀ y ∈ r, y ≤ 2 * b) :
    shiftL (Y.length * a * b) (lmulB a b X Y m) = lmul m (shiftL a X) (shiftL b Y) := by
  simp only [lmulB, lmul, shiftL, List.map_map]
  apply List.map_congr_left
  intro r hr
  obtain ⟨hrl, hrb⟩ := hX r hr
  have hne : r ≠ [] := by
    intro h; rw [h] at hrl; exact hk (List.length_eq_zero_iff.mp hrl.symm)
  simp only [Function.comp_def, List.map_map]
  apply List.map_congr_left
  intro j hj
  rw [List.mem_range] at hj
  have hcol := lcol_shiftL b Y j (fun r hr ↦ by rw [(hY r hr).1]; exact hj)
  rw [show (shiftL b Y) = Y.map (·.map (sh b)) from rfl] at hcol
  rw [hcol]
  have hw := radix_bound Y.length a b
  rw [← hrl] at hw ⊢
  have := entryB_eq (Nat.log2 (4 * r.length * r.length * a * b) + 1) a b r (colN Y j)
    (by rw [length_colN, hrl]) hne hrb (colN_le Y j _ (fun r hr ↦ (hY r hr).2)) hw
  simp only [sh] at this ⊢
  rw [this]
  push_cast
  ring

theorem lmulB_shape (a b m : ℕ) (X Y : List (List ℕ)) (hk : Y ≠ [])
    (hX : ∀ r ∈ X, r.length = Y.length ∧ ∀ x ∈ r, x ≤ 2 * a)
    (hY : ∀ r ∈ Y, r.length = m ∧ ∀ y ∈ r, y ≤ 2 * b) :
    ∀ r ∈ lmulB a b X Y m, r.length = m ∧ ∀ x ∈ r, x ≤ 2 * (Y.length * a * b) := by
  intro r hr
  simp only [lmulB, List.mem_map] at hr
  obtain ⟨r0, hr0, rfl⟩ := hr
  obtain ⟨hrl, hrb⟩ := hX r0 hr0
  have hne : r0 ≠ [] := by
    intro h; rw [h] at hrl; exact hk (List.length_eq_zero_iff.mp hrl.symm)
  refine ⟨by simp, ?_⟩
  intro x hx
  simp only [List.map_map, List.mem_map, List.mem_range, Function.comp] at hx
  obtain ⟨j, _, rfl⟩ := hx
  have hw := radix_bound Y.length a b
  rw [← hrl] at hw ⊢
  exact entryB_le _ a b r0 (colN Y j) (by rw [length_colN, hrl]) hne hrb
    (colN_le Y j _ (fun r hr ↦ (hY r hr).2)) hw

theorem length_lmulB (a b m : ℕ) (X Y : List (List ℕ)) : (lmulB a b X Y m).length = X.length := by
  simp [lmulB]

/-- The `n × m` integer matrix stored as `X` with offset `β`. -/
def toMatB (n m β : ℕ) (X : List (List ℕ)) : Matrix (Fin n) (Fin m) ℤ := toMat n m (shiftL β X)

theorem toMatB_lmulB {n k m : ℕ} (a b : ℕ) (X Y : List (List ℕ)) (hYk : Y.length = k)
    (hk : Y ≠ []) (hX : ∀ r ∈ X, r.length = k ∧ ∀ x ∈ r, x ≤ 2 * a)
    (hY : ∀ r ∈ Y, r.length = m ∧ ∀ y ∈ r, y ≤ 2 * b) :
    toMatB n m (k * a * b) (lmulB a b X Y m) = toMatB n k a X * toMatB k m b Y := by
  subst hYk
  unfold toMatB
  rw [lmulB_spec a b m X Y hk hX hY]
  apply toMat_lmul
  · intro r hr
    simp only [shiftL, List.mem_map] at hr
    obtain ⟨r0, hr0, rfl⟩ := hr
    simpa using (hX r0 hr0).1
  · simp [shiftL]

theorem getD_getD_eq (L : List (List ℕ)) (i j : ℕ) (hi : i < L.length) (hj : j < L[i].length) :
    (L.getD i []).getD j 0 = L[i][j] := by
  rw [List.getD_eq_getElem _ _ hi, List.getD_eq_getElem _ _ hj]

theorem lget_shiftL (β : ℕ) (L : List (List ℕ)) (i j : ℕ) (hi : i < L.length)
    (hj : j < L[i].length) : lget (shiftL β L) i j = sh β L[i][j] := by
  unfold lget shiftL
  have h1 : (L.map (·.map (sh β))).getD i [] = L[i].map (sh β) := by
    rw [List.getD_eq_getElem _ _ (by simpa using hi), List.getElem_map]
  rw [h1, List.getD_eq_getElem _ _ (by simpa using hj), List.getElem_map]

/-- `∑ᵢ (X Y)ᵢᵢ` as stored with offset `Y.length * a * b`, from the diagonal entries only. -/
def ltraceB (a b : ℕ) (X Y : List (List ℕ)) (n : ℕ) : ℕ :=
  (List.ofFn fun i : Fin n ↦ ((lmulB a b X Y n).getD i []).getD i 0).sum

theorem trace_toMatB_mul {n k : ℕ} (a b : ℕ) (X Y : List (List ℕ)) (hXn : X.length = n)
    (hYk : Y.length = k) (hk : Y ≠ []) (hX : ∀ r ∈ X, r.length = k ∧ ∀ x ∈ r, x ≤ 2 * a)
    (hY : ∀ r ∈ Y, r.length = n ∧ ∀ y ∈ r, y ≤ 2 * b) :
    Matrix.trace (toMatB n k a X * toMatB k n b Y) =
      (ltraceB a b X Y n : ℤ) - n * (k * a * b) := by
  rw [← toMatB_lmulB a b X Y hYk hk hX hY, Matrix.trace, ltraceB, List.sum_ofFn]
  have hshape := lmulB_shape a b n X Y hk (by rw [hYk]; exact hX) hY
  have hlen := length_lmulB a b n X Y
  push_cast
  have hc : (n : ℤ) * (k * a * b) = ∑ _i : Fin n, ((k * a * b : ℕ) : ℤ) := by simp
  rw [hc, ← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro i _
  have hi : (i : ℕ) < (lmulB a b X Y n).length := by rw [hlen, hXn]; exact i.isLt
  have hrow := (hshape _ (List.getElem_mem hi)).1
  have hj : (i : ℕ) < ((lmulB a b X Y n)[(i : ℕ)]).length := by rw [hrow]; exact i.isLt
  change lget (shiftL _ _) i i = _
  rw [lget_shiftL _ _ _ _ hi hj, getD_getD_eq _ _ _ hi hj, sh]

/-! ### Radix uniqueness -/

/-- Digits of absolute value below the radix are determined by the number: a vanishing radix
expansion has all digits zero. -/
theorem radix_eq_zero (R : ℤ) (hR : 0 < R) : ∀ (n : ℕ) (d : Fin n → ℤ), (∀ j, |d j| < R) →
    ∑ j, d j * R ^ (j : ℕ) = 0 → ∀ j, d j = 0
  | 0, _, _, _ => fun j ↦ j.elim0
  | n + 1, d, hd, hs => by
    rw [Fin.sum_univ_succ] at hs
    simp only [Fin.val_zero, pow_zero, mul_one, Fin.val_succ, pow_succ] at hs
    set S := ∑ j : Fin n, d j.succ * R ^ (j : ℕ) with hS
    have hsum : ∑ j : Fin n, d j.succ * (R ^ (j : ℕ) * R) = R * S := by
      rw [hS, Finset.mul_sum]; congr 1; ext j; ring
    rw [hsum] at hs
    have h0 : d 0 = 0 := Int.eq_zero_of_abs_lt_dvd ⟨-S, by linarith⟩ (hd 0)
    have hS0 : S = 0 := by
      rw [h0, zero_add] at hs
      exact (mul_eq_zero.mp hs).resolve_left hR.ne'
    intro j
    refine Fin.cases h0 (fun j ↦ ?_) j
    exact radix_eq_zero R hR n (fun j ↦ d j.succ) (fun j ↦ hd _) hS0 j

/-! ### Decoding packed data

Generated data are stored as single `Nat` literals (lanes of fixed width), which the elaborator
reads instantly and the kernel splits with a few big shifts per entry. -/

/-- The `n` lanes of width `w` of `x`, lane zero first. -/
def decodeN (w : ℕ) : ℕ → ℕ → List ℕ
  | 0, _ => []
  | n + 1, x => x % 2 ^ w :: decodeN w n (x >>> w)

/-- Split a list into `k` consecutive rows of length `m`. -/
def chunks (m : ℕ) : ℕ → List ℕ → List (List ℕ)
  | 0, _ => []
  | k + 1, xs => xs.take m :: chunks m k (xs.drop m)

/-- The `n × n` row-major matrix of an offset code `(β, w, x)`, still offset by `β`. -/
def decodeSquare (n : ℕ) (code : ℕ × ℕ × ℕ) : List (List ℕ) :=
  chunks n n (decodeN code.2.1 (n * n) code.2.2)

/-- A sparse grouped list: total number of entries, the group sizes `(width, code)`, the
indices `(width, code)`, and the values as an offset code `(β, width, code)`. -/
structure SparseCode where
  total : ℕ
  counts : ℕ × ℕ
  idx : ℕ × ℕ
  vals : ℕ × ℕ × ℕ

/-- Regroup a flat list by group sizes. -/
def regroup {α : Type} : List ℕ → List α → List (List α)
  | [], _ => []
  | c :: cs, xs => xs.take c :: regroup cs (xs.drop c)

/-- The `g` groups of `(index, value)` pairs of a sparse code. -/
def decodeSparse (g : ℕ) (s : SparseCode) : List (List (ℕ × ℤ)) :=
  regroup (decodeN s.counts.1 g s.counts.2)
    (List.zip (decodeN s.idx.1 s.total s.idx.2)
      ((decodeN s.vals.2.1 s.total s.vals.2.2).map (sh s.vals.1)))

end SeedKernel.Trial
