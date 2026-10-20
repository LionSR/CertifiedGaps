/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialLinear
import Mathlib.LinearAlgebra.Matrix.Kronecker

/-!
# Sparse matrices, entry bounds and Kronecker products

Generic facts used to check the intertwiner identities `T W = W B`, `O W = W C`, `V W = c 1`
by evaluating both sides on one radix vector:

* `sparseEntry g j`: entry `j` of a sparse row `g` of `(index, value)` pairs, and the sums
  `sum_sparseEntry_mul` that evaluate a sparse matrix on a vector without enumerating zeros.
* `EntryBound`: entrywise absolute bounds, closed under sums and products.
* `kronecker_mulVec`: `(A ⊗ B) vec M = vec (A M Bᵀ)`.
* `chunksZ`, `vadd`, `vsumN`: flat vectors and their rows, with their `getD` lemmas.
-/

namespace SeedKernel.Trial

open Matrix
open scoped Kronecker

/-! ### Sparse rows -/

/-- Entry `j` of the sparse row `g` (a list of `(index, value)` pairs; repeated indices add). -/
def sparseEntry (g : List (ℕ × ℤ)) (j : ℕ) : ℤ := (g.map fun t ↦ if t.1 = j then t.2 else 0).sum

/-- `∑ₜ vₜ F(jₜ)` over the pairs of a sparse row. -/
def sparseDot (g : List (ℕ × ℤ)) (F : ℕ → ℤ) : ℤ := (g.map fun t ↦ t.2 * F t.1).sum

/-- A sparse row times a vector, without enumerating zeros. -/
theorem sum_sparseEntry_mul {N : ℕ} (F : ℕ → ℤ) :
    ∀ g : List (ℕ × ℤ), (∀ t ∈ g, t.1 < N) →
      ∑ j : Fin N, sparseEntry g j * F j = sparseDot g F
  | [], _ => by simp [sparseEntry, sparseDot]
  | t :: g, h => by
    have ih := sum_sparseEntry_mul F g (fun u hu ↦ h u (by simp [hu]))
    have ht : t.1 < N := h t (by simp)
    simp only [sparseEntry, sparseDot, List.map_cons, List.sum_cons, add_mul,
      Finset.sum_add_distrib] at ih ⊢
    rw [ih]
    congr 1
    rw [Finset.sum_eq_single ⟨t.1, ht⟩]
    · simp
    · intro j _ hj
      have : t.1 ≠ (j : ℕ) := fun e ↦ hj (Fin.ext e.symm)
      simp [this]
    · simp

theorem abs_sparseEntry_le (g : List (ℕ × ℤ)) (j : ℕ) :
    |sparseEntry g j| ≤ (g.map fun t ↦ |t.2|).sum := by
  induction g with
  | nil => simp [sparseEntry]
  | cons t g ih =>
    simp only [sparseEntry, List.map_cons, List.sum_cons] at ih ⊢
    refine (abs_add_le _ _).trans (add_le_add ?_ ih)
    split_ifs <;> simp

/-- Dropping zero values does not change a sparse dot product. -/
theorem sparseDot_filter (g : List (ℕ × ℤ)) (F : ℕ → ℤ) :
    sparseDot (g.filter fun t ↦ t.2 ≠ 0) F = sparseDot g F := by
  induction g with
  | nil => rfl
  | cons t g ih =>
    simp only [sparseDot] at ih ⊢
    by_cases h : t.2 = 0
    · rw [List.filter_cons_of_neg (by simp [h]), ih, List.map_cons, List.sum_cons, h,
        zero_mul, zero_add]
    · rw [List.filter_cons_of_pos (by simp [h]), List.map_cons, List.map_cons, List.sum_cons,
        List.sum_cons, ih]

theorem sparseDot_ofFn' {n : ℕ} (f : Fin n → ℤ) (F : ℕ → ℤ) :
    sparseDot (List.ofFn fun j : Fin n ↦ ((j : ℕ), f j)) F = ∑ j : Fin n, f j * F j := by
  simp [sparseDot, List.map_ofFn, List.sum_ofFn, Function.comp]

theorem sparseDot_ofFn {n : ℕ} (f : ℕ → ℤ) (F : ℕ → ℤ) :
    sparseDot (List.ofFn fun j : Fin n ↦ ((j : ℕ), f j)) F = ∑ j : Fin n, f j * F j := by
  simp [sparseDot, List.map_ofFn, List.sum_ofFn, Function.comp]

/-! ### Entry bounds -/

/-- Every entry of `M` has absolute value at most `K`. -/
def EntryBound {α β : Type*} (M : Matrix α β ℤ) (K : ℤ) : Prop := ∀ i j, |M i j| ≤ K

section Bounds

variable {α β γ δ : Type*}

theorem EntryBound.add {M N : Matrix α β ℤ} {a b : ℤ} (hM : EntryBound M a)
    (hN : EntryBound N b) : EntryBound (M + N) (a + b) :=
  fun i j ↦ (abs_add_le _ _).trans (add_le_add (hM i j) (hN i j))

theorem EntryBound.sub {M N : Matrix α β ℤ} {a b : ℤ} (hM : EntryBound M a)
    (hN : EntryBound N b) : EntryBound (M - N) (a + b) :=
  fun i j ↦ (abs_sub _ _).trans (add_le_add (hM i j) (hN i j))

theorem EntryBound.mul [Fintype β] {M : Matrix α β ℤ} {N : Matrix β γ ℤ} {a b : ℤ}
    (hM : EntryBound M a) (hN : EntryBound N b) :
    EntryBound (M * N) (Fintype.card β * (a * b)) := by
  intro i j
  rw [mul_apply]
  calc |∑ k, M i k * N k j| ≤ ∑ k, |M i k * N k j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _k : β, a * b := by
        gcongr with k
        rw [abs_mul]
        exact mul_le_mul (hM i k) (hN k j) (abs_nonneg _) ((abs_nonneg _).trans (hM i k))
    _ = Fintype.card β * (a * b) := by simp

theorem EntryBound.kronecker {M : Matrix α β ℤ} {N : Matrix γ δ ℤ} {a b : ℤ}
    (hM : EntryBound M a) (hN : EntryBound N b) : EntryBound (M ⊗ₖ N) (a * b) := by
  intro i j
  rw [kroneckerMap_apply, abs_mul]
  exact mul_le_mul (hM _ _) (hN _ _) (abs_nonneg _) ((abs_nonneg _).trans (hM i.1 j.1))

theorem EntryBound.mono {M : Matrix α β ℤ} {a b : ℤ} (hM : EntryBound M a) (h : a ≤ b) :
    EntryBound M b := fun i j ↦ (hM i j).trans h

end Bounds

/-! ### Kronecker products on vectors -/

/-- The vector of a square matrix, indexed by pairs. -/
def vecM {n : Type*} (M : Matrix n n ℤ) : n × n → ℤ := fun p ↦ M p.1 p.2

theorem kronecker_mulVec {n : Type*} [Fintype n] (A B M : Matrix n n ℤ) :
    (A ⊗ₖ B) *ᵥ vecM M = vecM (A * M * Bᵀ) := by
  ext ⟨a, b⟩
  simp only [vecM, mulVec, dotProduct, kroneckerMap_apply, Fintype.sum_prod_type, mul_apply,
    transpose_apply, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun d _ ↦ Finset.sum_congr rfl fun c _ ↦ ?_
  ring

theorem vecM_add {n : Type*} (M N : Matrix n n ℤ) : vecM (M + N) = vecM M + vecM N := rfl

theorem vecM_sub {n : Type*} (M N : Matrix n n ℤ) : vecM (M - N) = vecM M - vecM N := rfl

/-! ### Flat vectors

A vector on the doubled index is computed as one flat list (one pass over each list of groups),
then cut into rows; random access into long lists is what the kernel evaluates slowly. -/

/-- `k` rows of length `m` of a flat list. -/
def chunksZ (m : ℕ) : ℕ → List ℤ → List (List ℤ)
  | 0, _ => []
  | k + 1, xs => xs.take m :: chunksZ m k (xs.drop m)

theorem lget_chunksZ (m : ℕ) : ∀ (k : ℕ) (xs : List ℤ) (a c : ℕ), a < k → c < m →
    lget (chunksZ m k xs) a c = xs.getD (m * a + c) 0
  | 0, _, _, _, ha, _ => absurd ha (Nat.not_lt_zero _)
  | k + 1, xs, 0, c, _, hc => by
    simp [lget, chunksZ, List.getD_eq_getElem?_getD, hc]
  | k + 1, xs, a + 1, c, ha, hc => by
    have ih := lget_chunksZ m k (xs.drop m) a c (by omega) hc
    simp only [lget, chunksZ, List.getD_cons_succ] at ih ⊢
    rw [ih, List.getD_eq_getElem?_getD, List.getElem?_drop, ← List.getD_eq_getElem?_getD]
    congr 1
    ring

/-- Entrywise sum of two lists. -/
def vadd (u v : List ℤ) : List ℤ := List.zipWith (· + ·) u v

/-- Entrywise negation. -/
def vneg (u : List ℤ) : List ℤ := u.map (- ·)

/-- Entrywise sum of lists of length `n`. -/
def vsumN (n : ℕ) (ls : List (List ℤ)) : List ℤ := ls.foldr vadd (List.replicate n 0)

theorem length_vsumN (n : ℕ) : ∀ ls : List (List ℤ), (∀ l ∈ ls, l.length = n) →
    (vsumN n ls).length = n
  | [], _ => by simp [vsumN]
  | l :: ls, h => by
    have ih := length_vsumN n ls (fun u hu ↦ h u (by simp [hu]))
    simp only [vsumN, List.foldr_cons, vadd, List.length_zipWith] at ih ⊢
    rw [ih, h l (by simp), Nat.min_self]

theorem getD_vsumN (n : ℕ) : ∀ ls : List (List ℤ), (∀ l ∈ ls, l.length = n) → ∀ p < n,
    (vsumN n ls).getD p 0 = (ls.map fun l ↦ l.getD p 0).sum
  | [], _, p, hp => by simp [vsumN, List.getD_eq_getElem?_getD, hp]
  | l :: ls, h, p, hp => by
    have hl := h l (by simp)
    have hls := length_vsumN n ls (fun u hu ↦ h u (by simp [hu]))
    have ih := getD_vsumN n ls (fun u hu ↦ h u (by simp [hu])) p hp
    simp only [List.map_cons, List.sum_cons]
    rw [← ih]
    change (List.zipWith (· + ·) l (vsumN n ls)).getD p 0 = _
    rw [List.getD_eq_getElem _ _ (by simp [hl, hls, hp]), List.getElem_zipWith,
      List.getD_eq_getElem _ _ (by omega), List.getD_eq_getElem _ _ (by omega)]

theorem length_decodeN (w : ℕ) : ∀ (n x : ℕ), (decodeN w n x).length = n
  | 0, _ => rfl
  | n + 1, x => by simp [decodeN, length_decodeN w n]

theorem length_regroup {α : Type} :
    ∀ (cs : List ℕ) (xs : List α), (regroup cs xs).length = cs.length
  | [], _ => rfl
  | _ :: cs, xs => by simp [regroup, length_regroup cs]

theorem length_decodeSparse (g : ℕ) (s : SparseCode) : (decodeSparse g s).length = g := by
  simp [decodeSparse, length_regroup, length_decodeN]

/-- `n` rows of length `n`. -/
def Sq (n : ℕ) (X : List (List ℤ)) : Prop := X.length = n ∧ ∀ r ∈ X, r.length = n

theorem getD_vadd {n : ℕ} (u v : List ℤ) (hu : u.length = n) (hv : v.length = n) (p : ℕ)
    (hp : p < n) : (vadd u v).getD p 0 = u.getD p 0 + v.getD p 0 := by
  unfold vadd
  rw [List.getD_eq_getElem _ _ (by simp [hu, hv, hp]), List.getElem_zipWith,
    List.getD_eq_getElem _ _ (by omega), List.getD_eq_getElem _ _ (by omega)]

theorem length_vadd (u v : List ℤ) : (vadd u v).length = min u.length v.length := by
  simp [vadd]

theorem getD_vneg (u : List ℤ) (p : ℕ) : (vneg u).getD p 0 = -u.getD p 0 := by
  rcases Nat.lt_or_ge p u.length with h | h
  · rw [List.getD_eq_getElem _ _ (by simpa [vneg] using h), List.getD_eq_getElem _ _ h]
    simp [vneg]
  · rw [List.getD_eq_default _ _ (by simpa [vneg] using h), List.getD_eq_default _ _ h]
    simp

theorem length_vneg (u : List ℤ) : (vneg u).length = u.length := by simp [vneg]

theorem sq_chunksZ (m : ℕ) : ∀ (k : ℕ) (xs : List ℤ), m * k ≤ xs.length →
    (chunksZ m k xs).length = k ∧ ∀ r ∈ chunksZ m k xs, r.length = m
  | 0, _, _ => by simp [chunksZ]
  | k + 1, xs, h => by
    obtain ⟨h1, h2⟩ := sq_chunksZ m k (xs.drop m) (by simp; rw [Nat.mul_succ] at h; omega)
    refine ⟨by simp [chunksZ, h1], ?_⟩
    intro r hr
    simp only [chunksZ, List.mem_cons] at hr
    rcases hr with rfl | hr
    · simp; rw [Nat.mul_succ] at h; omega
    · exact h2 r hr

end SeedKernel.Trial
