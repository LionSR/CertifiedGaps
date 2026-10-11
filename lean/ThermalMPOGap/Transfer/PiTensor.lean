/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Logic.Equiv.Fin.Rotate

/-!
# Tensor products of local matrices and ordered products

Generic tensor-network bookkeeping used for the transfer-matrix factorization of Trotterized
spin chains.

## Main definitions

* `piTensor A`: the tensor product `⨂ᵢ A i` of local matrices, as a matrix on configurations.

## Main results

* `list_prod_ofFn_sum`: an ordered product of sums is a sum of ordered products.
* `list_prod_ofFn_ite`: an ordered product of factors equal to `1` except at one position.
* `trace_pow_eq_sum_closedPath`: `Tr T ^ N` is the sum over closed paths of length `N`.
* `piTensor_one`, `piTensor_mul`, `list_prod_map_piTensor`, `list_prod_ofFn_piTensor`,
  `trace_piTensor`: the tensor
  product is multiplicative sitewise and its trace is the product of the local traces.
-/

open Matrix

namespace ThermalMPOGap

section Expand

/-- Splitting off the first entry of a tuple in a sum over tuples. -/
theorem sum_fin_succ_fun {R X : Type*} [AddCommMonoid R] [Fintype X] {n : ℕ}
    (g : (Fin (n + 1) → X) → R) :
    ∑ a, g a = ∑ x, ∑ a : Fin n → X, g (Fin.cons x a) := by
  rw [← (Fin.consEquiv fun _ ↦ X).sum_comp, Fintype.sum_prod_type]
  rfl

/-- An ordered product of sums expands into a sum of ordered products. -/
theorem list_prod_ofFn_sum {R X : Type*} [Semiring R] [Fintype X] :
    ∀ (n : ℕ) (f : Fin n → X → R),
      (List.ofFn fun t ↦ ∑ x, f t x).prod = ∑ a : Fin n → X, (List.ofFn fun t ↦ f t (a t)).prod
  | 0, f => by simp
  | n + 1, f => by
    rw [sum_fin_succ_fun, List.ofFn_succ, List.prod_cons, list_prod_ofFn_sum n,
      Finset.sum_mul_sum]
    simp [List.ofFn_succ]

/-- An ordered product of factors which are all `1` except at one position. -/
theorem list_prod_ofFn_eq_of_ne {M : Type*} [Monoid M] :
    ∀ (n : ℕ) (f : Fin n → M) (k₀ : Fin n), (∀ k, k ≠ k₀ → f k = 1) →
      (List.ofFn f).prod = f k₀
  | 0, _, k₀, _ => k₀.elim0
  | n + 1, f, k₀, h => by
    rw [List.ofFn_succ, List.prod_cons]
    induction k₀ using Fin.cases with
    | zero =>
      rw [List.prod_eq_one, mul_one]
      simp only [List.mem_ofFn, forall_exists_index]
      rintro _ k rfl
      exact h _ (Fin.succ_ne_zero k)
    | succ k₁ =>
      rw [h 0 (Fin.succ_ne_zero k₁).symm, one_mul,
        list_prod_ofFn_eq_of_ne n (fun k ↦ f k.succ) k₁ fun k hk ↦ h _ (by simpa using hk)]

/-- An ordered product of `if P k then A k else 1`, where `P` holds exactly at `k₀`. -/
theorem list_prod_ofFn_ite {M : Type*} [Monoid M] {n : ℕ} (P : Fin n → Prop) [DecidablePred P]
    (A : Fin n → M) (k₀ : Fin n) (hP : ∀ k, P k ↔ k = k₀) :
    (List.ofFn fun k ↦ if P k then A k else 1).prod = A k₀ := by
  rw [list_prod_ofFn_eq_of_ne n _ k₀ fun k hk ↦ by simp [(hP k).not.mpr hk]]
  simp [hP]

end Expand

section ClosedPath

variable {R D : Type*} [CommSemiring R] [Fintype D] [DecidableEq D]

/-- Matrix powers as sums over open paths. -/
theorem pow_succ_apply_eq_sum_path (T : Matrix D D R) :
    ∀ (n : ℕ) (x y : D), (T ^ (n + 1)) x y =
      ∑ q : Fin n → D, ∏ k : Fin (n + 1),
        T ((Fin.cons x q : Fin (n + 1) → D) k) ((Fin.snoc q y : Fin (n + 1) → D) k)
  | 0, x, y => by simp [Fin.snoc]
  | n + 1, x, y => by
    rw [pow_succ', mul_apply, sum_fin_succ_fun]
    refine Finset.sum_congr rfl fun z _ ↦ ?_
    rw [pow_succ_apply_eq_sum_path T n, Finset.mul_sum]
    refine Finset.sum_congr rfl fun q _ ↦ ?_
    rw [← Fin.cons_snoc_eq_snoc_cons, Fin.prod_univ_succ (n := n + 1)]
    simp only [Fin.cons_zero, Fin.cons_succ]

/-- The trace of `T ^ N` is the sum over closed paths of length `N` around a ring. -/
theorem trace_pow_eq_sum_closedPath (T : Matrix D D R) {N : ℕ} (hN : 1 ≤ N) :
    trace (T ^ N) = ∑ q : Fin N → D, ∏ k, T (q k) (q (finRotate N k)) := by
  obtain ⟨n, rfl⟩ : ∃ n, N = n + 1 := ⟨N - 1, by omega⟩
  rw [trace, sum_fin_succ_fun]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [diag_apply, pow_succ_apply_eq_sum_path]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  rw [Fin.snoc_eq_cons_rotate]

end ClosedPath

section PiTensor

variable {ι n R : Type*} [Fintype ι] [Fintype n] [CommSemiring R]

/-- The tensor product `⨂ᵢ A i` of a family of local matrices, as a matrix on configurations
`ι → n`: its entries are `∏ i, A i (σ i) (τ i)`. -/
def piTensor (A : ι → Matrix n n R) : Matrix (ι → n) (ι → n) R :=
  of fun σ τ ↦ ∏ i, A i (σ i) (τ i)

omit [Fintype n] in
/-- Entries of a tensor product. -/
theorem piTensor_apply (A : ι → Matrix n n R) (σ τ : ι → n) :
    piTensor A σ τ = ∏ i, A i (σ i) (τ i) := rfl

omit [Fintype n] in
/-- The tensor product of identities is the identity. -/
theorem piTensor_one [DecidableEq n] :
    piTensor (fun _ : ι ↦ (1 : Matrix n n R)) = 1 := by
  classical
  ext σ τ
  simp only [piTensor_apply, one_apply]
  rw [Finset.prod_ite_zero, Finset.prod_const_one]
  split_ifs with h₁ h₂ h₂ <;> simp_all [funext_iff]

/-- Tensor products multiply sitewise. -/
theorem piTensor_mul [DecidableEq ι] (A B : ι → Matrix n n R) :
    piTensor A * piTensor B = piTensor fun i ↦ A i * B i := by
  ext σ τ
  simp only [mul_apply, piTensor_apply, ← Finset.prod_mul_distrib]
  exact (Fintype.prod_sum fun i r ↦ A i (σ i) r * B i r (τ i)).symm

/-- An ordered product of tensor products is the tensor product of the sitewise ordered
products. -/
theorem list_prod_map_piTensor [DecidableEq ι] [DecidableEq n] {X : Type*}
    (A : X → ι → Matrix n n R) :
    ∀ l : List X, (l.map fun x ↦ piTensor (A x)).prod = piTensor fun i ↦ (l.map (A · i)).prod
  | [] => by simp [piTensor_one]
  | x :: l => by simp [list_prod_map_piTensor A l, piTensor_mul]

/-- An ordered product of tensor products, indexed by `Fin m`, is the tensor product of the
sitewise ordered products. -/
theorem list_prod_ofFn_piTensor [DecidableEq ι] [DecidableEq n] {m : ℕ}
    (A : Fin m → ι → Matrix n n R) :
    (List.ofFn fun k ↦ piTensor (A k)).prod = piTensor fun i ↦ (List.ofFn fun k ↦ A k i).prod := by
  simpa [List.ofFn_eq_map] using list_prod_map_piTensor A (List.finRange m)

/-- The trace of a tensor product is the product of the traces. -/
theorem trace_piTensor [DecidableEq ι] (A : ι → Matrix n n R) :
    trace (piTensor A) = ∏ i, trace (A i) := by
  simp only [trace, diag_apply, piTensor_apply]
  exact (Fintype.prod_sum fun i s ↦ A i s s).symm

end PiTensor

end ThermalMPOGap
