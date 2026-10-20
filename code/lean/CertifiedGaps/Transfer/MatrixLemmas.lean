/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Matrix.Spectrum
import CertifiedGaps.Transfer.PiTensor

/-!
# Matrix lemmas for transfer matrices: marked traces, blocks, eigenvalues

Generic bookkeeping for transfer matrices that are block diagonal by a symmetry character.
A matrix `T` on a finite index set `D` with a key `key : D → κ` is *block diagonal* if
`T x y = 0` whenever `key x ≠ key y`. Its block `blockOf T key c` keeps the entries between
indices of key `c` and sets all other entries to zero. A closed path of nonzero weight stays in
one block, so a trace of a power weighted by a function of the key splits into blocks:

  `Tr (diag(φ ∘ key) T^(n+1)) = ∑_c φ(c) Tr (blockOf T key c)^(n+1)`.

If a relabeling `σ` of `D` preserves `T` and moves keys by an injective map `π`, the blocks of
`c` and `π c` are conjugate, so they have equal traces of powers.

## Main results

* `trace_diagonal_mul_pow`: `Tr (diag f T^(n+1))` as a sum over closed paths.
* `trace_diagonal_mul_pow_eq_sum_blockOf`: the block decomposition above.
* `blockOf_reindex`, `trace_blockOf_pow_of_reindex`: relabeled blocks.
* `trace_list_prod_ofFn_map_ofReal`: complexifying commutes with ordered products and traces.
* `trace_reindex_self`: the trace is invariant under reindexing.
* `trace_pow_eq_sum_eigenvalues_real`: the trace of a power of a real symmetric matrix is the
  power sum of its eigenvalues (the complex version is
  `CertifiedGaps.trace_pow_eq_sum_eigenvalues` in `Transfer/MomentForm.lean`).
* `hasDerivAt_sum_prod_smul`: the derivative at `0` of `∑_e (∏ⱼ c_{e,j}(u)) • P_e`.
-/

open Matrix

namespace CertifiedGaps

/-! ### Closed paths with a marked point -/

/-- The trace of `diag f T^(n+1)` as a sum over closed paths of length `n + 1`, weighted by `f`
at the starting point. -/
theorem trace_diagonal_mul_pow {R D : Type*} [CommSemiring R] [Fintype D] [DecidableEq D]
    (T : Matrix D D R) (f : D → R) (n : ℕ) :
    trace (diagonal f * T ^ (n + 1)) =
      ∑ q : Fin (n + 1) → D, f (q 0) * ∏ k, T (q k) (q (finRotate (n + 1) k)) := by
  rw [trace, sum_fin_succ_fun]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [diag_apply, diagonal_mul, pow_succ_apply_eq_sum_path, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  rw [Fin.snoc_eq_cons_rotate, Fin.cons_zero]

/-- The marked trace with a relabeled start: `∑ₓ f x (T^(n+1)) (g x) x` is the sum over closed
paths `q` of length `n + 1` whose first step leaves from `g (q 0)` instead of `q 0`. With `g = id`
this is `trace_diagonal_mul_pow`; a twist of a ring by a relabeling of the transfer index takes
this form. -/
theorem sum_mul_pow_apply_eq_sum_path {R D : Type*} [CommSemiring R] [Fintype D] [DecidableEq D]
    (T : Matrix D D R) (f : D → R) (g : D → D) (n : ℕ) :
    ∑ x, f x * (T ^ (n + 1)) (g x) x =
      ∑ q : Fin (n + 1) → D, f (q 0) *
        ∏ k, T (Function.update q 0 (g (q 0)) k) (q (finRotate (n + 1) k)) := by
  rw [sum_fin_succ_fun]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [pow_succ_apply_eq_sum_path, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  rw [Fin.snoc_eq_cons_rotate, Fin.cons_zero, Fin.update_cons_zero]

/-! ### Complexification -/

/-- The trace of an ordered product of complexified real matrices is the complexified trace. -/
theorem trace_list_prod_ofFn_map_ofReal {d : Type*} [Fintype d] [DecidableEq d] {m : ℕ}
    (F : Fin m → Matrix d d ℝ) :
    trace (List.ofFn fun t ↦ (F t).map Complex.ofReal).prod =
      ((trace (List.ofFn F).prod : ℝ) : ℂ) := by
  change trace (List.ofFn (Complex.ofRealHom.mapMatrix ∘ F)).prod = _
  rw [← List.map_ofFn, ← map_list_prod, RingHom.mapMatrix_apply]
  simp [trace, Complex.ofReal_sum]

/-! ### Reindexing -/

/-- The trace is invariant under reindexing. -/
theorem trace_reindex_self {R n : Type*} [AddCommMonoid R] [Fintype n] (e : n ≃ n)
    (A : Matrix n n R) : trace (reindex e e A) = trace A := by
  simp only [trace, diag_apply, reindex_apply, submatrix_apply]
  exact Equiv.sum_comp e.symm fun i ↦ A i i

/-! ### Blocks of a block-diagonal matrix -/

section Blocks

variable {R D κ : Type*} [DecidableEq κ]

/-- The block of `T` on the key `c`: the entries between indices of key `c`, zero elsewhere. -/
def blockOf [Zero R] (T : Matrix D D R) (key : D → κ) (c : κ) : Matrix D D R :=
  of fun x y ↦ if key x = c ∧ key y = c then T x y else 0

/-- A block of a Hermitian matrix is Hermitian. -/
theorem blockOf_isHermitian [AddMonoid R] [StarAddMonoid R] {T : Matrix D D R}
    (hT : T.IsHermitian) (key : D → κ) (c : κ) : (blockOf T key c).IsHermitian := by
  ext x y
  simp only [conjTranspose_apply, blockOf, of_apply, and_comm]
  split_ifs
  · exact hT.apply x y
  · exact star_zero R

variable [CommSemiring R]

/-- The weight of a closed path in one block. -/
theorem prod_blockOf {n : ℕ} (T : Matrix D D R) (key : D → κ) (c : κ) (q : Fin (n + 1) → D) :
    ∏ k, blockOf T key c (q k) (q (finRotate (n + 1) k)) =
      if ∀ k, key (q k) = c then ∏ k, T (q k) (q (finRotate (n + 1) k)) else 0 := by
  split_ifs with h
  · exact Finset.prod_congr rfl fun k _ ↦ by simp [blockOf, h]
  · push Not at h
    obtain ⟨k, hk⟩ := h
    exact Finset.prod_eq_zero (Finset.mem_univ k) (by simp [blockOf, hk])

omit [DecidableEq κ] in
/-- In a block-diagonal matrix, a closed path of nonzero weight stays in one block. -/
theorem key_eq_of_prod_ne_zero {n : ℕ} {T : Matrix D D R} {key : D → κ}
    (hT : ∀ x y, key x ≠ key y → T x y = 0) {q : Fin (n + 1) → D}
    (h : ∏ k, T (q k) (q (finRotate (n + 1) k)) ≠ 0) :
    ∀ k, key (q k) = key (q 0) := by
  have hstep : ∀ k, key (q k) = key (q (finRotate (n + 1) k)) := fun k ↦ by
    by_contra hk
    exact h (Finset.prod_eq_zero (Finset.mem_univ k) (hT _ _ hk))
  have hi : ∀ i (hi : i < n + 1), key (q ⟨i, hi⟩) = key (q 0) := by
    intro i
    induction i with
    | zero => intro _; rfl
    | succ i ih =>
      intro hi
      have hrot : finRotate (n + 1) ⟨i, by omega⟩ = ⟨i + 1, hi⟩ := by
        ext; rw [coe_finRotate]; split_ifs with h <;> simp [Fin.ext_iff] at h ⊢; omega
      rw [← ih (by omega), hstep ⟨i, by omega⟩, hrot]
  exact fun k ↦ hi k.val k.isLt

variable [Fintype D] [DecidableEq D]

/-- **Block decomposition.** For a block-diagonal `T` and a weight `φ` of the key,
`Tr (diag(φ ∘ key) T^(n+1)) = ∑_c φ(c) Tr (blockOf T key c)^(n+1)`. -/
theorem trace_diagonal_mul_pow_eq_sum_blockOf [Fintype κ] (T : Matrix D D R) (key : D → κ)
    (hT : ∀ x y, key x ≠ key y → T x y = 0) (φ : κ → R) (n : ℕ) :
    trace (diagonal (fun x ↦ φ (key x)) * T ^ (n + 1)) =
      ∑ c, φ c * trace (blockOf T key c ^ (n + 1)) := by
  simp only [trace_pow_eq_sum_closedPath _ (Nat.succ_pos n), trace_diagonal_mul_pow,
    prod_blockOf, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  by_cases hP : ∏ k, T (q k) (q (finRotate (n + 1) k)) = 0
  · simp only [hP, ite_self, mul_zero, Finset.sum_const_zero]
  · have hiff : ∀ c, (∀ k, key (q k) = c) ↔ key (q 0) = c := fun c ↦
      ⟨fun h ↦ h 0, fun h k ↦ (key_eq_of_prod_ne_zero hT hP k).trans h⟩
    simp only [hiff, mul_ite, mul_zero]
    rw [Finset.sum_ite_eq]
    simp only [Finset.mem_univ, ite_true]

omit [Fintype D] [DecidableEq D] in
/-- A relabeling `σ` that preserves `T` and moves keys by an injective `π` maps the block of `c`
to the block of `π c`. -/
theorem blockOf_reindex (T : Matrix D D R) (key : D → κ) (σ : D ≃ D) {π : κ → κ}
    (hπ : Function.Injective π) (hT : ∀ x y, T (σ x) (σ y) = T x y)
    (hkey : ∀ x, key (σ x) = π (key x)) (c : κ) :
    blockOf T key (π c) = reindex σ σ (blockOf T key c) := by
  ext x y
  obtain ⟨x, rfl⟩ := σ.surjective x
  obtain ⟨y, rfl⟩ := σ.surjective y
  simp only [blockOf, reindex_apply, submatrix_apply, of_apply, Equiv.symm_apply_apply, hkey,
    hπ.eq_iff, hT]

/-- Under the hypotheses of `blockOf_reindex`, the blocks of `c` and `π c` have equal traces of
powers. -/
theorem trace_blockOf_pow_of_reindex (T : Matrix D D R) (key : D → κ) (σ : D ≃ D) {π : κ → κ}
    (hπ : Function.Injective π) (hT : ∀ x y, T (σ x) (σ y) = T x y)
    (hkey : ∀ x, key (σ x) = π (key x)) (c : κ) (N : ℕ) :
    trace (blockOf T key (π c) ^ N) = trace (blockOf T key c ^ N) := by
  have h := map_pow (reindexAlgEquiv R R σ) (blockOf T key c) N
  rw [coe_reindexAlgEquiv] at h
  rw [blockOf_reindex T key σ hπ hT hkey, ← h, trace_reindex_self]

omit [DecidableEq κ] in
/-- Powers of a block-diagonal matrix are block diagonal. -/
theorem pow_apply_eq_zero_of_key_ne {T : Matrix D D R} {key : D → κ}
    (hT : ∀ x y, key x ≠ key y → T x y = 0) (N : ℕ) {x y : D} (h : key x ≠ key y) :
    (T ^ N) x y = 0 := by
  induction N generalizing y with
  | zero => rw [pow_zero, one_apply_ne fun e ↦ h (congrArg key e)]
  | succ N ih =>
    rw [pow_succ, mul_apply]
    refine Finset.sum_eq_zero fun z _ ↦ ?_
    by_cases hz : key x = key z
    · rw [hT z y (hz ▸ h), mul_zero]
    · rw [ih hz, zero_mul]

/-- The positive powers of a block of a block-diagonal matrix are the blocks of its powers. -/
theorem blockOf_pow_succ_apply {T : Matrix D D R} {key : D → κ}
    (hT : ∀ x y, key x ≠ key y → T x y = 0) (c : κ) (N : ℕ) (x y : D) :
    (blockOf T key c ^ (N + 1)) x y =
      if key x = c ∧ key y = c then (T ^ (N + 1)) x y else 0 := by
  induction N generalizing y with
  | zero => simp [blockOf]
  | succ N ih =>
    rw [pow_succ, mul_apply, pow_succ T (N + 1), mul_apply]
    simp only [ih]
    simp only [blockOf, of_apply]
    split_ifs with hxy
    · refine Finset.sum_congr rfl fun z _ ↦ ?_
      by_cases hz : key z = c
      · simp [hz, hxy.1, hxy.2]
      · rw [hT z y (by rw [hxy.2]; exact hz)]
        simp
    · refine Finset.sum_eq_zero fun z _ ↦ ?_
      split_ifs <;> simp_all

/-- **Twisted traces see the fixed blocks.** If a relabeling `g` moves keys by `κ'`
(`key (g x) = κ' (key x)`) and `c` is the only key fixed by `κ'`, then in the marked trace
`∑ₓ f x (T^(N+1)) (g x) x` of a block-diagonal `T` only the block of `c` contributes. -/
theorem sum_mul_pow_apply_eq_blockOf {T : Matrix D D R} {key : D → κ}
    (hT : ∀ x y, key x ≠ key y → T x y = 0) {g : D → D} {κ' : κ → κ}
    (hg : ∀ x, key (g x) = κ' (key x)) {c : κ} (hc : ∀ c', κ' c' = c' ↔ c' = c) (f : D → R)
    (N : ℕ) :
    ∑ x, f x * (T ^ (N + 1)) (g x) x = ∑ x, f x * (blockOf T key c ^ (N + 1)) (g x) x := by
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [blockOf_pow_succ_apply hT, hg]
  by_cases hx : key x = c
  · have h1 : κ' c = c := (hc c).2 rfl
    simp [h1, hx]
  · have h1 : key (g x) ≠ key x := by rw [hg]; exact fun h ↦ hx ((hc _).1 h)
    simp [hx, pow_apply_eq_zero_of_key_ne hT _ h1]

end Blocks

/-! ### Eigenvalues -/

/-- The trace of a power of a real symmetric matrix is the power sum of its eigenvalues. -/
theorem trace_pow_eq_sum_eigenvalues_real {D : Type*} [Fintype D] [DecidableEq D]
    {T : Matrix D D ℝ} (hT : T.IsHermitian) (N : ℕ) :
    trace (T ^ N) = ∑ i, hT.eigenvalues i ^ N := by
  conv_lhs => rw [hT.spectral_theorem, ← map_pow, diagonal_pow, Unitary.conjStarAlgAut_apply]
  rw [trace_mul_cycle, Unitary.coe_star_mul_self, Matrix.one_mul, trace_diagonal]
  simp

/-! ### Derivatives of polynomial families -/

/-- Derivative at `0` of a finite sum of products of scalar functions times fixed elements of a
normed algebra. -/
theorem hasDerivAt_sum_prod_smul {𝔸 E J : Type*} [NormedRing 𝔸] [NormedAlgebra ℂ 𝔸]
    [Fintype E] [Fintype J] [DecidableEq J] (c : E → J → ℂ → ℂ) (c' : E → J → ℂ)
    (hc : ∀ e j, HasDerivAt (c e j) (c' e j) 0) (P : E → 𝔸) :
    HasDerivAt (fun u ↦ ∑ e, (∏ j, c e j u) • P e)
      (∑ e, (∑ j, (∏ k ∈ Finset.univ.erase j, c e k 0) * c' e j) • P e) 0 := by
  refine HasDerivAt.fun_sum fun e _ ↦ ?_
  have := (HasDerivAt.fun_finsetProd (u := Finset.univ) fun j _ ↦ hc e j).smul_const (P e)
  simpa [smul_eq_mul] using this

end CertifiedGaps
