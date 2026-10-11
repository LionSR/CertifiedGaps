/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import ThermalMPOGap.Transfer.Defs
import ThermalMPOGap.Transfer.PiTensor
import Mathlib.LinearAlgebra.Matrix.PosDef

/-!
# Transfer-matrix factorization of the Trotterized Heisenberg ring

For a ring of `L = 2 N` spin-one sites, the trace of `M` first-order Trotter slices
`Tr (trotterSlice (2 N) τ)^M` equals `Tr T^N` for a positive semidefinite matrix `T` that depends
only on `M` and `τ ≥ 0`, not on `N`.

## Proof outline

* Real form of a bond factor. With `C_none = 1` and `C_α = √τ ε_α` (`siteFactor`, real
  entries), `S^α ⊗ S^α = -ε_α ⊗ ε_α` gives for distinct sites `j ≠ k`
  `1 - τ Sⱼ · Sₖ = ∑_γ C_γ ⊗ C_γ` (`one_sub_smul_bondOp_exchange`), a sum of tensor products
  `piTensor`.
* One slice. Expanding the ordered product (even bonds left of odd bonds) over one label per
  bond writes the slice as `∑_{(b, c)} ⨂ᵢ wᵢ(b, c)` (`trotterSlice_eq_sum`), with local word
  `C_{b m} C_{c (m-1)}` at the even site `2 m` and `C_{b m} C_{c m}` at the odd site `2 m + 1`
  (`sliceSite_evenSite`, `sliceSite_oddSite`).
* `M` slices. Expanding again and taking the trace sitewise, the even site `2 m` contributes
  `F (p m) (q (m-1))` and the odd site `F (p m) (q m)`, with `F x y = Tr ∏ₜ C_{x t} C_{y t}`
  (`cellTransfer`) on the time-label space `Fin M → Option (Fin 3)`.
* Summing out the even-bond labels `p` leaves the closed-path sum of `Fᵀ F` around the ring of
  `N` unit cells, i.e. `Tr (Fᵀ F)^N` (`trace_trotterSlice_pow_eq_trace_transfer`). Since `F` has
  real entries, `Fᵀ F = Fᴴ F` is positive semidefinite.

## Main results

* `trace_trotterSlice_pow_eq_trace_transfer`: `Tr (trotterSlice (2 N) τ)^M = Tr (Fᵀ F)^N`.
* `exists_posSemidef_transfer`: the factorization through a positive semidefinite `T`.
-/

open Matrix ThermalMPOGap
open scoped Kronecker ComplexOrder

namespace SpinOneHeisenberg

/-! ### Local factors of a bond -/

/-- The real local factors of one Trotter bond factor: `C_none = 1` and
`C_α = √τ ε_α` with `(ε_α)_{s s'} = ε_{α s s'}`. -/
noncomputable def realSiteFactor (τ : ℝ) : Option (Fin 3) → Matrix (Fin 3) (Fin 3) ℝ
  | none => 1
  | some α => of fun s s' ↦ Real.sqrt τ * leviCivita α s s'

/-- The local factors of one Trotter bond factor, as complex matrices with real entries. -/
noncomputable def siteFactor (τ : ℝ) (γ : Option (Fin 3)) : Matrix (Fin 3) (Fin 3) ℂ :=
  (realSiteFactor τ γ).map (↑)

/-- The local factor with trivial label is the identity. -/
@[simp] theorem siteFactor_none (τ : ℝ) : siteFactor τ none = 1 := by
  simp [siteFactor, realSiteFactor]

/-- `C_α ⊗ C_α = -τ S^α ⊗ S^α`. -/
theorem kronecker_siteFactor_some {τ : ℝ} (hτ : 0 ≤ τ) (α : Fin 3) :
    kroneckerMap (· * ·) (siteFactor τ (some α)) (siteFactor τ (some α)) =
      (-(τ : ℂ)) • kroneckerMap (· * ·) (spinOne α) (spinOne α) := by
  ext ⟨a, b⟩ ⟨c, d⟩
  have h : (Real.sqrt τ : ℂ) * Real.sqrt τ = τ := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt hτ]
  simp only [kroneckerMap_apply, siteFactor, realSiteFactor, map_apply, of_apply, spinOne,
    Matrix.smul_apply, smul_eq_mul, Complex.ofReal_mul, Complex.ofReal_intCast]
  linear_combination (leviCivita α a c : ℂ) * leviCivita α b d * h +
    (τ : ℂ) * leviCivita α a c * leviCivita α b d * Complex.I_sq

/-- `bondOp` is additive in the two-site operator. -/
theorem bondOp_sum {L : ℕ} (j k : Fin L) {X : Type*} (s : Finset X)
    (B : X → Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ) :
    bondOp j k (∑ x ∈ s, B x) = ∑ x ∈ s, bondOp j k (B x) := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.sum_apply]
  split_ifs <;> simp

/-- `bondOp` is homogeneous in the two-site operator. -/
theorem bondOp_smul {L : ℕ} (j k : Fin L) (c : ℂ) (B : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ) :
    bondOp j k (c • B) = c • bondOp j k B := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.smul_apply, smul_eq_mul]
  split_ifs <;> simp

/-- A product operator `C ⊗ C` on two distinct sites is a tensor product over the chain. -/
theorem bondOp_kronecker_self {L : ℕ} {j k : Fin L} (hjk : j ≠ k) (C : Matrix (Fin 3) (Fin 3) ℂ) :
    bondOp j k (kroneckerMap (· * ·) C C) = piTensor fun i ↦ if i = j ∨ i = k then C else 1 := by
  ext σ τ
  have hk : k ∈ Finset.univ.erase j := Finset.mem_erase.2 ⟨hjk.symm, Finset.mem_univ _⟩
  rw [piTensor_apply, ← Finset.mul_prod_erase _ _ (Finset.mem_univ j),
    ← Finset.mul_prod_erase _ _ hk,
    Finset.prod_congr rfl (g := fun i ↦ if σ i = τ i then (1 : ℂ) else 0) fun i hi ↦ ?_,
    Finset.prod_ite_zero, Finset.prod_const_one]
  · simp only [bondOp, of_apply, kroneckerMap_apply, true_or, or_true, ite_true,
      Finset.mem_erase, Finset.mem_univ, and_true]
    split_ifs with h₁ h₂ h₂ <;> simp_all
  · simp only [Finset.mem_erase, Finset.mem_univ, and_true] at hi
    simp [hi.1, hi.2, one_apply]

/-- A first-order bond factor on two distinct sites as a sum of tensor products of the real
local factors: `1 - τ Sⱼ · Sₖ = ∑_γ C_γ ⊗ C_γ`. -/
theorem one_sub_smul_bondOp_exchange {L : ℕ} {j k : Fin L} (hjk : j ≠ k) {τ : ℝ} (hτ : 0 ≤ τ) :
    1 - (τ : ℂ) • bondOp j k exchange =
      ∑ γ, piTensor fun i ↦ if i = j ∨ i = k then siteFactor τ γ else 1 := by
  rw [Fintype.sum_option]
  simp only [siteFactor_none, ite_self, piTensor_one]
  rw [sub_eq_add_neg, exchange, bondOp_sum, Finset.smul_sum, ← Finset.sum_neg_distrib]
  congr 1
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [← bondOp_kronecker_self hjk, kronecker_siteFactor_some hτ, bondOp_smul, neg_smul]

/-! ### One Trotter slice as a sum of tensor products -/

/-- The cyclic successor in coordinates: `finRotate L j = j + 1`, wrapping `L - 1` to `0`. -/
theorem val_finRotate {L : ℕ} (j : Fin L) :
    (finRotate L j : ℕ) = if (j : ℕ) + 1 = L then 0 else (j : ℕ) + 1 := by
  cases L with
  | zero => exact j.elim0
  | succ n =>
    rw [coe_finRotate]
    split_ifs with h₁ h₂ h₂ <;> simp_all [Fin.ext_iff]

/-- On a ring of at least two sites, the two ends of a bond are distinct. -/
theorem ne_finRotate {L : ℕ} (hL : 2 ≤ L) (j : Fin L) : j ≠ finRotate L j := by
  intro h
  have h' := congrArg Fin.val h
  rw [val_finRotate] at h'
  split_ifs at h' <;> omega

/-- The tensor-product factor of label `γ` of the bond `j`: `C_γ` on the sites `j` and `j + 1`,
the identity elsewhere. -/
noncomputable def bondSite (L : ℕ) (τ : ℝ) (j : Fin L) (γ : Option (Fin 3)) :
    Fin L → Matrix (Fin 3) (Fin 3) ℂ :=
  fun i ↦ if i = j ∨ i = finRotate L j then siteFactor τ γ else 1

/-- A first-order bond factor as a sum over the bond label. -/
theorem one_sub_smul_bond {L : ℕ} (hL : 2 ≤ L) {τ : ℝ} (hτ : 0 ≤ τ) (j : Fin L) :
    1 - (τ : ℂ) • bond L j = ∑ γ, piTensor (bondSite L τ j γ) :=
  one_sub_smul_bondOp_exchange (ne_finRotate hL j) hτ

/-- The even site `2 k` of a ring of `2 N` sites. -/
def evenSite (N : ℕ) (k : Fin N) : Fin (2 * N) := ⟨2 * k, by omega⟩

/-- The odd site `2 k + 1` of a ring of `2 N` sites. -/
def oddSite (N : ℕ) (k : Fin N) : Fin (2 * N) := ⟨2 * k + 1, by omega⟩

/-- On a ring of even length the Trotter order lists the even bonds `2 k`, then the odd bonds
`2 k + 1`, each in increasing order. -/
theorem bondOrder_two_mul (N : ℕ) :
    bondOrder (2 * N) = List.ofFn (evenSite N) ++ List.ofFn (oddSite N) := by
  unfold bondOrder
  congr 1
  · refine List.SortedLT.eq_of_mem_iff ?_ ?_ fun j ↦ ?_
    · exact ((List.sortedLT_finRange _).pairwise.filter _).sortedLT
    · exact List.sortedLT_ofFn_iff.2 fun a b h ↦ by
        simp only [evenSite, Fin.lt_def] at h ⊢; omega
    · simp only [List.mem_filter, List.mem_finRange, true_and, decide_eq_true_eq, List.mem_ofFn]
      constructor
      · rintro ⟨r, hr⟩
        exact ⟨⟨r, by omega⟩, Fin.ext (by simp [evenSite]; omega)⟩
      · rintro ⟨k, rfl⟩
        exact ⟨k, by simp [evenSite]; omega⟩
  · refine List.SortedLT.eq_of_mem_iff ?_ ?_ fun j ↦ ?_
    · exact ((List.sortedLT_finRange _).pairwise.filter _).sortedLT
    · exact List.sortedLT_ofFn_iff.2 fun a b h ↦ by
        simp only [oddSite, Fin.lt_def] at h ⊢; omega
    · simp only [List.mem_filter, List.mem_finRange, true_and, Bool.not_eq_eq_eq_not,
        Bool.not_true, decide_eq_false_iff_not, Nat.not_even_iff_odd, List.mem_ofFn]
      constructor
      · rintro ⟨r, hr⟩
        exact ⟨⟨r, by omega⟩, Fin.ext (by simp [oddSite]; omega)⟩
      · rintro ⟨k, rfl⟩
        exact ⟨k, by simp [oddSite]⟩

/-- The local word at site `i` of one Trotter slice with even-bond labels `bc.1` and odd-bond
labels `bc.2`: the ordered product of the factors from the even bonds, then the odd bonds. -/
noncomputable def sliceSite (N : ℕ) (τ : ℝ)
    (bc : (Fin N → Option (Fin 3)) × (Fin N → Option (Fin 3)))
    (i : Fin (2 * N)) : Matrix (Fin 3) (Fin 3) ℂ :=
  (List.ofFn fun k ↦ bondSite (2 * N) τ (evenSite N k) (bc.1 k) i).prod *
    (List.ofFn fun k ↦ bondSite (2 * N) τ (oddSite N k) (bc.2 k) i).prod

/-- One Trotter slice on a ring of `2 N` sites is a sum of tensor products, one for each
labelling of the bonds. -/
theorem trotterSlice_eq_sum {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 ≤ τ) :
    trotterSlice (2 * N) τ = ∑ bc, piTensor (sliceSite N τ bc) := by
  rw [trotterSlice, bondOrder_two_mul, List.map_append, List.prod_append, List.map_ofFn,
    List.map_ofFn]
  simp only [Function.comp_def, one_sub_smul_bond (by omega : 2 ≤ 2 * N) hτ]
  rw [list_prod_ofFn_sum, list_prod_ofFn_sum, Finset.sum_mul_sum, ← Finset.sum_product']
  refine Finset.sum_congr rfl fun bc _ ↦ ?_
  rw [list_prod_ofFn_piTensor, list_prod_ofFn_piTensor, piTensor_mul]
  rfl

/-! ### The local words of one slice -/

/-- The even site `2 m` lies on the even bond `2 k` exactly when `k = m`. -/
theorem evenSite_mem_evenBond_iff {N : ℕ} (m k : Fin N) :
    evenSite N m = evenSite N k ∨ evenSite N m = finRotate (2 * N) (evenSite N k) ↔ k = m := by
  simp only [Fin.ext_iff, evenSite, val_finRotate]
  have := m.isLt
  have := k.isLt
  split_ifs <;> omega

/-- The odd site `2 m + 1` lies on the even bond `2 k` exactly when `k = m`. -/
theorem oddSite_mem_evenBond_iff {N : ℕ} (m k : Fin N) :
    oddSite N m = evenSite N k ∨ oddSite N m = finRotate (2 * N) (evenSite N k) ↔ k = m := by
  simp only [Fin.ext_iff, evenSite, oddSite, val_finRotate]
  have := m.isLt
  have := k.isLt
  split_ifs <;> omega

/-- The even site `2 m` lies on the odd bond `2 k + 1` exactly when `k = m - 1` cyclically. -/
theorem evenSite_mem_oddBond_iff {N : ℕ} (m k : Fin N) :
    evenSite N m = oddSite N k ∨ evenSite N m = finRotate (2 * N) (oddSite N k) ↔
      k = (finRotate N).symm m := by
  rw [Equiv.eq_symm_apply]
  simp only [Fin.ext_iff, evenSite, oddSite, val_finRotate]
  have := m.isLt
  have := k.isLt
  split_ifs <;> omega

/-- The odd site `2 m + 1` lies on the odd bond `2 k + 1` exactly when `k = m`. -/
theorem oddSite_mem_oddBond_iff {N : ℕ} (m k : Fin N) :
    oddSite N m = oddSite N k ∨ oddSite N m = finRotate (2 * N) (oddSite N k) ↔ k = m := by
  simp only [Fin.ext_iff, oddSite, val_finRotate]
  have := m.isLt
  have := k.isLt
  split <;> omega

/-- The word of a slice at the even site `2 m`: the even bond `2 m` acts first (leftmost), then
the odd bond `2 m - 1`. -/
theorem sliceSite_evenSite (N : ℕ) (τ : ℝ)
    (bc : (Fin N → Option (Fin 3)) × (Fin N → Option (Fin 3)))
    (m : Fin N) :
    sliceSite N τ bc (evenSite N m) =
      siteFactor τ (bc.1 m) * siteFactor τ (bc.2 ((finRotate N).symm m)) := by
  simp only [sliceSite, bondSite]
  rw [list_prod_ofFn_ite _ _ m (evenSite_mem_evenBond_iff m),
    list_prod_ofFn_ite _ _ _ (evenSite_mem_oddBond_iff m)]

/-- The word of a slice at the odd site `2 m + 1`: the even bond `2 m`, then the odd bond
`2 m + 1`. -/
theorem sliceSite_oddSite (N : ℕ) (τ : ℝ)
    (bc : (Fin N → Option (Fin 3)) × (Fin N → Option (Fin 3)))
    (m : Fin N) :
    sliceSite N τ bc (oddSite N m) = siteFactor τ (bc.1 m) * siteFactor τ (bc.2 m) := by
  simp only [sliceSite, bondSite]
  rw [list_prod_ofFn_ite _ _ m (oddSite_mem_evenBond_iff m),
    list_prod_ofFn_ite _ _ m (oddSite_mem_oddBond_iff m)]

/-- A product over the `2 N` sites of the ring, grouped into the unit cells `(2 m, 2 m + 1)`. -/
theorem prod_fin_two_mul {R : Type*} [CommMonoid R] (N : ℕ) (f : Fin (2 * N) → R) :
    ∏ i, f i = ∏ m, f (evenSite N m) * f (oddSite N m) := by
  rw [← (finProdFinEquiv.trans (finCongr (Nat.mul_comm N 2))).prod_comp, Fintype.prod_prod_type]
  refine Finset.prod_congr rfl fun m _ ↦ ?_
  rw [Fin.prod_univ_two]
  congr 2 <;> ext <;> simp [evenSite, oddSite]; ring

/-! ### Many slices and the transfer matrix -/

/-- The time-ordered word on one site over `M` slices, when the site receives the factors
`C_{p t}` and then `C_{q t}` in slice `t`. -/
noncomputable def timeWord (M : ℕ) (τ : ℝ) (p q : Fin M → Option (Fin 3)) :
    Matrix (Fin 3) (Fin 3) ℂ :=
  (List.ofFn fun t ↦ siteFactor τ (p t) * siteFactor τ (q t)).prod

/-- The one-site transfer matrix `F p q = Tr (timeWord p q)` on the time-label space
`Fin M → Option (Fin 3)`. The transfer matrix of a unit cell of two sites is `Fᵀ F`. -/
noncomputable def cellTransfer (M : ℕ) (τ : ℝ) :
    Matrix (Fin M → Option (Fin 3)) (Fin M → Option (Fin 3)) ℂ :=
  of fun p q ↦ trace (timeWord M τ p q)

/-- The trace of `M` Trotter slices as a sum over space-time bond labels of products of
one-site traces. -/
theorem trace_trotterSlice_pow_eq_sum {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 ≤ τ) (M : ℕ) :
    trace (trotterSlice (2 * N) τ ^ M) =
      ∑ a : Fin M → (Fin N → Option (Fin 3)) × (Fin N → Option (Fin 3)),
        ∏ m, cellTransfer M τ (fun t ↦ (a t).1 m) (fun t ↦ (a t).2 ((finRotate N).symm m)) *
          cellTransfer M τ (fun t ↦ (a t).1 m) (fun t ↦ (a t).2 m) := by
  rw [← List.prod_replicate, ← List.ofFn_const, trotterSlice_eq_sum hN hτ, list_prod_ofFn_sum,
    trace_sum]
  refine Finset.sum_congr rfl fun a _ ↦ ?_
  rw [list_prod_ofFn_piTensor, trace_piTensor, prod_fin_two_mul]
  simp only [sliceSite_evenSite, sliceSite_oddSite]
  rfl

/-- The trace of `M` Trotter slices on a ring of `2 N` sites is the trace of the `N`-th power
of the unit-cell transfer matrix `Fᵀ F`. -/
theorem trace_trotterSlice_pow_eq_trace_transfer {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 ≤ τ)
    (M : ℕ) :
    trace (trotterSlice (2 * N) τ ^ M) =
      trace (((cellTransfer M τ)ᵀ * cellTransfer M τ) ^ N) := by
  rw [trace_trotterSlice_pow_eq_sum hN hτ, trace_pow_eq_sum_closedPath _ hN]
  let F := cellTransfer M τ
  let e := (Equiv.arrowProdEquivProdArrow (Fin M) (fun _ ↦ Fin N → Option (Fin 3))
    (fun _ ↦ Fin N → Option (Fin 3))).trans
    (Equiv.prodCongr (Equiv.piComm _) (Equiv.piComm _))
  refine (Fintype.sum_equiv e _
    (fun pq ↦ ∏ m, F (pq.1 m) (pq.2 ((finRotate N).symm m)) * F (pq.1 m) (pq.2 m))
    fun a ↦ rfl).trans ?_
  rw [Fintype.sum_prod_type, Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  rw [← Fintype.prod_sum fun m x ↦ F x (q ((finRotate N).symm m)) * F x (q m),
    ← Equiv.prod_comp (finRotate N)]
  refine Finset.prod_congr rfl fun m _ ↦ ?_
  simp only [Equiv.symm_apply_apply, mul_apply, transpose_apply]
  rfl

/-- The one-site transfer matrix has real entries. -/
theorem star_cellTransfer_apply (M : ℕ) (τ : ℝ) (p q : Fin M → Option (Fin 3)) :
    star (cellTransfer M τ p q) = cellTransfer M τ p q := by
  have h : timeWord M τ p q = (List.ofFn fun t ↦
      realSiteFactor τ (p t) * realSiteFactor τ (q t)).prod.map (↑) := by
    change _ = (Complex.ofRealHom.mapMatrix : Matrix (Fin 3) (Fin 3) ℝ →+* _) _
    rw [map_list_prod, List.map_ofFn]
    simp only [timeWord, siteFactor, Function.comp_def, RingHom.mapMatrix_apply, Matrix.map_mul]
    rfl
  simp only [cellTransfer, of_apply, h, trace, diag_apply, map_apply, ← Complex.ofReal_sum]
  exact Complex.conj_ofReal _

/-- The conjugate transpose of the one-site transfer matrix is its transpose. -/
theorem conjTranspose_cellTransfer (M : ℕ) (τ : ℝ) :
    (cellTransfer M τ)ᴴ = (cellTransfer M τ)ᵀ := by
  ext p q
  exact star_cellTransfer_apply M τ q p

/-- **Transfer-matrix factorization of the Trotterized Heisenberg ring.** For every number of
slices `M` and every `τ ≥ 0` there is a positive semidefinite matrix `T`, independent of the
ring size, with `Tr (trotterSlice (2 N) τ)^M = Tr T^N` for every `N ≥ 1`. -/
theorem exists_posSemidef_transfer (M : ℕ) (τ : ℝ) (hτ : 0 ≤ τ) :
    ∃ (D : Type) (_ : Fintype D) (_ : DecidableEq D) (T : Matrix D D ℂ), T.PosSemidef ∧
      ∀ N : ℕ, 1 ≤ N → Matrix.trace (trotterSlice (2 * N) τ ^ M) = Matrix.trace (T ^ N) :=
  ⟨Fin M → Option (Fin 3), inferInstance, inferInstance,
    (cellTransfer M τ)ᴴ * cellTransfer M τ, posSemidef_conjTranspose_mul_self _,
    fun N hN ↦ by rw [conjTranspose_cellTransfer, trace_trotterSlice_pow_eq_trace_transfer hN hτ]⟩

end SpinOneHeisenberg
