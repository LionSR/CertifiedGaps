/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.Model
import CertifiedGaps.Ring.Transfer
import CertifiedGaps.Transfer.ProductLimit

/-!
# Transfer-matrix factorization of Trotterized open chains

Number the bonds of the open chain on the `n + 1` sites `0, …, n` by `k = 0, …, n + 1`, bond `k`
joining the sites `k - 1` and `k`: the bonds `0` and `n + 1` are the two end terms (the *caps*),
the others the bonds of the chain. One first-order Trotter slice of width `τ` is the ordered
product (odd bonds)(on-site factors `1 - τ oᵢ`)(even bonds), mirroring the ring's `trotterSlice`.
For a ring without internal reflection (`R.refl = 1`) and `S = Fin M → Option ι` the time labels
of one bond over `M` slices, with the one-cell transfer matrix `F = R.cellTransfer M τ`:

  `Tr (slice)^M = ⟨v, Xⁿ⁺¹ v⟩`,   `X = F J`,   `Tr (ring slice)^M = Tr X^{2N}` on `2 N` cells,

where `J` reverses time and `v` is the boundary vector of the end term. `X` is Hermitian, so the
open-chain and ring Trotter traces are the boundary and bulk moments of one Hermitian matrix.

## Proof outline

* A bond factor is `1 - τ h = ∑_γ C_γ ⊗ C_γ` with `C_none = 1`, `C_α = √τ C_α`, and a cap factor
  is `1 - τ b = ∑_γ ω(γ) C_γ` at one site with the scalar weights `ω(none) = 1`,
  `ω(α) = √τ w_α` (`one_sub_smul_openBondTerm`).
* The two stages are written as products over *all* bonds, a bond of the other stage acting at
  width `0`, where its expansion has only the label `none`. Reindexing the two stage labels by
  one label per bond expands the slice (`openSlice_eq_sum`): site `i` carries the word
  `C_p O C_q`, with `p` the label of its odd bond and `q` that of its even bond.
* Over `M` slices the site `i` contributes `F` (odd `i`) or `Fᵀ` (even `i`) and each cap the
  vector `c p = ∏ₜ ω(p t)`. An open-path expansion gives
  `Tr (slice)^M = c_Lᵀ Fᵀ F Fᵀ ⋯ c_R` (`trace_openSlice_pow_eq_dotProduct`).
* Every `C_γ` is symmetric up to the sign `ε_γ` and `O` is symmetric, so transposing the word
  gives `Fᵀ = J F J` with `J p p' = ε(p) [p' = p ∘ rev]` (`transpose_cellTransfer`). Hence
  `X = F J` is real symmetric, the alternating product is `J Xⁿ⁺¹ J` or `J Xⁿ⁺¹`, and `J c = c̄`.

## Main results

* `RPRing.openTransfer_isHermitian`: `X = F J` is Hermitian.
* `RPRing.trace_trotterSlice_pow_eq_trace_openTransfer`: `Tr (ring slice)^M = Tr X^{2N}`.
* `RPRing.trace_openSlice_pow`: `Tr (slice)^M = ⟨v, Xⁿ⁺¹ v⟩`.
* `RPRing.tendsto_openSliceTrace`: the Trotter traces converge to `Re Tr e^{-βH}`.
-/

open Filter Matrix
open scoped Kronecker

namespace CertifiedGaps

/-- The weights of an end term over one slice of width `τ`: `ω(none) = 1`, `ω(α) = √τ u_α`. -/
noncomputable def capWeight {ι : Type} (τ : ℝ) (u : ι → ℂ) : Option ι → ℂ
  | none => 1
  | some α => Real.sqrt τ * u α

/-- The boundary vector of an end term over `M` slices: `c p = ∏ₜ ω(p t)`. -/
noncomputable def capVector {ι : Type} (M : ℕ) (τ : ℝ) (u : ι → ℂ) : (Fin M → Option ι) → ℂ :=
  fun p ↦ ∏ t, capWeight τ u (p t)

/-- The width of bond `k` in the first stage of a slice: the odd bonds act there. -/
def firstWidth (τ : ℝ) (k : ℕ) : ℝ := if k % 2 = 1 then τ else 0

/-- The width of bond `k` in the last stage of a slice: the even bonds act there. -/
def lastWidth (τ : ℝ) (k : ℕ) : ℝ := if k % 2 = 1 then 0 else τ

/-- The odd one of the two bonds `i, i + 1` at the site `i` of a chain of `n + 1` sites. -/
def oddBond (n : ℕ) (i : Fin (n + 1)) : Fin (n + 2) :=
  if (i : ℕ) % 2 = 1 then i.castSucc else i.succ

/-- The even one of the two bonds `i, i + 1` at the site `i` of a chain of `n + 1` sites. -/
def evenBond (n : ℕ) (i : Fin (n + 1)) : Fin (n + 2) :=
  if (i : ℕ) % 2 = 1 then i.succ else i.castSucc

namespace RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-! ### Signs of the bond factors -/

namespace Boundary

variable {R}

open Classical in
/-- The sign `ε_α = ±1` of a bond factor: `(C_α)ᵀ = ε_α C_α` and `w̄_α = ε_α w_α`. When `C_α` is
both symmetric and antisymmetric (that is, zero) the weight decides. -/
noncomputable def sign (b : R.Boundary) (α : ι) : ℝ :=
  if (R.C α)ᵀ = R.C α ∧ star (b.w α) = b.w α then 1 else -1

/-- The two defining properties of the sign. -/
theorem sign_spec (b : R.Boundary) (α : ι) :
    (R.C α)ᵀ = b.sign α • R.C α ∧ star (b.w α) = (b.sign α : ℂ) * b.w α := by
  unfold sign
  split_ifs with h
  · simpa using h
  · rcases b.w_symm α with h' | ⟨hC, hw⟩
    · exact (h h').elim
    · simp [hC, hw]

/-- The sign squares to one. -/
theorem sign_mul_self (b : R.Boundary) (α : ι) : b.sign α * b.sign α = 1 := by
  unfold sign; split_ifs <;> norm_num

/-- The sign of a bond label, with `ε(none) = 1`. -/
noncomputable def labelSign (b : R.Boundary) : Option ι → ℝ
  | none => 1
  | some α => b.sign α

/-- The label sign squares to one. -/
theorem labelSign_mul_self (b : R.Boundary) (γ : Option ι) :
    b.labelSign γ * b.labelSign γ = 1 := by
  cases γ with
  | none => simp [labelSign]
  | some α => exact b.sign_mul_self α

/-- The weights of the end term at the last site of a chain of `n` sites: `w` if `n` is odd,
`ε w` if `n` is even. -/
noncomputable def endW (b : R.Boundary) (n : ℕ) : ι → ℂ :=
  if Even n then fun α ↦ b.sign α * b.w α else b.w

/-- The transposed end term `bᵀ = -∑_α ε_α w_α C_α`, so both end terms have the form
`-∑_α u_α C_α`. -/
theorem endMat_eq (b : R.Boundary) (n : ℕ) : b.endMat n = -∑ α, b.endW n α • R.Cc α := by
  unfold endMat endW
  split_ifs
  · rw [mat, transpose_neg, transpose_sum]
    congr 1
    refine Finset.sum_congr rfl fun α _ ↦ ?_
    ext a c
    have h := congrArg (fun A ↦ A c a) (b.sign_spec α).1
    simp only [transpose_apply, Matrix.smul_apply, smul_eq_mul] at h
    have hs : (b.sign α : ℂ) * b.sign α = 1 := by exact_mod_cast b.sign_mul_self α
    simp only [Cc, map_apply, transpose_apply, Matrix.smul_apply, smul_eq_mul, h]
    push_cast
    linear_combination (-(b.w α * R.C α c a : ℂ)) * hs
  · rfl

end Boundary

/-! ### Local factors -/

/-- Without an internal reflection the right factors are the left factors. -/
theorem rightFactor_eq_of_refl (hrefl : R.refl = 1) (τ : ℝ) (γ : Option ι) :
    R.rightFactor τ γ = R.siteFactor τ γ := by
  have h : R.Rc = 1 := by simp [Rc, hrefl, Matrix.map_one]
  rw [rightFactor_eq, h, Matrix.one_mul, Matrix.mul_one]

/-- The scaled bond factor `C_α = √τ C_α` as a complex matrix. -/
theorem siteFactor_some (τ : ℝ) (α : ι) :
    R.siteFactor τ (some α) = (Real.sqrt τ : ℂ) • R.Cc α := by
  ext a c; simp [siteFactor, realSiteFactor, Cc]

/-- A cap factor as a sum over labels: `1 - τ (-∑_α u_α C_α) = ∑_γ ω(γ) C_γ`. -/
theorem one_sub_smul_cap {τ : ℝ} (hτ : 0 ≤ τ) (u : ι → ℂ) :
    1 - (τ : ℂ) • (-∑ α, u α • R.Cc α) = ∑ γ, capWeight τ u γ • R.siteFactor τ γ := by
  have h : (Real.sqrt τ : ℂ) * Real.sqrt τ = τ := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt hτ]
  rw [Fintype.sum_option, smul_neg, sub_neg_eq_add, Finset.smul_sum]
  simp only [capWeight, siteFactor_none, one_smul, siteFactor_some, smul_smul]
  congr 1
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [mul_right_comm, h, mul_comm]

/-- Each scaled bond factor is symmetric up to its sign. -/
theorem transpose_siteFactor (b : R.Boundary) (τ : ℝ) (γ : Option ι) :
    (R.siteFactor τ γ)ᵀ = (b.labelSign γ : ℂ) • R.siteFactor τ γ := by
  cases γ with
  | none => simp [Boundary.labelSign]
  | some α =>
    ext a c
    have h := congrArg (fun A ↦ A c a) (b.sign_spec α).1
    simp only [transpose_apply, Matrix.smul_apply, smul_eq_mul] at h
    have hs : (b.sign α : ℂ) * b.sign α = 1 := by exact_mod_cast b.sign_mul_self α
    simp [siteFactor, realSiteFactor, Boundary.labelSign, h]
    linear_combination (-(Real.sqrt τ * R.C α c a : ℂ)) * hs

/-- The on-site factor is symmetric. -/
theorem transpose_onsite (τ : ℝ) : (R.onsite τ)ᵀ = R.onsite τ := by
  ext a c
  simp only [transpose_apply, onsite, realOnsite, map_apply, Matrix.sub_apply, Matrix.smul_apply,
    smul_eq_mul]
  rw [← transpose_apply R.o, R.o_symm, one_apply, one_apply]
  simp [eq_comm]

/-! ### Bonds of the open chain -/

/-- The term of bond `k` of the chain on the sites `0, …, n`: the end term `b` at `k = 0`, the
end term `b.endMat (n + 1)` at `k = n + 1`, and the bond `(k - 1, k)` otherwise. -/
noncomputable def openBondTerm (b : R.Boundary) (n : ℕ) (k : Fin (n + 2)) :
    Matrix (Fin (n + 1) → d) (Fin (n + 1) → d) ℂ :=
  if (k : ℕ) = 0 then siteOp 0 b.mat
  else if h : (k : ℕ) = n + 1 then siteOp (Fin.last n) (b.endMat (n + 1))
  else bondOp ⟨k - 1, by omega⟩ ⟨k, by have := k.isLt; omega⟩ R.bondMatrix

/-- The scalar weight of label `γ` on bond `k`: the cap weights on the two caps, `1` on the
bonds of the chain. -/
noncomputable def openBondWeight (b : R.Boundary) (n : ℕ) (τ : ℝ) (k : Fin (n + 2))
    (γ : Option ι) : ℂ :=
  if (k : ℕ) = 0 then capWeight τ b.w γ
  else if (k : ℕ) = n + 1 then capWeight τ (b.endW (n + 1)) γ else 1

/-- The tensor factor of label `γ` on bond `k`: `C_γ` on the sites `k - 1` and `k` that the bond
touches. -/
noncomputable def openBondFactor (n : ℕ) (τ : ℝ) (k : Fin (n + 2)) (γ : Option ι) :
    Fin (n + 1) → Matrix d d ℂ :=
  fun i ↦ if (k : ℕ) = i ∨ (k : ℕ) = i + 1 then R.siteFactor τ γ else 1

/-- The label `none` has weight `1`. -/
@[simp] theorem openBondWeight_none (b : R.Boundary) (n : ℕ) (τ : ℝ) (k : Fin (n + 2)) :
    R.openBondWeight b n τ k none = 1 := by
  unfold openBondWeight; split_ifs <;> rfl

/-- The label `none` carries the identity. -/
@[simp] theorem openBondFactor_none (n : ℕ) (τ : ℝ) (k : Fin (n + 2)) :
    R.openBondFactor n τ k none = fun _ ↦ 1 := by
  funext i; simp [openBondFactor]

/-- At width `0` a label other than `none` carries the zero matrix. -/
theorem piTensor_openBondFactor_zero (n : ℕ) (k : Fin (n + 2)) (α : ι) :
    piTensor (R.openBondFactor n 0 k (some α)) = 0 := by
  refine piTensor_eq_zero (i₀ := ⟨min (k : ℕ) n, by omega⟩) ?_
  have : (k : ℕ) = min (k : ℕ) n ∨ (k : ℕ) = min (k : ℕ) n + 1 := by have := k.isLt; omega
  simp only [openBondFactor]
  split_ifs
  simp [siteFactor_some]

/-- A first-order bond factor as a sum over the bond label: `1 - τ h_k = ∑_γ ω_k(γ) C_γ ⊗ C_γ`
on the sites of bond `k`. -/
theorem one_sub_smul_openBondTerm (hrefl : R.refl = 1) (b : R.Boundary) (n : ℕ) {τ : ℝ}
    (hτ : 0 ≤ τ) (k : Fin (n + 2)) :
    1 - (τ : ℂ) • R.openBondTerm b n k =
      ∑ γ, R.openBondWeight b n τ k γ • piTensor (R.openBondFactor n τ k γ) := by
  unfold openBondTerm openBondWeight
  split_ifs with h0 hl
  · have hf : ∀ γ, R.openBondFactor n τ k γ = fun i ↦ if i = 0 then R.siteFactor τ γ else 1 :=
      fun γ ↦ funext fun i ↦ by
        simp only [openBondFactor, h0, Fin.ext_iff, Fin.val_zero]
        congr 1
        apply propext; omega
    simp only [hf]
    rw [one_sub_smul_siteOp, sum_smul_piTensor_single, Boundary.mat, R.one_sub_smul_cap hτ]
  · have hf : ∀ γ, R.openBondFactor n τ k γ =
        fun i ↦ if i = Fin.last n then R.siteFactor τ γ else 1 :=
      fun γ ↦ funext fun i ↦ by
        have := i.isLt
        simp only [openBondFactor, hl, Fin.ext_iff, Fin.val_last]
        congr 1
        apply propext; omega
    simp only [hf]
    rw [one_sub_smul_siteOp, sum_smul_piTensor_single, Boundary.endMat_eq, R.one_sub_smul_cap hτ]
  · have hlr : (⟨k - 1, by omega⟩ : Fin (n + 1)) ≠ ⟨k, by have := k.isLt; omega⟩ := by
      simp only [ne_eq, Fin.mk.injEq]; omega
    have hf : ∀ γ, R.openBondFactor n τ k γ = fun i ↦
        if i = ⟨k - 1, by omega⟩ then R.siteFactor τ γ
        else if i = ⟨k, by have := k.isLt; omega⟩ then R.siteFactor τ γ else 1 :=
      fun γ ↦ funext fun i ↦ by
        simp only [openBondFactor, Fin.ext_iff]
        split_ifs <;> first | rfl | omega
    simp only [one_smul, hf]
    rw [Fintype.sum_option, bondMatrix, bondOp_neg, bondOp_sum, smul_neg, sub_neg_eq_add,
      Finset.smul_sum]
    simp only [siteFactor_none, ite_self, piTensor_one]
    congr 1
    refine Finset.sum_congr rfl fun α _ ↦ ?_
    have hk := kronecker_siteFactor_some R hτ α
    rw [R.rightFactor_eq_of_refl hrefl, R.Dc_eq_of_refl hrefl] at hk
    rw [← bondOp_kronecker hlr, hk, bondOp_smul, R.Dc_eq_of_refl hrefl]

/-! ### Trotter slices of the open chain -/

/-- One stage of a slice: the factors `1 - w_k h_k` of all bonds in increasing order. A bond of
the other stage has width `0` and acts as the identity; keeping it makes every stage a product
over all bonds, which expands over one label per bond. -/
noncomputable def openStage (b : R.Boundary) (n : ℕ) (w : ℕ → ℝ) :
    Matrix (Fin (n + 1) → d) (Fin (n + 1) → d) ℂ :=
  (List.ofFn fun k : Fin (n + 2) ↦ 1 - (w k : ℂ) • R.openBondTerm b n k).prod

/-- One first-order Trotter slice of the open chain on `n + 1` sites: odd bonds, then all
on-site terms, then even bonds. -/
noncomputable def openSlice (b : R.Boundary) (n : ℕ) (τ : ℝ) :
    Matrix (Fin (n + 1) → d) (Fin (n + 1) → d) ℂ :=
  R.openStage b n (firstWidth τ) * (List.ofFn fun i ↦ 1 - (τ : ℂ) • R.site (n + 1) i).prod *
    R.openStage b n (lastWidth τ)

/-- A stage as a sum over one label per bond. -/
theorem openStage_eq_sum (hrefl : R.refl = 1) (b : R.Boundary) (n : ℕ) (w : ℕ → ℝ)
    (hw : ∀ k, 0 ≤ w k) :
    R.openStage b n w = ∑ ℓ : Fin (n + 2) → Option ι,
      (List.ofFn fun k ↦ R.openBondWeight b n (w k) k (ℓ k) •
        piTensor (R.openBondFactor n (w k) k (ℓ k))).prod := by
  unfold openStage
  simp only [R.one_sub_smul_openBondTerm hrefl b n (hw _)]
  exact list_prod_ofFn_sum _ _

/-- One term of a stage as a scaled tensor product. -/
theorem list_prod_openBond_smul (b : R.Boundary) (n : ℕ) (w : ℕ → ℝ)
    (ℓ : Fin (n + 2) → Option ι) :
    (List.ofFn fun k ↦ R.openBondWeight b n (w k) k (ℓ k) •
        piTensor (R.openBondFactor n (w k) k (ℓ k))).prod =
      (∏ k : Fin (n + 2), R.openBondWeight b n (w k) k (ℓ k)) •
        piTensor fun i ↦ (List.ofFn fun k ↦ R.openBondFactor n (w k) k (ℓ k) i).prod := by
  rw [list_prod_ofFn_smul, list_prod_ofFn_piTensor]

/-- The first stage at site `i`, with labels on the odd bonds only, is the factor of the odd
bond of `i`. -/
theorem list_prod_openBondFactor_first (n : ℕ) (τ : ℝ) (ℓ : Fin (n + 2) → Option ι)
    (i : Fin (n + 1)) :
    (List.ofFn fun k ↦ R.openBondFactor n (firstWidth τ k) k
        (if (k : ℕ) % 2 = 1 then ℓ k else none) i).prod = R.siteFactor τ (ℓ (oddBond n i)) := by
  rw [list_prod_ofFn_eq_of_ne _ _ (oddBond n i) fun k hk ↦ ?_]
  · rcases Nat.mod_two_eq_zero_or_one (i : ℕ) with hi | hi
    · have h1 : ((i : ℕ) + 1) % 2 = 1 := by omega
      have ho : oddBond n i = i.succ := by simp [oddBond, hi]
      rw [ho]
      simp [openBondFactor, firstWidth, h1]
    · have ho : oddBond n i = i.castSucc := by simp [oddBond, hi]
      rw [ho]
      simp [openBondFactor, firstWidth, hi]
  · by_cases hc : (k : ℕ) = i ∨ (k : ℕ) = i + 1
    · have hk2 : ¬ (k : ℕ) % 2 = 1 := fun h ↦ hk (by
        apply Fin.ext; unfold oddBond; split_ifs <;> simp <;> omega)
      simp [openBondFactor, hk2]
    · simp [openBondFactor, hc]

/-- The last stage at site `i`, with labels on the even bonds only, is the factor of the even
bond of `i`. -/
theorem list_prod_openBondFactor_last (n : ℕ) (τ : ℝ) (ℓ : Fin (n + 2) → Option ι)
    (i : Fin (n + 1)) :
    (List.ofFn fun k ↦ R.openBondFactor n (lastWidth τ k) k
        (if (k : ℕ) % 2 = 1 then none else ℓ k) i).prod = R.siteFactor τ (ℓ (evenBond n i)) := by
  rw [list_prod_ofFn_eq_of_ne _ _ (evenBond n i) fun k hk ↦ ?_]
  · rcases Nat.mod_two_eq_zero_or_one (i : ℕ) with hi | hi
    · have he : evenBond n i = i.castSucc := by simp [evenBond, hi]
      rw [he]
      simp [openBondFactor, lastWidth, hi]
    · have h1 : ((i : ℕ) + 1) % 2 = 0 := by omega
      have he : evenBond n i = i.succ := by simp [evenBond, hi]
      rw [he]
      simp [openBondFactor, lastWidth, h1]
  · by_cases hc : (k : ℕ) = i ∨ (k : ℕ) = i + 1
    · have hk2 : (k : ℕ) % 2 = 1 := by
        by_contra h
        exact hk (by apply Fin.ext; unfold evenBond; split_ifs <;> simp <;> omega)
      simp [openBondFactor, hk2]
    · simp [openBondFactor, hc]

/-- The widths of both stages are nonnegative. -/
theorem firstWidth_nonneg {τ : ℝ} (hτ : 0 ≤ τ) (k : ℕ) : 0 ≤ firstWidth τ k := by
  unfold firstWidth; split_ifs <;> simp [hτ]

/-- The widths of both stages are nonnegative. -/
theorem lastWidth_nonneg {τ : ℝ} (hτ : 0 ≤ τ) (k : ℕ) : 0 ≤ lastWidth τ k := by
  unfold lastWidth; split_ifs <;> simp [hτ]

/-- One slice of the open chain is a sum of tensor products: with one label `ℓ k` per bond, the
site `i` carries the word `C_{ℓ(odd bond)} O C_{ℓ(even bond)}` and the caps their weights. -/
theorem openSlice_eq_sum (hrefl : R.refl = 1) (b : R.Boundary) (n : ℕ) {τ : ℝ} (hτ : 0 ≤ τ) :
    R.openSlice b n τ = ∑ ℓ : Fin (n + 2) → Option ι, (∏ k, R.openBondWeight b n τ k (ℓ k)) •
      piTensor fun i ↦
        R.siteFactor τ (ℓ (oddBond n i)) * R.onsite τ * R.siteFactor τ (ℓ (evenBond n i)) := by
  have hsites : (List.ofFn fun i ↦ 1 - (τ : ℂ) • R.site (n + 1) i).prod =
      piTensor fun _ ↦ R.onsite τ := by
    simp only [one_sub_smul_site, list_prod_ofFn_piTensor]
    congr 1
    ext1 k
    exact list_prod_ofFn_ite (fun i ↦ k = i) _ k fun i ↦ eq_comm
  rw [openSlice, R.openStage_eq_sum hrefl b n _ (firstWidth_nonneg hτ),
    R.openStage_eq_sum hrefl b n _ (lastWidth_nonneg hτ), hsites, Finset.sum_mul,
    Finset.sum_mul_sum, ← Fintype.sum_prod_type']
  -- Only the pairs of stage labels supported on the bonds of their own stage contribute; they
  -- are the images of single labels `ℓ`.
  symm
  refine Fintype.sum_of_injective (fun ℓ : Fin (n + 2) → Option ι ↦
    (fun k : Fin (n + 2) ↦ if (k : ℕ) % 2 = 1 then ℓ k else none,
      fun k : Fin (n + 2) ↦ if (k : ℕ) % 2 = 1 then none else ℓ k)) ?_ _ _ ?_ ?_
  · intro ℓ ℓ' h
    simp only [Prod.mk.injEq] at h
    funext k
    have h1 := congrFun h.1 k
    have h2 := congrFun h.2 k
    split_ifs at h1 h2 <;> assumption
  · rintro ⟨ℓ₁, ℓ₃⟩ hne
    have : ∃ k : Fin (n + 2), (¬ (k : ℕ) % 2 = 1 ∧ ℓ₁ k ≠ none) ∨
        ((k : ℕ) % 2 = 1 ∧ ℓ₃ k ≠ none) := by
      by_contra hc
      push Not at hc
      refine hne ⟨fun k ↦ if (k : ℕ) % 2 = 1 then ℓ₁ k else ℓ₃ k, ?_⟩
      simp only [Prod.mk.injEq]
      constructor <;> funext k <;> have := hc k <;> split_ifs <;> simp_all
    obtain ⟨k, ⟨hk, hα⟩ | ⟨hk, hα⟩⟩ := this
    · obtain ⟨α, hα⟩ := Option.ne_none_iff_exists'.1 hα
      have h0 : R.openBondWeight b n (firstWidth τ k) k (ℓ₁ k) •
          piTensor (R.openBondFactor n (firstWidth τ k) k (ℓ₁ k)) = 0 := by
        simp only [hα, firstWidth, hk, ite_false, piTensor_openBondFactor_zero, smul_zero]
      simp only
      rw [List.prod_eq_zero (List.mem_ofFn.2 ⟨k, h0⟩), zero_mul, zero_mul]
    · obtain ⟨α, hα⟩ := Option.ne_none_iff_exists'.1 hα
      have h0 : R.openBondWeight b n (lastWidth τ k) k (ℓ₃ k) •
          piTensor (R.openBondFactor n (lastWidth τ k) k (ℓ₃ k)) = 0 := by
        simp only [hα, lastWidth, hk, ite_true, piTensor_openBondFactor_zero, smul_zero]
      simp only
      have hz : (List.ofFn fun k : Fin (n + 2) ↦ R.openBondWeight b n (lastWidth τ k) k (ℓ₃ k) •
          piTensor (R.openBondFactor n (lastWidth τ k) k (ℓ₃ k))).prod = 0 :=
        List.prod_eq_zero (List.mem_ofFn.2 ⟨k, h0⟩)
      rw [hz, mul_zero]
  · intro ℓ
    simp only [list_prod_openBond_smul, smul_mul_assoc, mul_smul_comm, smul_smul, piTensor_mul]
    congr 1
    · rw [← Finset.prod_mul_distrib]
      refine Finset.prod_congr rfl fun k _ ↦ ?_
      simp only [firstWidth, lastWidth]
      split_ifs <;> simp
    · congr 1
      funext i
      rw [R.list_prod_openBondFactor_first, R.list_prod_openBondFactor_last]

/-! ### Many slices -/

/-- The trace of `M` slices of the open chain as a sum over the time labels `q k` of the bonds:
the caps contribute their weights and site `i` the one-cell transfer matrix between the labels
of its odd and its even bond. -/
theorem trace_openSlice_pow_eq_sum (hrefl : R.refl = 1) (b : R.Boundary) (n : ℕ) {τ : ℝ}
    (hτ : 0 ≤ τ) (M : ℕ) :
    trace (R.openSlice b n τ ^ M) = ∑ q : Fin (n + 2) → Fin M → Option ι,
      (∏ k, ∏ t, R.openBondWeight b n τ k (q k t)) *
        ∏ i, R.cellTransfer M τ (q (oddBond n i)) (q (evenBond n i)) := by
  rw [← List.prod_replicate, ← List.ofFn_const, R.openSlice_eq_sum hrefl b n hτ,
    list_prod_ofFn_sum, trace_sum]
  refine Fintype.sum_equiv (Equiv.piComm fun _ _ ↦ Option ι) _ _ fun a ↦ ?_
  rw [list_prod_ofFn_smul, trace_smul, list_prod_ofFn_piTensor, trace_piTensor, smul_eq_mul,
    Finset.prod_comm]
  congr 1
  refine Finset.prod_congr rfl fun i _ ↦ ?_
  simp only [cellTransfer, timeWord, of_apply, R.rightFactor_eq_of_refl hrefl]
  rfl

/-- The trace of `M` slices as a matrix element of the alternating product `Fᵀ F Fᵀ ⋯` of `n + 1`
one-cell transfer matrices between the boundary vectors of the two caps. -/
theorem trace_openSlice_pow_eq_dotProduct [DecidableEq ι] (hrefl : R.refl = 1)
    (b : R.Boundary) (n : ℕ) {τ : ℝ} (hτ : 0 ≤ τ) (M : ℕ) :
    trace (R.openSlice b n τ ^ M) = capVector M τ b.w ⬝ᵥ
      ((List.ofFn fun i : Fin (n + 1) ↦
          if (i : ℕ) % 2 = 1 then R.cellTransfer M τ else (R.cellTransfer M τ)ᵀ).prod *ᵥ
        capVector M τ (b.endW (n + 1))) := by
  rw [R.trace_openSlice_pow_eq_sum hrefl b n hτ, dotProduct_list_prod_mulVec]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  have hin : ∀ j : Fin n,
      ∏ t, R.openBondWeight b n τ j.castSucc.succ (q j.castSucc.succ t) = 1 := fun j ↦
    Finset.prod_eq_one fun t _ ↦ by
      have : (j : ℕ) ≠ n := j.isLt.ne
      simp [openBondWeight, this]
  rw [Fin.prod_univ_succ, Fin.prod_univ_castSucc, Finset.prod_eq_one fun j _ ↦ hin j, one_mul,
    Fin.succ_last]
  have hF : ∏ i, R.cellTransfer M τ (q (oddBond n i)) (q (evenBond n i)) =
      ∏ i : Fin (n + 1), (if (i : ℕ) % 2 = 1 then R.cellTransfer M τ else (R.cellTransfer M τ)ᵀ)
        (q i.castSucc) (q i.succ) :=
    Finset.prod_congr rfl fun i _ ↦ by
      unfold oddBond evenBond
      split_ifs <;> rfl
  simp only [hF, openBondWeight, capVector, Fin.val_zero, Fin.val_last, ite_true,
    Nat.succ_ne_zero, ite_false]
  ring

/-! ### Time reversal -/

namespace Boundary

variable {R}

/-- The sign `s(p) = ∏ₜ ε(p t)` of a sequence of time labels. -/
noncomputable def timeSign (b : R.Boundary) {M : ℕ} (p : Fin M → Option ι) : ℝ :=
  ∏ t, b.labelSign (p t)

/-- The sign of a sequence is invariant under time reversal. -/
theorem timeSign_rev (b : R.Boundary) {M : ℕ} (p : Fin M → Option ι) :
    b.timeSign (fun t ↦ p t.rev) = b.timeSign p :=
  Fintype.prod_equiv Fin.revPerm _ _ fun _ ↦ rfl

/-- The sign of a sequence squares to one. -/
theorem timeSign_mul_self (b : R.Boundary) {M : ℕ} (p : Fin M → Option ι) :
    b.timeSign p * b.timeSign p = 1 := by
  rw [timeSign, ← Finset.prod_mul_distrib]
  exact Finset.prod_eq_one fun t _ ↦ b.labelSign_mul_self _

variable [DecidableEq ι]

/-- The time reversal `J p p' = s(p) [p' = p ∘ rev]` on the time labels of one bond. -/
noncomputable def timeReversal (b : R.Boundary) (M : ℕ) :
    Matrix (Fin M → Option ι) (Fin M → Option ι) ℂ :=
  of fun p p' ↦ if p' = (fun t ↦ p t.rev) then (b.timeSign p : ℂ) else 0

/-- `J` is symmetric. -/
theorem transpose_timeReversal (b : R.Boundary) (M : ℕ) :
    (b.timeReversal M)ᵀ = b.timeReversal M := by
  ext p q
  simp only [transpose_apply, timeReversal, of_apply]
  by_cases h : q = fun t ↦ p t.rev
  · subst h
    simp [timeSign_rev]
  · have h' : ¬ p = fun t ↦ q t.rev := fun h' ↦ h (by subst h'; funext t; simp)
    simp [h, h']

/-- `J` is real and symmetric, hence Hermitian. -/
theorem conjTranspose_timeReversal (b : R.Boundary) (M : ℕ) :
    (b.timeReversal M)ᴴ = b.timeReversal M := by
  ext p q
  rw [conjTranspose_apply, ← transpose_apply (b.timeReversal M), transpose_timeReversal]
  simp only [timeReversal, of_apply]
  split_ifs <;> simp

/-- Left multiplication by `J` reverses the row labels. -/
theorem timeReversal_mul_apply (b : R.Boundary) {M : ℕ} (A : Matrix (Fin M → Option ι)
    (Fin M → Option ι) ℂ) (p q : Fin M → Option ι) :
    (b.timeReversal M * A) p q = b.timeSign p * A (fun t ↦ p t.rev) q := by
  simp [mul_apply, timeReversal, ite_mul]

/-- Right multiplication by `J` reverses the column labels. -/
theorem mul_timeReversal_apply (b : R.Boundary) {M : ℕ} (A : Matrix (Fin M → Option ι)
    (Fin M → Option ι) ℂ) (p q : Fin M → Option ι) :
    (A * b.timeReversal M) p q = A p (fun t ↦ q t.rev) * b.timeSign q := by
  have hJ : ∀ y, b.timeReversal M y q =
      if y = (fun t ↦ q t.rev) then (b.timeSign q : ℂ) else 0 := fun y ↦ by
    rw [← transpose_apply (b.timeReversal M) q y, transpose_timeReversal]; rfl
  simp only [mul_apply, hJ, mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

/-- `J` is an involution. -/
theorem timeReversal_mul_self (b : R.Boundary) (M : ℕ) :
    b.timeReversal M * b.timeReversal M = 1 := by
  ext p r
  rw [timeReversal_mul_apply, one_apply]
  simp only [timeReversal, of_apply, Fin.rev_rev, timeSign_rev]
  by_cases h : p = r
  · subst h
    simp only [ite_true]
    exact_mod_cast b.timeSign_mul_self p
  · simp [Ne.symm h, h]

/-- The action of `J` on a boundary vector multiplies each weight by its sign. -/
theorem mulVec_capVector (b : R.Boundary) (M : ℕ) (τ : ℝ) (u : ι → ℂ) :
    b.timeReversal M *ᵥ capVector M τ u = capVector M τ fun α ↦ b.sign α * u α := by
  funext p
  simp only [mulVec, dotProduct, timeReversal, of_apply, ite_mul, zero_mul,
    Finset.sum_ite_eq', Finset.mem_univ, ite_true, capVector, timeSign]
  rw [Fintype.prod_equiv Fin.revPerm (fun t ↦ capWeight τ u (p t.rev))
    (fun t ↦ capWeight τ u (p t)) fun _ ↦ rfl, Complex.ofReal_prod, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun t _ ↦ ?_
  cases p t with
  | none => simp [labelSign, capWeight]
  | some α => simp only [labelSign, capWeight]; ring

omit [DecidableEq ι] in
/-- The conjugate of the boundary vector of `b` is its time reversal. -/
theorem star_capVector (b : R.Boundary) (M : ℕ) (τ : ℝ) :
    star (capVector M τ b.w) = capVector M τ fun α ↦ b.sign α * b.w α := by
  funext p
  simp only [Pi.star_apply, capVector, star_prod]
  refine Finset.prod_congr rfl fun t _ ↦ ?_
  cases p t with
  | none => simp [capWeight]
  | some α => simp [capWeight, (b.sign_spec α).2]

end Boundary

variable [DecidableEq ι]

omit [DecidableEq ι] in
/-- Each word of the transfer matrix transposes to the reversed word, up to the signs. -/
theorem transpose_word (b : R.Boundary) (τ : ℝ) (γ γ' : Option ι) :
    (R.siteFactor τ γ * R.onsite τ * R.siteFactor τ γ')ᵀ =
      ((b.labelSign γ' * b.labelSign γ : ℝ) : ℂ) •
        (R.siteFactor τ γ' * R.onsite τ * R.siteFactor τ γ) := by
  rw [transpose_mul, transpose_mul, R.transpose_siteFactor b, R.transpose_siteFactor b,
    transpose_onsite]
  simp only [Matrix.mul_smul, Matrix.smul_mul, smul_smul, Matrix.mul_assoc, Complex.ofReal_mul]
  congr 1
  ring

omit [DecidableEq ι] in
/-- Reversing time in both labels: `F q p = s(p) s(q) F (p ∘ rev) (q ∘ rev)`. -/
theorem cellTransfer_swap (hrefl : R.refl = 1) (b : R.Boundary) (M : ℕ) (τ : ℝ)
    (p q : Fin M → Option ι) :
    R.cellTransfer M τ q p = b.timeSign p * b.timeSign q *
      R.cellTransfer M τ (fun t ↦ p t.rev) (fun t ↦ q t.rev) := by
  simp only [cellTransfer, timeWord, of_apply, R.rightFactor_eq_of_refl hrefl]
  rw [← trace_transpose, transpose_list_prod, List.map_ofFn, list_ofFn_reverse]
  simp only [Function.comp_apply, R.transpose_word b]
  rw [list_prod_ofFn_smul, trace_smul, smul_eq_mul]
  congr 1
  rw [← Complex.ofReal_mul, ← b.timeSign_rev p, ← b.timeSign_rev q, Boundary.timeSign,
    Boundary.timeSign, ← Finset.prod_mul_distrib, Complex.ofReal_prod]

/-- **Time-reversal symmetry.** `Fᵀ = J F J`. -/
theorem transpose_cellTransfer (hrefl : R.refl = 1) (b : R.Boundary) (M : ℕ) (τ : ℝ) :
    (R.cellTransfer M τ)ᵀ = b.timeReversal M * R.cellTransfer M τ * b.timeReversal M := by
  ext p q
  rw [transpose_apply, Boundary.mul_timeReversal_apply, Boundary.timeReversal_mul_apply,
    R.cellTransfer_swap hrefl b]
  ring

/-! ### The Hermitian transfer matrix -/

/-- The transfer matrix `X = F J` of the open chain and of the ring. -/
noncomputable def openTransfer (b : R.Boundary) (M : ℕ) (τ : ℝ) :
    Matrix (Fin M → Option ι) (Fin M → Option ι) ℂ :=
  R.cellTransfer M τ * b.timeReversal M

/-- `X = F J` is Hermitian: `F` is real and `Fᵀ = J F J`. -/
theorem openTransfer_isHermitian (hrefl : R.refl = 1) (b : R.Boundary) (M : ℕ) (τ : ℝ) :
    (R.openTransfer b M τ).IsHermitian := by
  change (R.cellTransfer M τ * b.timeReversal M)ᴴ = R.cellTransfer M τ * b.timeReversal M
  rw [conjTranspose_mul, Boundary.conjTranspose_timeReversal, conjTranspose_cellTransfer,
    R.transpose_cellTransfer hrefl b, ← Matrix.mul_assoc, ← Matrix.mul_assoc,
    Boundary.timeReversal_mul_self, Matrix.one_mul]

/-- The ring traces are the even power traces of `X`: `Tr (slice)^M = Tr X^{2N}` on `2 N`
cells. -/
theorem trace_trotterSlice_pow_eq_trace_openTransfer (hrefl : R.refl = 1) (b : R.Boundary)
    {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 ≤ τ) (M : ℕ) :
    trace (R.trotterSlice (2 * N) τ ^ M) = trace (R.openTransfer b M τ ^ (2 * N)) := by
  have h : (R.cellTransfer M τ)ᵀ * R.cellTransfer M τ =
      b.timeReversal M * R.openTransfer b M τ ^ 2 * b.timeReversal M := by
    rw [R.transpose_cellTransfer hrefl b, openTransfer, pow_two]
    simp only [Matrix.mul_assoc, Boundary.timeReversal_mul_self, Matrix.mul_one]
  rw [R.trace_trotterSlice_pow_eq_trace_transfer hN hτ, h,
    pow_conj_of_mul_self (b.timeReversal_mul_self M), trace_mul_cycle,
    b.timeReversal_mul_self, Matrix.one_mul, pow_mul]

/-- The alternating product `Fᵀ F Fᵀ ⋯` of `m` factors is `J Xᵐ J` for even `m` and `J Xᵐ` for
odd `m`. -/
theorem list_prod_alternating (hrefl : R.refl = 1) (b : R.Boundary) (M : ℕ) (τ : ℝ) :
    ∀ m : ℕ, (List.ofFn fun i : Fin m ↦
        if (i : ℕ) % 2 = 1 then R.cellTransfer M τ else (R.cellTransfer M τ)ᵀ).prod =
      b.timeReversal M * R.openTransfer b M τ ^ m *
        (if m % 2 = 0 then b.timeReversal M else 1)
  | 0 => by simp [b.timeReversal_mul_self]
  | m + 1 => by
    rw [List.ofFn_succ', List.prod_concat]
    simp only [Fin.val_castSucc, Fin.val_last]
    rw [list_prod_alternating hrefl b M τ m, pow_succ, openTransfer]
    rcases Nat.mod_two_eq_zero_or_one m with h | h
    · have h' : ¬ (m + 1) % 2 = 0 := by omega
      have h'' : ¬ m % 2 = 1 := by omega
      simp only [h, h', zero_ne_one, ite_true, ite_false, R.transpose_cellTransfer hrefl b,
        Matrix.mul_assoc, Matrix.mul_one]
      rw [← Matrix.mul_assoc (b.timeReversal M) (b.timeReversal M),
        b.timeReversal_mul_self, Matrix.one_mul]
    · have h' : (m + 1) % 2 = 0 := by omega
      have h'' : ¬ m % 2 = 0 := by omega
      simp only [h, h', one_ne_zero, ite_true, ite_false, Matrix.mul_assoc, Matrix.mul_one,
        b.timeReversal_mul_self]

/-- **Open-chain factorization.** `Tr (slice)^M = ⟨v, Xⁿ⁺¹ v⟩` on the chain of `n + 1` sites,
with `v` the boundary vector of `b`. -/
theorem trace_openSlice_pow (hrefl : R.refl = 1) (b : R.Boundary) (n : ℕ) {τ : ℝ} (hτ : 0 ≤ τ)
    (M : ℕ) :
    trace (R.openSlice b n τ ^ M) =
      star (capVector M τ b.w) ⬝ᵥ (R.openTransfer b M τ ^ (n + 1)) *ᵥ capVector M τ b.w := by
  have hE : (if (n + 1) % 2 = 0 then b.timeReversal M else 1) *ᵥ
      capVector M τ (b.endW (n + 1)) = capVector M τ b.w := by
    unfold Boundary.endW
    rcases Nat.mod_two_eq_zero_or_one (n + 1) with h | h
    · have he : Even (n + 1) := Nat.even_iff.2 h
      rw [ite_eq_left h, ite_eq_left he, b.mulVec_capVector]
      congr 1
      funext α
      rw [← mul_assoc, ← Complex.ofReal_mul, b.sign_mul_self, Complex.ofReal_one, one_mul]
    · have he : ¬ Even (n + 1) := by rw [Nat.not_even_iff]; exact h
      rw [ite_eq_right (by omega), ite_eq_right he, one_mulVec]
  rw [R.trace_openSlice_pow_eq_dotProduct hrefl b n hτ, R.list_prod_alternating hrefl b,
    ← mulVec_mulVec, hE, ← mulVec_mulVec, dotProduct_mulVec, ← mulVec_transpose,
    Boundary.transpose_timeReversal, b.mulVec_capVector, b.star_capVector]

/-! ### Convergence -/

omit [DecidableEq ι] in
/-- The terms of all bonds are the bonds of the chain and the two end terms. -/
theorem sum_openBondTerm (b : R.Boundary) (n : ℕ) :
    ∑ k, R.openBondTerm b n k = ∑ j : Fin n, bondOp j.castSucc j.succ R.bondMatrix +
      siteOp 0 b.mat + siteOp (Fin.last n) (b.endMat (n + 1)) := by
  have hj : ∀ j : Fin n,
      R.openBondTerm b n j.castSucc.succ = bondOp j.castSucc j.succ R.bondMatrix := by
    intro j
    have : (j : ℕ) ≠ n := j.isLt.ne
    unfold openBondTerm
    rw [ite_eq_right (by simp), dite_eq_right (by simp [this])]
    congr 1
  rw [Fin.sum_univ_succ, Fin.sum_univ_castSucc, Finset.sum_congr rfl fun j _ ↦ hj j]
  have h0 : R.openBondTerm b n 0 = siteOp 0 b.mat := by simp [openBondTerm]
  have hl : R.openBondTerm b n (Fin.last n).succ = siteOp (Fin.last n) (b.endMat (n + 1)) := by
    simp [openBondTerm]
  rw [h0, hl]
  abel

omit [DecidableEq ι] in
/-- The Trotter traces of the open chain on `n + 1` sites converge to its partition function. -/
theorem tendsto_openSliceTrace (b : R.Boundary) (β : ℝ) (n : ℕ) :
    Tendsto (fun M : ℕ ↦ (trace (R.openSlice b n (β / M) ^ M)).re) atTop
      (nhds (R.openPartitionFunction b (n + 1) β)) := by
  let g := R.openBondTerm b n
  let hs : List (Matrix (Fin (n + 1) → d) (Fin (n + 1) → d) ℂ) :=
    List.ofFn (fun k : Fin (n + 2) ↦ if (k : ℕ) % 2 = 1 then -g k else 0) ++
      List.ofFn (fun i ↦ -R.site (n + 1) i) ++
      List.ofFn (fun k : Fin (n + 2) ↦ if (k : ℕ) % 2 = 1 then 0 else -g k)
  have hslice : ∀ τ : ℝ, R.openSlice b n τ = (hs.map fun h ↦ 1 + (τ : ℂ) • h).prod := by
    intro τ
    have hk1 : ∀ k : Fin (n + 2), (1 + (τ : ℂ) • if (k : ℕ) % 2 = 1 then -g k else 0) =
        1 - (firstWidth τ k : ℂ) • R.openBondTerm b n k := fun k ↦ by
      simp only [firstWidth, g]; split_ifs <;> simp [sub_eq_add_neg]
    have hk3 : ∀ k : Fin (n + 2), (1 + (τ : ℂ) • if (k : ℕ) % 2 = 1 then 0 else -g k) =
        1 - (lastWidth τ k : ℂ) • R.openBondTerm b n k := fun k ↦ by
      simp only [lastWidth, g]; split_ifs <;> simp [sub_eq_add_neg]
    simp only [hs, openSlice, openStage, List.map_append, List.prod_append, List.map_ofFn,
      Function.comp_def, hk1, hk3, smul_neg, ← sub_eq_add_neg]
  have hsum : hs.sum = -R.openHamiltonian b (n + 1) := by
    have hb : ∀ k : Fin (n + 2), ((if (k : ℕ) % 2 = 1 then -g k else 0) +
        if (k : ℕ) % 2 = 1 then 0 else -g k) = -g k := fun k ↦ by split_ifs <;> simp
    simp only [hs, List.sum_append, List.sum_ofFn]
    rw [add_right_comm, ← Finset.sum_add_distrib, Finset.sum_congr rfl fun k _ ↦ hb k,
      Finset.sum_neg_distrib, Finset.sum_neg_distrib, R.sum_openBondTerm]
    simp only [openHamiltonian, site]
    abel
  have h := tendsto_list_prod_one_add_smul_pow hs β
  rw [hsum, smul_neg, ← neg_smul] at h
  have hc : Continuous fun A : Matrix (Fin (n + 1) → d) (Fin (n + 1) → d) ℂ ↦ (trace A).re :=
    Complex.continuous_re.comp (continuous_id.matrix_trace)
  have := (hc.tendsto _).comp h
  simpa [Function.comp_def, hslice, openPartitionFunction] using this

end RPRing

end CertifiedGaps
