/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Model
import Mathlib.LinearAlgebra.Matrix.PosDef

/-!
# Transfer-matrix factorization of Trotterized reflection-positive rings

One first-order Trotter slice of width `τ` is the ordered product
`∏_{even j} (1 - τ hⱼ) · ∏ᵢ (1 - τ oᵢ) · ∏_{odd j} (1 - τ hⱼ)`. On a ring of `L = 2 N` cells,
`Tr (slice)^M = Tr T^N` for a positive semidefinite matrix `T` that depends only on `M` and `τ`.

## Proof outline

* With `C_none = 1` and `C_α = √τ C_α`, a bond factor is `1 - τ h = ∑_γ C_γ ⊗ C_γ`, a sum of
  tensor products (`one_sub_smul_bond`), and an on-site factor is the tensor product with
  `O = 1 - τ o` at one site (`one_sub_smul_site`).
* Expanding one slice over one label per bond, the cell `i` carries the word
  `C_{b} O C_{c}`, with `b` the label of its even bond and `c` that of its odd bond
  (`sliceSite_evenSite`, `sliceSite_oddSite`).
* Over `M` slices the cell `i` contributes `F p q = Tr ∏ₜ C_{p t} O C_{q t}` (`cellTransfer`),
  the same function for even and odd cells. Summing the even-bond labels gives `Tr (Fᵀ F)^N`.
* `F` has real entries because `C_α` and `o` are real, so `Fᵀ F = Fᴴ F ⪰ 0`. This is where
  reflection positivity enters.

## Main results

* `RPRing.trace_trotterSlice_pow_eq_trace_transfer`: `Tr (slice)^M = Tr (Fᵀ F)^N`.
* `RPRing.exists_posSemidef_transfer`: the factorization through a positive semidefinite `T`.
-/

open Matrix
open scoped Kronecker ComplexOrder

namespace CertifiedGaps

variable {d : Type} [Fintype d] [DecidableEq d]

/-! ### Ring geometry -/

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

/-- The even site `2 k` of a ring of `2 N` sites. -/
def evenSite (N : ℕ) (k : Fin N) : Fin (2 * N) := ⟨2 * k, by omega⟩

/-- The odd site `2 k + 1` of a ring of `2 N` sites. -/
def oddSite (N : ℕ) (k : Fin N) : Fin (2 * N) := ⟨2 * k + 1, by omega⟩

/-- The even bonds of a ring, in increasing order. -/
def evenBonds (L : ℕ) : List (Fin L) := (List.finRange L).filter fun j ↦ decide (Even j.val)

/-- The odd bonds of a ring, in increasing order. -/
def oddBonds (L : ℕ) : List (Fin L) := (List.finRange L).filter fun j ↦ !decide (Even j.val)

theorem evenBonds_two_mul (N : ℕ) : evenBonds (2 * N) = List.ofFn (evenSite N) := by
  refine List.SortedLT.eq_of_mem_iff ?_ ?_ fun j ↦ ?_
  · exact ((List.sortedLT_finRange _).pairwise.filter _).sortedLT
  · exact List.sortedLT_ofFn_iff.2 fun a b h ↦ by
      simp only [evenSite, Fin.lt_def] at h ⊢; omega
  · simp only [evenBonds, List.mem_filter, List.mem_finRange, true_and, decide_eq_true_eq,
      List.mem_ofFn]
    constructor
    · rintro ⟨r, hr⟩
      exact ⟨⟨r, by omega⟩, Fin.ext (by simp [evenSite]; omega)⟩
    · rintro ⟨k, rfl⟩
      exact ⟨k, by simp [evenSite]; omega⟩

theorem oddBonds_two_mul (N : ℕ) : oddBonds (2 * N) = List.ofFn (oddSite N) := by
  refine List.SortedLT.eq_of_mem_iff ?_ ?_ fun j ↦ ?_
  · exact ((List.sortedLT_finRange _).pairwise.filter _).sortedLT
  · exact List.sortedLT_ofFn_iff.2 fun a b h ↦ by
      simp only [oddSite, Fin.lt_def] at h ⊢; omega
  · simp only [oddBonds, List.mem_filter, List.mem_finRange, true_and, Bool.not_eq_eq_eq_not,
      Bool.not_true, decide_eq_false_iff_not, Nat.not_even_iff_odd, List.mem_ofFn]
    constructor
    · rintro ⟨r, hr⟩
      exact ⟨⟨r, by omega⟩, Fin.ext (by simp [oddSite]; omega)⟩
    · rintro ⟨k, rfl⟩
      exact ⟨k, by simp [oddSite]⟩

/-- A product over the `2 N` sites of the ring, grouped into the pairs `(2 m, 2 m + 1)`. -/
theorem prod_fin_two_mul {R : Type*} [CommMonoid R] (N : ℕ) (f : Fin (2 * N) → R) :
    ∏ i, f i = ∏ m, f (evenSite N m) * f (oddSite N m) := by
  rw [← (finProdFinEquiv.trans (finCongr (Nat.mul_comm N 2))).prod_comp, Fintype.prod_prod_type]
  refine Finset.prod_congr rfl fun m _ ↦ ?_
  rw [Fin.prod_univ_two]
  congr 2 <;> ext <;> simp [evenSite, oddSite]; ring

/-! ### Local factors -/

omit [Fintype d] in
/-- `1 - c A` at one site, as a tensor product. -/
theorem one_sub_smul_siteOp {L : ℕ} (j : Fin L) (A : Matrix d d ℂ) (c : ℂ) :
    1 - c • siteOp j A = piTensor fun i ↦ if i = j then 1 - c • A else 1 := by
  ext σ τ
  have key : ∀ f : Fin L → Matrix d d ℂ, (∀ i, i ≠ j → f i = 1) →
      ∏ i, f i (σ i) (τ i) =
        f j (σ j) (τ j) * ∏ i ∈ Finset.univ.erase j, (1 : Matrix d d ℂ) (σ i) (τ i) := by
    intro f hf
    rw [← Finset.mul_prod_erase _ _ (Finset.mem_univ j)]
    congr 1
    exact Finset.prod_congr rfl fun i hi ↦ by rw [hf i (Finset.ne_of_mem_erase hi)]
  have h1 : (1 : Matrix (Fin L → d) (Fin L → d) ℂ) σ τ =
      (1 : Matrix d d ℂ) (σ j) (τ j) * ∏ i ∈ Finset.univ.erase j, (1 : Matrix d d ℂ) (σ i) (τ i) :=
    by rw [← piTensor_one, piTensor_apply]; exact key _ fun _ _ ↦ rfl
  rw [Matrix.sub_apply, Matrix.smul_apply, h1, siteOp, piTensor_apply, piTensor_apply,
    key _ fun i hi ↦ by simp [hi], key _ fun i hi ↦ by simp [hi]]
  simp only [ite_true, Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul]
  ring

/-! ### Which end of which bond a site is -/

theorem evenSite_eq_evenSite_iff {N : ℕ} (m k : Fin N) : evenSite N m = evenSite N k ↔ k = m := by
  simp only [Fin.ext_iff, evenSite]; omega

theorem evenSite_ne_rot_evenSite {N : ℕ} (m k : Fin N) :
    evenSite N m ≠ finRotate (2 * N) (evenSite N k) := by
  simp only [ne_eq, Fin.ext_iff, evenSite, val_finRotate]
  have := m.isLt; have := k.isLt
  split_ifs <;> omega

theorem evenSite_ne_oddSite {N : ℕ} (m k : Fin N) : evenSite N m ≠ oddSite N k := by
  simp only [ne_eq, Fin.ext_iff, evenSite, oddSite]; omega

theorem evenSite_eq_rot_oddSite_iff {N : ℕ} (m k : Fin N) :
    evenSite N m = finRotate (2 * N) (oddSite N k) ↔ k = (finRotate N).symm m := by
  rw [Equiv.eq_symm_apply]
  simp only [Fin.ext_iff, evenSite, oddSite, val_finRotate]
  have := m.isLt; have := k.isLt
  split_ifs <;> omega

theorem oddSite_ne_evenSite {N : ℕ} (m k : Fin N) : oddSite N m ≠ evenSite N k := by
  simp only [ne_eq, Fin.ext_iff, evenSite, oddSite]; omega

theorem oddSite_eq_rot_evenSite_iff {N : ℕ} (m k : Fin N) :
    oddSite N m = finRotate (2 * N) (evenSite N k) ↔ k = m := by
  simp only [Fin.ext_iff, evenSite, oddSite, val_finRotate]
  have := m.isLt; have := k.isLt
  split_ifs <;> omega

theorem oddSite_eq_oddSite_iff {N : ℕ} (m k : Fin N) : oddSite N m = oddSite N k ↔ k = m := by
  simp only [Fin.ext_iff, oddSite]; omega

theorem oddSite_ne_rot_oddSite {N : ℕ} (m k : Fin N) :
    oddSite N m ≠ finRotate (2 * N) (oddSite N k) := by
  simp only [ne_eq, Fin.ext_iff, oddSite, val_finRotate]
  have := m.isLt; have := k.isLt
  split_ifs <;> first | exact not_false | omega

/-- A product of conjugates by an involution is the conjugate of the product. -/
theorem list_prod_map_conj {M X : Type*} [Monoid M] {r : M} (hr : r * r = 1) (f : X → M) :
    ∀ l : List X, (l.map fun x ↦ r * f x * r).prod = r * (l.map f).prod * r
  | [] => by simp [hr]
  | x :: l => by
    rw [List.map_cons, List.prod_cons, list_prod_map_conj hr f l, List.map_cons, List.prod_cons]
    simp only [mul_assoc]
    rw [← mul_assoc r r, hr, one_mul]

namespace RPRing

variable {ι : Type} [Fintype ι] (R : RPRing d ι)

/-- The real local factors of one Trotter bond factor at the left end of the bond:
`C_none = 1`, `C_α = √τ C_α`. -/
noncomputable def realSiteFactor (τ : ℝ) : Option ι → Matrix d d ℝ
  | none => 1
  | some α => Real.sqrt τ • R.C α

/-- The left local factors as complex matrices with real entries. -/
noncomputable def siteFactor (τ : ℝ) (γ : Option ι) : Matrix d d ℂ :=
  (R.realSiteFactor τ γ).map (↑)

/-- The real local factors at the right end of a bond: `refl C_γ refl`. -/
noncomputable def realRightFactor (τ : ℝ) (γ : Option ι) : Matrix d d ℝ :=
  R.refl * R.realSiteFactor τ γ * R.refl

/-- The right local factors as complex matrices with real entries. -/
noncomputable def rightFactor (τ : ℝ) (γ : Option ι) : Matrix d d ℂ :=
  (R.realRightFactor τ γ).map (↑)

/-- The real on-site factor `O = 1 - τ o`. -/
def realOnsite (τ : ℝ) : Matrix d d ℝ := 1 - τ • R.o

/-- The on-site factor `O = 1 - τ o`, as a complex matrix with real entries. -/
def onsite (τ : ℝ) : Matrix d d ℂ := (R.realOnsite τ).map (↑)

theorem onsite_eq (τ : ℝ) : R.onsite τ = 1 - (τ : ℂ) • R.oc := by
  ext a b
  simp [onsite, realOnsite, oc, one_apply, apply_ite]

theorem Rc_mul_Rc : R.Rc * R.Rc = 1 := by
  change Complex.ofRealHom.mapMatrix R.refl * Complex.ofRealHom.mapMatrix R.refl = 1
  rw [← map_mul, R.refl_sq, map_one]

theorem rightFactor_eq (τ : ℝ) (γ : Option ι) :
    R.rightFactor τ γ = R.Rc * R.siteFactor τ γ * R.Rc := by
  change Complex.ofRealHom.mapMatrix (R.refl * R.realSiteFactor τ γ * R.refl) = _
  rw [map_mul, map_mul]
  rfl

@[simp] theorem siteFactor_none (τ : ℝ) : R.siteFactor τ none = 1 := by
  simp [siteFactor, realSiteFactor]

@[simp] theorem rightFactor_none (τ : ℝ) : R.rightFactor τ none = 1 := by
  rw [rightFactor_eq, siteFactor_none, mul_one, Rc_mul_Rc]

/-- `C_α ⊗ (refl C_α refl) = τ (C_α ⊗ D_α)` for the scaled factors. -/
theorem kronecker_siteFactor_some {τ : ℝ} (hτ : 0 ≤ τ) (α : ι) :
    R.siteFactor τ (some α) ⊗ₖ R.rightFactor τ (some α) = (τ : ℂ) • (R.Cc α ⊗ₖ R.Dc α) := by
  have hs : R.siteFactor τ (some α) = (Real.sqrt τ : ℂ) • R.Cc α := by
    ext a b; simp [siteFactor, realSiteFactor, Cc]
  have h : (Real.sqrt τ : ℂ) * Real.sqrt τ = τ := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt hτ]
  rw [rightFactor_eq, hs, Dc, Matrix.mul_smul, Matrix.smul_mul, smul_kronecker, kronecker_smul,
    smul_smul, h]

/-- The tensor-product factor of label `γ` of the bond `j`: `C_γ` at `j`, `refl C_γ refl` at
`j + 1`. -/
noncomputable def bondSite (L : ℕ) (τ : ℝ) (j : Fin L) (γ : Option ι) :
    Fin L → Matrix d d ℂ :=
  fun i ↦ if i = j then R.siteFactor τ γ else if i = finRotate L j then R.rightFactor τ γ else 1

/-- A first-order bond factor as a sum over the bond label: `1 - τ h = ∑_γ C_γ ⊗ D_γ`. -/
theorem one_sub_smul_bond {L : ℕ} (hL : 2 ≤ L) {τ : ℝ} (hτ : 0 ≤ τ) (j : Fin L) :
    1 - (τ : ℂ) • R.bond L j = ∑ γ, piTensor (R.bondSite L τ j γ) := by
  have hjk := ne_finRotate hL j
  have h1 : R.bondSite L τ j none = fun _ ↦ 1 := funext fun i ↦ by
    unfold bondSite; rw [siteFactor_none, rightFactor_none, ite_self, ite_self]
  rw [Fintype.sum_option, h1, piTensor_one]
  rw [bond, bondMatrix, bondOp_neg, bondOp_sum, smul_neg, sub_neg_eq_add, Finset.smul_sum]
  congr 1
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  change _ = piTensor fun i ↦ if i = j then R.siteFactor τ (some α) else
    if i = finRotate L j then R.rightFactor τ (some α) else 1
  rw [← bondOp_kronecker hjk, kronecker_siteFactor_some R hτ, bondOp_smul]

/-- A first-order on-site factor as a tensor product. -/
theorem one_sub_smul_site (L : ℕ) (τ : ℝ) (i : Fin L) :
    1 - (τ : ℂ) • R.site L i = piTensor fun k ↦ if k = i then R.onsite τ else 1 := by
  rw [site, one_sub_smul_siteOp, onsite_eq]

/-! ### Trotter slices -/

/-- One first-order Trotter slice: even bonds, then all on-site terms, then odd bonds. -/
def trotterSlice (L : ℕ) (τ : ℝ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  ((evenBonds L).map fun j ↦ 1 - (τ : ℂ) • R.bond L j).prod *
    (List.ofFn fun i ↦ 1 - (τ : ℂ) • R.site L i).prod *
    ((oddBonds L).map fun j ↦ 1 - (τ : ℂ) • R.bond L j).prod

/-- The Trotter approximant `Re Tr (trotterSlice L (β/M))^M` of `Z_L(β)`. -/
noncomputable def trotterTrace (M : ℕ) (β : ℝ) (L : ℕ) : ℝ :=
  (trace (R.trotterSlice L (β / M) ^ M)).re

/-- The local word at cell `i` of one slice with even-bond labels `bc.1` and odd-bond labels
`bc.2`. -/
noncomputable def sliceSite (N : ℕ) (τ : ℝ) (bc : (Fin N → Option ι) × (Fin N → Option ι))
    (i : Fin (2 * N)) : Matrix d d ℂ :=
  (List.ofFn fun k ↦ R.bondSite (2 * N) τ (evenSite N k) (bc.1 k) i).prod * R.onsite τ *
    (List.ofFn fun k ↦ R.bondSite (2 * N) τ (oddSite N k) (bc.2 k) i).prod

/-- One Trotter slice on a ring of `2 N` cells is a sum of tensor products. -/
theorem trotterSlice_eq_sum {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 ≤ τ) :
    R.trotterSlice (2 * N) τ = ∑ bc, piTensor (R.sliceSite N τ bc) := by
  have hsites : (List.ofFn fun i ↦ 1 - (τ : ℂ) • R.site (2 * N) i).prod =
      piTensor fun _ ↦ R.onsite τ := by
    simp only [one_sub_smul_site, list_prod_ofFn_piTensor]
    congr 1
    ext1 k
    exact list_prod_ofFn_ite (fun i ↦ k = i) _ k fun i ↦ eq_comm
  rw [trotterSlice, evenBonds_two_mul, oddBonds_two_mul, List.map_ofFn, List.map_ofFn, hsites]
  simp only [Function.comp_def, R.one_sub_smul_bond (by omega : 2 ≤ 2 * N) hτ]
  rw [list_prod_ofFn_sum, list_prod_ofFn_sum, Finset.sum_mul, Finset.sum_mul_sum,
    ← Finset.sum_product']
  refine Finset.sum_congr rfl fun bc _ ↦ ?_
  rw [list_prod_ofFn_piTensor, list_prod_ofFn_piTensor, piTensor_mul, piTensor_mul]
  rfl

/-- The word of a slice at the even cell `2 m`: the left factor of its even bond, the on-site
factor, the right factor of its odd bond. -/
theorem sliceSite_evenSite (N : ℕ) (τ : ℝ) (bc : (Fin N → Option ι) × (Fin N → Option ι))
    (m : Fin N) :
    R.sliceSite N τ bc (evenSite N m) =
      R.siteFactor τ (bc.1 m) * R.onsite τ * R.rightFactor τ (bc.2 ((finRotate N).symm m)) := by
  simp only [sliceSite, bondSite, evenSite_ne_rot_evenSite, evenSite_ne_oddSite, ite_false]
  rw [list_prod_ofFn_ite _ _ m (evenSite_eq_evenSite_iff m),
    list_prod_ofFn_ite _ _ _ (evenSite_eq_rot_oddSite_iff m)]

/-- The word of a slice at the odd cell `2 m + 1`: the right factor of its even bond, the on-site
factor, the left factor of its odd bond. -/
theorem sliceSite_oddSite (N : ℕ) (τ : ℝ) (bc : (Fin N → Option ι) × (Fin N → Option ι))
    (m : Fin N) :
    R.sliceSite N τ bc (oddSite N m) =
      R.rightFactor τ (bc.1 m) * R.onsite τ * R.siteFactor τ (bc.2 m) := by
  simp only [sliceSite, bondSite, oddSite_ne_evenSite, oddSite_ne_rot_oddSite, ite_false]
  rw [list_prod_ofFn_ite _ _ m (oddSite_eq_rot_evenSite_iff m),
    list_prod_ofFn_ite _ _ m (oddSite_eq_oddSite_iff m)]

/-! ### Many slices and the transfer matrix -/

/-- The time-ordered word on an even cell over `M` slices, `∏ₜ C_{p t} O D_{q t}`. -/
noncomputable def timeWord (M : ℕ) (τ : ℝ) (p q : Fin M → Option ι) : Matrix d d ℂ :=
  (List.ofFn fun t ↦ R.siteFactor τ (p t) * R.onsite τ * R.rightFactor τ (q t)).prod

/-- The one-cell transfer matrix `F p q = Tr (timeWord p q)` on time labels. -/
noncomputable def cellTransfer (M : ℕ) (τ : ℝ) :
    Matrix (Fin M → Option ι) (Fin M → Option ι) ℂ :=
  of fun p q ↦ trace (R.timeWord M τ p q)

/-- The word on an odd cell is the reflection of the word on an even cell, so it has the same
trace: `Tr ∏ₜ D_{p t} O C_{q t} = F p q`. -/
theorem trace_oddWord (M : ℕ) (τ : ℝ) (p q : Fin M → Option ι) :
    trace (List.ofFn fun t ↦ R.rightFactor τ (p t) * R.onsite τ * R.siteFactor τ (q t)).prod =
      R.cellTransfer M τ p q := by
  have hreal : R.refl * R.o = R.o * R.refl := by
    have h := congrArg (· * R.refl) R.refl_o
    simp only [Matrix.mul_assoc, R.refl_sq, Matrix.mul_one] at h
    exact h
  have hO : R.Rc * R.onsite τ = R.onsite τ * R.Rc := by
    change Complex.ofRealHom.mapMatrix R.refl * Complex.ofRealHom.mapMatrix (1 - τ • R.o) =
      Complex.ofRealHom.mapMatrix (1 - τ • R.o) * Complex.ofRealHom.mapMatrix R.refl
    rw [← map_mul, ← map_mul, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, Matrix.one_mul,
      Matrix.mul_smul, Matrix.smul_mul, hreal]
  have hword : ∀ t, R.rightFactor τ (p t) * R.onsite τ * R.siteFactor τ (q t) =
      R.Rc * (R.siteFactor τ (p t) * R.onsite τ * R.rightFactor τ (q t)) * R.Rc := by
    intro t
    rw [rightFactor_eq, rightFactor_eq]
    calc R.Rc * R.siteFactor τ (p t) * R.Rc * R.onsite τ * R.siteFactor τ (q t)
        = R.Rc * R.siteFactor τ (p t) * (R.Rc * R.onsite τ) * R.siteFactor τ (q t) := by
          simp only [Matrix.mul_assoc]
      _ = R.Rc * R.siteFactor τ (p t) * (R.onsite τ * R.Rc) * R.siteFactor τ (q t) := by rw [hO]
      _ = R.Rc * (R.siteFactor τ (p t) * R.onsite τ * (R.Rc * R.siteFactor τ (q t) * R.Rc)) *
            R.Rc := by
          simp only [Matrix.mul_assoc, R.Rc_mul_Rc, Matrix.mul_one]
  simp only [hword, List.ofFn_eq_map]
  rw [list_prod_map_conj R.Rc_mul_Rc, trace_mul_cycle, R.Rc_mul_Rc, Matrix.one_mul]
  simp only [cellTransfer, of_apply, timeWord, List.ofFn_eq_map]

/-- The trace of `M` Trotter slices as a sum over space-time bond labels. -/
theorem trace_trotterSlice_pow_eq_sum {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 ≤ τ) (M : ℕ) :
    trace (R.trotterSlice (2 * N) τ ^ M) =
      ∑ a : Fin M → (Fin N → Option ι) × (Fin N → Option ι),
        ∏ m, R.cellTransfer M τ (fun t ↦ (a t).1 m) (fun t ↦ (a t).2 ((finRotate N).symm m)) *
          R.cellTransfer M τ (fun t ↦ (a t).1 m) (fun t ↦ (a t).2 m) := by
  rw [← List.prod_replicate, ← List.ofFn_const, R.trotterSlice_eq_sum hN hτ, list_prod_ofFn_sum,
    trace_sum]
  refine Finset.sum_congr rfl fun a _ ↦ ?_
  rw [list_prod_ofFn_piTensor, trace_piTensor, prod_fin_two_mul]
  refine Finset.prod_congr rfl fun m _ ↦ ?_
  simp only [sliceSite_evenSite, sliceSite_oddSite]
  rw [R.trace_oddWord]
  rfl
/-- **Transfer factorization.** `Tr (slice)^M = Tr (Fᵀ F)^N` on a ring of `2 N` cells. -/
theorem trace_trotterSlice_pow_eq_trace_transfer [DecidableEq ι] {N : ℕ} (hN : 1 ≤ N) {τ : ℝ}
    (hτ : 0 ≤ τ) (M : ℕ) :
    trace (R.trotterSlice (2 * N) τ ^ M) =
      trace (((R.cellTransfer M τ)ᵀ * R.cellTransfer M τ) ^ N) := by
  rw [R.trace_trotterSlice_pow_eq_sum hN hτ, trace_pow_eq_sum_closedPath _ hN]
  let F := R.cellTransfer M τ
  let e := (Equiv.arrowProdEquivProdArrow (Fin M) (fun _ ↦ Fin N → Option ι)
    (fun _ ↦ Fin N → Option ι)).trans
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

/-- The one-cell transfer matrix has real entries. -/
theorem star_cellTransfer_apply (M : ℕ) (τ : ℝ) (p q : Fin M → Option ι) :
    star (R.cellTransfer M τ p q) = R.cellTransfer M τ p q := by
  have h : R.timeWord M τ p q = (List.ofFn fun t ↦ R.realSiteFactor τ (p t) *
      R.realOnsite τ * R.realRightFactor τ (q t)).prod.map (↑) := by
    change _ = (Complex.ofRealHom.mapMatrix : Matrix d d ℝ →+* _) _
    rw [map_list_prod, List.map_ofFn]
    simp only [timeWord, siteFactor, rightFactor, onsite, Function.comp_def,
      RingHom.mapMatrix_apply, Matrix.map_mul]
    rfl
  simp only [cellTransfer, of_apply, h, trace, diag_apply, map_apply, ← Complex.ofReal_sum]
  exact Complex.conj_ofReal _

/-- The conjugate transpose of the one-cell transfer matrix is its transpose. -/
theorem conjTranspose_cellTransfer (M : ℕ) (τ : ℝ) :
    (R.cellTransfer M τ)ᴴ = (R.cellTransfer M τ)ᵀ := by
  ext p q
  exact R.star_cellTransfer_apply M τ q p

/-- **Transfer-matrix factorization of a Trotterized reflection-positive ring.** For every number
of slices `M` and every `τ ≥ 0` there is a positive semidefinite `T`, independent of the ring
size, with `Tr (slice)^M = Tr T^N` on every ring of `2 N ≥ 2` cells. -/
theorem exists_posSemidef_transfer (M : ℕ) (τ : ℝ) (hτ : 0 ≤ τ) :
    ∃ (D : Type) (_ : Fintype D) (_ : DecidableEq D) (T : Matrix D D ℂ), T.PosSemidef ∧
      ∀ N : ℕ, 1 ≤ N → trace (R.trotterSlice (2 * N) τ ^ M) = trace (T ^ N) := by
  classical
  exact ⟨Fin M → Option ι, inferInstance, inferInstance,
      (R.cellTransfer M τ)ᴴ * R.cellTransfer M τ, posSemidef_conjTranspose_mul_self _,
      fun N hN ↦ by
        rw [conjTranspose_cellTransfer, R.trace_trotterSlice_pow_eq_trace_transfer hN hτ]⟩

end RPRing

end CertifiedGaps
