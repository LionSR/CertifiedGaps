/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.Logic.Equiv.Fin.Rotate
import CertifiedGaps.Transfer.PiTensor

/-!
# Reflection-positive rings

A translation-invariant ring of `L` cells with local space `d` and Hamiltonian

  `H_L = ∑ⱼ hⱼ,ⱼ₊₁ + ∑ⱼ oⱼ`,  `-h = ∑_α C_α ⊗ C_α`,

where the bond factors `C_α` and the on-site term `o` are real matrices. This is the reflection
positivity hypothesis of the thermal-MPO method in its real form: in the chosen local basis the
bond is minus a sum of squares of real tensor products. Each `C_α` is required to be symmetric or
antisymmetric and `o` symmetric, which makes `H_L` Hermitian.

A model with a multi-site unit cell is a ring of cells: `d` is the cell space, the bond couples
neighbouring cells and the intra-cell terms are part of `o`.

## Definitions

* `RPRing d ι`: the data `C : ι → Matrix d d ℝ` and `o : Matrix d d ℝ` with the symmetry
  conditions.
* `bondOp j k B`, `siteOp j A`: two-site and one-site operators on the ring.
* `RPRing.hamiltonian L`, `RPRing.spectrum L`, `RPRing.partitionFunction L β`.

## Main results

* `RPRing.hamiltonian_isHermitian`.
* `RPRing.partitionFunction_eq_sum`: `Re Tr e^{-βH_L} = ∑ₘ e^{-β Eₘ}`.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace CertifiedGaps

variable {d : Type} [Fintype d] [DecidableEq d]

/-! ### Operators on the ring -/

/-- The operator `B` on two sites `j, k` of a ring of `L` sites, acting as the identity on all
other sites. -/
def bondOp {L : ℕ} (j k : Fin L) (B : Matrix (d × d) (d × d) ℂ) :
    Matrix (Fin L → d) (Fin L → d) ℂ :=
  of fun σ τ ↦ if ∀ i, i ≠ j → i ≠ k → σ i = τ i then B (σ j, σ k) (τ j, τ k) else 0

/-- The operator `A` on the site `j` of a ring of `L` sites, acting as the identity elsewhere. -/
def siteOp {L : ℕ} (j : Fin L) (A : Matrix d d ℂ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  piTensor fun i ↦ if i = j then A else 1

omit [Fintype d] in
/-- A Hermitian two-site operator embeds as a Hermitian operator on the ring. -/
theorem bondOp_isHermitian {L : ℕ} (j k : Fin L) {B : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) : (bondOp j k B).IsHermitian := by
  ext σ τ
  simp only [conjTranspose_apply, bondOp, of_apply]
  split_ifs with h₁ h₂ h₂
  · exact hB.apply _ _
  · exact (h₂ fun i hj hk ↦ (h₁ i hj hk).symm).elim
  · exact (h₁ fun i hj hk ↦ (h₂ i hj hk).symm).elim
  · exact star_zero _

omit [Fintype d] [DecidableEq d] in
/-- The conjugate transpose of a tensor product is the tensor product of the conjugate
transposes. -/
theorem conjTranspose_piTensor {ι : Type*} [Fintype ι] (A : ι → Matrix d d ℂ) :
    (piTensor A)ᴴ = piTensor fun i ↦ (A i)ᴴ := by
  ext σ τ
  simp [piTensor_apply, conjTranspose_apply, star_prod]

omit [Fintype d] in
/-- A Hermitian one-site operator embeds as a Hermitian operator on the ring. -/
theorem siteOp_isHermitian {L : ℕ} (j : Fin L) {A : Matrix d d ℂ} (hA : A.IsHermitian) :
    (siteOp j A).IsHermitian := by
  change (siteOp j A)ᴴ = siteOp j A
  rw [siteOp, conjTranspose_piTensor]
  congr 1
  ext1 i
  split_ifs
  · exact hA
  · exact conjTranspose_one

omit [Fintype d] in
/-- `bondOp` is additive in the two-site operator. -/
theorem bondOp_sum {L : ℕ} (j k : Fin L) {X : Type*} (s : Finset X)
    (B : X → Matrix (d × d) (d × d) ℂ) :
    bondOp j k (∑ x ∈ s, B x) = ∑ x ∈ s, bondOp j k (B x) := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.sum_apply]
  split_ifs <;> simp

omit [Fintype d] in
/-- `bondOp` is homogeneous in the two-site operator. -/
theorem bondOp_smul {L : ℕ} (j k : Fin L) (c : ℂ) (B : Matrix (d × d) (d × d) ℂ) :
    bondOp j k (c • B) = c • bondOp j k B := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.smul_apply, smul_eq_mul]
  split_ifs <;> simp

omit [Fintype d] in
/-- `bondOp` commutes with negation. -/
theorem bondOp_neg {L : ℕ} (j k : Fin L) (B : Matrix (d × d) (d × d) ℂ) :
    bondOp j k (-B) = -bondOp j k B := by
  simpa using bondOp_smul j k (-1) B

omit [Fintype d] in
/-- A product operator `A ⊗ B` on two distinct sites is a tensor product over the ring. -/
theorem bondOp_kronecker {L : ℕ} {j k : Fin L} (hjk : j ≠ k) (A B : Matrix d d ℂ) :
    bondOp j k (A ⊗ₖ B) = piTensor fun i ↦ if i = j then A else if i = k then B else 1 := by
  ext σ τ
  have hk : k ∈ Finset.univ.erase j := Finset.mem_erase.2 ⟨hjk.symm, Finset.mem_univ _⟩
  rw [piTensor_apply, ← Finset.mul_prod_erase _ _ (Finset.mem_univ j),
    ← Finset.mul_prod_erase _ _ hk,
    Finset.prod_congr rfl (g := fun i ↦ if σ i = τ i then (1 : ℂ) else 0) fun i hi ↦ ?_,
    Finset.prod_ite_zero, Finset.prod_const_one]
  · simp only [bondOp, of_apply, kroneckerMap_apply, ite_true, Ne.symm hjk, ite_false,
      Finset.mem_erase, Finset.mem_univ, and_true]
    split_ifs with h₁ h₂ h₂ <;> simp_all
  · simp only [Finset.mem_erase, Finset.mem_univ, and_true] at hi
    simp [hi.1, hi.2, one_apply]

/-! ### The trace of a Gibbs operator -/

/-- For a Hermitian matrix `A` and a real scalar `c`, `Tr exp (c A) = ∑ᵢ exp (c λᵢ)`. -/
theorem trace_exp_smul_eq_sum_of_isHermitian {n : Type*} [Fintype n] [DecidableEq n]
    {A : Matrix n n ℂ} (hA : A.IsHermitian) (c : ℝ) :
    trace (exp ((c : ℂ) • A)) = ∑ i, (Real.exp (c * hA.eigenvalues i) : ℂ) := by
  set U : Matrix n n ℂ := (hA.eigenvectorUnitary : Matrix n n ℂ)
  have hU : star U * U = 1 := Unitary.coe_star_mul_self _
  have hUinv : U⁻¹ = star U := inv_eq_left_inv hU
  have hUunit : IsUnit U := (isUnit_iff_isUnit_det U).mpr
    (isUnit_det_of_left_inverse hU)
  have hA' : (c : ℂ) • A = U * diagonal (fun i ↦ ((c * hA.eigenvalues i : ℝ) : ℂ)) * U⁻¹ := by
    conv_lhs => rw [hA.spectral_theorem]
    rw [Unitary.conjStarAlgAut_apply, hUinv, ← smul_mul_assoc, ← mul_smul_comm, ← diagonal_smul]
    congr
    ext i
    simp
  rw [hA', exp_conj _ _ hUunit, exp_diagonal, trace_mul_cycle, hUinv, hU, one_mul,
    trace_diagonal]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [Pi.coe_exp, ← Complex.exp_eq_exp_ℂ, Complex.ofReal_exp]

/-! ### Reflection-positive rings -/

/-- A reflection-positive ring: real bond factors `C_α`, each symmetric or antisymmetric, a real
symmetric on-site term `o`, and an internal reflection `refl` of the cell (a real symmetric
involution commuting with `o`, the identity unless the cell has several sites). The bond is
`h = -∑_α C_α ⊗ (refl C_α refl)`. -/
structure RPRing (d ι : Type) [Fintype d] [DecidableEq d] [Fintype ι] where
  /-- The real bond factors. -/
  C : ι → Matrix d d ℝ
  /-- The real on-site term. -/
  o : Matrix d d ℝ
  /-- Each bond factor is symmetric or antisymmetric. -/
  C_symm : ∀ α, (C α)ᵀ = C α ∨ (C α)ᵀ = -C α
  /-- The on-site term is symmetric. -/
  o_symm : oᵀ = o
  /-- The internal reflection of a cell. -/
  refl : Matrix d d ℝ := 1
  /-- The reflection is symmetric. -/
  refl_symm : reflᵀ = refl := by simp
  /-- The reflection is an involution. -/
  refl_sq : refl * refl = 1 := by simp
  /-- The reflection commutes with the on-site term. -/
  refl_o : refl * o * refl = o := by simp

namespace RPRing

variable {ι : Type} [Fintype ι] (R : RPRing d ι)

/-- The bond factors as complex matrices. -/
def Cc (α : ι) : Matrix d d ℂ := (R.C α).map (↑)

/-- The on-site term as a complex matrix. -/
def oc : Matrix d d ℂ := R.o.map (↑)

/-- The reflection as a complex matrix. -/
def Rc : Matrix d d ℂ := R.refl.map (↑)

/-- The reflected bond factors `refl C_α refl` as complex matrices. -/
def Dc (α : ι) : Matrix d d ℂ := R.Rc * R.Cc α * R.Rc

/-- The two-site bond `h = -∑_α C_α ⊗ (refl C_α refl)`. -/
def bondMatrix : Matrix (d × d) (d × d) ℂ := -∑ α, R.Cc α ⊗ₖ R.Dc α

/-- Without an internal reflection the reflected factors are the factors. -/
theorem Dc_eq_of_refl (h : R.refl = 1) (α : ι) : R.Dc α = R.Cc α := by
  simp [Dc, Rc, h, Matrix.map_one]

/-- The bond term `hⱼ,ⱼ₊₁` of the ring. -/
def bond (L : ℕ) (j : Fin L) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  bondOp j (finRotate L j) R.bondMatrix

/-- The on-site term `oⱼ` of the ring. -/
def site (L : ℕ) (j : Fin L) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  siteOp j R.oc

/-- The ring Hamiltonian `H_L = ∑ⱼ hⱼ,ⱼ₊₁ + ∑ⱼ oⱼ`. -/
def hamiltonian (L : ℕ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  ∑ j, R.bond L j + ∑ j, R.site L j

theorem conjTranspose_Cc (α : ι) : (R.Cc α)ᴴ = R.Cc α ∨ (R.Cc α)ᴴ = -R.Cc α := by
  have h : (R.Cc α)ᴴ = ((R.C α)ᵀ).map (↑) := by
    ext a b; simp [Cc, conjTranspose_apply]
  rcases R.C_symm α with hC | hC
  · exact Or.inl (by rw [h, hC]; rfl)
  · exact Or.inr (by rw [h, hC]; ext; simp [Cc])

theorem conjTranspose_Rc : R.Rcᴴ = R.Rc := by
  ext a b
  simp only [conjTranspose_apply, Rc, map_apply, Complex.star_def, Complex.conj_ofReal]
  rw [← transpose_apply R.refl, R.refl_symm]

theorem bondMatrix_isHermitian : R.bondMatrix.IsHermitian := by
  refine IsHermitian.neg (isSelfAdjoint_sum _ fun α _ ↦ ?_)
  change (R.Cc α ⊗ₖ R.Dc α)ᴴ = R.Cc α ⊗ₖ R.Dc α
  rw [conjTranspose_kronecker, Dc, conjTranspose_mul, conjTranspose_mul, conjTranspose_Rc,
    ← Matrix.mul_assoc]
  rcases R.conjTranspose_Cc α with h | h <;> rw [h]
  ext ⟨a, b⟩ ⟨c, e⟩
  simp [kroneckerMap_apply]

theorem oc_isHermitian : R.oc.IsHermitian := by
  ext a b
  simp only [conjTranspose_apply, oc, map_apply, Complex.star_def, Complex.conj_ofReal]
  rw [← transpose_apply R.o, R.o_symm]

/-- The ring Hamiltonian is Hermitian. -/
theorem hamiltonian_isHermitian (L : ℕ) : (R.hamiltonian L).IsHermitian :=
  (isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian j _ R.bondMatrix_isHermitian).add
    (isSelfAdjoint_sum _ fun j _ ↦ siteOp_isHermitian j R.oc_isHermitian)

/-- The spectrum of `H_L` with multiplicity, indexed by configurations. -/
noncomputable def spectrum (L : ℕ) : (Fin L → d) → ℝ :=
  (R.hamiltonian_isHermitian L).eigenvalues

/-- The partition function `Z_L(β) = Re Tr e^{-βH_L}`. -/
noncomputable def partitionFunction (L : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • R.hamiltonian L))).re

/-- `Re Tr e^{-βH_L} = ∑ₘ e^{-β Eₘ}`, the partition function of the abstract criterion. -/
theorem partitionFunction_eq_sum (L : ℕ) (β : ℝ) :
    R.partitionFunction L β = ∑ m, Real.exp (-(β * R.spectrum L m)) := by
  have h := trace_exp_smul_eq_sum_of_isHermitian (R.hamiltonian_isHermitian L) (-β)
  rw [Complex.ofReal_neg] at h
  rw [partitionFunction, h, ← Complex.ofReal_sum, Complex.ofReal_re]
  simp [spectrum, neg_mul]

end RPRing

end CertifiedGaps
