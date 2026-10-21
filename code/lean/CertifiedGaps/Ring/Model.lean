/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Hamiltonian

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
* `RPRing.hamiltonian L`, `RPRing.spectrum L`, `RPRing.partitionFunction L β`, built from the
  generic `bondOp`, `siteOp` and `gibbsTrace` of `Ring/Hamiltonian.lean`.

## Main results

* `RPRing.hamiltonian_isHermitian`.
* `RPRing.hamiltonian_eq_ringHamiltonian`: `H_L = ringHamiltonian h L + siteSum o L`.
* `RPRing.partitionFunction_eq_gibbsTrace`: `Z_L(β) = gibbsTrace H_L β`.
* `RPRing.partitionFunction_eq_sum`: `Re Tr e^{-βH_L} = ∑ₘ e^{-β Eₘ}`.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace CertifiedGaps

variable {d : Type} [Fintype d] [DecidableEq d]

/-! ### Symmetric and antisymmetric real factors -/

/-- A Kronecker product of symmetric or antisymmetric matrices is symmetric or antisymmetric. -/
theorem transpose_kronecker_eq_or {m n : Type*} {A : Matrix m m ℝ} {B : Matrix n n ℝ}
    (hA : Aᵀ = A ∨ Aᵀ = -A) (hB : Bᵀ = B ∨ Bᵀ = -B) :
    (A ⊗ₖ B)ᵀ = A ⊗ₖ B ∨ (A ⊗ₖ B)ᵀ = -(A ⊗ₖ B) := by
  rw [← kroneckerMap_transpose]
  rcases hA with hA | hA <;> rcases hB with hB | hB <;> rw [hA, hB]
  · exact .inl rfl
  · exact .inr (by ext; simp [kroneckerMap_apply])
  · exact .inr (by ext; simp [kroneckerMap_apply])
  · exact .inl (by ext; simp [kroneckerMap_apply])

/-- A real multiple of a symmetric or antisymmetric matrix is symmetric or antisymmetric. -/
theorem transpose_smul_eq_or {m : Type*} (c : ℝ) {A : Matrix m m ℝ} (hA : Aᵀ = A ∨ Aᵀ = -A) :
    (c • A)ᵀ = c • A ∨ (c • A)ᵀ = -(c • A) := by
  rw [transpose_smul]
  rcases hA with hA | hA <;> rw [hA]
  · exact .inl rfl
  · exact .inr (smul_neg c A)

/-- An exchange `∑ₖ Aₖ ⊗ Aₖ` of symmetric or antisymmetric matrices is symmetric. -/
theorem transpose_sum_kronecker_self {m κ : Type*} [Fintype κ] (A : κ → Matrix m m ℝ)
    (hA : ∀ k, (A k)ᵀ = A k ∨ (A k)ᵀ = -A k) : (∑ k, A k ⊗ₖ A k)ᵀ = ∑ k, A k ⊗ₖ A k := by
  rw [transpose_sum]
  refine Finset.sum_congr rfl fun k _ ↦ ?_
  rw [← kroneckerMap_transpose]
  rcases hA k with h | h <;> rw [h]
  ext; simp [kroneckerMap_apply]

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

/-- The complex bond factors are Hermitian or anti-Hermitian. -/
theorem conjTranspose_Cc (α : ι) : (R.Cc α)ᴴ = R.Cc α ∨ (R.Cc α)ᴴ = -R.Cc α := by
  have h : (R.Cc α)ᴴ = ((R.C α)ᵀ).map (↑) := by
    ext a b; simp [Cc, conjTranspose_apply]
  rcases R.C_symm α with hC | hC
  · exact Or.inl (by rw [h, hC]; rfl)
  · exact Or.inr (by rw [h, hC]; ext; simp [Cc])

/-- The complex reflection is Hermitian. -/
theorem conjTranspose_Rc : R.Rcᴴ = R.Rc := by
  ext a b
  simp only [conjTranspose_apply, Rc, map_apply, Complex.star_def, Complex.conj_ofReal]
  rw [← transpose_apply R.refl, R.refl_symm]

/-- The bond is Hermitian. -/
theorem bondMatrix_isHermitian : R.bondMatrix.IsHermitian := by
  refine IsHermitian.neg (isSelfAdjoint_sum _ fun α _ ↦ ?_)
  change (R.Cc α ⊗ₖ R.Dc α)ᴴ = R.Cc α ⊗ₖ R.Dc α
  rw [conjTranspose_kronecker, Dc, conjTranspose_mul, conjTranspose_mul, conjTranspose_Rc,
    ← Matrix.mul_assoc]
  rcases R.conjTranspose_Cc α with h | h <;> rw [h]
  ext ⟨a, b⟩ ⟨c, e⟩
  simp [kroneckerMap_apply]

/-- The complex on-site term is Hermitian. -/
theorem oc_isHermitian : R.oc.IsHermitian := by
  ext a b
  simp only [conjTranspose_apply, oc, map_apply, Complex.star_def, Complex.conj_ofReal]
  rw [← transpose_apply R.o, R.o_symm]

/-- The ring Hamiltonian is the generic ring of the bond `h` plus the site sum of `o`. -/
theorem hamiltonian_eq_ringHamiltonian (L : ℕ) :
    R.hamiltonian L = ringHamiltonian R.bondMatrix L + siteSum R.oc L := rfl

/-- The ring Hamiltonian is Hermitian. -/
theorem hamiltonian_isHermitian (L : ℕ) : (R.hamiltonian L).IsHermitian := by
  rw [hamiltonian_eq_ringHamiltonian]
  exact (ringHamiltonian_isHermitian R.bondMatrix_isHermitian L).add
    (siteSum_isHermitian R.oc_isHermitian L)

/-- The spectrum of `H_L` with multiplicity, indexed by configurations. -/
noncomputable def spectrum (L : ℕ) : (Fin L → d) → ℝ :=
  (R.hamiltonian_isHermitian L).eigenvalues

/-- The partition function `Z_L(β) = Re Tr e^{-βH_L}`. -/
noncomputable def partitionFunction (L : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • R.hamiltonian L))).re

/-- The partition function is the Gibbs trace of the ring Hamiltonian. -/
theorem partitionFunction_eq_gibbsTrace (L : ℕ) (β : ℝ) :
    R.partitionFunction L β = gibbsTrace (R.hamiltonian L) β := rfl

/-- `Re Tr e^{-βH_L} = ∑ₘ e^{-β Eₘ}`, the partition function of the abstract criterion. -/
theorem partitionFunction_eq_sum (L : ℕ) (β : ℝ) :
    R.partitionFunction L β = ∑ m, Real.exp (-(β * R.spectrum L m)) := by
  rw [partitionFunction_eq_gibbsTrace, gibbsTrace_eq_sum (R.hamiltonian_isHermitian L)]
  rfl

/-- The partition function of a ring with a nonempty local space is positive. -/
theorem partitionFunction_pos [Nonempty d] (L : ℕ) (β : ℝ) : 0 < R.partitionFunction L β :=
  gibbsTrace_pos (R.hamiltonian_isHermitian L) β

end RPRing

end CertifiedGaps
