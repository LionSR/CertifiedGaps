/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Model

/-!
# Open chains of a reflection-positive ring

The open chain of `n` sites cut from a reflection-positive ring (`RPRing`) has the Hamiltonian

  `H_n = ∑_{j=0}^{n-2} hⱼ,ⱼ₊₁ + ∑ᵢ oᵢ + b₀ + b'_{n-1}`,

with the ring's bond `h` and on-site term `o`, and an end term at each end. The end terms are
built from the bond factors, `b = -∑_α w_α C_α`, with weights `w_α` real or imaginary according
as `C_α` is symmetric or antisymmetric, so that `b` is Hermitian. The right end term `b'` is `b`
on chains of odd length and `bᵀ` on chains of even length: this is the choice for which the
open-chain trace is a matrix element `⟨v, Xⁿ v⟩` of the ring's transfer operator `X` with one
boundary vector `v` (`CertifiedGaps.Open.Transfer`). For the spin-1 Heisenberg chain with
`w_z = -i h`, `b = -h S^z`, and the chain carries equal end fields at odd length and opposite
ones at even length.

## Definitions

* `RPRing.Boundary R`: the boundary weights `w`.
* `RPRing.Boundary.mat`, `RPRing.Boundary.endMat`: the end terms `b` and `b'`.
* `RPRing.openHamiltonian b n`, `RPRing.openSpectrum b n`, `RPRing.openPartitionFunction b n β`.

## Main results

* `RPRing.openHamiltonian_isHermitian`.
* `RPRing.openPartitionFunction_eq_sum`: `Re Tr e^{-βH_n} = ∑ₘ e^{-β Eₘ}`.
-/

open Matrix NormedSpace

namespace CertifiedGaps

namespace RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- Boundary weights `w` of an end term `b = -∑_α w_α C_α`; each weight is real or imaginary
according as `C_α` is symmetric or antisymmetric, which makes `b` Hermitian. -/
structure Boundary (R : RPRing d ι) where
  /-- The weight of each bond factor in the end term. -/
  w : ι → ℂ
  /-- A weight is real on a symmetric factor and imaginary on an antisymmetric one. -/
  w_symm : ∀ α, ((R.C α)ᵀ = R.C α ∧ star (w α) = w α) ∨
    ((R.C α)ᵀ = -R.C α ∧ star (w α) = -w α)

namespace Boundary

variable {R}

/-- The end term `b = -∑_α w_α C_α`. -/
def mat (b : R.Boundary) : Matrix d d ℂ := -∑ α, b.w α • R.Cc α

/-- The end term is Hermitian. -/
theorem mat_isHermitian (b : R.Boundary) : b.mat.IsHermitian := by
  refine IsHermitian.neg (isSelfAdjoint_sum _ fun α _ ↦ ?_)
  change (b.w α • R.Cc α)ᴴ = b.w α • R.Cc α
  have h : (R.Cc α)ᴴ = ((R.C α)ᵀ).map (↑) := by
    ext a c; simp [Cc, conjTranspose_apply]
  rw [conjTranspose_smul, h]
  rcases b.w_symm α with ⟨hC, hw⟩ | ⟨hC, hw⟩ <;> rw [hC, hw]
  · rfl
  · ext a c; simp [Cc]

/-- The transposed end term is Hermitian. -/
theorem transpose_mat_isHermitian (b : R.Boundary) : b.matᵀ.IsHermitian :=
  b.mat_isHermitian.transpose

/-- The end term at the last site of a chain of `n` sites: `b` if `n` is odd, `bᵀ` if `n` is
even. -/
def endMat (b : R.Boundary) (n : ℕ) : Matrix d d ℂ := if Even n then b.matᵀ else b.mat

/-- The end term at the last site is Hermitian. -/
theorem endMat_isHermitian (b : R.Boundary) (n : ℕ) : (b.endMat n).IsHermitian := by
  unfold endMat
  split_ifs
  · exact b.transpose_mat_isHermitian
  · exact b.mat_isHermitian

/-- If `bᵀ = b`, the end term at the last site is `b` at every length. -/
theorem endMat_of_transpose_eq (b : R.Boundary) (h : b.matᵀ = b.mat) (n : ℕ) :
    b.endMat n = b.mat := by
  unfold endMat; split_ifs <;> simp [h]

/-- Free ends, `b = 0`: the zero weights, for a ring whose bond factors are each symmetric or
antisymmetric. Its open chains are the chains with open boundary conditions. -/
def free (hC : ∀ α, (R.C α)ᵀ = R.C α ∨ (R.C α)ᵀ = -R.C α) : R.Boundary where
  w _ := 0
  w_symm α := (hC α).imp (fun h ↦ ⟨h, by simp⟩) (fun h ↦ ⟨h, by simp⟩)

/-- The free end term is zero. -/
theorem free_mat (hC : ∀ α, (R.C α)ᵀ = R.C α ∨ (R.C α)ᵀ = -R.C α) : (free hC).mat = 0 := by
  simp [mat, free]

/-- The free end terms vanish at every length. -/
theorem free_endMat (hC : ∀ α, (R.C α)ᵀ = R.C α ∨ (R.C α)ᵀ = -R.C α) (n : ℕ) :
    (free hC).endMat n = 0 := by
  rw [endMat_of_transpose_eq _ (by rw [free_mat, transpose_zero]), free_mat]

end Boundary

/-- The open chain Hamiltonian on `n` sites `0, …, n - 1`: the bonds `(j, j + 1)`, the on-site
terms, the end term `b` at site `0`, and the end term `b.endMat n` at site `n - 1`. The empty
chain has `H_0 = 0`. -/
def openHamiltonian (b : R.Boundary) : (n : ℕ) → Matrix (Fin n → d) (Fin n → d) ℂ
  | 0 => 0
  | n + 1 => ∑ j : Fin n, bondOp j.castSucc j.succ R.bondMatrix + ∑ i, siteOp i R.oc +
      siteOp 0 b.mat + siteOp (Fin.last n) (b.endMat (n + 1))

/-- The open chain Hamiltonian is Hermitian. -/
theorem openHamiltonian_isHermitian (b : R.Boundary) (n : ℕ) :
    (R.openHamiltonian b n).IsHermitian := by
  cases n with
  | zero => exact isHermitian_zero
  | succ n =>
    exact (((isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian _ _ R.bondMatrix_isHermitian).add
      (isSelfAdjoint_sum _ fun i _ ↦ siteOp_isHermitian i R.oc_isHermitian)).add
      (siteOp_isHermitian _ b.mat_isHermitian)).add
      (siteOp_isHermitian _ (b.endMat_isHermitian _))

/-- With free ends the open chain is the chain with open boundary conditions: the bonds and the
on-site terms only. -/
theorem openHamiltonian_free (hC : ∀ α, (R.C α)ᵀ = R.C α ∨ (R.C α)ᵀ = -R.C α) (n : ℕ) :
    R.openHamiltonian (Boundary.free hC) (n + 1) =
      ∑ j : Fin n, bondOp j.castSucc j.succ R.bondMatrix + ∑ i, siteOp i R.oc := by
  simp only [openHamiltonian, Boundary.free_mat, Boundary.free_endMat, siteOp_zero, add_zero]

/-- The spectrum of `H_n` with multiplicity, indexed by configurations. -/
noncomputable def openSpectrum (b : R.Boundary) (n : ℕ) : (Fin n → d) → ℝ :=
  (R.openHamiltonian_isHermitian b n).eigenvalues

/-- The open-chain partition function `G_n(β) = Re Tr e^{-βH_n}`. -/
noncomputable def openPartitionFunction (b : R.Boundary) (n : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • R.openHamiltonian b n))).re

/-- `Re Tr e^{-βH_n} = ∑ₘ e^{-β Eₘ}` for the open chain. -/
theorem openPartitionFunction_eq_sum (b : R.Boundary) (n : ℕ) (β : ℝ) :
    R.openPartitionFunction b n β = ∑ m, Real.exp (-(β * R.openSpectrum b n m)) := by
  have h := trace_exp_smul_eq_sum_of_isHermitian (R.openHamiltonian_isHermitian b n) (-β)
  rw [Complex.ofReal_neg] at h
  rw [openPartitionFunction, h, ← Complex.ofReal_sum, Complex.ofReal_re]
  simp [openSpectrum, neg_mul]

end RPRing

end CertifiedGaps
