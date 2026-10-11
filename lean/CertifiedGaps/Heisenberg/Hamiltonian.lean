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

/-!
# The periodic spin-one Heisenberg chain and its partition function

This file defines the antiferromagnetic spin-one Heisenberg Hamiltonian
`H_L = ∑ⱼ Sⱼ · Sⱼ₊₁` on a ring of `L` sites, as an explicit complex matrix indexed by
configurations `Fin L → Fin 3`, and identifies the trace of its Gibbs operator with the
eigenvalue sum `∑ₘ exp (-β Eₘ(L))`. This is the concrete Hamiltonian-side input for the
abstract two-purity criterion of `ThermalMPOGap.PurityCriterion`, whose partition function is
exactly such an eigenvalue sum.

## Definitions

* `leviCivita`: the Levi-Civita symbol `ε_{abc}` on `Fin 3`.
* `spinOne α`: the Cartesian spin-one matrices `(S^α)_{s s'} = -i ε_{α s s'}`.
* `bondOp j k B`: a two-site operator `B` acting on sites `j` and `k` of the chain, and as the
  identity elsewhere.
* `exchange`: the two-site exchange `S · S = ∑_α S^α ⊗ S^α`.
* `hamiltonian L`: the periodic chain `∑ⱼ bondOp j (j + 1) exchange`, where the cyclic successor
  is `finRotate L`. For `L = 2` both bonds join sites `0` and `1`, so the Hamiltonian is the
  double-bond ring `2 S₁ · S₂`; for `L = 1` it is the single self-bond `S₀ · S₀ = 2`.
* `spectrum L`: the eigenvalues of `hamiltonian L` with multiplicity, indexed by
  `Fin L → Fin 3`.
* `partitionFunction L β`: the real part of `Tr exp (-β H_L)`.

## Main results

* `spinOne_isHermitian`, `sum_spinOne_mul_self`: the spin matrices are Hermitian and
  `∑_α (S^α)² = 2`, the Casimir value `s (s + 1)` for `s = 1`.
* `bondOp_isHermitian`, `exchange_isHermitian`, `hamiltonian_isHermitian`.
* `trace_exp_smul_eq_sum_of_isHermitian`: for a Hermitian matrix `A` and real `c`,
  `Tr exp (c A) = ∑ᵢ exp (c λᵢ)`.
* `trace_exp_hamiltonian`: `Tr exp (-β H_L) = ∑ₘ exp (-(β Eₘ(L)))` as complex numbers.
* `partitionFunction_eq_sum`: `partitionFunction L β = ∑ₘ exp (-(β Eₘ(L)))`, the form of
  `ThermalPurityGap.partitionFunction spectrum L β`.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace SpinOneHeisenberg

/-! ### Spin-one matrices -/

/-- The Levi-Civita symbol on `Fin 3`: `ε₀₁₂ = 1`, totally antisymmetric. -/
def leviCivita : Fin 3 → Fin 3 → Fin 3 → ℤ :=
  ![![![0, 0, 0], ![0, 0, 1], ![0, -1, 0]],
    ![![0, 0, -1], ![0, 0, 0], ![1, 0, 0]],
    ![![0, 1, 0], ![-1, 0, 0], ![0, 0, 0]]]

/-- The Levi-Civita symbol is antisymmetric in its last two indices. -/
theorem leviCivita_swap (a b c : Fin 3) : leviCivita a c b = -leviCivita a b c := by
  fin_cases a <;> fin_cases b <;> fin_cases c <;> rfl

/-- The Cartesian spin-one matrices, `(S^α)_{s s'} = -i ε_{α s s'}`. -/
def spinOne (α : Fin 3) : Matrix (Fin 3) (Fin 3) ℂ :=
  of fun s s' ↦ -Complex.I * (leviCivita α s s' : ℂ)

/-- Each spin-one matrix is Hermitian. -/
theorem spinOne_isHermitian (α : Fin 3) : (spinOne α).IsHermitian := by
  ext s s'
  simp only [conjTranspose_apply, spinOne, of_apply, leviCivita_swap α s s', star_mul',
    star_neg, Complex.star_def, Complex.conj_I, Int.cast_neg, ← Complex.ofReal_intCast,
    Complex.conj_ofReal]
  push_cast
  ring

/-- The spin-one Casimir: `∑_α (S^α)² = s (s + 1) = 2`. -/
theorem sum_spinOne_mul_self :
    ∑ α, spinOne α * spinOne α = (2 : ℕ) • (1 : Matrix (Fin 3) (Fin 3) ℂ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Fin.sum_univ_three, mul_apply, spinOne, leviCivita, one_apply] <;> norm_num

/-! ### Two-site operators on the chain -/

/-- The operator `B` on two sites `j, k` of a chain of `L` spin-one sites, acting as the identity
on all other sites: its matrix element between configurations `σ` and `τ` is
`B (σ j, σ k) (τ j, τ k)` when `σ` and `τ` agree away from `j` and `k`, and `0` otherwise. -/
def bondOp {L : ℕ} (j k : Fin L) (B : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ) :
    Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  of fun σ τ ↦ if ∀ i, i ≠ j → i ≠ k → σ i = τ i then B (σ j, σ k) (τ j, τ k) else 0

/-- A Hermitian two-site operator embeds as a Hermitian operator on the chain, for any pair of
sites `j, k` (including `j = k`). -/
theorem bondOp_isHermitian {L : ℕ} (j k : Fin L) {B : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ}
    (hB : B.IsHermitian) : (bondOp j k B).IsHermitian := by
  ext σ τ
  simp only [conjTranspose_apply, bondOp, of_apply]
  split_ifs with h₁ h₂ h₂
  · exact hB.apply _ _
  · exact (h₂ fun i hj hk ↦ (h₁ i hj hk).symm).elim
  · exact (h₁ fun i hj hk ↦ (h₂ i hj hk).symm).elim
  · exact star_zero _

/-- The two-site exchange interaction `S · S = ∑_α S^α ⊗ S^α`. -/
def exchange : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ :=
  ∑ α, kroneckerMap (· * ·) (spinOne α) (spinOne α)

/-- The exchange interaction is Hermitian. -/
theorem exchange_isHermitian : exchange.IsHermitian :=
  isSelfAdjoint_sum _ fun α _ ↦ by
    change (spinOne α ⊗ₖ spinOne α)ᴴ = spinOne α ⊗ₖ spinOne α
    rw [conjTranspose_kronecker, (spinOne_isHermitian α).eq]

/-! ### The periodic Hamiltonian -/

/-- The periodic spin-one Heisenberg Hamiltonian `H_L = ∑ⱼ Sⱼ · Sⱼ₊₁` on `L` sites, with the
cyclic successor `finRotate L j = j + 1` (`finRotate_succ_apply`). For `L = 2` it is the
double-bond ring `2 S₀ · S₁`. -/
def hamiltonian (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  ∑ j : Fin L, bondOp j (finRotate L j) exchange

/-- The periodic Heisenberg Hamiltonian is Hermitian. -/
theorem hamiltonian_isHermitian (L : ℕ) : (hamiltonian L).IsHermitian :=
  isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian j (finRotate L j) exchange_isHermitian

/-- The spectrum of `hamiltonian L`: its eigenvalues with multiplicity, indexed by
configurations. -/
noncomputable def spectrum (L : ℕ) : (Fin L → Fin 3) → ℝ :=
  (hamiltonian_isHermitian L).eigenvalues

/-- The partition function `Z_L(β) = Tr exp (-β H_L)`, taken as the real part of the complex
trace. -/
noncomputable def partitionFunction (L : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • hamiltonian L))).re

/-! ### The trace of the Gibbs operator -/

/-- For a Hermitian matrix `A` and a real scalar `c`, `Tr exp (c A) = ∑ᵢ exp (c λᵢ)`, where `λᵢ`
are the eigenvalues of `A`. -/
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

/-- The trace of the Gibbs operator of the Heisenberg chain is the eigenvalue sum
`∑ₘ exp (-(β Eₘ(L)))`. -/
theorem trace_exp_hamiltonian (L : ℕ) (β : ℝ) :
    trace (exp ((-(β : ℂ)) • hamiltonian L)) =
      ((∑ m, Real.exp (-(β * spectrum L m)) : ℝ) : ℂ) := by
  have h := trace_exp_smul_eq_sum_of_isHermitian (hamiltonian_isHermitian L) (-β)
  rw [Complex.ofReal_neg] at h
  rw [h, Complex.ofReal_sum]
  simp [spectrum, neg_mul]

/-- The partition function of the Heisenberg chain is the eigenvalue sum
`∑ₘ exp (-(β Eₘ(L)))`, i.e. `ThermalPurityGap.partitionFunction spectrum L β`. -/
theorem partitionFunction_eq_sum (L : ℕ) (β : ℝ) :
    partitionFunction L β = ∑ m, Real.exp (-(β * spectrum L m)) := by
  rw [partitionFunction, trace_exp_hamiltonian, Complex.ofReal_re]

end SpinOneHeisenberg
