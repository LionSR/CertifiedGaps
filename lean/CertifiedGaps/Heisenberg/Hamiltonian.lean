/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Model

/-!
# The periodic spin-one Heisenberg chain and its partition function

This file defines the antiferromagnetic spin-one Heisenberg Hamiltonian
`H_L = ∑ⱼ Sⱼ · Sⱼ₊₁` on a ring of `L` sites, as an explicit complex matrix indexed by
configurations `Fin L → Fin 3`, and identifies the trace of its Gibbs operator with the
eigenvalue sum `∑ₘ exp (-β Eₘ(L))`. This is the concrete Hamiltonian-side input for the
abstract two-purity criterion of `CertifiedGaps.PurityCriterion`, whose partition function is
exactly such an eigenvalue sum.

## Definitions

* `leviCivita`: the Levi-Civita symbol `ε_{abc}` on `Fin 3`.
* `spinOne α`: the Cartesian spin-one matrices `(S^α)_{s s'} = -i ε_{α s s'}`.
* `exchange`: the two-site exchange `S · S = ∑_α S^α ⊗ S^α`.
* `hamiltonian L`: the periodic chain `∑ⱼ bondOp j (j + 1) exchange`, where the cyclic successor
  is `finRotate L`. For `L = 2` both bonds join sites `0` and `1`, so the Hamiltonian is the
  double-bond ring `2 S₀ · S₁`.
* `spectrum L`: the eigenvalues of `hamiltonian L` with multiplicity, indexed by
  `Fin L → Fin 3`.
* `partitionFunction L β`: the real part of `Tr exp (-β H_L)`.

## Main results

* `spinOne_isHermitian`: the spin matrices are Hermitian.
* `exchange_isHermitian`, `hamiltonian_isHermitian`.
* `trace_exp_hamiltonian`: `Tr exp (-β H_L) = ∑ₘ exp (-(β Eₘ(L)))` as complex numbers.
* `partitionFunction_eq_sum`: `partitionFunction L β = ∑ₘ exp (-(β Eₘ(L)))`, the form of
  `CertifiedGaps.partitionFunction spectrum L β`.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

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
`∑ₘ exp (-(β Eₘ(L)))`, i.e. `CertifiedGaps.partitionFunction spectrum L β`. -/
theorem partitionFunction_eq_sum (L : ℕ) (β : ℝ) :
    partitionFunction L β = ∑ m, Real.exp (-(β * spectrum L m)) := by
  rw [partitionFunction, trace_exp_hamiltonian, Complex.ofReal_re]

end CertifiedGaps.SpinOneHeisenberg
