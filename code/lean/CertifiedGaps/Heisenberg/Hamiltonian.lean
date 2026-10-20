/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Local.SpinOne
import CertifiedGaps.Ring.Hamiltonian

/-!
# The periodic spin-one Heisenberg chain and its partition function

This file defines the antiferromagnetic spin-one Heisenberg Hamiltonian
`H_L = ∑ⱼ Sⱼ · Sⱼ₊₁` on a ring of `L` sites, as the generic ring Hamiltonian of the exchange
`SpinOne.exchange` on configurations `Fin L → Fin 3`, and identifies its Gibbs trace with the
eigenvalue sum `∑ₘ exp (-β Eₘ(L))`. This is the concrete Hamiltonian-side input for the
abstract two-purity criterion of `CertifiedGaps.PurityCriterion`, whose partition function is
exactly such an eigenvalue sum.

## Definitions

* `hamiltonian L = ringHamiltonian SpinOne.exchange L`, the periodic chain
  `∑ⱼ bondOp j (j + 1) (S · S)`, where the cyclic successor is `finRotate L`. For `L = 2` both
  bonds join sites `0` and `1`, so the Hamiltonian is the double-bond ring `2 S₀ · S₁`.
* `spectrum L`: the eigenvalues of `hamiltonian L` with multiplicity, indexed by
  `Fin L → Fin 3`.
* `partitionFunction L β`: the real part of `Tr exp (-β H_L)`, the Gibbs trace of `H_L`
  (`partitionFunction_eq_gibbsTrace`).
* `twistedHamiltonian L`, `twistedPartitionFunction L β`: the ring whose closing bond `L - 1 → 0`
  is twisted by the `π` rotation about `z`, `-SˣSˣ - SʸSʸ + SᶻSᶻ` (the twisted traces `Z_L(β, P)`
  of the seed tables).

## Main results

* `hamiltonian_isHermitian`.
* `trace_exp_hamiltonian`: `Tr exp (-β H_L) = ∑ₘ exp (-(β Eₘ(L)))` as complex numbers.
* `partitionFunction_eq_sum`: `partitionFunction L β = ∑ₘ exp (-(β Eₘ(L)))`, the form of
  `CertifiedGaps.partitionFunction spectrum L β`.
-/

open Matrix NormedSpace

namespace CertifiedGaps.SpinOneHeisenberg

/-! ### The periodic Hamiltonian -/

/-- The periodic spin-one Heisenberg Hamiltonian `H_L = ∑ⱼ Sⱼ · Sⱼ₊₁` on `L` sites, with the
cyclic successor `finRotate L j = j + 1` (`finRotate_succ_apply`). For `L = 2` it is the
double-bond ring `2 S₀ · S₁`. -/
def hamiltonian (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  ringHamiltonian SpinOne.exchange L

/-- The periodic Heisenberg Hamiltonian is Hermitian. -/
theorem hamiltonian_isHermitian (L : ℕ) : (hamiltonian L).IsHermitian :=
  ringHamiltonian_isHermitian SpinOne.exchange_isHermitian L

/-- The spectrum of `hamiltonian L`: its eigenvalues with multiplicity, indexed by
configurations. -/
noncomputable def spectrum (L : ℕ) : (Fin L → Fin 3) → ℝ :=
  (hamiltonian_isHermitian L).eigenvalues

/-- The partition function `Z_L(β) = Tr exp (-β H_L)`, taken as the real part of the complex
trace. -/
noncomputable def partitionFunction (L : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • hamiltonian L))).re

/-- The partition function is the Gibbs trace of the ring Hamiltonian. -/
theorem partitionFunction_eq_gibbsTrace (L : ℕ) (β : ℝ) :
    partitionFunction L β = gibbsTrace (hamiltonian L) β := rfl

/-! ### The trace of the Gibbs operator -/

/-- The partition function of the Heisenberg chain is the eigenvalue sum
`∑ₘ exp (-(β Eₘ(L)))`, i.e. `CertifiedGaps.partitionFunction spectrum L β`. -/
theorem partitionFunction_eq_sum (L : ℕ) (β : ℝ) :
    partitionFunction L β = ∑ m, Real.exp (-(β * spectrum L m)) :=
  gibbsTrace_eq_sum (hamiltonian_isHermitian L) β

/-- The trace of the Gibbs operator of the Heisenberg chain is the eigenvalue sum
`∑ₘ exp (-(β Eₘ(L)))`. -/
theorem trace_exp_hamiltonian (L : ℕ) (β : ℝ) :
    trace (exp ((-(β : ℂ)) • hamiltonian L)) =
      ((∑ m, Real.exp (-(β * spectrum L m)) : ℝ) : ℂ) := by
  have h := trace_exp_smul_eq_sum_of_isHermitian (hamiltonian_isHermitian L) (-β)
  rw [Complex.ofReal_neg] at h
  rw [h, Complex.ofReal_sum]
  simp [spectrum, neg_mul]

/-! ### The twisted ring -/

/-- The signs of the `π` rotation about `z` on the Cartesian spin components. -/
def twistSign : Fin 3 → ℂ := ![-1, -1, 1]

/-- The exchange across the twisted bond, `-SˣSˣ - SʸSʸ + SᶻSᶻ`. -/
def twistedExchange : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ :=
  ∑ α, twistSign α • kroneckerMap (· * ·) (SpinOne.S α) (SpinOne.S α)

/-- The twisted exchange is Hermitian. -/
theorem twistedExchange_isHermitian : twistedExchange.IsHermitian :=
  isSelfAdjoint_sum _ fun α _ ↦ by
    change (twistSign α • kroneckerMap (· * ·) (SpinOne.S α) (SpinOne.S α))ᴴ = _
    rw [conjTranspose_smul, conjTranspose_kronecker, (SpinOne.S_isHermitian α).eq]
    congr 1
    fin_cases α <;> simp [twistSign]

/-- The Heisenberg ring with its closing bond `L - 1 → 0` twisted by the `π` rotation about
`z`. -/
def twistedHamiltonian (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  ∑ j : Fin L, bondOp j (finRotate L j)
    (if j.val + 1 = L then twistedExchange else SpinOne.exchange)

/-- The twisted partition function `Z_L(β, P) = Re Tr exp (-β H_L(P))`. -/
noncomputable def twistedPartitionFunction (L : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • twistedHamiltonian L))).re


end CertifiedGaps.SpinOneHeisenberg
