/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Hamiltonian
import CertifiedGaps.Ring.Conjugation

/-!
# From the Cartesian basis to the `S^z` basis

The Heisenberg Hamiltonian of `CertifiedGaps.SpinOneHeisenberg.hamiltonian` is written in the
Cartesian basis `|x⟩, |y⟩, |z⟩`, where the reflection-positivity factorization is real. The
seed computations (trace table, trial state) work in the `S^z` basis, where the exchange is the
integer matrix `exchangeZ`. This file proves that the two are unitarily equivalent site by site
(the one-site unitary `SpinOne.szBasis`), so both bases give the same partition functions, with
and without the twist.

Both rings are instances of the generic `twistedRingHamiltonian` (`twistedHamiltonian_eq`,
`twistedHamiltonianZ`), so the ring-level step is the single generic lemma
`conj_twistedRingHamiltonian`; this file supplies only the two-site identities.

## Main definitions

* `twistedExchangeZ`: the twisted exchange `-SˣSˣ - SʸSʸ + SᶻSᶻ` in the `S^z` basis, the integer
  matrix `SᶻSᶻ - (S⁺S⁻ + S⁻S⁺)/2`; `twistedHamiltonianZ L`.
* `szFrame L = szBasis^{⊗L}` on the ring.

## Main results

* `conj_exchange`: `(U ⊗ U)ᴴ (S · S) (U ⊗ U) = exchangeZ`, and `conj_twistedExchange`.
* `conj_hamiltonian`, `conj_twistedHamiltonian`: the ring Hamiltonians, for `L ≠ 1`.
* `partitionFunction_eq_partitionFunctionZ`: `Z_L(β)` in either basis.
* `twistedPartitionFunction_eq_gibbsTrace_twistedHamiltonianZ`: `Z_L(β, P)` in the `S^z`
  basis.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

open SpinOne

/-! ### The `S^z` basis -/

/-- The magnetic quantum number `m = i - 1` of the basis label `i : Fin 3`. -/
def mz (i : Fin 3) : ℤ := (i.val : ℤ) - 1

/-- The exchange `S · S = SᶻSᶻ + (S⁺S⁻ + S⁻S⁺)/2` in the `S^z` basis, with rows the output
pair and columns the input pair. All entries are integers. -/
def exchangeZ : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℤ := fun p q ↦
  (if p = q then mz p.1 * mz p.2 else 0) +
    (if p.1.val = q.1.val + 1 ∧ p.2.val + 1 = q.2.val then 1 else 0) +
    (if p.1.val + 1 = q.1.val ∧ p.2.val = q.2.val + 1 then 1 else 0)

/-- The Heisenberg ring in the `S^z` basis. -/
def hamiltonianZ (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  ∑ j : Fin L, bondOp j (finRotate L j) (exchangeZ.map (↑))

/-- The partition function `Z_L(β) = Re Tr exp (-β H_L)` in the `S^z` basis. -/
noncomputable def partitionFunctionZ (L : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • hamiltonianZ L))).re

/-! ### The conventions as generic rings -/

/-- The twist signs of `twistedExchange` are those of the `π` rotation about `z`. -/
theorem twistSign_eq (α : Fin 3) : twistSign α = (rotZ α : ℂ) := by
  fin_cases α <;> simp [twistSign, rotZ]

/-- The twisted Heisenberg ring is the generic twisted ring. -/
theorem twistedHamiltonian_eq (L : ℕ) :
    twistedHamiltonian L = twistedRingHamiltonian exchange twistedExchange L := rfl

/-- The twisted partition function is a Gibbs trace. -/
theorem twistedPartitionFunction_eq (L : ℕ) (β : ℝ) :
    twistedPartitionFunction L β = gibbsTrace (twistedHamiltonian L) β := rfl

/-! ### The two-site identities -/

/-- The exchange built from `spinZ` is the integer matrix `exchangeZ` (81 entries). -/
theorem sum_spinZ_kronecker : ∑ α, spinZ α ⊗ₖ spinZ α = exchangeZ.map (↑) := by
  ext ⟨a, b⟩ ⟨c, d⟩
  simp only [Fin.sum_univ_three, Matrix.add_apply, kroneckerMap_apply, map_apply]
  fin_cases a <;> fin_cases b <;> fin_cases c <;> fin_cases d <;>
    simp [spinZ, exchangeZ, mz] <;> ring_nf <;> simp [invSqrtTwo_sq] <;> norm_num

/-- The two-site identity `(U ⊗ U)ᴴ (S · S) (U ⊗ U) = exchangeZ`. -/
theorem conj_exchange :
    (szBasis ⊗ₖ szBasis)ᴴ * exchange * (szBasis ⊗ₖ szBasis) = exchangeZ.map (↑) := by
  rw [exchange, Finset.mul_sum, Finset.sum_mul, ← sum_spinZ_kronecker]
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [conjTranspose_kronecker, ← mul_kronecker_mul, ← mul_kronecker_mul, conj_S]

/-- The twisted exchange `-SˣSˣ - SʸSʸ + SᶻSᶻ = SᶻSᶻ - (S⁺S⁻ + S⁻S⁺)/2` in the `S^z` basis:
`exchangeZ` with the off-diagonal (spin-flip) entries negated. -/
def twistedExchangeZ : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℤ := fun p q ↦
  (if p = q then mz p.1 * mz p.2 else 0) -
    (if p.1.val = q.1.val + 1 ∧ p.2.val + 1 = q.2.val then 1 else 0) -
    (if p.1.val + 1 = q.1.val ∧ p.2.val = q.2.val + 1 then 1 else 0)

/-- The twisted exchange built from `spinZ` is the integer matrix `twistedExchangeZ`. -/
theorem sum_twistSign_spinZ_kronecker :
    ∑ α, twistSign α • spinZ α ⊗ₖ spinZ α = twistedExchangeZ.map (↑) := by
  ext ⟨a, b⟩ ⟨c, d⟩
  simp only [Fin.sum_univ_three, Matrix.add_apply, Matrix.smul_apply, kroneckerMap_apply,
    map_apply]
  fin_cases a <;> fin_cases b <;> fin_cases c <;> fin_cases d <;>
    simp [spinZ, twistSign, twistedExchangeZ, mz] <;> ring_nf <;> simp [invSqrtTwo_sq] <;>
    norm_num

/-- The two-site identity for the twisted bond. -/
theorem conj_twistedExchange :
    (szBasis ⊗ₖ szBasis)ᴴ * twistedExchange * (szBasis ⊗ₖ szBasis) =
      twistedExchangeZ.map (↑) := by
  rw [twistedExchange, Finset.mul_sum, Finset.sum_mul, ← sum_twistSign_spinZ_kronecker]
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [Matrix.mul_smul, Matrix.smul_mul, conjTranspose_kronecker, ← mul_kronecker_mul,
    ← mul_kronecker_mul, conj_S]

/-! ### The ring -/

/-- The twisted Heisenberg ring in the `S^z` basis: closing bond `L - 1 → 0` carries
`twistedExchangeZ`. -/
def twistedHamiltonianZ (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  twistedRingHamiltonian (exchangeZ.map (↑)) (twistedExchangeZ.map (↑)) L

/-- The change of basis `szBasis^{⊗L}` on a ring of `L` sites. -/
noncomputable def szFrame (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  piTensor fun _ ↦ szBasis

/-- `szFrame L` is unitary. -/
theorem szFrame_conjTranspose_mul_self (L : ℕ) : (szFrame L)ᴴ * szFrame L = 1 :=
  piTensor_const_conjTranspose_mul_self _ szBasis_conjTranspose_mul_self

/-- The Heisenberg ring in the `S^z` basis is `hamiltonianZ`. A ring of one site is excluded
because its bond joins the site to itself. -/
theorem conj_hamiltonian {L : ℕ} (hL : L ≠ 1) :
    (szFrame L)ᴴ * hamiltonian L * szFrame L = hamiltonianZ L := by
  have h := conj_ringHamiltonian szBasis szBasis_conjTranspose_mul_self hL exchange
  rwa [conj_exchange] at h

/-- The twisted Heisenberg ring in the `S^z` basis is `twistedHamiltonianZ`. -/
theorem conj_twistedHamiltonian {L : ℕ} (hL : L ≠ 1) :
    (szFrame L)ᴴ * twistedHamiltonian L * szFrame L = twistedHamiltonianZ L :=
  conj_twistedRingHamiltonian szBasis szBasis_conjTranspose_mul_self hL conj_exchange
    conj_twistedExchange

/-- The partition function is the same in the Cartesian and the `S^z` basis. -/
theorem partitionFunction_eq_partitionFunctionZ {L : ℕ} (hL : L ≠ 1) (β : ℝ) :
    partitionFunction L β = partitionFunctionZ L β :=
  (gibbsTrace_ringHamiltonian_conj szBasis szBasis_conjTranspose_mul_self hL conj_exchange
    β).symm

/-- The twisted partition function `Z_L(β, P)` computed in the `S^z` basis. -/
theorem twistedPartitionFunction_eq_gibbsTrace_twistedHamiltonianZ {L : ℕ} (hL : L ≠ 1)
    (β : ℝ) : twistedPartitionFunction L β = gibbsTrace (twistedHamiltonianZ L) β :=
  (gibbsTrace_twistedRingHamiltonian_conj szBasis szBasis_conjTranspose_mul_self hL
    conj_exchange conj_twistedExchange β).symm

end CertifiedGaps.SpinOneHeisenberg
