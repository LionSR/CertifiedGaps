/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.MatrixProduct
import CertifiedGaps.Heisenberg.SzBridge

/-!
# The matrix-product trial bound on the spin-one partition function

A matrix-product trial state `ψ(σ) = Tr (A σ₀ ⋯ A σ_{N-1})` on the spin-one Heisenberg ring
bounds the partition function from below: if its energy density lies below `-κ`, then
`Z_N(β) > e^{βκN}`. This file instantiates the generic bound of
`CertifiedGaps/Transfer/MatrixProduct.lean` with the `S^z`-basis exchange `exchangeZ`: the
spin-one definitions `mpsState`, `mpsTransfer`, `mpsInsertion` (defined here) are the generic
ones by `rfl`
(`mpsTransfer_eq`, `mpsInsertion_eq`, `hamiltonianZ_eq`, `partitionFunctionZ_eq`), and the only
model-specific fact is that `exchangeZ` is symmetric.

## Main results

* `hamiltonianZ_isHermitian`.
* `sum_mpsState_sq`: `∑_σ ψ(σ)² = Tr Tᴺ`.
* `dotProduct_hamiltonianZ_mulVec`: `⟨ψ, H_N ψ⟩ = N Tr (T^(N-2) O)` for `N ≥ 2`.
* `exp_lt_partitionFunctionZ_of_mps`: `Tr (T^(N-2) O) < -κ Tr Tᴺ → e^{βκN} < Z_N(β)`, with the
  integer forms `exp_lt_partitionFunctionZ_of_mps_int` and
  `exp_lt_partitionFunctionZ_of_mps_700741` (`κ = 700741/500000`).
-/

open Matrix

namespace CertifiedGaps.SpinOneHeisenberg

section Definitions

open scoped Kronecker

/-! ### Matrix-product trial states -/

variable {D : Type*} [Fintype D] [DecidableEq D]

/-- The periodic matrix-product state `ψ(σ) = Tr (A σ₀ ⋯ A σ_{N-1})` in the `S^z` basis. -/
def mpsState {N : ℕ} (A : Fin 3 → Matrix D D ℤ) (σ : Fin N → Fin 3) : ℤ :=
  trace (List.ofFn fun j ↦ A (σ j)).prod

/-- The transfer matrix `T = ∑ₛ Aₛ ⊗ Aₛ` of a real matrix-product state. -/
def mpsTransfer (A : Fin 3 → Matrix D D ℤ) : Matrix (D × D) (D × D) ℤ :=
  ∑ s, A s ⊗ₖ A s

/-- The two-site energy insertion `O = ∑ h(st; s't') (Aₛ Aₜ) ⊗ (Aₛ' Aₜ')`, with `h` the
`S^z`-basis exchange `exchangeZ`. -/
def mpsInsertion (A : Fin 3 → Matrix D D ℤ) : Matrix (D × D) (D × D) ℤ :=
  ∑ p : Fin 3 × Fin 3, ∑ q : Fin 3 × Fin 3,
    exchangeZ p q • ((A p.1 * A p.2) ⊗ₖ (A q.1 * A q.2))


end Definitions

/-! ### The spin-one definitions as instances of the generic ones -/

section Bridge

variable {D : Type*} [Fintype D] [DecidableEq D] (A : Fin 3 → Matrix D D ℤ)

/-- The spin-one trial state is the generic matrix-product state. -/
theorem mpsState_eq {N : ℕ} (σ : Fin N → Fin 3) : mpsState A σ = MPS.state A σ := rfl

omit [Fintype D] [DecidableEq D] in
/-- The spin-one transfer matrix is the generic one. -/
theorem mpsTransfer_eq : mpsTransfer A = MPS.transfer A := rfl

omit [DecidableEq D] in
/-- The spin-one energy insertion is the generic insertion of the bond `exchangeZ`. -/
theorem mpsInsertion_eq : mpsInsertion A = MPS.insertion exchangeZ A := rfl

/-- The trial state as a complex vector on configurations. -/
def mpsVec (N : ℕ) : (Fin N → Fin 3) → ℂ := MPS.vec A N

end Bridge

/-- The Heisenberg ring in the `S^z` basis is the ring Hamiltonian of `exchangeZ`. -/
theorem hamiltonianZ_eq (L : ℕ) : hamiltonianZ L = ringHamiltonian (exchangeZ.map (↑)) L := rfl

/-- The `S^z`-basis partition function is the Gibbs trace of the ring of `exchangeZ`. -/
theorem partitionFunctionZ_eq (L : ℕ) (β : ℝ) :
    partitionFunctionZ L β = gibbsTrace (ringHamiltonian (exchangeZ.map (↑)) L) β := rfl

/-! ### The Hamiltonian in the `S^z` basis -/

/-- The `S^z`-basis exchange is a symmetric integer matrix. -/
theorem exchangeZ_transpose : exchangeZᵀ = exchangeZ := by
  have h : ∀ p q : Fin 3 × Fin 3, exchangeZ q p = exchangeZ p q := by
    intro ⟨a, b⟩ ⟨c, d⟩
    fin_cases a <;> fin_cases b <;> fin_cases c <;> fin_cases d <;> rfl
  ext p q
  exact h p q

/-- The `S^z`-basis exchange is Hermitian. -/
theorem exchangeZ_isHermitian : (exchangeZ.map (Int.cast : ℤ → ℂ)).IsHermitian :=
  isHermitian_map_intCast exchangeZ_transpose

/-- The Heisenberg ring in the `S^z` basis is Hermitian. -/
theorem hamiltonianZ_isHermitian (L : ℕ) : (hamiltonianZ L).IsHermitian :=
  ringHamiltonian_isHermitian exchangeZ_isHermitian L

/-! ### The norm and the energy of the trial state -/

section TrialState

variable {D : Type*} [Fintype D] [DecidableEq D] (A : Fin 3 → Matrix D D ℤ)

/-- **Norm of the trial state.** `∑_σ ψ(σ)² = Tr Tᴺ`. -/
theorem sum_mpsState_sq (N : ℕ) :
    ∑ σ : Fin N → Fin 3, mpsState A σ ^ 2 = trace (mpsTransfer A ^ N) :=
  MPS.sum_state_sq A N

/-- **Energy of the trial state.** On a ring of `N ≥ 2` sites,
`⟨ψ, H ψ⟩ = N Tr (T^(N-2) O)`: every bond contributes the same by cyclicity of the trace. -/
theorem dotProduct_hamiltonianZ_mulVec {N : ℕ} (hN : 2 ≤ N) :
    star (mpsVec A N) ⬝ᵥ (hamiltonianZ N *ᵥ mpsVec A N) =
      (N : ℂ) * ((trace (mpsTransfer A ^ (N - 2) * mpsInsertion A) : ℤ) : ℂ) :=
  MPS.dotProduct_ringHamiltonian_mulVec A exchangeZ hN

end TrialState

/-! ### The trial-state bound on the partition function -/

/-- **Trial-state bound.** If the matrix-product trial state on a ring of `N ≥ 3` sites has
energy density below `-κ`, i.e. `Tr (T^(N-2) O) < -κ Tr Tᴺ`, then `Z_N(β) > e^{βκN}` for every
`β > 0`. -/
theorem exp_lt_partitionFunctionZ_of_mps {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ) {N : ℕ} (hN : 3 ≤ N) (κ : ℚ)
    (hU : ((trace (mpsTransfer A ^ (N - 2) * mpsInsertion A) : ℤ) : ℚ) <
      -κ * ((trace (mpsTransfer A ^ N) : ℤ) : ℚ))
    {β : ℝ} (hβ : 0 < β) : Real.exp (β * κ * N) < partitionFunctionZ N β :=
  MPS.exp_lt_gibbsTrace A exchangeZ exchangeZ_transpose (by omega) κ hU hβ

/-- The trial-state bound with an integer certificate `a Tr (T^(N-2) O) < -b Tr Tᴺ`, `a > 0`:
then `Z_N(β) > e^{β (b/a) N}`. -/
theorem exp_lt_partitionFunctionZ_of_mps_int {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ) {N : ℕ} (hN : 3 ≤ N) {a b : ℤ} (ha : 0 < a)
    (hU : a * trace (mpsTransfer A ^ (N - 2) * mpsInsertion A) < -b * trace (mpsTransfer A ^ N))
    {β : ℝ} (hβ : 0 < β) : Real.exp (β * (b / a) * N) < partitionFunctionZ N β :=
  MPS.exp_lt_gibbsTrace_int A exchangeZ exchangeZ_transpose (by omega) ha hU hβ

/-- The trial-state bound with the spin-one certificate constant `κ = 700741/500000`. -/
theorem exp_lt_partitionFunctionZ_of_mps_700741 {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ) {N : ℕ} (hN : 3 ≤ N)
    (hU : 500000 * trace (mpsTransfer A ^ (N - 2) * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ N))
    {β : ℝ} (hβ : 0 < β) : Real.exp (β * (700741 / 500000) * N) < partitionFunctionZ N β := by
  simpa using exp_lt_partitionFunctionZ_of_mps_int A hN (a := 500000) (b := 700741)
    (by norm_num) hU hβ

/-- **The trial state bounds the ground-state energy of the Heisenberg ring.** If an integer
matrix-product state has `500000 Tr (T^(N-2) O) < -700741 Tr Tᴺ` on `N ≥ 2` sites, the ring has an
eigenvalue below `-κN`, `κ = 700741/500000`. The `S^z`-basis bound transfers to the Cartesian
ring because the two are unitarily conjugate (`conj_hamiltonian`). -/
theorem exists_spectrum_lt_of_mps {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ) {N : ℕ} (hN : 2 ≤ N)
    (hU : 500000 * trace (mpsTransfer A ^ (N - 2) * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ N)) :
    ∃ m, spectrum N m < -(700741 / 500000 * N) := by
  obtain ⟨i, hi⟩ := MPS.exists_eigenvalue_lt_int A exchangeZ exchangeZ_transpose hN
    (a := 500000) (b := 700741) (by norm_num) hU
  have hU1 : (szFrame N)ᴴᴴ * (szFrame N)ᴴ = 1 := by
    rw [conjTranspose_conjTranspose]; exact mul_eq_one_comm.mp (szFrame_conjTranspose_mul_self N)
  have heq := eigenvalues_eq_of_conj (hamiltonian_isHermitian N)
    (ringHamiltonian_isHermitian (isHermitian_map_intCast exchangeZ_transpose) N) hU1
    (by rw [conjTranspose_conjTranspose, conj_hamiltonian (by omega)]; rfl)
  refine ⟨i, ?_⟩
  rw [spectrum, ← heq]
  push_cast at hi
  linarith

end CertifiedGaps.SpinOneHeisenberg
