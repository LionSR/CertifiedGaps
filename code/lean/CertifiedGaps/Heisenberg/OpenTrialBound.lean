/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.OpenMatrixProduct
import CertifiedGaps.Heisenberg.OpenChain
import CertifiedGaps.Heisenberg.TrialBound

/-!
# The matrix-product trial bound on the opposite-field chain

An open matrix-product trial vector `ψ(σ) = bᵀ A_{σ₁} ⋯ A_{σₘ} c` with integer data bounds the
ground-state energy of the open spin-one chain with end fields,

  `H = ∑ⱼ Sⱼ · Sⱼ₊₁ - h₀ S₁ᶻ - h₁ Sₘᶻ`   (`endFieldHamiltonian h₀ h₁ m`),

from above. This file instantiates the generic bound of
`CertifiedGaps/Transfer/OpenMatrixProduct.lean` in the `S^z` basis: the chain is unitarily
conjugate (`szFrame`) to the open chain of the integer exchange `exchangeZ` with end terms
`-h₀ Sᶻ` and `-h₁ Sᶻ`, `Sᶻ = diag(-1, 0, 1)` (`szDiag`), so the norm and the energy of `ψ` are
the integer contractions `OpenMPS.norm`, `OpenMPS.bondEnergy`, `OpenMPS.leftField`,
`OpenMPS.rightField`.

The instance is the opposite-field chain `H_m^{opp}` (`h₀ = 3/5`, `h₁ = -3/5`) of OpenAI's
*A boundary-field gap for the spin-one Heisenberg chain* (Proposition `apptr:bound`), with the
integer certificate in the paper's scaling: the energy `500000 ⟨ψ, (H + a m) ψ⟩`,
`a = 700741/500000`, below `202100 ⟨ψ, ψ⟩` gives an eigenvalue below `-a m + α`,
`α = 2021/5000`.

## Main results

* `conj_endFieldHamiltonian`: the chain with end fields in the `S^z` basis.
* `exists_endFieldSpectrum_lt_of_openMPS`: the variational bound for any end fields.
* `exists_oppositeFieldSpectrum_lt_of_openMPS`: the opposite-field chain with the paper's
  integer certificate.
-/

open Matrix

namespace CertifiedGaps.SpinOneHeisenberg

/-- `Sᶻ = diag(-1, 0, 1)` in the `S^z` basis, as an integer matrix. -/
def szDiag : Matrix (Fin 3) (Fin 3) ℤ := diagonal mz

/-- `Sᶻ` is symmetric. -/
theorem szDiag_transpose : szDiagᵀ = szDiag := diagonal_transpose _

/-- `Sᶻ` of `SpinOne.spinZ` is the cast of `szDiag`. -/
theorem spinZ_two_eq : SpinOne.spinZ 2 = szDiag.map (↑) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [SpinOne.spinZ, szDiag, mz]

/-- The chain with end fields is the generic open chain of the exchange with end terms
`-h₀ Sᶻ`, `-h₁ Sᶻ`. -/
theorem endFieldHamiltonian_eq (h₀ h₁ : ℝ) (n : ℕ) :
    endFieldHamiltonian h₀ h₁ n =
      openChainHamiltonian SpinOne.exchange (-(h₀ : ℂ) • SpinOne.S 2) (-(h₁ : ℂ) • SpinOne.S 2)
        n := by
  cases n <;> rfl

/-- **The chain with end fields in the `S^z` basis**: the open chain of `exchangeZ` with end
terms `-h₀ Sᶻ` and `-h₁ Sᶻ`. -/
theorem conj_endFieldHamiltonian (h₀ h₁ : ℝ) (n : ℕ) :
    (szFrame n)ᴴ * endFieldHamiltonian h₀ h₁ n * szFrame n =
      openChainHamiltonian (exchangeZ.map (↑)) (((-h₀ : ℝ) : ℂ) • szDiag.map (↑))
        (((-h₁ : ℝ) : ℂ) • szDiag.map (↑)) n := by
  have hS (h : ℝ) : SpinOne.szBasisᴴ * (-(h : ℂ) • SpinOne.S 2) * SpinOne.szBasis =
      ((-h : ℝ) : ℂ) • szDiag.map (↑) := by
    rw [Matrix.mul_smul, Matrix.smul_mul, SpinOne.conj_S, spinZ_two_eq, Complex.ofReal_neg]
  rw [endFieldHamiltonian_eq, szFrame,
    conj_openChainHamiltonian _ SpinOne.szBasis_conjTranspose_mul_self, conj_exchange, hS, hS]

/-- **The trial vector bounds the ground-state energy of the chain with end fields.** If an
integer open matrix-product state on `n + 1` sites has Rayleigh quotient below `E` in the
`S^z`-basis chain, `∑ⱼ L₀ Tʲ O T^(n-1-j) u - h₀ L₀ G Tⁿ u - h₁ L₀ Tⁿ G u < E L₀ T^(n+1) u`,
then the chain `endFieldHamiltonian h₀ h₁ (n + 1)` has an eigenvalue below `E`. -/
theorem exists_endFieldSpectrum_lt_of_openMPS {D : Type*} [Fintype D] [DecidableEq D]
    (b c : D → ℤ) (A : Fin 3 → Matrix D D ℤ) (h₀ h₁ : ℝ) (n : ℕ) {E : ℝ}
    (hE : (OpenMPS.bondEnergy exchangeZ b c A n : ℝ) + -h₀ * OpenMPS.leftField szDiag b c A n +
      -h₁ * OpenMPS.rightField szDiag b c A n < E * OpenMPS.norm b c A (n + 1)) :
    ∃ i, (endFieldHamiltonian_isHermitian h₀ h₁ (n + 1)).eigenvalues i < E := by
  obtain ⟨i, hi⟩ := OpenMPS.exists_eigenvalue_lt b c A exchangeZ exchangeZ_transpose
    szDiag_transpose szDiag_transpose (-h₀) (-h₁) n hE
  have hU1 : (szFrame (n + 1))ᴴᴴ * (szFrame (n + 1))ᴴ = 1 := by
    rw [conjTranspose_conjTranspose]
    exact mul_eq_one_comm.mp (szFrame_conjTranspose_mul_self (n + 1))
  have heq := eigenvalues_eq_of_conj (endFieldHamiltonian_isHermitian h₀ h₁ (n + 1))
    (openChainHamiltonian_isHermitian (isHermitian_map_intCast exchangeZ_transpose)
      (OpenMPS.isHermitian_smul_map_intCast szDiag_transpose (-h₀))
      (OpenMPS.isHermitian_smul_map_intCast szDiag_transpose (-h₁)) (n + 1)) hU1
    (by rw [conjTranspose_conjTranspose, conj_endFieldHamiltonian])
  exact ⟨i, heq ▸ hi⟩

/-- **The opposite-field trial bound with an integer certificate.** For the opposite-field
chain `H^{opp} = ∑ⱼ Sⱼ · Sⱼ₊₁ - (3/5) S₁ᶻ + (3/5) Sₘᶻ` on `m = n + 1` sites: if an integer open
matrix-product state satisfies the paper's comparison

  `500000 ∑ⱼ L₀ Tʲ O T^(n-1-j) u - 300000 L₀ G Tⁿ u + 300000 L₀ Tⁿ G u + 700741 m L₀ Tᵐ u
    < 202100 L₀ Tᵐ u`,

i.e. `500000 ⟨ψ, (H^{opp} + a m) ψ⟩ < 202100 ⟨ψ, ψ⟩`, then `H^{opp}` has an eigenvalue below
`-a m + α` with `a = 700741/500000`, `α = 2021/5000`. -/
theorem exists_oppositeFieldSpectrum_lt_of_openMPS {D : Type*} [Fintype D] [DecidableEq D]
    (b c : D → ℤ) (A : Fin 3 → Matrix D D ℤ) (n : ℕ)
    (hU : 500000 * OpenMPS.bondEnergy exchangeZ b c A n -
        300000 * OpenMPS.leftField szDiag b c A n + 300000 * OpenMPS.rightField szDiag b c A n +
        700741 * (n + 1) * OpenMPS.norm b c A (n + 1) < 202100 * OpenMPS.norm b c A (n + 1)) :
    ∃ i, oppositeFieldSpectrum (3 / 5) (n + 1) i <
      -(700741 / 500000 : ℝ) * (n + 1 : ℕ) + 2021 / 5000 := by
  refine exists_endFieldSpectrum_lt_of_openMPS b c A (3 / 5) (-(3 / 5)) n ?_
  have h := (Int.cast_lt (R := ℝ)).2 hU
  push_cast at h ⊢
  linarith

end CertifiedGaps.SpinOneHeisenberg
