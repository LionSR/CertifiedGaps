/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.DimerizedHeisenberg
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Local.SpinHalfRotation

/-!
# The dimerized spin-one-half chain in the physical frame

`physicalHamiltonian δ L` is the dimerized chain
`∑ⱼ (1 + δ) S₂ⱼ · S₂ⱼ₊₁ + (1 - δ) S₂ⱼ₊₁ · S₂ⱼ₊₂` on `2L` spins, written as a ring of `L` strong
pairs: the on-site term is `(1 + δ) S · S` within a pair and the bond is
`(1 - δ) ∑_α (1 ⊗ S_α) ⊗ (S_α ⊗ 1)` from the second spin of a pair to the first spin of the next.
For `δ ≤ 1` it is the rotated chain `model δ` conjugated by `e^{iπSʸ}` on the second spin of every
pair (`physicalHamiltonian_eq`), so the gap theorem holds for the physical chain (`physical_gap`).
-/

open Matrix Real
open scoped Kronecker

namespace CertifiedGaps.DimerizedHeisenberg

open SpinHalf

/-- The physical dimerized chain as a ring of strong pairs. -/
noncomputable def physicalHamiltonian (δ : ℝ) (L : ℕ) :
    Matrix (Fin L → Fin 2 × Fin 2) (Fin L → Fin 2 × Fin 2) ℂ :=
  ∑ j : Fin L, ∑ α, bondOp j (finRotate L j)
      (((1 - δ : ℝ) : ℂ) • ((1 ⊗ₖ spinHalf α) ⊗ₖ (spinHalf α ⊗ₖ 1))) +
    ∑ j : Fin L, siteOp j (((1 + δ : ℝ) : ℂ) • ∑ α, spinHalf α ⊗ₖ spinHalf α)

/-- The physical dimerized chain is Hermitian. -/
theorem physicalHamiltonian_isHermitian (δ : ℝ) (L : ℕ) :
    (physicalHamiltonian δ L).IsHermitian := by
  have hS := spinHalf_isHermitian
  have hreal : ∀ r : ℝ, star (r : ℂ) = r := fun r ↦ Complex.conj_ofReal r
  refine (isSelfAdjoint_sum _ fun j _ ↦ isSelfAdjoint_sum _ fun α _ ↦
    bondOp_isHermitian j _ ?_).add (isSelfAdjoint_sum _ fun j _ ↦ siteOp_isHermitian j ?_)
  · change (_ • _)ᴴ = _
    rw [conjTranspose_smul, conjTranspose_kronecker, conjTranspose_kronecker,
      conjTranspose_kronecker, (hS α).eq, conjTranspose_one, hreal]
  · change (_ • _)ᴴ = _
    rw [conjTranspose_smul, conjTranspose_sum, hreal]
    simp only [conjTranspose_kronecker, (hS _).eq]

/-- The spectrum of the physical dimerized chain with multiplicity. -/
noncomputable def physicalSpectrum (δ : ℝ) (L : ℕ) : (Fin L → Fin 2 × Fin 2) → ℝ :=
  (physicalHamiltonian_isHermitian δ L).eigenvalues

/-- `e^{iπSʸ}` on the second spin of a pair. -/
def cellRot : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ := 1 ⊗ₖ rot

/-- The rotation of a pair is unitary. -/
theorem cellRot_mul_conjTranspose : cellRot * cellRotᴴ = 1 := by
  rw [cellRot, conjTranspose_kronecker, ← mul_kronecker_mul, conjTranspose_one, Matrix.mul_one,
    rot_mul_conjTranspose, one_kronecker_one]

/-- The rotation of a pair acts on its second spin. -/
theorem cellRot_conj (A B : Matrix (Fin 2) (Fin 2) ℂ) :
    cellRot * (A ⊗ₖ B) * cellRotᴴ = A ⊗ₖ (rot * B * rotᴴ) := by
  rw [cellRot, kronecker_conj, Matrix.one_mul, conjTranspose_one, Matrix.mul_one]

/-- The bond factors of the rotated chain as complex matrices. -/
theorem Cc_eq (δ : ℝ) (α : Fin 3) :
    (model δ).Cc α = ((Real.sqrt (1 - δ) : ℝ) : ℂ) • (1 ⊗ₖ Pc α) := by
  ext ⟨a, b⟩ ⟨c', e⟩
  by_cases h : a = c' <;> simp [RPRing.Cc, model, Pc, kroneckerMap_apply, one_apply, h]

/-- The reflected bond factors of the rotated chain as complex matrices. -/
theorem Dc_eq (δ : ℝ) (α : Fin 3) :
    (model δ).Dc α = ((Real.sqrt (1 - δ) : ℝ) : ℂ) • (Pc α ⊗ₖ 1) := by
  have h : (model δ).Dc α =
      Complex.ofRealHom.mapMatrix ((model δ).refl * (model δ).C α * (model δ).refl) := by
    rw [map_mul, map_mul]; rfl
  rw [h]
  change Complex.ofRealHom.mapMatrix (cellSwap * (Real.sqrt (1 - δ) • (1 ⊗ₖ P α)) * cellSwap) = _
  rw [Matrix.mul_smul, Matrix.smul_mul, cellSwap_mul_kronecker_mul_cellSwap]
  ext ⟨a, b⟩ ⟨c', e⟩
  by_cases h : b = e <;> simp [Pc, kroneckerMap_apply, one_apply, h]

/-- The on-site term of the rotated chain as a complex matrix. -/
theorem oc_eq (δ : ℝ) :
    (model δ).oc = ((-(1 + δ) : ℝ) : ℂ) • ∑ α, Pc α ⊗ₖ Pc α := by
  ext ⟨a, b⟩ ⟨c', e⟩
  simp [RPRing.oc, model, Pc, kroneckerMap_apply, Matrix.sum_apply, Finset.mul_sum]

/-- **The rotation.** For `δ ≤ 1` the physical dimerized chain is the rotated chain conjugated
by `e^{iπSʸ}` on the second spin of every pair. -/
theorem physicalHamiltonian_eq {δ : ℝ} (hδ : δ ≤ 1) {L : ℕ} (hL : 2 ≤ L) :
    physicalHamiltonian δ L = piTensor (fun _ : Fin L ↦ cellRot) * (model δ).hamiltonian L *
      piTensor (fun _ ↦ cellRotᴴ) := by
  have ha : ((Real.sqrt (1 - δ) : ℝ) : ℂ) * ((Real.sqrt (1 - δ) : ℝ) : ℂ) = ((1 - δ : ℝ) : ℂ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by linarith)]
  rw [(model δ).conj_hamiltonian _ (fun _ ↦ cellRot_mul_conjTranspose) (by omega),
    physicalHamiltonian]
  congr 1
  · refine Finset.sum_congr rfl fun j _ ↦ Finset.sum_congr rfl fun α _ ↦ ?_
    rw [Cc_eq, Dc_eq]
    congr 1
    rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_smul, Matrix.smul_mul, cellRot_conj,
      cellRot_conj, rot_conj_Pc, Matrix.mul_one, rot_mul_conjTranspose, Pc_eq]
    simp only [smul_kronecker, kronecker_smul, smul_smul, ← neg_smul]
    congr 1
    linear_combination
      -(((Real.sqrt (1 - δ) : ℝ) : ℂ) * ((Real.sqrt (1 - δ) : ℝ) : ℂ)) * neg_s_mul_c_sq α - ha
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [oc_eq]
    congr 1
    rw [Matrix.mul_smul, Matrix.smul_mul, Finset.mul_sum, Finset.sum_mul, Finset.smul_sum,
      Finset.smul_sum]
    refine Finset.sum_congr rfl fun α _ ↦ ?_
    rw [cellRot_conj, rot_conj_Pc, Pc_eq]
    simp only [smul_kronecker, kronecker_smul, smul_smul]
    congr 1
    push_cast
    linear_combination (-(1 + δ : ℂ)) * neg_s_mul_c_sq α

/-- For `δ ≤ 1` the physical and the rotated chain have the same spectrum. -/
theorem physicalSpectrum_eq {δ : ℝ} (hδ : δ ≤ 1) {L : ℕ} (hL : 2 ≤ L) :
    physicalSpectrum δ L = (model δ).spectrum L :=
  (model δ).eigenvalues_eq_spectrum _ (fun _ ↦ cellRot_mul_conjTranspose) _
    (physicalHamiltonian_eq hδ hL)

/-- **The physical dimerized chain.** At every point `p ∈ points`, every even ring of `L ≥ 4`
pairs has a unique ground state and gap `Δ_L > log r / β₀` (`≈ 0.5320, 0.1723, 0.0572` at
`δ = 1/2, 3/10, 1/4`). -/
theorem physical_gap {p : Point} (hp : p ∈ points) {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum p.δ L) (log (p.cert.r : ℝ) / p.β₀) := by
  rw [physicalSpectrum_eq (by exact_mod_cast points_δ_le_one p hp) (by omega)]
  exact gap hp hL hLe

end CertifiedGaps.DimerizedHeisenberg
