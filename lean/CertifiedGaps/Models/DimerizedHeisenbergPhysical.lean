/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.DimerizedHeisenberg
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Models.SpinHalfRotation

/-!
# The dimerized spin-one-half chain in the physical frame

`physicalHamiltonian δ L` is the dimerized chain
`∑ⱼ (1 + δ) S₂ⱼ · S₂ⱼ₊₁ + (1 - δ) S₂ⱼ₊₁ · S₂ⱼ₊₂` on `2L` spins, written as a ring of `L` strong
pairs: the on-site term is `(1 + δ) S · S` within a pair and the bond is
`(1 - δ) ∑_α (1 ⊗ S_α) ⊗ (S_α ⊗ 1)` from the second spin of a pair to the first spin of the next.
For `δ ≤ 1` it is the rotated chain `model δ` conjugated by `e^{iπSʸ}` on the second spin of every
pair (`physicalHamiltonian_eq`), so the gap theorems hold for the physical chain.
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

theorem cellRot_mul_conjTranspose : cellRot * cellRotᴴ = 1 := by
  rw [cellRot, conjTranspose_kronecker, ← mul_kronecker_mul, conjTranspose_one, Matrix.mul_one,
    rot_mul_conjTranspose, one_kronecker_one]

theorem conjTranspose_cellRot_mul : cellRotᴴ * cellRot = 1 := by
  rw [cellRot, conjTranspose_kronecker, ← mul_kronecker_mul, conjTranspose_one, Matrix.mul_one,
    conjTranspose_rot_mul, one_kronecker_one]

theorem cellRot_conj (A B : Matrix (Fin 2) (Fin 2) ℂ) :
    cellRot * (A ⊗ₖ B) * cellRotᴴ = A ⊗ₖ (rot * B * rotᴴ) := by
  rw [cellRot, conjTranspose_kronecker, ← mul_kronecker_mul, ← mul_kronecker_mul,
    conjTranspose_one, Matrix.one_mul, Matrix.mul_one]

theorem Cc_eq (δ : ℝ) (α : Fin 3) :
    (model δ).Cc α = ((Real.sqrt (1 - δ) : ℝ) : ℂ) • (1 ⊗ₖ Pc α) := by
  ext ⟨a, b⟩ ⟨c', e⟩
  by_cases h : a = c' <;> simp [RPRing.Cc, model, Pc, kroneckerMap_apply, one_apply, h]

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

theorem oc_eq (δ : ℝ) :
    (model δ).oc = ((-(1 + δ) : ℝ) : ℂ) • ∑ α, Pc α ⊗ₖ Pc α := by
  ext ⟨a, b⟩ ⟨c', e⟩
  simp [RPRing.oc, model, Pc, kroneckerMap_apply, Matrix.sum_apply, Finset.mul_sum]

/-- **The rotation.** For `δ ≤ 1` the physical dimerized chain is the rotated chain conjugated
by `e^{iπSʸ}` on the second spin of every pair. -/
theorem physicalHamiltonian_eq {δ : ℝ} (hδ : δ ≤ 1) {L : ℕ} (hL : 2 ≤ L) :
    physicalHamiltonian δ L = piTensor (fun _ : Fin L ↦ cellRot) * (model δ).hamiltonian L *
      piTensor (fun _ ↦ cellRotᴴ) := by
  have hu : ∀ _ : Fin L, cellRot * cellRotᴴ = 1 := fun _ ↦ cellRot_mul_conjTranspose
  have ha : ((Real.sqrt (1 - δ) : ℝ) : ℂ) * ((Real.sqrt (1 - δ) : ℝ) : ℂ) = ((1 - δ : ℝ) : ℂ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by linarith)]
  rw [RPRing.hamiltonian, Matrix.mul_add, Matrix.add_mul, Finset.mul_sum, Finset.sum_mul,
    Finset.mul_sum, Finset.sum_mul, physicalHamiltonian]
  congr 1
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [RPRing.bond, RPRing.bondMatrix, bondOp_neg, bondOp_sum, Matrix.mul_neg, Matrix.neg_mul,
      Finset.mul_sum, Finset.sum_mul, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun α _ ↦ ?_
    rw [conj_bondOp_kronecker _ hu (ne_finRotate hL j), ← bondOp_neg, Cc_eq, Dc_eq]
    congr 1
    rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_smul, Matrix.smul_mul, cellRot_conj,
      cellRot_conj, rot_conj_Pc, Matrix.mul_one, rot_mul_conjTranspose, Pc_eq]
    simp only [smul_kronecker, kronecker_smul, smul_smul, ← neg_smul]
    congr 1
    linear_combination
      -(((Real.sqrt (1 - δ) : ℝ) : ℂ) * ((Real.sqrt (1 - δ) : ℝ) : ℂ)) * neg_s_mul_c_sq α - ha
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [RPRing.site, conj_siteOp _ hu, oc_eq]
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
    physicalSpectrum δ L = (model δ).spectrum L := by
  refine eigenvalues_eq_of_conj _ _ (U := piTensor (fun _ : Fin L ↦ cellRot)) ?_ ?_
  · rw [conjTranspose_piTensor]
    exact piTensor_conjTranspose_mul_self _ fun _ ↦ conjTranspose_cellRot_mul
  · rw [conjTranspose_piTensor]; exact physicalHamiltonian_eq hδ hL

/-- **The physical dimerized chain at `δ = 1/2`.** Every even ring of `L ≥ 4` pairs has a unique
ground state and gap `Δ_L > log 24.33 / 6 ≈ 0.5320`. -/
theorem physical_gap_d12 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (1 / 2) L) (log 24.33 / 6) := by
  rw [physicalSpectrum_eq (by norm_num) (by omega)]; exact gap_d12 hL hLe

/-- The physical dimerized chain at `δ = 3/10`: `Δ_L > log 3.065 / (13/2) ≈ 0.1723`. -/
theorem physical_gap_d310 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (3 / 10) L) (log 3.065 / (13 / 2)) := by
  rw [physicalSpectrum_eq (by norm_num) (by omega)]; exact gap_d310 hL hLe

/-- The physical dimerized chain at `δ = 1/4`: `Δ_L > log 1.514 / (29/4) ≈ 0.0572`. -/
theorem physical_gap_d14 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (1 / 4) L) (log 1.514 / (29 / 4)) := by
  rw [physicalSpectrum_eq (by norm_num) (by omega)]; exact gap_d14 hL hLe

end CertifiedGaps.DimerizedHeisenberg
