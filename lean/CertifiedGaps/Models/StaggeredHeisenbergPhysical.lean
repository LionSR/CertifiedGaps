/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.StaggeredHeisenberg
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Models.SpinHalfRotation

/-!
# The staggered-field chain in the physical frame

`physicalHamiltonian h L = ∑ⱼ Sⱼ · Sⱼ₊₁ - h ∑ⱼ (-1)ʲ Sᶻⱼ` with the complex spin one-half matrices.
On an even ring it is the rotated chain `model h` conjugated by `e^{iπSʸ}` on the odd sites
(`physicalHamiltonian_eq`), so the two have the same spectrum (`physicalSpectrum_eq`) and the gap
theorems of `StaggeredHeisenberg` hold for the physical chain (`physical_gap_h1`, …).
-/

open Matrix Real
open scoped Kronecker

namespace CertifiedGaps.StaggeredHeisenberg

open SpinHalf

/-- The physical chain `∑ⱼ Sⱼ · Sⱼ₊₁ - h ∑ⱼ (-1)ʲ Sᶻⱼ`. -/
noncomputable def physicalHamiltonian (h : ℝ) (L : ℕ) :
    Matrix (Fin L → Fin 2) (Fin L → Fin 2) ℂ :=
  ∑ j : Fin L, ∑ α, bondOp j (finRotate L j) (spinHalf α ⊗ₖ spinHalf α) +
    ∑ j : Fin L, siteOp j (((-h * (-1) ^ (j : ℕ) : ℝ) : ℂ) • spinHalf 2)

/-- The physical chain is Hermitian. -/
theorem physicalHamiltonian_isHermitian (h : ℝ) (L : ℕ) :
    (physicalHamiltonian h L).IsHermitian := by
  refine (isSelfAdjoint_sum _ fun j _ ↦ isSelfAdjoint_sum _ fun α _ ↦
    bondOp_isHermitian j _ ?_).add (isSelfAdjoint_sum _ fun j _ ↦ siteOp_isHermitian j ?_)
  · change (spinHalf α ⊗ₖ spinHalf α)ᴴ = _
    rw [conjTranspose_kronecker, (spinHalf_isHermitian α).eq]
  · change (_ • spinHalf 2)ᴴ = _
    rw [conjTranspose_smul, (spinHalf_isHermitian 2).eq, Complex.star_def, Complex.conj_ofReal]

/-- The spectrum of the physical chain with multiplicity. -/
noncomputable def physicalSpectrum (h : ℝ) (L : ℕ) : (Fin L → Fin 2) → ℝ :=
  (physicalHamiltonian_isHermitian h L).eigenvalues

/-- `e^{iπSʸ}` on the odd sites. -/
noncomputable def rotation (L : ℕ) (i : Fin L) : Matrix (Fin 2) (Fin 2) ℂ :=
  if Odd (i : ℕ) then rot else 1

/-- The rotation is unitary. -/
theorem rotation_mul_conjTranspose (L : ℕ) (i : Fin L) : rotation L i * (rotation L i)ᴴ = 1 := by
  unfold rotation; split_ifs
  · exact rot_mul_conjTranspose
  · simp

/-- The rotated bond in either orientation is the exchange `S_α ⊗ S_α`. -/
theorem rotated_bond (h : ℝ) (α : Fin 3) {x y : Matrix (Fin 2) (Fin 2) ℂ}
    (hxy : (x = 1 ∧ y = rot) ∨ (x = rot ∧ y = 1)) :
    -((x * (model h).Cc α * xᴴ) ⊗ₖ (y * (model h).Cc α * yᴴ)) = spinHalf α ⊗ₖ spinHalf α := by
  rcases hxy with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
  · simp only [Matrix.one_mul, Matrix.mul_one, conjTranspose_one, rot_conj]
    fin_cases α <;> ext ⟨a, b⟩ ⟨c, e⟩ <;> fin_cases a <;> fin_cases b <;> fin_cases c <;>
      fin_cases e <;>
      simp [RPRing.Cc, model, P, spinHalf, kroneckerMap_apply] <;> ring_nf <;>
      simp [Complex.ext_iff] <;> norm_num

/-- The rotated on-site term at a site of parity `j` is `-h (-1)ʲ Sᶻ`. -/
theorem rotated_site (h : ℝ) (L : ℕ) (j : Fin L) :
    rotation L j * (model h).oc * (rotation L j)ᴴ =
      (((-h * (-1) ^ (j : ℕ) : ℝ) : ℂ)) • spinHalf 2 := by
  unfold rotation
  split_ifs with hj
  · rw [hj.neg_one_pow, rot_conj]
    ext a b; fin_cases a <;> fin_cases b <;>
      simp [RPRing.oc, model, P, spinHalf] <;> norm_num
  · rw [(Nat.not_odd_iff_even.1 hj).neg_one_pow]
    ext a b; fin_cases a <;> fin_cases b <;>
      simp [RPRing.oc, model, P, spinHalf]

/-- **The rotation.** On an even ring the physical chain is the rotated chain conjugated by
`e^{iπSʸ}` on the odd sites. -/
theorem physicalHamiltonian_eq (h : ℝ) {L : ℕ} (hL : Even L) :
    physicalHamiltonian h L = piTensor (rotation L) * (model h).hamiltonian L *
      piTensor (fun i ↦ (rotation L i)ᴴ) := by
  rw [(model h).conj_hamiltonian _ (rotation_mul_conjTranspose L)
    (by rintro rfl; exact Nat.not_even_one hL), physicalHamiltonian]
  congr 1
  · refine Finset.sum_congr rfl fun j _ ↦ Finset.sum_congr rfl fun α _ ↦ ?_
    rw [(model h).Dc_eq_of_refl rfl, rotated_bond]
    unfold rotation
    have := odd_finRotate_iff hL j
    by_cases hj : Odd (j : ℕ)
    · right; exact ⟨ite_eq_left hj, ite_eq_right (by rw [this]; exact not_not.mpr hj)⟩
    · left; exact ⟨ite_eq_right hj, ite_eq_left (this.mpr hj)⟩
  · exact Finset.sum_congr rfl fun j _ ↦ by rw [rotated_site]

/-- On an even ring the physical and the rotated chain have the same spectrum. -/
theorem physicalSpectrum_eq (h : ℝ) {L : ℕ} (hL : Even L) :
    physicalSpectrum h L = (model h).spectrum L :=
  (model h).eigenvalues_eq_spectrum _ (rotation_mul_conjTranspose L) _ (physicalHamiltonian_eq h hL)

/-- **The physical staggered-field chain at `h = 1`.** Every even ring of `L ≥ 8` sites has a
unique ground state and gap `Δ_L > log 413.4 / (21/4) ≈ 1.1475`. -/
theorem physical_gap_h1 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum 1 L) (log 413.4 / (21 / 4)) := by
  rw [physicalSpectrum_eq 1 hLe]; exact gap_h1 hL hLe

/-- The physical staggered-field chain at `h = 3/4`: `Δ_L > log 106.9 / (21/4) ≈ 0.8899`. -/
theorem physical_gap_h34 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (3 / 4) L) (log 106.9 / (21 / 4)) := by
  rw [physicalSpectrum_eq _ hLe]; exact gap_h34 hL hLe

/-- The physical staggered-field chain at `h = 1/2`: `Δ_L > log 26.23 / (11/2) ≈ 0.5940`. -/
theorem physical_gap_h12 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (1 / 2) L) (log 26.23 / (11 / 2)) := by
  rw [physicalSpectrum_eq _ hLe]; exact gap_h12 hL hLe

/-- The physical staggered-field chain at `h = 7/20`: `Δ_L > log 9.035 / (23/4) ≈ 0.3828`. -/
theorem physical_gap_h720 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (7 / 20) L) (log 9.035 / (23 / 4)) := by
  rw [physicalSpectrum_eq _ hLe]; exact gap_h720 hL hLe

/-- The physical staggered-field chain at `h = 1/4`: `Δ_L > log 3.806 / (25/4) ≈ 0.2139`. -/
theorem physical_gap_h14 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum (1 / 4) L) (log 3.806 / (25 / 4)) := by
  rw [physicalSpectrum_eq _ hLe]; exact gap_h14 hL hLe

end CertifiedGaps.StaggeredHeisenberg
