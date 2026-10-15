/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.TwoLegLadder
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Models.SpinHalfRotation

/-!
# The two-leg ladder in the physical frame

`physicalHamiltonian J⊥ L = ∑_{legs} Sᵢ · Sⱼ + J⊥ ∑_{rungs} Sᵢ · Sⱼ` on a ring of `L` rungs. On an
even ring it is the rotated ladder `model J⊥` conjugated by `e^{iπSʸ}` on the checkerboard
sublattice (leg `b` of the even rungs and leg `a` of the odd rungs; `physicalHamiltonian_eq`),
so the gap theorem holds for the physical ladder (`physical_gap_J4`).
-/

open Matrix Real
open scoped Kronecker

namespace CertifiedGaps.TwoLegLadder

open SpinHalf

/-- The spin matrix `S_α` on leg `α.1` of a rung. -/
noncomputable def legS (α : Fin 2 × Fin 3) : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  if α.1 = 0 then spinHalf α.2 ⊗ₖ 1 else 1 ⊗ₖ spinHalf α.2

/-- The physical ladder: leg exchanges `∑_α legS α ⊗ legS α` and rung exchange `J⊥ S · S`. -/
noncomputable def physicalHamiltonian (Jperp : ℝ) (L : ℕ) :
    Matrix (Fin L → Fin 2 × Fin 2) (Fin L → Fin 2 × Fin 2) ℂ :=
  ∑ j : Fin L, ∑ α, bondOp j (finRotate L j) (legS α ⊗ₖ legS α) +
    ∑ j : Fin L, siteOp j ((Jperp : ℂ) • ∑ k, spinHalf k ⊗ₖ spinHalf k)

/-- The leg spin matrices are Hermitian. -/
theorem legS_isHermitian (α : Fin 2 × Fin 3) : (legS α).IsHermitian := by
  unfold legS
  split_ifs
  · change (_ ⊗ₖ _)ᴴ = _; rw [conjTranspose_kronecker, (spinHalf_isHermitian _).eq,
      conjTranspose_one]
  · change (_ ⊗ₖ _)ᴴ = _; rw [conjTranspose_kronecker, (spinHalf_isHermitian _).eq,
      conjTranspose_one]

/-- The physical ladder is Hermitian. -/
theorem physicalHamiltonian_isHermitian (Jperp : ℝ) (L : ℕ) :
    (physicalHamiltonian Jperp L).IsHermitian := by
  refine (isSelfAdjoint_sum _ fun j _ ↦ isSelfAdjoint_sum _ fun α _ ↦
    bondOp_isHermitian j _ ?_).add (isSelfAdjoint_sum _ fun j _ ↦ siteOp_isHermitian j ?_)
  · change (legS α ⊗ₖ legS α)ᴴ = _
    rw [conjTranspose_kronecker, (legS_isHermitian α).eq]
  · change (_ • _)ᴴ = _
    rw [conjTranspose_smul, conjTranspose_sum, Complex.star_def, Complex.conj_ofReal]
    simp only [conjTranspose_kronecker, (spinHalf_isHermitian _).eq]

/-- The spectrum of the physical ladder with multiplicity. -/
noncomputable def physicalSpectrum (Jperp : ℝ) (L : ℕ) : (Fin L → Fin 2 × Fin 2) → ℝ :=
  (physicalHamiltonian_isHermitian Jperp L).eigenvalues

/-- The checkerboard rotation on rung `j`: leg `a` if `j` is odd, leg `b` if `j` is even. -/
noncomputable def rungRot (L : ℕ) (j : Fin L) : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  if Odd (j : ℕ) then rot ⊗ₖ 1 else 1 ⊗ₖ rot

/-- The checkerboard rotation is unitary. -/
theorem rungRot_mul_conjTranspose (L : ℕ) (j : Fin L) : rungRot L j * (rungRot L j)ᴴ = 1 := by
  unfold rungRot
  split_ifs <;>
    rw [conjTranspose_kronecker, ← mul_kronecker_mul, conjTranspose_one] <;>
    simp [rot_mul_conjTranspose]

/-- The sign that the rotation of rung `j` puts on the bond factor `α`: `s` on the rotated
leg, `1` on the other. -/
noncomputable def sign (L : ℕ) (j : Fin L) (α : Fin 2 × Fin 3) : ℂ :=
  if (α.1 = 0 ↔ Odd (j : ℕ)) then s α.2 else 1

/-- The bond factors of the rotated ladder are `c_α S_α` on one leg. -/
theorem Cc_eq (Jperp : ℝ) (α : Fin 2 × Fin 3) :
    (model Jperp).Cc α = c α.2 • legS α := by
  have h : (model Jperp).Cc α = if α.1 = 0 then Pc α.2 ⊗ₖ 1 else 1 ⊗ₖ Pc α.2 := by
    unfold RPRing.Cc
    simp only [model]
    split_ifs <;> ext ⟨a, b⟩ ⟨c', e⟩ <;> by_cases hbe : b = e <;> by_cases hac : a = c' <;>
      simp [Pc, kroneckerMap_apply, one_apply, hbe, hac]
  rw [h, legS]
  split_ifs <;> rw [Pc_eq] <;> simp [smul_kronecker, kronecker_smul]

/-- The checkerboard rotation of rung `j` multiplies a bond factor by `sign L j α`. -/
theorem conj_legS (L : ℕ) (j : Fin L) (α : Fin 2 × Fin 3) :
    rungRot L j * (c α.2 • legS α) * (rungRot L j)ᴴ = (sign L j α * c α.2) • legS α := by
  rw [Matrix.mul_smul, Matrix.smul_mul, ← smul_smul, smul_comm (sign L j α)]
  congr 1
  unfold rungRot sign legS
  by_cases hj : Odd (j : ℕ) <;> by_cases hl : α.1 = 0 <;>
    simp only [hj, hl, iff_true, iff_false, not_true_eq_false,
      ite_true, ite_false, kronecker_conj, rot_conj_spinHalf, Matrix.one_mul, Matrix.mul_one,
      conjTranspose_one, rot_mul_conjTranspose, one_smul, smul_kronecker, kronecker_smul]

/-- On an even ring the rotation flips the sign of exactly one end of every leg bond. -/
theorem sign_mul_sign {L : ℕ} (hL : Even L) (j : Fin L) (α : Fin 2 × Fin 3) :
    sign L j α * sign L (finRotate L j) α = s α.2 := by
  have := odd_finRotate_iff hL j
  unfold sign
  by_cases hl : α.1 = 0 <;> by_cases hj : Odd (j : ℕ) <;>
    simp only [hl, hj, this, iff_true, iff_false, not_true_eq_false,
      not_false_eq_true, ite_true, ite_false, one_mul, mul_one]

/-- The on-site term of the rotated ladder as a complex matrix. -/
theorem oc_eq (Jperp : ℝ) : (model Jperp).oc = ((-Jperp : ℝ) : ℂ) • ∑ k, Pc k ⊗ₖ Pc k := by
  ext ⟨a, b⟩ ⟨c', e⟩
  simp [RPRing.oc, model, rungExchange, Pc, kroneckerMap_apply, Matrix.sum_apply, Finset.mul_sum]

/-- **The rotation.** On an even ring the physical ladder is the rotated ladder conjugated by
`e^{iπSʸ}` on the checkerboard sublattice. -/
theorem physicalHamiltonian_eq (Jperp : ℝ) {L : ℕ} (hL : Even L) :
    physicalHamiltonian Jperp L = piTensor (rungRot L) * (model Jperp).hamiltonian L *
      piTensor (fun i ↦ (rungRot L i)ᴴ) := by
  rw [(model Jperp).conj_hamiltonian _ (rungRot_mul_conjTranspose L)
    (by rintro rfl; exact Nat.not_even_one hL), physicalHamiltonian]
  congr 1
  · refine Finset.sum_congr rfl fun j _ ↦ Finset.sum_congr rfl fun α _ ↦ ?_
    rw [(model Jperp).Dc_eq_of_refl rfl, Cc_eq, conj_legS, conj_legS]
    congr 1
    rw [smul_kronecker, kronecker_smul, smul_smul, ← neg_smul]
    have h1 : -(sign L j α * c α.2 * (sign L (finRotate L j) α * c α.2)) = 1 := by
      linear_combination (-(c α.2 ^ 2)) * sign_mul_sign hL j α + neg_s_mul_c_sq α.2
    rw [h1, one_smul]
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [oc_eq]
    congr 1
    rw [Matrix.mul_smul, Matrix.smul_mul, Finset.mul_sum, Finset.sum_mul, Finset.smul_sum,
      Finset.smul_sum]
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    have hk : rungRot L j * (Pc k ⊗ₖ Pc k) * (rungRot L j)ᴴ = s k • (Pc k ⊗ₖ Pc k) := by
      unfold rungRot
      split_ifs <;> rw [kronecker_conj, rot_conj_Pc] <;>
        simp [smul_kronecker, kronecker_smul]
    rw [hk, Pc_eq]
    simp only [smul_kronecker, kronecker_smul, smul_smul]
    congr 1
    push_cast
    linear_combination (-(Jperp : ℂ)) * neg_s_mul_c_sq k

/-- On an even ring the physical and the rotated ladder have the same spectrum. -/
theorem physicalSpectrum_eq (Jperp : ℝ) {L : ℕ} (hL : Even L) :
    physicalSpectrum Jperp L = (model Jperp).spectrum L :=
  (model Jperp).eigenvalues_eq_spectrum _ (rungRot_mul_conjTranspose L) _
    (physicalHamiltonian_eq Jperp hL)

/-- **The physical two-leg ladder at `J⊥ = 4`.** Every even ring of `L ≥ 4` rungs has a unique
ground state and gap `Δ_L > log 9.033 / 2 ≈ 1.1004`. -/
theorem physical_gap_J4 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT (physicalSpectrum 4 L) (log 9.033 / 2) := by
  rw [physicalSpectrum_eq 4 hLe]; exact gap_J4 hL hLe

end CertifiedGaps.TwoLegLadder
