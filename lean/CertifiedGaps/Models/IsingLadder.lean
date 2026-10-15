/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Gap

/-!
# The two-leg transverse-field Ising ladder

`H_L = -∑_{legs} ZᵢZⱼ - ∑_{rungs} ZᵢZⱼ - g ∑ᵢ Xᵢ` (Pauli matrices) on a periodic ring of `L`
rungs. It is a reflection-positive ring of rungs (`model g`) with cell space `Fin 2 × Fin 2`
(leg `a` ⊗ leg `b`), bond factors `Z ⊗ 1` and `1 ⊗ Z` (so `-h = ∑_{legs} Z ⊗ Z`), and the rung
and field as the on-site term `o = -Z ⊗ Z - g (X ⊗ 1 + 1 ⊗ X)`. No basis change is needed.

The seeds are certified in `models/two_leg_tfim` (complete spectra of the rings of `4` and `8`
rungs with rigorous eigenvalue enclosures); they are the only `sorry`s.
-/

open Matrix Real
open scoped Kronecker

namespace CertifiedGaps.IsingLadder

/-- The Pauli matrix `Z`. -/
def pauliZ : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, -1]

/-- The Pauli matrix `X`. -/
def pauliX : Matrix (Fin 2) (Fin 2) ℝ := !![0, 1; 1, 0]

theorem pauliZ_transpose : pauliZᵀ = pauliZ := by
  ext i j; fin_cases i <;> fin_cases j <;> rfl

theorem pauliX_transpose : pauliXᵀ = pauliX := by
  ext i j; fin_cases i <;> fin_cases j <;> rfl

/-- The two-leg Ising ladder in the transverse field `g` as a reflection-positive ring of rungs. -/
noncomputable def model (g : ℝ) : RPRing (Fin 2 × Fin 2) (Fin 2) where
  C := ![pauliZ ⊗ₖ 1, 1 ⊗ₖ pauliZ]
  o := -(pauliZ ⊗ₖ pauliZ) - g • (pauliX ⊗ₖ 1 + 1 ⊗ₖ pauliX)
  C_symm α := by
    left
    fin_cases α <;>
      simp [← kroneckerMap_transpose, pauliZ_transpose]
  o_symm := by
    simp [transpose_sub, transpose_add, ← kroneckerMap_transpose, pauliZ_transpose,
      pauliX_transpose]

/-- The rational chain certificate for `g = 3`, seed `(n₀, β₀) = (4, 3/2)`. -/
def certG3 : ChainCertificate where
  s := [
    1125e-4,
    5736e-5,
    1035e-5,
    2422e-7,
    1186e-10,
    2814e-17,
    1584e-30,
    5019e-57,
    5039e-110]
  t := [
    5162e-5,
    1611e-5,
    9601e-7,
    2346e-9,
    1158e-14,
    2689e-25,
    1447e-46,
    4188e-89,
    3508e-174]
  r := 2585e-3
  ε := 5039e-110

theorem certG3_valid : certG3.Valid := by decide +kernel

/-- Seed: `P_s(4, 3/2) ≥ 1 - 0.1125`. Certified in `models/two_leg_tfim`
(`εₛ ≈ 0.11240836635`). -/
theorem seed_spatial_G3 : 1 - (1125e-4 : ℝ) ≤ (model 3).spatialPurity 4 (3 / 2) := by
  sorry

/-- Seed: `P_τ(8, 3/2) ≥ 1 - 0.05162`. Certified in `models/two_leg_tfim`
(`εₜ ≈ 0.05161523886`). -/
theorem seed_thermal_G3 : 1 - (5162e-5 : ℝ) ≤ (model 3).thermalPurity 8 (3 / 2) := by
  sorry

/-- **The two-leg Ising ladder at `g = 3`.** Every even ring of `L ≥ 4` rungs has a unique ground
state and gap `Δ_L > log 2.585 / (3/2) ≈ 0.6332`. -/
theorem gap_G3 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT ((model 3).spectrum L) (log 2.585 / (3 / 2)) := by
  have h := (model 3).gapGT_of_chainCertificate (n₀ := 4) (β₀ := 3 / 2) (by norm_num)
    (by decide) (by norm_num) certG3 certG3_valid
    (by simpa [certG3, ChainCertificate.sAt] using seed_spatial_G3)
    (by simpa [certG3, ChainCertificate.tAt] using seed_thermal_G3) hL hLe
  convert h using 3
  simp [certG3]

end CertifiedGaps.IsingLadder
