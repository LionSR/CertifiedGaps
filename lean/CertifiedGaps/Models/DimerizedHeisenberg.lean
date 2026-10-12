/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Gap
import CertifiedGaps.Models.SpinHalf
import CertifiedGaps.Models.CellSwap

/-!
# The dimerized spin-one-half Heisenberg chain

`H = ∑ⱼ (1 + δ) S₂ⱼ · S₂ⱼ₊₁ + (1 - δ) S₂ⱼ₊₁ · S₂ⱼ₊₂` on a periodic ring of `2L` spins, as a ring
of `L` cells (the strong pairs). Rotating the second spin of every cell by `e^{iπSʸ}` (a product
of one-site unitaries, so the spectrum is unchanged) turns every bond into
`-SˣSˣ + SʸSʸ - SᶻSᶻ`. In the rotated frame (`model δ`):

* the bond from the second spin of a cell to the first spin of the next is
  `-h = (1 - δ) ∑_α (1 ⊗ Pα) ⊗ (Pα ⊗ 1)` with `Pα = Sˣ, iSʸ, Sᶻ`, that is
  `-h = ∑_α Cα ⊗ (F Cα F)` with `Cα = √(1 - δ) (1 ⊗ Pα)` and `F` the swap of the cell;
* the strong pair is the on-site term `o = -(1 + δ) ∑_α Pα ⊗ Pα`, which commutes with `F`.

So it is a reflection-positive ring with internal reflection `F`. The seeds are certified in
`models/dimerized_heisenberg` (complete spectra of the rings of `8` and `16` spins); they are the
only `sorry`s.
-/

open Matrix Real ThermalPurityGap
open scoped Kronecker

namespace CertifiedGaps.DimerizedHeisenberg

open SpinHalf

/-- The rotated dimerized chain as a reflection-positive ring of strong pairs. -/
noncomputable def model (δ : ℝ) : RPRing (Fin 2 × Fin 2) (Fin 3) where
  C α := Real.sqrt (1 - δ) • (1 ⊗ₖ P α)
  o := -(1 + δ) • ∑ α, P α ⊗ₖ P α
  C_symm α := by
    rw [transpose_smul, ← kroneckerMap_transpose, transpose_one]
    rcases P_symm α with h | h <;> rw [h]
    · exact Or.inl rfl
    · right; ext; simp [kroneckerMap_apply]
  o_symm := by
    rw [transpose_smul, transpose_sum]
    congr 1
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [← kroneckerMap_transpose]
    rcases P_symm k with h | h <;> rw [h]
    ext; simp [kroneckerMap_apply]
  refl := cellSwap
  refl_symm := cellSwap_transpose
  refl_sq := cellSwap_mul_cellSwap
  refl_o := by
    rw [Matrix.mul_smul, Matrix.smul_mul, Finset.mul_sum, Finset.sum_mul]
    simp_rw [cellSwap_mul_kronecker_mul_cellSwap]

/-! ### `δ = 1 / 2` -/

/-- The rational chain certificate for `δ = 1 / 2`, seed `(n₀, β₀) = (4, 6)`. -/
def cert_d12 : ChainCertificate where
  s := [
    1133e-5,
    6383e-7,
    129e-8,
    3908e-15,
    3128e-26,
    1958e-48,
    7668e-93,
    1176e-181,
    2766e-359]
  t := [
    1224e-5,
    3275e-7,
    2156e-10,
    9297e-17,
    1729e-29,
    5979e-55,
    715e-105,
    1023e-207,
    2094e-411]
  r := 2433e-2
  ε := 2766e-359

theorem cert_d12_valid : cert_d12.Valid := by decide +kernel

/-- Seed: `P_s(4, 6) ≥ 1 - 1133e-5`. Certified in `models/dimerized_heisenberg`. -/
theorem seed_spatial_d12 : 1 - (1133e-5 : ℝ) ≤ (model (1 / 2)).spatialPurity 4 (6) := by
  sorry

/-- Seed: `P_τ(8, 6) ≥ 1 - 1224e-5`. Certified in `models/dimerized_heisenberg`. -/
theorem seed_thermal_d12 : 1 - (1224e-5 : ℝ) ≤ (model (1 / 2)).thermalPurity 8 (6) := by
  sorry

/-- **The dimerized chain at `δ = 1 / 2`.** Every even ring of `L ≥ 4` cells has a unique ground
state and gap `Δ_L > log 24.33 / (6) ≈ 0.5320`. -/
theorem gap_d12 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT ((model (1 / 2)).spectrum L) (log 24.33 / (6)) := by
  have h := (model (1 / 2)).gapGT_of_chainCertificate (n₀ := 4) (β₀ := 6) (by norm_num)
    (by decide) (by norm_num) cert_d12 cert_d12_valid
    (by simpa [cert_d12, ChainCertificate.sAt] using seed_spatial_d12)
    (by simpa [cert_d12, ChainCertificate.tAt] using seed_thermal_d12) hL hLe
  convert h using 3
  simp [cert_d12]

/-! ### `δ = 3 / 10` -/

/-- The rational chain certificate for `δ = 3 / 10`, seed `(n₀, β₀) = (4, 13 / 2)`. -/
def cert_d310 : ChainCertificate where
  s := [
    9516e-5,
    4151e-5,
    5272e-6,
    6201e-8,
    7763e-12,
    1206e-19,
    2909e-35,
    1693e-66,
    5733e-129]
  t := [
    518e-4,
    1287e-5,
    502e-6,
    5691e-10,
    6567e-16,
    8627e-28,
    1489e-51,
    4435e-99,
    3934e-194]
  r := 3065e-3
  ε := 5733e-129

theorem cert_d310_valid : cert_d310.Valid := by decide +kernel

/-- Seed: `P_s(4, 13 / 2) ≥ 1 - 9516e-5`. Certified in `models/dimerized_heisenberg`. -/
theorem seed_spatial_d310 : 1 - (9516e-5 : ℝ) ≤ (model (3 / 10)).spatialPurity 4 (13 / 2) := by
  sorry

/-- Seed: `P_τ(8, 13 / 2) ≥ 1 - 518e-4`. Certified in `models/dimerized_heisenberg`. -/
theorem seed_thermal_d310 : 1 - (518e-4 : ℝ) ≤ (model (3 / 10)).thermalPurity 8 (13 / 2) := by
  sorry

/-- **The dimerized chain at `δ = 3 / 10`.** Every even ring of `L ≥ 4` cells has a unique ground
state and gap `Δ_L > log 3.065 / (13 / 2) ≈ 0.1723`. -/
theorem gap_d310 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT ((model (3 / 10)).spectrum L) (log 3.065 / (13 / 2)) := by
  have h := (model (3 / 10)).gapGT_of_chainCertificate (n₀ := 4) (β₀ := 13 / 2) (by norm_num)
    (by decide) (by norm_num) cert_d310 cert_d310_valid
    (by simpa [cert_d310, ChainCertificate.sAt] using seed_spatial_d310)
    (by simpa [cert_d310, ChainCertificate.tAt] using seed_thermal_d310) hL hLe
  convert h using 3
  simp [cert_d310]

/-! ### `δ = 1 / 4` -/

/-- The rational chain certificate for `δ = 1 / 4`, seed `(n₀, β₀) = (4, 29 / 4)`. -/
def cert_d14 : ChainCertificate where
  s := [
    1658e-4,
    1325e-4,
    7345e-5,
    1653e-5,
    6135e-7,
    759e-9,
    1153e-15,
    2659e-27,
    1415e-50]
  t := [
    5123e-5,
    3937e-5,
    1438e-5,
    1091e-6,
    3923e-9,
    3703e-14,
    2829e-24,
    1603e-44,
    514e-84]
  r := 1514e-3
  ε := 1415e-50

theorem cert_d14_valid : cert_d14.Valid := by decide +kernel

/-- Seed: `P_s(4, 29 / 4) ≥ 1 - 1658e-4`. Certified in `models/dimerized_heisenberg`. -/
theorem seed_spatial_d14 : 1 - (1658e-4 : ℝ) ≤ (model (1 / 4)).spatialPurity 4 (29 / 4) := by
  sorry

/-- Seed: `P_τ(8, 29 / 4) ≥ 1 - 5123e-5`. Certified in `models/dimerized_heisenberg`. -/
theorem seed_thermal_d14 : 1 - (5123e-5 : ℝ) ≤ (model (1 / 4)).thermalPurity 8 (29 / 4) := by
  sorry

/-- **The dimerized chain at `δ = 1 / 4`.** Every even ring of `L ≥ 4` cells has a unique ground
state and gap `Δ_L > log 1.514 / (29 / 4) ≈ 0.0572`. -/
theorem gap_d14 {L : ℕ} (hL : 4 ≤ L) (hLe : Even L) :
    GapGT ((model (1 / 4)).spectrum L) (log 1.514 / (29 / 4)) := by
  have h := (model (1 / 4)).gapGT_of_chainCertificate (n₀ := 4) (β₀ := 29 / 4) (by norm_num)
    (by decide) (by norm_num) cert_d14 cert_d14_valid
    (by simpa [cert_d14, ChainCertificate.sAt] using seed_spatial_d14)
    (by simpa [cert_d14, ChainCertificate.tAt] using seed_thermal_d14) hL hLe
  convert h using 3
  simp [cert_d14]

end CertifiedGaps.DimerizedHeisenberg
