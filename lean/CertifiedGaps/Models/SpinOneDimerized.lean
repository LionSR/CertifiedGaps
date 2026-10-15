/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.SpinOneLargeD
import CertifiedGaps.Models.CellSwap

/-!
# The dimerized spin-one Heisenberg chain

`H = ∑ⱼ (1 - δ) S₂ⱼ · S₂ⱼ₊₁ + (1 + δ) S₂ⱼ₊₁ · S₂ⱼ₊₂` on a periodic ring of `2L` spin-one sites,
as a ring of `L` cells (the weak pairs), in the Cartesian basis. With the real antisymmetric
`Gα = ε_α = i Sα` (`SpinOneLargeD.G`), `-S · S = ∑_α Gα ⊗ Gα`, so (`model δ`):

* the strong bond between cells is `-h = ∑_α Cα ⊗ (F Cα F)` with `Cα = √(1 + δ) (1 ⊗ Gα)` and
  `F` the swap of the cell;
* the weak pair is the on-site term `o = -(1 - δ) ∑_α Gα ⊗ Gα`, which commutes with `F`.

So it is a reflection-positive ring with internal reflection `F`, with no basis change. The
seed is certified in `models/spin1_dimerized` (complete spectra of the rings of `4` and `8`
spins); it is the only `sorry`. The `δ = 7/10` certificate of that directory uses an odd seed
(`n₀ = 3` cells) and is not covered here.
-/

open Matrix Real
open scoped Kronecker

namespace CertifiedGaps.SpinOneDimerized

open SpinOneLargeD

/-- The dimerized spin-one chain as a reflection-positive ring of weak pairs. -/
noncomputable def model (δ : ℝ) : RPRing (Fin 3 × Fin 3) (Fin 3) where
  C α := Real.sqrt (1 + δ) • (1 ⊗ₖ G α)
  o := -(1 - δ) • ∑ α, G α ⊗ₖ G α
  C_symm α := transpose_smul_eq_or _ (transpose_kronecker_eq_or (.inl transpose_one)
    (.inr (G_transpose α)))
  o_symm := by rw [transpose_smul, transpose_sum_kronecker_self G fun α ↦ .inr (G_transpose α)]
  refl := cellSwap
  refl_symm := cellSwap_transpose
  refl_sq := cellSwap_mul_cellSwap
  refl_o := cellSwap_conj_smul_sum_kronecker_self _ _

/-- The rational chain certificate for `δ = 9/10`, seed `(n₀, β₀) = (2, 7/2)`. -/
def cert_d910 : ChainCertificate where
  s := [4972e-5, 1038e-5, 294e-6, 1857e-10, 6931e-17, 9608e-30, 1847e-55, 6823e-107, 9311e-210]
  t := [3206e-5, 3084e-6, 2108e-8, 8967e-13, 1609e-21, 5178e-39, 5363e-74, 5753e-144, 662e-283]
  r := 6339e-3
  ε := 9311e-210

/-- The certificate is valid, checked by the kernel. -/
theorem cert_d910_valid : cert_d910.Valid := by decide +kernel

/-- Seed: `P_s(2, 7/2) ≥ 1 - 0.04972`. Certified in `models/spin1_dimerized`. -/
theorem seed_spatial_d910 : 1 - (4972e-5 : ℝ) ≤ (model (9 / 10)).spatialPurity 2 (7 / 2) := by
  sorry

/-- Seed: `P_τ(4, 7/2) ≥ 1 - 0.03206`. Certified in `models/spin1_dimerized`. -/
theorem seed_thermal_d910 : 1 - (3206e-5 : ℝ) ≤ (model (9 / 10)).thermalPurity 4 (7 / 2) := by
  sorry

/-- **The dimerized spin-one chain at `δ = 9/10`.** Every even ring of `L ≥ 2` cells has a unique
ground state and gap `Δ_L > log 6.339 / (7/2) ≈ 0.5276`. -/
theorem gap_d910 {L : ℕ} (hL : 2 ≤ L) (hLe : Even L) :
    GapGT ((model (9 / 10)).spectrum L) (log 6.339 / (7 / 2)) :=
  (model (9 / 10)).gapGT_of_chainCertificate le_rfl (by decide) (by norm_num)
    cert_d910_valid (by norm_num [cert_d910, ChainCertificate.sAt, ChainCertificate.tAt])
    seed_spatial_d910 seed_thermal_d910 hL hLe

end CertifiedGaps.SpinOneDimerized
