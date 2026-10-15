/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Gap
import CertifiedGaps.Models.SpinHalf

/-!
# The spin-one-half Heisenberg chain in a staggered field

`H_L = ∑ⱼ Sⱼ · Sⱼ₊₁ - h ∑ⱼ (-1)ʲ Sᶻⱼ` on a periodic ring. The rotation `e^{iπSʸ}` on every other
site maps every bond to `-SˣSˣ + SʸSʸ - SᶻSᶻ` and the staggered field to the uniform field
`-h Sᶻ` (since `e^{iπSʸ} Sᶻ e^{-iπSʸ} = -Sᶻ`); on an even ring it is a unitary, so the spectrum
is unchanged. This file states the result for the rotated chain `model h`: a reflection-positive
ring with bond factors `Sˣ, iSʸ, Sᶻ` (`SpinHalf.P`) and on-site term `-h Sᶻ`.

The seeds are certified in `models/heisenberg_staggered_field` (complete spectra of the rings of
`8` and `16` sites with rigorous eigenvalue enclosures); they are the only `sorry`s.
-/

open Matrix Real

namespace CertifiedGaps.StaggeredHeisenberg

open SpinHalf

/-- The rotated chain in the field `h` as a reflection-positive ring. -/
noncomputable def model (h : ℝ) : RPRing (Fin 2) (Fin 3) where
  C := P
  o := -h • P 2
  C_symm := P_symm
  o_symm := by
    rw [transpose_smul]
    congr 1
    ext i j; fin_cases i <;> fin_cases j <;> norm_num [P]

/-! ### `h = 1` -/

/-- The rational chain certificate for `h = 1`, seed `(n₀, β₀) = (8, 21 / 4)`. -/
def cert_h1 : ChainCertificate where
  s := [6937e-7, 2254e-9, 1562e-14, 5636e-25, 6477e-46, 8394e-88, 141e-170, 3977e-339, 3164e-674]
  t := [7327e-7, 108e-8, 2333e-15, 1089e-26, 2372e-49, 1126e-94, 2536e-185, 1287e-366, 3313e-729]
  r := 4134e-1
  ε := 3164e-674

/-- The certificate is valid, checked by the kernel. -/
theorem cert_h1_valid : cert_h1.Valid := by decide +kernel

/-- Seed: `P_s(8, 21 / 4) ≥ 1 - 6937e-7`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_spatial_h1 : 1 - (6937e-7 : ℝ) ≤ (model 1).spatialPurity 8 (21 / 4) := by
  sorry

/-- Seed: `P_τ(16, 21 / 4) ≥ 1 - 7327e-7`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_thermal_h1 : 1 - (7327e-7 : ℝ) ≤ (model 1).thermalPurity 16 (21 / 4) := by
  sorry

/-- **The staggered-field chain at `h = 1`.** Every even ring of `L ≥ 8` sites has a unique
ground state and gap `Δ_L > log 413.4 / (21/4) ≈ 1.1475`. -/
theorem gap_h1 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT ((model 1).spectrum L) (log 413.4 / (21 / 4)) :=
  (model 1).gapGT_of_chainCertificate (by norm_num) (by decide) (by norm_num)
    cert_h1_valid (by norm_num [cert_h1, ChainCertificate.sAt, ChainCertificate.tAt])
    seed_spatial_h1 seed_thermal_h1 hL hLe

/-! ### `h = 3/4` -/

/-- The rational chain certificate for `h = 3/4`, seed `(n₀, β₀) = (8, 21 / 4)`. -/
def cert_h34 : ChainCertificate where
  s := [2343e-6, 3082e-8, 3324e-12, 2768e-20, 1603e-36, 515e-68, 5305e-134, 5629e-264, 6338e-524]
  t := [3123e-6, 1989e-8, 7915e-13, 1253e-21, 3141e-39, 1974e-74, 7794e-145, 1215e-285, 2953e-567]
  r := 1069e-1
  ε := 6338e-524

/-- The certificate is valid, checked by the kernel. -/
theorem cert_h34_valid : cert_h34.Valid := by decide +kernel

/-- Seed: `P_s(8, 21 / 4) ≥ 1 - 2343e-6`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_spatial_h34 : 1 - (2343e-6 : ℝ) ≤ (model (3 / 4)).spatialPurity 8 (21 / 4) := by
  sorry

/-- Seed: `P_τ(16, 21 / 4) ≥ 1 - 3123e-6`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_thermal_h34 : 1 - (3123e-6 : ℝ) ≤ (model (3 / 4)).thermalPurity 16 (21 / 4) := by
  sorry

/-- **The staggered-field chain at `h = 3/4`.** Every even ring of `L ≥ 8` sites has a unique
ground state and gap `Δ_L > log 106.9 / (21/4) ≈ 0.8899`. -/
theorem gap_h34 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT ((model (3 / 4)).spectrum L) (log 106.9 / (21 / 4)) :=
  (model (3 / 4)).gapGT_of_chainCertificate (by norm_num) (by decide) (by norm_num)
    cert_h34_valid (by norm_num [cert_h34, ChainCertificate.sAt, ChainCertificate.tAt])
    seed_spatial_h34 seed_thermal_h34 hL hLe

/-! ### `h = 1/2` -/

/-- The rational chain certificate for `h = 1/2`, seed `(n₀, β₀) = (8, 11 / 2)`. -/
def cert_h12 : ChainCertificate where
  s := [1075e-5, 5573e-7, 9615e-10, 2144e-15, 9381e-27, 1761e-49, 6203e-95, 7696e-186, 1185e-367]
  t := [1116e-5, 2708e-7, 1474e-10, 4346e-17, 3778e-30, 2855e-56, 1631e-108, 5321e-213, 5663e-422]
  r := 2623e-2
  ε := 1185e-367

/-- The certificate is valid, checked by the kernel. -/
theorem cert_h12_valid : cert_h12.Valid := by decide +kernel

/-- Seed: `P_s(8, 11 / 2) ≥ 1 - 1075e-5`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_spatial_h12 : 1 - (1075e-5 : ℝ) ≤ (model (1 / 2)).spatialPurity 8 (11 / 2) := by
  sorry

/-- Seed: `P_τ(16, 11 / 2) ≥ 1 - 1116e-5`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_thermal_h12 : 1 - (1116e-5 : ℝ) ≤ (model (1 / 2)).thermalPurity 16 (11 / 2) := by
  sorry

/-- **The staggered-field chain at `h = 1/2`.** Every even ring of `L ≥ 8` sites has a unique
ground state and gap `Δ_L > log 26.23 / (11/2) ≈ 0.5940`. -/
theorem gap_h12 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT ((model (1 / 2)).spectrum L) (log 26.23 / (11 / 2)) :=
  (model (1 / 2)).gapGT_of_chainCertificate (by norm_num) (by decide) (by norm_num)
    cert_h12_valid (by norm_num [cert_h12, ChainCertificate.sAt, ChainCertificate.tAt])
    seed_spatial_h12 seed_thermal_h12 hL hLe

/-! ### `h = 7/20` -/

/-- The rational chain certificate for `h = 7/20`, seed `(n₀, β₀) = (8, 23 / 4)`. -/
def cert_h720 : ChainCertificate where
  s := [3266e-5, 4883e-6, 6973e-8, 1087e-11, 2391e-19, 1144e-34, 2618e-65, 1371e-126, 376e-248]
  t := [2744e-5, 1948e-6, 791e-8, 1254e-13, 3146e-23, 198e-41, 7841e-81, 123e-156, 3026e-311]
  r := 9035e-3
  ε := 376e-248

/-- The certificate is valid, checked by the kernel. -/
theorem cert_h720_valid : cert_h720.Valid := by decide +kernel

/-- Seed: `P_s(8, 23 / 4) ≥ 1 - 3266e-5`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_spatial_h720 : 1 - (3266e-5 : ℝ) ≤ (model (7 / 20)).spatialPurity 8 (23 / 4) := by
  sorry

/-- Seed: `P_τ(16, 23 / 4) ≥ 1 - 2744e-5`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_thermal_h720 : 1 - (2744e-5 : ℝ) ≤ (model (7 / 20)).thermalPurity 16 (23 / 4) := by
  sorry

/-- **The staggered-field chain at `h = 7/20`.** Every even ring of `L ≥ 8` sites has a unique
ground state and gap `Δ_L > log 9.035 / (23/4) ≈ 0.3828`. -/
theorem gap_h720 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT ((model (7 / 20)).spectrum L) (log 9.035 / (23 / 4)) :=
  (model (7 / 20)).gapGT_of_chainCertificate (by norm_num) (by decide) (by norm_num)
    cert_h720_valid (by norm_num [cert_h720, ChainCertificate.sAt, ChainCertificate.tAt])
    seed_spatial_h720 seed_thermal_h720 hL hLe

/-! ### `h = 1/4` -/

/-- The rational chain certificate for `h = 1/4`, seed `(n₀, β₀) = (8, 25 / 4)`. -/
def cert_h14 : ChainCertificate where
  s := [7987e-5, 2797e-5, 2248e-6, 1099e-8, 2431e-13, 1182e-22, 2795e-41, 1563e-78, 4886e-153]
  t := [4477e-5, 8125e-6, 1757e-7, 657e-10, 8665e-18, 1502e-31, 4513e-59, 4074e-114, 332e-223]
  r := 3806e-3
  ε := 4886e-153

/-- The certificate is valid, checked by the kernel. -/
theorem cert_h14_valid : cert_h14.Valid := by decide +kernel

/-- Seed: `P_s(8, 25 / 4) ≥ 1 - 7987e-5`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_spatial_h14 : 1 - (7987e-5 : ℝ) ≤ (model (1 / 4)).spatialPurity 8 (25 / 4) := by
  sorry

/-- Seed: `P_τ(16, 25 / 4) ≥ 1 - 4477e-5`. Certified in `models/heisenberg_staggered_field`. -/
theorem seed_thermal_h14 : 1 - (4477e-5 : ℝ) ≤ (model (1 / 4)).thermalPurity 16 (25 / 4) := by
  sorry

/-- **The staggered-field chain at `h = 1/4`.** Every even ring of `L ≥ 8` sites has a unique
ground state and gap `Δ_L > log 3.806 / (25/4) ≈ 0.2139`. -/
theorem gap_h14 {L : ℕ} (hL : 8 ≤ L) (hLe : Even L) :
    GapGT ((model (1 / 4)).spectrum L) (log 3.806 / (25 / 4)) :=
  (model (1 / 4)).gapGT_of_chainCertificate (by norm_num) (by decide) (by norm_num)
    cert_h14_valid (by norm_num [cert_h14, ChainCertificate.sAt, ChainCertificate.tAt])
    seed_spatial_h14 seed_thermal_h14 hL hLe

end CertifiedGaps.StaggeredHeisenberg
