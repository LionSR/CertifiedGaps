/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.StaggeredHeisenbergPhysical

/-!
# Records for the spin-½ Heisenberg ring `∑ Sⱼ·Sⱼ₊₁ - h ∑ (-1)ʲ Sᶻⱼ`

The hypothesis of each record is the seed inequality of `CertifiedGaps.StaggeredHeisenberg.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/heisenberg_staggered_field`.

The records are about the physical chain (`physicalSpectrum h`); the seeds are stated for the
rotated reflection-positive frame `model h`, which has the same spectrum on even rings.
-/

namespace Leaderboard.Records.StaggeredHeisenberg

open CertifiedGaps CertifiedGaps.StaggeredHeisenberg

/-- The spin-½ Heisenberg ring `∑ Sⱼ·Sⱼ₊₁ - h ∑ (-1)ʲ Sᶻⱼ`, `h = 1`.
Gap `Δ_L > log 4134e-1 / (21/4)` for every even `L ≥ 8`, given the seeds
`P_s(8, 21/4) ≥ 1 - 6937e-7` and `P_τ(16, 21/4) ≥ 1 - 7327e-7`. -/
@[gap_record] theorem h_1 : Record.Conditional
    (1 - (6937e-7 : ℝ) ≤ (model 1).spatialPurity 8 (21 / 4) ∧
      1 - (7327e-7 : ℝ) ≤ (model 1).thermalPurity 16 (21 / 4))
    (physicalSpectrum 1) 8 4134e-1 (21 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[0].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg ring `∑ Sⱼ·Sⱼ₊₁ - h ∑ (-1)ʲ Sᶻⱼ`, `h = 3/4`.
Gap `Δ_L > log 1069e-1 / (21/4)` for every even `L ≥ 8`, given the seeds
`P_s(8, 21/4) ≥ 1 - 2343e-6` and `P_τ(16, 21/4) ≥ 1 - 3123e-6`. -/
@[gap_record] theorem h_3_4 : Record.Conditional
    (1 - (2343e-6 : ℝ) ≤ (model (3 / 4)).spatialPurity 8 (21 / 4) ∧
      1 - (3123e-6 : ℝ) ≤ (model (3 / 4)).thermalPurity 16 (21 / 4))
    (physicalSpectrum (3 / 4)) 8 1069e-1 (21 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[1].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg ring `∑ Sⱼ·Sⱼ₊₁ - h ∑ (-1)ʲ Sᶻⱼ`, `h = 1/2`.
Gap `Δ_L > log 2623e-2 / (11/2)` for every even `L ≥ 8`, given the seeds
`P_s(8, 11/2) ≥ 1 - 1075e-5` and `P_τ(16, 11/2) ≥ 1 - 1116e-5`. -/
@[gap_record] theorem h_1_2 : Record.Conditional
    (1 - (1075e-5 : ℝ) ≤ (model (1 / 2)).spatialPurity 8 (11 / 2) ∧
      1 - (1116e-5 : ℝ) ≤ (model (1 / 2)).thermalPurity 16 (11 / 2))
    (physicalSpectrum (1 / 2)) 8 2623e-2 (11 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[2].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg ring `∑ Sⱼ·Sⱼ₊₁ - h ∑ (-1)ʲ Sᶻⱼ`, `h = 7/20`.
Gap `Δ_L > log 9035e-3 / (23/4)` for every even `L ≥ 8`, given the seeds
`P_s(8, 23/4) ≥ 1 - 3266e-5` and `P_τ(16, 23/4) ≥ 1 - 2744e-5`. -/
@[gap_record] theorem h_7_20 : Record.Conditional
    (1 - (3266e-5 : ℝ) ≤ (model (7 / 20)).spatialPurity 8 (23 / 4) ∧
      1 - (2744e-5 : ℝ) ≤ (model (7 / 20)).thermalPurity 16 (23 / 4))
    (physicalSpectrum (7 / 20)) 8 9035e-3 (23 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[3].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg ring `∑ Sⱼ·Sⱼ₊₁ - h ∑ (-1)ʲ Sᶻⱼ`, `h = 1/4`.
Gap `Δ_L > log 3806e-3 / (25/4)` for every even `L ≥ 8`, given the seeds
`P_s(8, 25/4) ≥ 1 - 7987e-5` and `P_τ(16, 25/4) ≥ 1 - 4477e-5`. -/
@[gap_record] theorem h_1_4 : Record.Conditional
    (1 - (7987e-5 : ℝ) ≤ (model (1 / 4)).spatialPurity 8 (25 / 4) ∧
      1 - (4477e-5 : ℝ) ≤ (model (1 / 4)).thermalPurity 16 (25 / 4))
    (physicalSpectrum (1 / 4)) 8 3806e-3 (25 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[4].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

end Leaderboard.Records.StaggeredHeisenberg
