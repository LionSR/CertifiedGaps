/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.DimerizedHeisenbergPhysical

/-!
# Records for the dimerized spin-½ Heisenberg chain (rings of strong pairs)

The hypothesis of each record is the seed inequality of `CertifiedGaps.DimerizedHeisenberg.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/dimerized_heisenberg`.

Ring lengths count 2-site cells. The records are about the physical chain
(`physicalSpectrum δ`); the seeds are stated for the rotated frame `model δ`, which has the same
spectrum for `δ ≤ 1` and `L ≥ 2`.
-/

namespace Leaderboard.Records.DimerizedHeisenberg

open CertifiedGaps CertifiedGaps.DimerizedHeisenberg

/-- The dimerized spin-½ Heisenberg chain (rings of strong pairs), `δ = 1/2`.
Gap `Δ_L > log 2433e-2 / 6` for every even `L ≥ 4`, given the seeds
`P_s(4, 6) ≥ 1 - 1133e-5` and `P_τ(8, 6) ≥ 1 - 1224e-5`. -/
@[gap_record] theorem delta_1_2 : Record.Conditional
    (1 - (1133e-5 : ℝ) ≤ (model (1 / 2)).spatialPurity 4 6 ∧
      1 - (1224e-5 : ℝ) ≤ (model (1 / 2)).thermalPurity 8 6)
    (physicalSpectrum (1 / 2)) 4 2433e-2 6 :=
  RPRing.conditional _ (by norm_num) (by decide) points[0].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

/-- The dimerized spin-½ Heisenberg chain (rings of strong pairs), `δ = 3/10`.
Gap `Δ_L > log 3065e-3 / (13/2)` for every even `L ≥ 4`, given the seeds
`P_s(4, 13/2) ≥ 1 - 9516e-5` and `P_τ(8, 13/2) ≥ 1 - 518e-4`. -/
@[gap_record] theorem delta_3_10 : Record.Conditional
    (1 - (9516e-5 : ℝ) ≤ (model (3 / 10)).spatialPurity 4 (13 / 2) ∧
      1 - (518e-4 : ℝ) ≤ (model (3 / 10)).thermalPurity 8 (13 / 2))
    (physicalSpectrum (3 / 10)) 4 3065e-3 (13 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[1].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

/-- The dimerized spin-½ Heisenberg chain (rings of strong pairs), `δ = 1/4`.
Gap `Δ_L > log 1514e-3 / (29/4)` for every even `L ≥ 4`, given the seeds
`P_s(4, 29/4) ≥ 1 - 1658e-4` and `P_τ(8, 29/4) ≥ 1 - 5123e-5`. -/
@[gap_record] theorem delta_1_4 : Record.Conditional
    (1 - (1658e-4 : ℝ) ≤ (model (1 / 4)).spatialPurity 4 (29 / 4) ∧
      1 - (5123e-5 : ℝ) ≤ (model (1 / 4)).thermalPurity 8 (29 / 4))
    (physicalSpectrum (1 / 4)) 4 1514e-3 (29 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[2].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

end Leaderboard.Records.DimerizedHeisenberg
