/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.SpinOneLargeD

/-!
# Records for the spin-1 Heisenberg ring with single-ion anisotropy `D ∑ (Sᶻⱼ)²`

The hypothesis of each record is the seed inequality of `CertifiedGaps.SpinOneLargeD.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/spin1_largeD`.

The Lean points are `D = 4` with seed `(6, 9/4)` and `D = 1603/500` with seed `(6, 381/125)`.
-/

namespace Leaderboard.Records.SpinOneLargeD

open CertifiedGaps CertifiedGaps.SpinOneLargeD

/-- The spin-1 Heisenberg ring with single-ion anisotropy `D ∑ (Sᶻⱼ)²`, `D = 4`.
Gap `Δ_L > log 4755e-3 / (9/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 9/4) ≥ 1 - 6411e-5` and `P_τ(12, 9/4) ≥ 1 - 4053e-5`. -/
@[gap_record] theorem D_4 : Record.Conditional
    (1 - (6411e-5 : ℝ) ≤ (model 4).spatialPurity 6 (9 / 4) ∧
      1 - (4053e-5 : ℝ) ≤ (model 4).thermalPurity 12 (9 / 4))
    (model 4).spectrum 6 4755e-3 (9 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[0].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The spin-1 Heisenberg ring with single-ion anisotropy `D ∑ (Sᶻⱼ)²`, `D = 1603/500`.
Gap `Δ_L > log 10065429e-7 / (381/125)` for every even `L ≥ 6`, given the seeds
`P_s(6, 381/125) ≥ 1 - 19164415e-8` and `P_τ(12, 381/125) ≥ 1 - 46462453e-9`. -/
@[gap_record] theorem D_1603_500 : Record.Conditional
    (1 - (19164415e-8 : ℝ) ≤ (model (1603 / 500)).spatialPurity 6 (381 / 125) ∧
      1 - (46462453e-9 : ℝ) ≤ (model (1603 / 500)).thermalPurity 12 (381 / 125))
    (model (1603 / 500)).spectrum 6 10065429e-7 (381 / 125) :=
  RPRing.conditional _ (by norm_num) (by decide) points[1].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

end Leaderboard.Records.SpinOneLargeD
