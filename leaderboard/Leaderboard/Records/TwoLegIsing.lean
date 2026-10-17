/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.TwoLegIsing

/-!
# Records for the two-leg transverse-field Ising ladder (rings of rungs)

The hypothesis of each record is the seed inequality of `CertifiedGaps.TwoLegIsing.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/two_leg_tfim`.

Ring lengths count rungs.
-/

namespace Leaderboard.Records.TwoLegIsing

open CertifiedGaps CertifiedGaps.TwoLegIsing

/-- The two-leg transverse-field Ising ladder (rings of rungs), `g = 3`.
Gap `Δ_L > log 2585e-3 / (3/2)` for every even `L ≥ 4`, given the seeds
`P_s(4, 3/2) ≥ 1 - 1125e-4` and `P_τ(8, 3/2) ≥ 1 - 5162e-5`. -/
@[gap_record] theorem g_3 : Record.Conditional
    (1 - (1125e-4 : ℝ) ≤ (model 3).spatialPurity 4 (3 / 2) ∧
      1 - (5162e-5 : ℝ) ≤ (model 3).thermalPurity 8 (3 / 2))
    (model 3).spectrum 4 2585e-3 (3 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[0].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

end Leaderboard.Records.TwoLegIsing
