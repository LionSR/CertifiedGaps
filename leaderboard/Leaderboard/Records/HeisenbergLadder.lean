/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.HeisenbergLadderPhysical

/-!
# Records for the two-leg spin-½ Heisenberg ladder (rings of rungs)

The hypothesis of each record is the seed inequality of `CertifiedGaps.HeisenbergLadder.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/heisenberg_ladder`.

Ring lengths count rungs. The records are about the physical ladder (`physicalSpectrum J⊥`); the
seeds are stated for the rotated frame `model J⊥`, which has the same spectrum on even rings.
Only `J⊥ = 4` has a Lean chain certificate; the points `J⊥ = 7/2, 3, 5/2` of
`models/heisenberg_ladder` have none yet, so they have no record.
-/

namespace Leaderboard.Records.HeisenbergLadder

open CertifiedGaps CertifiedGaps.HeisenbergLadder

/-- The two-leg spin-½ Heisenberg ladder (rings of rungs), `J⊥ = 4`.
Gap `Δ_L > log 9033e-3 / 2` for every even `L ≥ 4`, given the seeds
`P_s(4, 2) ≥ 1 - 3509e-5` and `P_τ(8, 2) ≥ 1 - 2463e-5`. -/
@[gap_record] theorem Jperp_4 : Record.Conditional
    (1 - (3509e-5 : ℝ) ≤ (model 4).spatialPurity 4 2 ∧
      1 - (2463e-5 : ℝ) ≤ (model 4).thermalPurity 8 2)
    (physicalSpectrum 4) 4 9033e-3 2 :=
  RPRing.conditional _ (by norm_num) (by decide) points[0].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)
    |>.congr fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

end Leaderboard.Records.HeisenbergLadder
