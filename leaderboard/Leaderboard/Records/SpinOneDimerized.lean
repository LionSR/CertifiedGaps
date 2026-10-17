/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.SpinOneDimerized

/-!
# Records for the dimerized spin-1 Heisenberg chain (rings of weak pairs)

The hypothesis of each record is the seed inequality of `CertifiedGaps.SpinOneDimerized.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/spin1_dimerized`.

Ring lengths count 2-site cells. The `δ = 7/10` certificate of `models/spin1_dimerized` uses an
odd seed (`n₀ = 3` cells), which the Lean model does not cover, so it has no record.
-/

namespace Leaderboard.Records.SpinOneDimerized

open CertifiedGaps CertifiedGaps.SpinOneDimerized

/-- The dimerized spin-1 Heisenberg chain (rings of weak pairs), `δ = 9/10`.
Gap `Δ_L > log 6339e-3 / (7/2)` for every even `L ≥ 2`, given the seeds
`P_s(2, 7/2) ≥ 1 - 4972e-5` and `P_τ(4, 7/2) ≥ 1 - 3206e-5`. -/
@[gap_record] theorem delta_9_10 : Record.Conditional
    (1 - (4972e-5 : ℝ) ≤ (model (9 / 10)).spatialPurity 2 (7 / 2) ∧
      1 - (3206e-5 : ℝ) ≤ (model (9 / 10)).thermalPurity 4 (7 / 2))
    (model (9 / 10)).spectrum 2 6339e-3 (7 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[0].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

end Leaderboard.Records.SpinOneDimerized
