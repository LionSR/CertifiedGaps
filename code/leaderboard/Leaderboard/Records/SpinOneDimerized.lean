/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.OddSeed.Record
import Leaderboard.Chains.SpinOneDimerized
import CertifiedGaps.Models.SpinOneDimerized

/-! # Records: the dimerized spin-1 Heisenberg chain (2-site cells) -/

namespace Leaderboard.Records.SpinOneDimerized

open CertifiedGaps CertifiedGaps.SpinOneDimerized

/-- The dimerized spin-1 Heisenberg chain, `δ = 9/10`. -/
theorem delta_9_10 :
    Record.Seeded (model (9 / 10)).spectrum 2 63645965e-7 (7 / 2) 2 49715666e-9 3205449e-8 :=
  RPRing.seeded (model (9 / 10)) _ Chains.SpinOneDimerized.delta_9_10_valid

/-- The dimerized spin-1 Heisenberg chain, `δ = 7/10`, with the odd seed `n₀ = 3` cells. The
transfer operator is positive semidefinite (`X_β = K K* ⪰ 0`,
`Leaderboard.OddSeed.Dimerized.hasSpatialMomentBoundsAll`), so the two-purity criterion applies at
the odd seed. -/
theorem delta_7_10 :
    Record.Seeded (model (7 / 10)).spectrum 4 25340426e-7 (17 / 4) 3 11719569e-8 48363911e-9 :=
  RPRing.seededAll (model (7 / 10)) _ Chains.SpinOneDimerized.delta_7_10_valid
    (hasSpatialMomentBoundsAll (7 / 10) (by norm_num) (by norm_num))

end Leaderboard.Records.SpinOneDimerized
