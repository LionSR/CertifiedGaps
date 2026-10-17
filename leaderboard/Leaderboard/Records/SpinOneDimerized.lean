/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import CertifiedGaps.Models.SpinOneDimerized

/-! # Records: the dimerized spin-1 Heisenberg chain (2-site cells) -/

namespace Leaderboard.Records.SpinOneDimerized

open CertifiedGaps CertifiedGaps.SpinOneDimerized

/-- The dimerized spin-1 Heisenberg chain, `δ = 9/10`. -/
theorem delta_9_10 :
    Record.Seeded (model (9 / 10)).spectrum 2 6339e-3 (7 / 2) 2 4972e-5 3206e-5 :=
  RPRing.seeded (model (9 / 10)) points[0].cert (points_valid _ (List.getElem_mem _))

end Leaderboard.Records.SpinOneDimerized
