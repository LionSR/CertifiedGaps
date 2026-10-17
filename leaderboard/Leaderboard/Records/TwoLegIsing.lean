/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.RPRing
import CertifiedGaps.Models.TwoLegIsing

/-! # Records: the two-leg transverse-field Ising ladder (rung cells) -/

namespace Leaderboard.Records.TwoLegIsing

open CertifiedGaps CertifiedGaps.TwoLegIsing

/-- The two-leg transverse-field Ising ladder, `g = 3`. -/
theorem g_3 :
    Record.Seeded (model 3).spectrum 4 2585e-3 (3 / 2) 4 1125e-4 5162e-5 :=
  RPRing.seeded (model 3) points[0].cert (points_valid _ (List.getElem_mem _))

end Leaderboard.Records.TwoLegIsing
