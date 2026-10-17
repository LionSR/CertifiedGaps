/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.TwoLegIsing
import CertifiedGaps.Models.TwoLegIsing

/-! # Records: the two-leg transverse-field Ising ladder (rung cells) -/

namespace Leaderboard.Records.TwoLegIsing

open CertifiedGaps CertifiedGaps.TwoLegIsing

/-- The two-leg transverse-field Ising ladder, `g = 3`. -/
theorem g_3 :
    Record.Seeded (model 3).spectrum 4 25974347e-7 (3 / 2) 4 11240837e-8 51615239e-9 :=
  RPRing.seeded (model 3) _ Chains.TwoLegIsing.g_3_valid

end Leaderboard.Records.TwoLegIsing
