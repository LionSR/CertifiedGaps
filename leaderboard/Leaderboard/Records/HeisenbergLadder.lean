/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.RPRing
import CertifiedGaps.Models.HeisenbergLadderPhysical

/-! # Records: the two-leg spin-½ Heisenberg ladder (physical frame, rung cells) -/

namespace Leaderboard.Records.HeisenbergLadder

open CertifiedGaps CertifiedGaps.HeisenbergLadder

/-- The two-leg spin-½ Heisenberg ladder, `J⊥ = 4`. -/
theorem Jperp_4 :
    Record.Seeded (physicalSpectrum 4) 4 9033e-3 2 4 3509e-5 2463e-5 :=
  RPRing.seeded (model 4) points[0].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

end Leaderboard.Records.HeisenbergLadder
