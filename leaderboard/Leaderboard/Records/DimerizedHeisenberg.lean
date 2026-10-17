/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import CertifiedGaps.Models.DimerizedHeisenbergPhysical

/-! # Records: the dimerized spin-½ Heisenberg chain (physical frame, 2-site cells) -/

namespace Leaderboard.Records.DimerizedHeisenberg

open CertifiedGaps CertifiedGaps.DimerizedHeisenberg

/-- The dimerized spin-½ Heisenberg chain, `δ = 1/2`. -/
theorem delta_1_2 :
    Record.Seeded (physicalSpectrum (1 / 2)) 4 2433e-2 6 4 1133e-5 1224e-5 :=
  RPRing.seeded (model (1 / 2)) points[0].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

/-- The dimerized spin-½ Heisenberg chain, `δ = 3/10`. -/
theorem delta_3_10 :
    Record.Seeded (physicalSpectrum (3 / 10)) 4 3065e-3 (13 / 2) 4 9516e-5 518e-4 :=
  RPRing.seeded (model (3 / 10)) points[1].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

/-- The dimerized spin-½ Heisenberg chain, `δ = 1/4`. -/
theorem delta_1_4 :
    Record.Seeded (physicalSpectrum (1 / 4)) 4 1514e-3 (29 / 4) 4 1658e-4 5123e-5 :=
  RPRing.seeded (model (1 / 4)) points[2].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

end Leaderboard.Records.DimerizedHeisenberg
