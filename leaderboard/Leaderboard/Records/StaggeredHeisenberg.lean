/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.RPRing
import CertifiedGaps.Models.StaggeredHeisenbergPhysical

/-! # Records: the spin-½ Heisenberg chain in a staggered field (physical frame) -/

namespace Leaderboard.Records.StaggeredHeisenberg

open CertifiedGaps CertifiedGaps.StaggeredHeisenberg

/-- The spin-½ Heisenberg chain in a staggered field, `h = 1`. -/
theorem h_1 :
    Record.Seeded (physicalSpectrum 1) 8 4134e-1 (21 / 4) 8 6937e-7 7327e-7 :=
  RPRing.seeded (model 1) points[0].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 3/4`. -/
theorem h_3_4 :
    Record.Seeded (physicalSpectrum (3 / 4)) 8 1069e-1 (21 / 4) 8 2343e-6 3123e-6 :=
  RPRing.seeded (model (3 / 4)) points[1].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 1/2`. -/
theorem h_1_2 :
    Record.Seeded (physicalSpectrum (1 / 2)) 8 2623e-2 (11 / 2) 8 1075e-5 1116e-5 :=
  RPRing.seeded (model (1 / 2)) points[2].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 7/20`. -/
theorem h_7_20 :
    Record.Seeded (physicalSpectrum (7 / 20)) 8 9035e-3 (23 / 4) 8 3266e-5 2744e-5 :=
  RPRing.seeded (model (7 / 20)) points[3].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 1/4`. -/
theorem h_1_4 :
    Record.Seeded (physicalSpectrum (1 / 4)) 8 3806e-3 (25 / 4) 8 7987e-5 4477e-5 :=
  RPRing.seeded (model (1 / 4)) points[4].cert (points_valid _ (List.getElem_mem _))
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

end Leaderboard.Records.StaggeredHeisenberg
