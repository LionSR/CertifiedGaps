/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.StaggeredHeisenberg
import CertifiedGaps.Models.StaggeredHeisenbergPhysical

/-! # Records: the spin-½ Heisenberg chain in a staggered field (physical frame) -/

namespace Leaderboard.Records.StaggeredHeisenberg

open CertifiedGaps CertifiedGaps.StaggeredHeisenberg

/-- The spin-½ Heisenberg chain in a staggered field, `h = 1`. -/
theorem h_1 :
    Record.Seeded (physicalSpectrum 1) 8 41501591e-5 (21 / 4) 8 69368678e-11 73264843e-11 :=
  RPRing.seeded (model 1) _ Chains.StaggeredHeisenberg.h_1_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 3/4`. -/
theorem h_3_4 :
    Record.Seeded (physicalSpectrum (3 / 4)) 8 10739897e-5 (21 / 4) 8 23425197e-10 31223351e-10 :=
  RPRing.seeded (model (3 / 4)) _ Chains.StaggeredHeisenberg.h_3_4_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 1/2`. -/
theorem h_1_2 :
    Record.Seeded (physicalSpectrum (1 / 2)) 8 26344587e-6 (11 / 2) 8 10740847e-9 11157006e-9 :=
  RPRing.seeded (model (1 / 2)) _ Chains.StaggeredHeisenberg.h_1_2_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 7/20`. -/
theorem h_7_20 :
    Record.Seeded (physicalSpectrum (7 / 20)) 8 90713109e-7 (23 / 4) 8 32653342e-9 2743591e-8 :=
  RPRing.seeded (model (7 / 20)) _ Chains.StaggeredHeisenberg.h_7_20_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The spin-½ Heisenberg chain in a staggered field, `h = 1/4`. -/
theorem h_1_4 :
    Record.Seeded (physicalSpectrum (1 / 4)) 8 38208854e-7 (25 / 4) 8 79866451e-9 44763872e-9 :=
  RPRing.seeded (model (1 / 4)) _ Chains.StaggeredHeisenberg.h_1_4_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

end Leaderboard.Records.StaggeredHeisenberg
