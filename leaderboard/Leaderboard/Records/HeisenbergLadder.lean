/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.HeisenbergLadder
import CertifiedGaps.Models.HeisenbergLadderPhysical

/-! # Records: the two-leg spin-½ Heisenberg ladder (physical frame, rung cells) -/

namespace Leaderboard.Records.HeisenbergLadder

open CertifiedGaps CertifiedGaps.HeisenbergLadder

/-- The two-leg spin-½ Heisenberg ladder, `J⊥ = 4`. -/
theorem Jperp_4 :
    Record.Seeded (physicalSpectrum 4) 4 90679864e-7 2 4 3508321e-8 246242e-7 :=
  RPRing.seeded (model 4) _ Chains.HeisenbergLadder.Jperp_4_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The two-leg spin-½ Heisenberg ladder, `J⊥ = 7/2`. -/
theorem Jperp_7_2 :
    Record.Seeded (physicalSpectrum (7 / 2)) 4 5765561e-6 (9 / 4) 4 5806853e-8 2987064e-8 :=
  RPRing.seeded (model (7 / 2)) _ Chains.HeisenbergLadder.Jperp_7_2_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The two-leg spin-½ Heisenberg ladder, `J⊥ = 3`. -/
theorem Jperp_3 :
    Record.Seeded (physicalSpectrum 3) 4 31124723e-7 (5 / 2) 4 9896358e-8 449773e-7 :=
  RPRing.seeded (model 3) _ Chains.HeisenbergLadder.Jperp_3_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The two-leg spin-½ Heisenberg ladder, `J⊥ = 5/2`. -/
theorem Jperp_5_2 :
    Record.Seeded (physicalSpectrum (5 / 2)) 4 11356953e-7 3 4 18558434e-8 5114736e-8 :=
  RPRing.seeded (model (5 / 2)) _ Chains.HeisenbergLadder.Jperp_5_2_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

/-- The two-leg spin-½ Heisenberg ladder, `J⊥ = 163069/65536`. -/
theorem Jperp_163069_65536 :
    Record.Seeded (physicalSpectrum (163069 / 65536)) 4 1004547e-6 3 4 18772796e-8 5269516e-8 :=
  RPRing.seeded (model (163069 / 65536)) _ Chains.HeisenbergLadder.Jperp_163069_65536_valid
    |>.congr (by decide) fun _ _ hLe ↦ (physicalSpectrum_eq _ hLe).symm

end Leaderboard.Records.HeisenbergLadder
