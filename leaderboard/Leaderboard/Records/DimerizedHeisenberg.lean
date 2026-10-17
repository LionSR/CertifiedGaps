/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.DimerizedHeisenberg
import CertifiedGaps.Models.DimerizedHeisenbergPhysical

/-! # Records: the dimerized spin-½ Heisenberg chain (physical frame, 2-site cells) -/

namespace Leaderboard.Records.DimerizedHeisenberg

open CertifiedGaps CertifiedGaps.DimerizedHeisenberg

/-- The dimerized spin-½ Heisenberg chain, `δ = 1/2`. -/
theorem delta_1_2 :
    Record.Seeded (physicalSpectrum (1 / 2)) 4 24741506e-6 6 4 11267e-6 12008e-6 :=
  RPRing.seeded (model (1 / 2)) _ Chains.DimerizedHeisenberg.delta_1_2_valid
    |>.congr (by decide) fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

/-- The dimerized spin-½ Heisenberg chain, `δ = 3/10`. -/
theorem delta_3_10 :
    Record.Seeded (physicalSpectrum (3 / 10)) 4 30824555e-7 (13 / 2) 4 95116e-6 51608e-6 :=
  RPRing.seeded (model (3 / 10)) _ Chains.DimerizedHeisenberg.delta_3_10_valid
    |>.congr (by decide) fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

/-- The dimerized spin-½ Heisenberg chain, `δ = 1/4`. -/
theorem delta_1_4 :
    Record.Seeded (physicalSpectrum (1 / 4)) 4 15251257e-7 (29 / 4) 4 165679e-6 51038e-6 :=
  RPRing.seeded (model (1 / 4)) _ Chains.DimerizedHeisenberg.delta_1_4_valid
    |>.congr (by decide) fun _ hL _ ↦ (physicalSpectrum_eq (by norm_num) (by omega)).symm

end Leaderboard.Records.DimerizedHeisenberg
