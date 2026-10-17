/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.SpinOneLargeD
import CertifiedGaps.Models.SpinOneLargeD

/-! # Records: the spin-1 Heisenberg chain with single-ion anisotropy `D` -/

namespace Leaderboard.Records.SpinOneLargeD

open CertifiedGaps CertifiedGaps.SpinOneLargeD

/-- The spin-1 Heisenberg chain with single-ion anisotropy, `D = 4`. -/
theorem D_4 :
    Record.Seeded (model 4).spectrum 6 47731477e-7 (9 / 4) 6 6410493e-8 40529394e-9 :=
  RPRing.seeded (model 4) _ Chains.SpinOneLargeD.D_4_valid

/-- The spin-1 Heisenberg chain with single-ion anisotropy, `D = 1603/500`. -/
theorem D_1603_500 :
    Record.Seeded (model (1603 / 500)).spectrum 6 10067118e-7 (381 / 125) 6
      19164415e-8 46462453e-9 :=
  RPRing.seeded (model (1603 / 500)) _ Chains.SpinOneLargeD.D_1603_500_valid

/-- The spin-1 Heisenberg chain with single-ion anisotropy, `D ∈ [3.2059998, 3.206]`. The
seeds hold on the whole interval by a perturbation bound from `D = 3.206`
(`models/spin1_largeD`). -/
theorem D_interval_3206 :
    Record.SeededOn (fun D ↦ (model D).spectrum) (16029999 / 5000000) (1603 / 500) 6
      10029978e-7 (381 / 125) 6 19165007e-8 46476403e-9 :=
  fun D _ _ ↦ RPRing.seeded (model D) _ Chains.SpinOneLargeD.D_interval_3206_valid

/-- The spin-1 Heisenberg chain with single-ion anisotropy, `D ∈ [3.9979773, 4]`. The seeds
hold on the whole interval by a perturbation bound from `D = 4` (`models/spin1_largeD`). -/
theorem D_interval_4 :
    Record.SeededOn (fun D ↦ (model D).spectrum) (39979773 / 10000000) 4 6 10032634e-7
      (23 / 10) 6 11589869e-8 13768987e-8 :=
  fun D _ _ ↦ RPRing.seeded (model D) _ Chains.SpinOneLargeD.D_interval_4_valid

end Leaderboard.Records.SpinOneLargeD
