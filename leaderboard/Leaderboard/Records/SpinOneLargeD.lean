/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.RPRing
import CertifiedGaps.Models.SpinOneLargeD

/-! # Records: the spin-1 Heisenberg chain with single-ion anisotropy `D` -/

namespace Leaderboard.Records.SpinOneLargeD

open CertifiedGaps CertifiedGaps.SpinOneLargeD

/-- The spin-1 Heisenberg chain with single-ion anisotropy, `D = 4`. -/
theorem D_4 :
    Record.Seeded (model 4).spectrum 6 4755e-3 (9 / 4) 6 6411e-5 4053e-5 :=
  RPRing.seeded (model 4) points[0].cert (points_valid _ (List.getElem_mem _))

/-- The spin-1 Heisenberg chain with single-ion anisotropy, `D = 1603/500`. -/
theorem D_1603_500 :
    Record.Seeded (model (1603 / 500)).spectrum 6 10065429e-7 (381 / 125) 6 19164415e-8
      46462453e-9 :=
  RPRing.seeded (model (1603 / 500)) points[1].cert (points_valid _ (List.getElem_mem _))

end Leaderboard.Records.SpinOneLargeD
