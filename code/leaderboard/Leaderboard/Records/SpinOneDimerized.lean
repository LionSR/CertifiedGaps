/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.SpinOneDimerized
import CertifiedGaps.Models.SpinOneDimerized

/-! # Records: the dimerized spin-1 Heisenberg chain (2-site cells) -/

namespace Leaderboard.Records.SpinOneDimerized

open CertifiedGaps CertifiedGaps.SpinOneDimerized

/-- The dimerized spin-1 Heisenberg chain, `δ = 9/10`. -/
theorem delta_9_10 :
    Record.Seeded (model (9 / 10)).spectrum 2 63645965e-7 (7 / 2) 2 49715666e-9 3205449e-8 :=
  RPRing.seeded (model (9 / 10)) _ Chains.SpinOneDimerized.delta_9_10_valid

end Leaderboard.Records.SpinOneDimerized
