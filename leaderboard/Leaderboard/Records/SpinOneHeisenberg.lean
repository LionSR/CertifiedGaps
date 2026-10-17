/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Record
import CertifiedGaps.Heisenberg.Main

/-!
# Records: the spin-1 Heisenberg chain

The seeds `P_s(2304, 784), P_τ(4608, 784) ≥ 0.99` are the `sorry`s `seed_spatial` and
`seed_thermal` of `CertifiedGaps.Heisenberg.Main`, stated for the spectrum. The proof uses the
proved moment inequalities and `gapGT_of_seeds`, never the seeds.
-/

namespace Leaderboard.Records.SpinOneHeisenberg

open CertifiedGaps CertifiedGaps.SpinOneHeisenberg

/-- The spin-1 Heisenberg chain. -/
theorem heisenberg :
    Record.Seeded spectrum 2304 20 784 2304 (1 / 100) (1 / 100) := by
  refine ⟨⟨by norm_num, by norm_num⟩, fun ⟨hs, ht⟩ L hL hLe ↦ ?_⟩
  have hZ : HasSpatialMomentBounds spectrum 2304 := fun β hβ ↦ by
    rw [partitionFunction_spectrum]
    exact (spatialMomentBounds β hβ).mono (by norm_num)
  have hs' : 1 - (1 / 100 : ℝ) ≤ spatialPurity spectrum 2304 784 := by
    push_cast at hs; exact hs
  have ht' : 1 - (1 / 100 : ℝ) ≤ thermalPurity spectrum (2 * 2304) 784 := by
    push_cast at ht; exact ht
  have h := gapGT_of_seeds spectrum (by norm_num) (by decide) (by norm_num) hZ (by norm_num)
    le_rfl hs' ht' hL hLe
  norm_num at h ⊢
  exact h

end Leaderboard.Records.SpinOneHeisenberg
