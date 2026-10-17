/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.SpinOneHeisenberg
import CertifiedGaps.Heisenberg.Main

/-!
# Records: the spin-1 Heisenberg chain

The record starts from the ring of 2304 sites at `β₀ = 784`, with the seed defects
`εₛ = 84513/10⁷` and `εₜ = 27493/(5 · 10⁶)` of `models/spin1/results/scalars.json` (six exact
doubling steps from the ring of 36 sites). They are smaller than the `1/100` of
`CertifiedGaps.Heisenberg.Main`, whose `sorry`s `seed_spatial` and `seed_thermal` state that
weaker seed. The proof uses the proved moment inequalities and the chain, never the seeds.
-/

namespace Leaderboard.Records.SpinOneHeisenberg

open CertifiedGaps CertifiedGaps.SpinOneHeisenberg

/-- The spin-1 Heisenberg chain. -/
theorem heisenberg :
    Record.Seeded spectrum 2304 41196029e-6 784 2304 84513e-7 54986e-7 :=
  Record.Seeded.of_chain spectrum (fun β hβ ↦ by
    rw [partitionFunction_spectrum]
    exact (spatialMomentBounds β hβ).mono (by norm_num)) _
    Chains.SpinOneHeisenberg.heisenberg_valid

end Leaderboard.Records.SpinOneHeisenberg
