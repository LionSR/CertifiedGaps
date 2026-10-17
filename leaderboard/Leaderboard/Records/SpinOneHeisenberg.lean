/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.TwoSeeded
import Leaderboard.Chains.SpinOneHeisenberg
import CertifiedGaps.Heisenberg.Main

/-!
# Records: the spin-1 Heisenberg chain

`heisenberg` starts from the ring of 2304 sites at `β₀ = 784`, with the seed defects
`εₛ = 84513/10⁷` and `εₜ = 27493/(5 · 10⁶)` of `models/spin1/results/scalars.json` (six exact
doubling steps from the ring of 36 sites). They are smaller than the `1/100` of
`CertifiedGaps.Heisenberg.Main`, whose `sorry`s `seed_spatial` and `seed_thermal` state that
weaker seed. The proof uses the proved moment inequalities and the chain, never the seeds.

`temperature_lifting` starts from the ring of 60 sites at `β₀ = 833/32`, whose seed defects
`models/spin1/experiments/temperature_lift.jl` bounds from the same finite inputs by raising the
temperature, and covers the rings of 36 to 58 sites from the seed at 36 sites, `β = 49/4`, in one
step. Its rate `r^{8/17}` at the bridge is checked as `r⁸ z¹⁷ < (1 - z)¹⁷`, since
`(49/4) / (833/32) = 8/17`.
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

/-- The spin-1 Heisenberg chain, from the temperature-lifted seed at 60 sites, with the seed at
36 sites for the rings of 36 to 58 sites. -/
theorem temperature_lifting :
    Record.TwoSeeded spectrum 36 66756076e-7 (833 / 32) 60 30530206e-9 47955678e-9 36 (49 / 4)
      4496169e-8 17686162e-8 :=
  Record.TwoSeeded.of_chain spectrum (fun β hβ ↦ by
    rw [partitionFunction_spectrum]
    exact (spatialMomentBounds β hβ).mono (by norm_num)) _
    Chains.SpinOneHeisenberg.temperature_lifting_valid 8 17

end Leaderboard.Records.SpinOneHeisenberg
