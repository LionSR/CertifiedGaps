/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.Record
import SeedKernel.Ising.Seed

/-!
# The unconditional Ising record

The longitudinal-field Ising chain at `g = 1/4`, `h = 1/20` is the first model whose seed is
proved in Lean (`SeedKernel.Ising.first_point_seed`, from kernel-checked integer arithmetic on the
rings of 6 and 12 sites), so its gap theorem `SeedKernel.Ising.ising_gap` has no hypothesis.
-/

namespace Leaderboard.Records.IsingLongitudinal

open CertifiedGaps

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/4`, `h = 1/20`, unconditionally.
Gap `Δ_L > log 13480 / (97/4) ≈ 0.392` for every even `L ≥ 6`
(`SeedKernel.Ising.ising_gap`). -/
@[gap_record] theorem g_1_4_h_1_20 :
    Record.Unconditional (IsingLongitudinal.model (1 / 4) (1 / 20)).spectrum 6 13480 (97 / 4) :=
  ⟨⟨by norm_num, by norm_num⟩, fun _ hL hLe ↦ by
    push_cast
    exact SeedKernel.Ising.ising_gap hL hLe⟩

end Leaderboard.Records.IsingLongitudinal
