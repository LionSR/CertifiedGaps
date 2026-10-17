/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.Record
import CertifiedGaps.Heisenberg.Main

/-!
# The spin-1 Heisenberg record

Theorem 1 of the thermal-MPO note as a conditional record. The hypothesis is the pair of seed
inequalities `CertifiedGaps.SpinOneHeisenberg.seed_spatial` and `seed_thermal` (the two `sorry`s
of `CertifiedGaps.Heisenberg.Main`), stated verbatim. The proof uses the proved moment
inequalities `SpinOneHeisenberg.spatialMomentBounds` and the generic `gapGT_log_twenty`, never
the `sorry`-closed `thermalMPOInputs`, so the record carries no `sorryAx`.
-/

namespace Leaderboard.Records.SpinOneHeisenberg

open CertifiedGaps CertifiedGaps.SpinOneHeisenberg

/-- The spin-1 Heisenberg ring `∑ Sⱼ·Sⱼ₊₁`. Gap `Δ_L > log 20 / 784 ≈ 0.00382` for every even
`L ≥ 2304`, given the seeds `P_s(2304, 784) ≥ 0.99` and `P_τ(4608, 784) ≥ 0.99` (Proposition 5.5
of OpenAI's *The periodic spin-one Haldane gap*, reproduced in `models/spin1`). -/
@[gap_record] theorem heisenberg : Record.Conditional
    (99 / 100 ≤ spatialPurity 2304 784 ∧ 99 / 100 ≤ thermalPurity 4608 784)
    spectrum 2304 20 784 :=
  ⟨⟨by norm_num, by norm_num⟩, fun ⟨hs, ht⟩ _ hL hLe ↦ by
    have h : ThermalMPOInputs 2304 784 (1 / 100) :=
      ⟨fun β hβ ↦ (spatialMomentBounds β hβ).mono (by norm_num), by linarith,
        by norm_num; linarith⟩
    push_cast
    exact gapGT_log_twenty h hL hLe⟩

end Leaderboard.Records.SpinOneHeisenberg
