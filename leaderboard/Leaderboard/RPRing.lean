/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Record
import CertifiedGaps.Ring.Gap

/-!
# Records of reflection-positive rings from chain certificates

`RPRing.conditional` turns a kernel-checked chain certificate of a reflection-positive ring into
a conditional record whose hypothesis is exactly the pair of seed inequalities of the repository's
model files (`1 - εₛ ≤ P_s(n₀, β₀)` and `1 - εₜ ≤ P_τ(2n₀, β₀)`), so that a later proof of a seed
turns the record unconditional with `Record.Conditional.unconditional`.

The defects, the rate and the inverse temperature are passed with equations to the literals of
the record statement, so that the statement shows the numbers instead of `points[i].cert.r`.
-/

namespace Leaderboard

open CertifiedGaps

/-- **Records from chain certificates.** A valid chain certificate `c` at the even seed length
`n₀ ≥ 2` gives the conditional record `Δ_L > log r / β` for every even `L ≥ n₀`, with hypothesis
the two seed inequalities `1 - s ≤ P_s(n₀, β)` and `1 - t ≤ P_τ(2n₀, β)`. Here `r = c.r`, and
`s`, `t`, `β'` are the real forms of the certified defects and of `β`. -/
theorem RPRing.conditional {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] [Nonempty d]
    (R : RPRing d ι) {n₀ : ℕ} (hn₀ : 2 ≤ n₀) (hn₀e : Even n₀) (c : ChainCertificate)
    (hc : c.Valid) {r β : ℚ} (hr : c.r = r) (hβ : 0 < β) (hr₁ : 1 < r) {s t β' : ℝ}
    (hs : (c.sAt 0 : ℝ) = s) (ht : (c.tAt 0 : ℝ) = t) (hβ' : (β : ℝ) = β') :
    Record.Conditional (1 - s ≤ R.spatialPurity n₀ β' ∧ 1 - t ≤ R.thermalPurity (2 * n₀) β')
      R.spectrum n₀ r β := by
  refine ⟨⟨hβ, hr₁⟩, fun ⟨hsp, hth⟩ L hL hLe ↦ ?_⟩
  subst hβ'
  exact R.gapGT_of_chainCertificate hn₀ hn₀e (by exact_mod_cast hβ) hc
    ⟨hs, ht, by rw [hr]⟩ hsp hth hL hLe

end Leaderboard
