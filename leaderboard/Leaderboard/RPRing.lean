/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Record
import CertifiedGaps.Ring.Gap

/-!
# Seeded records of reflection-positive rings

Every model in `lean/CertifiedGaps/Models` is a reflection-positive ring (`RPRing`) with a list
`points` of parameter points, each carrying a kernel-checked chain certificate `cert`. The
model's single `sorry`, `seeds`, states the two seed inequalities at every point. A record of a
point is `Record.Seeded` with the certificate's defects `s = cert.sAt 0`, `t = cert.tAt 0`, its
rate `r = cert.r` and the seed inverse temperature `β₀`. It does not use `seeds`, so it carries
no `sorryAx`. A proof of the seeds of a point turns its record unconditional
(`Record.Seeded.unconditional`).

`RPRing.seeded` proves such a record in one line: the equations between the record's literals
and the certificate's entries hold by `rfl`, and the positivity of `n₀`, `β`, `r - 1` by
`norm_num`. Records of a physical frame transfer with `Record.Seeded.congr`.
-/

namespace Leaderboard

open CertifiedGaps

/-- **Seeded records from chain certificates.** A valid chain certificate `c` at the even seed
length `n₀ ≥ 2` gives the seeded record `Δ_L > log r / β` for every even `L ≥ n₀`, with seeds
`P_s(n₀, β) ≥ 1 - s` and `P_τ(2n₀, β) ≥ 1 - t`, where `r`, `s`, `t` are the rate and the first
defects of `c`. -/
theorem RPRing.seeded {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] [Nonempty d]
    (R : RPRing d ι) {n₀ : ℕ} (c : ChainCertificate) (hc : c.Valid) {r β s t : ℚ}
    (hn₀ : 2 ≤ n₀ := by norm_num) (hn₀e : Even n₀ := by decide) (hr : c.r = r := by rfl)
    (hs : c.sAt 0 = s := by rfl) (ht : c.tAt 0 = t := by rfl) (hβ : 0 < β := by norm_num)
    (hr₁ : 1 < r := by norm_num) :
    Record.Seeded R.spectrum n₀ r β n₀ s t := by
  refine ⟨⟨hβ, hr₁⟩, fun ⟨hsp, hth⟩ L hL hLe ↦ ?_⟩
  subst hr hs ht
  exact gapGT_of_chainCertificate R.spectrum (by omega) hn₀e (by exact_mod_cast hβ)
    (fun β hβ ↦ (R.hasSpatialMomentBounds β hβ).mono hn₀) c hc hsp hth hL hLe

end Leaderboard
