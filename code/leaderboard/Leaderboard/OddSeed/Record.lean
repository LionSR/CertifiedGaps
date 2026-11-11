/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.OddSeed.Criterion
import Leaderboard.OddSeed.Dimerized
import Leaderboard.RPRing

/-!
# Seeded records at an odd seed length

`RPRing.seededAll` is `RPRing.seeded` for a reflection-positive ring whose transfer operator is
positive semidefinite (`HasSpatialMomentBoundsAll`, `Leaderboard.OddSeed.Dimerized`): the seed
length `n₀` may be odd, as for the dimerized spin-1 chain at `δ = 7/10` (`n₀ = 3` cells). The chain
certificate is checked by the kernel as usual; the seed inequalities at `(n₀, β)` remain the
hypotheses, certified by computer in `code/models/<slug>/`.

## Main results

* `RPRing.seededAll`: a seeded record from an all-lengths moment form and a valid chain certificate.
-/

namespace Leaderboard

open CertifiedGaps

/-- **Seeded records of reflection-positive rings with a positive transfer operator.** `n₀` may be
odd: the moment inequalities hold at every length, so `gapGT_of_chainCertificate_all` applies. -/
theorem RPRing.seededAll {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] [Nonempty d]
    (R : RPRing d ι) {n₀ : ℕ} (c : ChainCertificate) (hc : c.Valid)
    (hZ : HasSpatialMomentBoundsAll R.spectrum 2) (hn₀ : 2 ≤ n₀ := by norm_num) {r β s t : ℚ}
    (hr : c.r = r := by rfl) (hs : c.sAt 0 = s := by rfl) (ht : c.tAt 0 = t := by rfl)
    (hβ : 0 < β := by norm_num) (hr₁ : 1 < r := by norm_num) {L₀ : ℕ} (hL₀ : n₀ ≤ L₀ := by omega) :
    Record.Seeded R.spectrum L₀ r β n₀ s t := by
  refine ⟨⟨hβ, hr₁⟩, fun hseed L hL hLe ↦ ?_⟩
  obtain ⟨hsp, hth⟩ := hseed
  subst hr hs ht
  exact gapGT_of_chainCertificate_all R.spectrum (by omega) (by exact_mod_cast hβ)
    (fun β' hβ' ↦ (hZ β' hβ').mono hn₀) c hc hsp hth (by omega) hLe

end Leaderboard
