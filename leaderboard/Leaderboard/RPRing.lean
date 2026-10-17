/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.Record
import CertifiedGaps.Ring.Gap

/-!
# Seeded records from chain certificates

A seeded record follows from a valid `ChainCertificate` `c` (the finite part of the criterion,
checked by `decide +kernel`) and the moment inequalities of the spectrum: its defects are
`s = c.sAt 0`, `t = c.tAt 0` and its rate is `r = c.r`. `Record.Seeded.of_chain` proves this for
any spectrum with the moment inequalities, and `RPRing.seeded` for a reflection-positive ring,
whose moment inequalities are proved once in `CertifiedGaps.Ring`. The chains of the records
are generated in `Leaderboard.Chains` by `tools/leaderboard/lean_chains.jl` from the seed
defects of the computer certificates.

In both lemmas the equations between the record's literals and the certificate's entries hold
by `rfl`, and the positivity of `n₀`, `β`, `r - 1` by `norm_num`. Records of a physical frame
transfer with `Record.Seeded.congr`; a proof of the seeds makes a record unconditional
(`Record.Seeded.unconditional`).
-/

namespace Leaderboard

open CertifiedGaps

/-- **Seeded records from chain certificates.** Under the moment inequalities above the even
seed length `n₀ ≥ 2`, a valid chain certificate `c` gives the seeded record `Δ_L > log r / β` for
every even `L ≥ n₀`, with seeds `P_s(n₀, β) ≥ 1 - s` and `P_τ(2n₀, β) ≥ 1 - t`, where `r`, `s`,
`t` are the rate and the first defects of `c`. -/
theorem Record.Seeded.of_chain {d : Type} [Fintype d] [DecidableEq d] [Nonempty d]
    (spectrum : (L : ℕ) → (Fin L → d) → ℝ) {n₀ : ℕ} (hZ : HasSpatialMomentBounds spectrum n₀)
    (c : ChainCertificate) (hc : c.Valid) {r β s t : ℚ} (hn₀ : 2 ≤ n₀ := by norm_num)
    (hn₀e : Even n₀ := by decide) (hr : c.r = r := by rfl) (hs : c.sAt 0 = s := by rfl)
    (ht : c.tAt 0 = t := by rfl) (hβ : 0 < β := by norm_num) (hr₁ : 1 < r := by norm_num) :
    Record.Seeded spectrum n₀ r β n₀ s t := by
  refine ⟨⟨hβ, hr₁⟩, fun ⟨hsp, hth⟩ L hL hLe ↦ ?_⟩
  subst hr hs ht
  exact gapGT_of_chainCertificate spectrum (by omega) hn₀e (by exact_mod_cast hβ) hZ c hc hsp
    hth hL hLe

/-- **Seeded records of reflection-positive rings.** `Record.Seeded.of_chain` for the spectrum
of a reflection-positive ring `R`. -/
theorem RPRing.seeded {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] [Nonempty d]
    (R : RPRing d ι) {n₀ : ℕ} (c : ChainCertificate) (hc : c.Valid) {r β s t : ℚ}
    (hn₀ : 2 ≤ n₀ := by norm_num) (hn₀e : Even n₀ := by decide) (hr : c.r = r := by rfl)
    (hs : c.sAt 0 = s := by rfl) (ht : c.tAt 0 = t := by rfl) (hβ : 0 < β := by norm_num)
    (hr₁ : 1 < r := by norm_num) :
    Record.Seeded R.spectrum n₀ r β n₀ s t :=
  Record.Seeded.of_chain R.spectrum (fun β hβ ↦ (R.hasSpatialMomentBounds β hβ).mono hn₀) c hc
    hn₀ hn₀e hr hs ht hβ hr₁

/-- Seed inequalities stated with the purities `RPRing.spatialPurity`, `RPRing.thermalPurity` of
a reflection-positive ring are its `SeedBounds`. -/
theorem RPRing.seedBounds {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)
    {n : ℕ} {β s t : ℚ} (hs : 1 - (s : ℝ) ≤ R.spatialPurity n β)
    (ht : 1 - (t : ℝ) ≤ R.thermalPurity (2 * n) β) : SeedBounds R.spectrum n β s t :=
  ⟨by simpa only [CertifiedGaps.spatialPurity, CertifiedGaps.partitionFunction,
      RPRing.spatialPurity, RPRing.partitionFunction_eq_sum] using hs,
    by simpa only [CertifiedGaps.thermalPurity, CertifiedGaps.partitionFunction,
      RPRing.thermalPurity, RPRing.partitionFunction_eq_sum] using ht⟩

end Leaderboard
