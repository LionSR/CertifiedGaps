/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.PurityCriterion

/-!
# Leaderboard records

A leaderboard record states a uniform spectral gap of a family of periodic rings. The family is
given by its spectra `spectrum L : (Fin L → d) → ℝ`, one eigenvalue per basis configuration of a
ring of `L` unit cells with cell space `d` (as `CertifiedGaps.RPRing.spectrum`). A record with
length `L₀` and rationals `r`, `β` states, for every even `L ≥ L₀`, a unique ground state and

  `Δ_L > log r / β`

(`CertifiedGaps.GapGT`), together with `0 < β` and `1 < r`, so that the bound is positive.

Every chain certificate of the repository gives a bound of this form (`log r / β₀`, with `r` the
rate of the certificate and `β₀` the seed inverse temperature), and any positive gap bound can be
weakened to it, so the rational pair `(r, β)` loses nothing and keeps the statement free of real
literals that the audit would have to evaluate.

* `Record.Unconditional spectrum L₀ r β`: the gap holds outright.
* `Record.Conditional H spectrum L₀ r β`: the gap holds under the hypothesis `H`, in practice the
  two seed inequalities that the computer certificates of `models/` establish outside Lean.

Only even ring lengths are covered, since the spacetime bootstrap proves gaps on even rings. Odd
rings or open chains would get their own record heads (with their own length predicate), added
here and to the audit; the existing records would not change.
-/

namespace Leaderboard

open CertifiedGaps

/-- The bound `log r / β` is positive: `0 < β` and `1 < r`. -/
def PositiveBound (r β : ℚ) : Prop := 0 < β ∧ 1 < r

/-- For every even `L ≥ L₀`, the spectrum `spectrum L` has a nondegenerate lowest level and every
other level lies more than `log r / β` above it. -/
def UniformGap {d : Type*} (spectrum : (L : ℕ) → (Fin L → d) → ℝ) (L₀ : ℕ) (r β : ℚ) : Prop :=
  ∀ L, L₀ ≤ L → Even L → GapGT (spectrum L) (Real.log r / β)

namespace Record

/-- **An unconditional record.** `Δ_L > log r / β > 0` for every even ring length `L ≥ L₀`. -/
def Unconditional {d : Type*} (spectrum : (L : ℕ) → (Fin L → d) → ℝ) (L₀ : ℕ) (r β : ℚ) :
    Prop :=
  PositiveBound r β ∧ UniformGap spectrum L₀ r β

/-- **A conditional record.** `log r / β > 0`, and if `H` holds then `Δ_L > log r / β` for every
even ring length `L ≥ L₀`. The positivity of the bound does not depend on `H`. -/
def Conditional (H : Prop) {d : Type*} (spectrum : (L : ℕ) → (Fin L → d) → ℝ) (L₀ : ℕ)
    (r β : ℚ) : Prop :=
  PositiveBound r β ∧ (H → UniformGap spectrum L₀ r β)

variable {H : Prop} {d : Type*} {spectrum : (L : ℕ) → (Fin L → d) → ℝ} {L₀ : ℕ} {r β : ℚ}

/-- An unconditional record is conditional on any hypothesis. -/
theorem Unconditional.conditional (h : Unconditional spectrum L₀ r β) :
    Conditional H spectrum L₀ r β :=
  ⟨h.1, fun _ ↦ h.2⟩

/-- A conditional record whose hypothesis is proved is unconditional. -/
theorem Conditional.unconditional (h : Conditional H spectrum L₀ r β) (hH : H) :
    Unconditional spectrum L₀ r β :=
  ⟨h.1, h.2 hH⟩

/-- A conditional record transfers to a hypothesis that implies the old one. -/
theorem Conditional.mono {H' : Prop} (h : Conditional H spectrum L₀ r β) (hH : H' → H) :
    Conditional H' spectrum L₀ r β :=
  ⟨h.1, fun h' ↦ h.2 (hH h')⟩

/-- Records transfer between spectra that agree on the even lengths `L ≥ L₀`, for example a
physical Hamiltonian and its unitarily rotated reflection-positive frame. -/
theorem Conditional.congr {spectrum' : (L : ℕ) → (Fin L → d) → ℝ}
    (h : Conditional H spectrum L₀ r β) (he : ∀ L, L₀ ≤ L → Even L → spectrum L = spectrum' L) :
    Conditional H spectrum' L₀ r β :=
  ⟨h.1, fun hH L hL hLe ↦ he L hL hLe ▸ h.2 hH L hL hLe⟩

end Record

end Leaderboard
