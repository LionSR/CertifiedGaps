/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import CertifiedGaps.PurityCriterion

/-!
# Leaderboard records

A record states, for every even ring length `L ≥ L₀`, a unique ground state and
`Δ_L > log r / β` (`CertifiedGaps.GapGT`), with rational `r > 1` and `β > 0`. The ring is given by
its spectra `spectrum L : (Fin L → d) → ℝ`, as `CertifiedGaps.RPRing.spectrum`; `L` counts unit
cells.

* `Record.Unconditional spectrum L₀ r β`: the gap holds outright.
* `Record.Seeded spectrum L₀ r β n₀ s t`: the gap holds if the two seed inequalities
  `P_s(n₀, β) ≥ 1 - s` and `P_τ(2n₀, β) ≥ 1 - t` hold for the same spectrum. The definition fixes
  the hypothesis, so the rationals `n₀`, `s`, `t` determine it.
* `Record.SeededOn family a b L₀ r β n₀ s t`: `Record.Seeded (family x) L₀ r β n₀ s t` for every
  real parameter `x` in `[a, b]`. A computer certificate that bounds the two defects on the whole
  interval, for example by a perturbation bound, proves the seeds at every `x`.

Only even rings are covered. Odd rings or open chains would need new record forms.
-/

namespace Leaderboard

open CertifiedGaps

/-- The bound `log r / β` is positive: `0 < β` and `1 < r`. -/
def PositiveBound (r β : ℚ) : Prop := 0 < β ∧ 1 < r

/-- For every even `L ≥ L₀`, `spectrum L` has a nondegenerate lowest level and every other level
lies more than `log r / β` above it. -/
def UniformGap {d : Type*} (spectrum : (L : ℕ) → (Fin L → d) → ℝ) (L₀ : ℕ) (r β : ℚ) : Prop :=
  ∀ L, L₀ ≤ L → Even L → GapGT (spectrum L) (Real.log r / β)

/-- The seed inequalities `P_s(n₀, β) ≥ 1 - s` and `P_τ(2n₀, β) ≥ 1 - t` of `spectrum`. -/
def SeedBounds {d : Type*} [Fintype d] [DecidableEq d] (spectrum : (L : ℕ) → (Fin L → d) → ℝ)
    (n₀ : ℕ) (β s t : ℚ) : Prop :=
  1 - (s : ℝ) ≤ spatialPurity spectrum n₀ β ∧ 1 - (t : ℝ) ≤ thermalPurity spectrum (2 * n₀) β

namespace Record

/-- **An unconditional record.** `Δ_L > log r / β > 0` for every even `L ≥ L₀`. -/
def Unconditional {d : Type*} (spectrum : (L : ℕ) → (Fin L → d) → ℝ) (L₀ : ℕ) (r β : ℚ) :
    Prop :=
  PositiveBound r β ∧ UniformGap spectrum L₀ r β

/-- **A seeded record.** `log r / β > 0`, and if `P_s(n₀, β) ≥ 1 - s` and `P_τ(2n₀, β) ≥ 1 - t`,
then `Δ_L > log r / β` for every even `L ≥ L₀`. -/
def Seeded {d : Type*} [Fintype d] [DecidableEq d] (spectrum : (L : ℕ) → (Fin L → d) → ℝ)
    (L₀ : ℕ) (r β : ℚ) (n₀ : ℕ) (s t : ℚ) : Prop :=
  PositiveBound r β ∧ (SeedBounds spectrum n₀ β s t → UniformGap spectrum L₀ r β)

/-- **A seeded record on a parameter interval.** For every real `x` with `a ≤ x ≤ b`, the seeded
record `Seeded (family x) L₀ r β n₀ s t`, whose seeds are those of `family x`. -/
def SeededOn {d : Type*} [Fintype d] [DecidableEq d]
    (family : ℝ → (L : ℕ) → (Fin L → d) → ℝ) (a b : ℚ) (L₀ : ℕ) (r β : ℚ) (n₀ : ℕ) (s t : ℚ) :
    Prop :=
  ∀ x : ℝ, (a : ℝ) ≤ x → x ≤ b → Seeded (family x) L₀ r β n₀ s t

variable {d : Type*} [Fintype d] [DecidableEq d] {spectrum spectrum' : (L : ℕ) → (Fin L → d) → ℝ}
  {L₀ n₀ : ℕ} {r β s t : ℚ}

/-- A seeded record whose seeds are proved is unconditional. -/
theorem Seeded.unconditional (h : Seeded spectrum L₀ r β n₀ s t)
    (hs : SeedBounds spectrum n₀ β s t) : Unconditional spectrum L₀ r β :=
  ⟨h.1, h.2 hs⟩

/-- Seeded records transfer between spectra that agree at the even lengths `≥ min L₀ n₀`, for
example a physical Hamiltonian and its rotated reflection-positive frame. The seed lengths `n₀`
and `2n₀` are among them when `n₀` is even. -/
theorem Seeded.congr (h : Seeded spectrum L₀ r β n₀ s t) (hn₀ : Even n₀)
    (he : ∀ L, min L₀ n₀ ≤ L → Even L → spectrum L = spectrum' L) :
    Seeded spectrum' L₀ r β n₀ s t := by
  have hZ : ∀ L, n₀ ≤ L → Even L → ∀ β : ℝ,
      partitionFunction spectrum L β = partitionFunction spectrum' L β := fun L hL hLe β ↦ by
    simp only [partitionFunction, he L (le_trans (min_le_right _ _) hL) hLe]
  have h2 : n₀ ≤ 2 * n₀ := by omega
  have he2 : Even (2 * n₀) := even_two_mul n₀
  refine ⟨h.1, fun ⟨hs, ht⟩ L hL hLe ↦ ?_⟩
  rw [← he L (le_trans (min_le_left _ _) hL) hLe]
  refine h.2 ⟨?_, ?_⟩ L hL hLe
  · simpa only [spatialPurity, hZ _ le_rfl hn₀, hZ _ h2 he2] using hs
  · simpa only [thermalPurity, hZ _ h2 he2] using ht

end Record

end Leaderboard
