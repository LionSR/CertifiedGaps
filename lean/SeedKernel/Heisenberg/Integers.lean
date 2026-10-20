/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryInitTable
import CertifiedGaps.Heisenberg.OpenChain
import CertifiedGaps.Heisenberg.ThermalTable60
import SeedKernel.Heisenberg.BoundaryTraces
import SeedKernel.Heisenberg.OpenTrial
import SeedKernel.Heisenberg.Reseed120
import SeedKernel.Heisenberg.Seed36

/-!
# The spin-one gaps from the papers' printed integers

This file collects the headline theorems of both OpenAI papers in the form the formalization now
reaches: every statement below holds without `sorry`, and its hypotheses are exactly the printed
integers that the Lean kernel cannot recompute in reasonable memory, plus (for the boundary-field
gap) the spectral witness at `b₀ = 49/4`. Every trial state is kernel-checked.

* **Periodic paper, seed at `(36, 49/4)`.** `HornerTotals` says the floored integer recurrence of
  `ThermalTable` reproduces the eighteen printed totals `W_{n,g}` (`n = 4, …, 12`, `g = 0, 1`).
  `gap_of_hornerTotals`: every even ring `L ≥ 36` has gap `Δ_L > log(1.05882)/(49/4)`, and
  `gap_log_twenty_of_hornerTotals` is OpenAI's eq. (1.4).
* **Periodic paper, seed at `(60, 105/4)`.** `HornerTotals60` says the same for the nine totals of
  the `b = 21/2` table. `gap_sixty_of_hornerTotals60`: `Δ_L > log(1.16321)/(105/4)` for even
  `L ≥ 60`.
* **Boundary-field paper.** `BoundaryTotalsLarge` says the boundary recurrence of
  `BoundaryThermal` reproduces the six printed integers `k_{n,g}` with `n = 5, 6` (the nine with
  `n ≤ 4` are kernel-checked). `BoundaryWitness` is the spectral part of the initialization at
  `b₀`: a witness of the boundary moment inequalities with `λ₀ ≥ 9/10`, `ρ₀ ≤ .874` and
  `c₀ P_*(λ₀)² ≤ .003745`. From these and `HornerTotals`, `boundaryFieldInputs_of_integers`
  proves the four inputs of `BoundaryFieldInputs (3/5)` that `OpenChain.lean` leaves as `sorry`, so
  `spin_one_boundary_field_gap_of_integers` (odd `n ≥ 1921`) and
  `spin_one_opposite_field_gap_of_integers` (even `n ≥ 1920`) hold with gap
  `log(10⁶/2241)/392 ≈ 0.01556`.

`OpenChain.lean` cannot import this file (it sits below the kernel library), so its `sorry`
seeds stay in place; the theorems here are the `sorry`-free versions.
-/

open Real

namespace CertifiedGaps.SpinOneHeisenberg

/-! ### The periodic gaps -/

/-- **The periodic gap from the eighteen totals of the `b = 49/4` table**: every even ring of
`L ≥ 36` sites has a unique ground state and `Δ_L > log(1.05882)/(49/4) ≈ 0.00467`. -/
theorem gap_of_hornerTotals (h : HornerTotals) {L : ℕ} (hL : 36 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log (105882 / 100000) / (49 / 4)) :=
  gap_of_table (traceTable_of_hornerTotals h) hL hLe

/-- **OpenAI's eq. (1.4) from the eighteen totals**: `Δ_L > (log 20)/784` for every even
`L ≥ 2304`. -/
theorem gap_log_twenty_of_hornerTotals (h : HornerTotals) {L : ℕ} (hL : 2304 ≤ L)
    (hLe : Even L) : GapGT (spectrum L) (log 20 / 784) :=
  gap_log_twenty_of_table (traceTable_of_hornerTotals h) hL hLe

/-- **The periodic gap from the nine totals of the `b = 21/2` table**: every even ring of
`L ≥ 60` sites has a unique ground state and `Δ_L > log(1.16321)/(105/4) ≈ 0.00576`. -/
theorem gap_sixty_of_hornerTotals60 (h : HornerTotals60) {L : ℕ} (hL : 60 ≤ L)
    (hLe : Even L) : GapGT (spectrum L) (log (116321 / 100000) / (105 / 4)) :=
  gap_sixty_of_hornerTotals60_mps h trialA trial_check_120.2 hL hLe

/-! ### The boundary-field gap -/

namespace BoundaryInit

/-- **The spectral part of the initialization at `b₀ = 49/4`**: a witness `(λ₀, c₀, ρ₀, s₀)` of
the shifted boundary moment inequalities with `λ₀ ≥ 9/10`, `ρ₀ ≤ .874` and
`c₀ P_*(λ₀)² ≤ .003745` (the `witness` field of `ThermalCore`). -/
def BoundaryWitness : Prop :=
  ∃ μ c ρ s : ℝ, OpenMomentWitness (shiftedRingZ (49 / 4)) (shiftedOpenZ (49 / 4)) μ c ρ s ∧
    9 / 10 ≤ μ ∧ ρ ≤ 874 / 1000 ∧ c * (μ ^ 3 - 7 / 50 * μ ^ 2 - 9 / 50 * μ) ^ 2 ≤ 3745 / 10 ^ 6

/-- `ThermalCore` from the witness and the six boundary integers (`Y₆(b₀) < .027`). -/
theorem thermalCore_of_integers (hB : BoundaryThermal.BoundaryTotalsLarge)
    (hw : BoundaryWitness) : ThermalCore :=
  ⟨hw, (BoundaryThermal.sixth_of_large hB).le⟩

/-- A kernel-checked trial bound `E + am < α` in the form `TrialBelow m`. -/
private theorem trialBelow_of_lt {m : ℕ}
    (h : ∃ i, oppositeFieldSpectrum (3 / 5) m i < -(700741 / 500000 : ℝ) * m + 2021 / 5000) :
    TrialBelow m := by
  obtain ⟨i, hi⟩ := h
  exact ⟨i, by linarith⟩

/-- **H2 in the kernel**: the four open-chain trial bounds from the matrix-product trial vectors
of `SeedKernel.OpenTrial`. -/
theorem trialData : TrialData where
  trial34 := trialBelow_of_lt (by exact_mod_cast SeedKernel.OpenTrial.oppositeField_trial_34)
  trial78 := trialBelow_of_lt (by exact_mod_cast SeedKernel.OpenTrial.oppositeField_trial_78)
  trial142 := trialBelow_of_lt (by exact_mod_cast SeedKernel.OpenTrial.oppositeField_trial_142)
  trial240 := trialBelow_of_lt (by exact_mod_cast SeedKernel.OpenTrial.oppositeField_trial_240)

end BoundaryInit

open BoundaryInit in
/-- **The four inputs of the boundary-field gap from the integers**: the ring inputs from the
`b = 49/4` totals and the kernel-checked periodic trial states, the open-chain inputs from the
same totals, the six boundary integers, the boundary witness and the kernel-checked open-chain
trial states. -/
theorem boundaryFieldInputs_of_integers (hW : HornerTotals)
    (hB : BoundaryThermal.BoundaryTotalsLarge) (hw : BoundaryWitness) :
    BoundaryFieldInputs (3 / 5) :=
  have htab := traceTable_of_hornerTotals hW
  have hcore := thermalCore_of_integers hB hw
  ⟨seed_ring_spatial_of_table htab, seed_ring_thermal_of_table htab,
    boundaryDefect_le_of_table hcore trialData htab trialA trial_check_72.2
      exists_spectrum_lt_onetwenty,
    oppositeFieldThermalPurity_ge_of_table hcore trialData htab⟩

/-- **The boundary-field gap from the integers.** For every odd `n ≥ 1921` the open spin-one
Heisenberg chain with equal end fields `h = 3/5` has a unique ground state and gap
`Δ_n > log(10⁶/2241)/392 ≈ 0.01556`, given the printed integers of both papers that the kernel
does not recompute and the boundary witness at `b₀`. -/
theorem spin_one_boundary_field_gap_of_integers (hW : HornerTotals)
    (hB : BoundaryThermal.BoundaryTotalsLarge) (hw : BoundaryInit.BoundaryWitness) {n : ℕ}
    (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (equalFieldSpectrum (3 / 5) n) (log (1000000 / 2241) / 392) :=
  gapGT_of_boundaryFieldInputs_sharp (boundaryFieldInputs_of_integers hW hB hw) hn hno

/-- **The opposite-field gap from the integers**: `Δ_n > log(10⁶/2241)/392` for every even
`n ≥ 1920`, under the same hypotheses. -/
theorem spin_one_opposite_field_gap_of_integers (hW : HornerTotals)
    (hB : BoundaryThermal.BoundaryTotalsLarge) (hw : BoundaryInit.BoundaryWitness) {n : ℕ}
    (hn : 1920 ≤ n) (hne : Even n) :
    GapGT (oppositeFieldSpectrum (3 / 5) n) (log (1000000 / 2241) / 392) :=
  oppositeFieldGapGT_of_boundaryFieldInputs_sharp (boundaryFieldInputs_of_integers hW hB hw) hn
    hne

end CertifiedGaps.SpinOneHeisenberg
