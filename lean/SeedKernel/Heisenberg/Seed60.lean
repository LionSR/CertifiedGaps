/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Main
import CertifiedGaps.Heisenberg.SectorTransferC
import CertifiedGaps.Heisenberg.Seed60
import CertifiedGaps.Heisenberg.TrialBound
import SeedKernel.Heisenberg.Trial

/-!
# The spin-one gap for `L ≥ 60` from the source-2 trace table alone

OpenAI's Proposition 5.4 (`seed60_of_inputs`) needs three inputs besides its trace table at
`b = 21/2`: the `I/O/N` sector representation of the untwisted, `π`-twisted and `2π/3`-twisted
traces (`exists_ion_approx`), a Boltzmann form of the `C`-twisted trace (it is the Gibbs trace of
a Hermitian Hamiltonian), and an eigenvalue of the length-120 ring below `-120 κ`, which the
kernel-checked trial state gives (`trial_check_120`, `exists_spectrum_lt_of_mps`). So only the
table remains (`TraceTable60`, OpenAI's Lemma "finite thermal inputs" at `b = 21/2`):

* `seed60_of_table`: `P_s(60, 105/4), P_τ(120, 105/4) ≥ 7/8`;
* `gap_sixty_of_table`: OpenAI's eq. (1.3), `Δ_L > (4/105) log(80/79)` for every even `L ≥ 60`;
* `gap_sixty_chain_of_table`: the sharper `Δ_L > log(1.16321)/(105/4) ≈ 0.00576`.
-/

open Real

namespace CertifiedGaps.SpinOneHeisenberg

/-- The trial state gives an eigenvalue of the length-120 ring below `-120 κ`. -/
theorem exists_spectrum_lt_onetwenty : ∃ m, spectrum 120 m < -(700741 / 500000 * 120) := by
  obtain ⟨m, hm⟩ := exists_spectrum_lt_of_mps trialA (N := 120) (by norm_num)
    trial_check_120.2
  exact ⟨m, by push_cast at hm; exact hm⟩

/-- **OpenAI's Proposition 5.4 from the trace table.** -/
theorem seed60_of_table (htab : TraceTable60 cTwistedPartitionFunction) :
    7 / 8 ≤ spatialPurity 60 (105 / 4) ∧ 7 / 8 ≤ thermalPurity 120 (105 / 4) :=
  seed60_of_inputs cTwistedPartitionFunction (fun _ hβ ns hns ↦ exists_ion_approx hβ ns hns)
    (boltzmannRep_gibbsTrace (cTwistedHamiltonian_isHermitian 30)) htab
    exists_spectrum_lt_onetwenty

/-- **OpenAI's eq. (1.3) from the trace table**: every even `L ≥ 60` has a unique ground state
and `Δ_L > (4/105) log(80/79)`. -/
theorem gap_sixty_of_table (htab : TraceTable60 cTwistedPartitionFunction) {L : ℕ}
    (hL : 60 ≤ L) (hLe : Even L) : GapGT (spectrum L) (4 / 105 * log (80 / 79)) :=
  gapGT_sixty_of_seeds (seed60_of_table htab).1 (seed60_of_table htab).2 hL hLe

/-- **The sharper gap for `L ≥ 60` from the trace table**: `Δ_L > log(1.16321)/(105/4)`. -/
theorem gap_sixty_chain_of_table (htab : TraceTable60 cTwistedPartitionFunction) {L : ℕ}
    (hL : 60 ≤ L) (hLe : Even L) : GapGT (spectrum L) (log (116321 / 100000) / (105 / 4)) :=
  gapGT_sixty_chain_of_seeds (seed60_of_table htab).1 (seed60_of_table htab).2 hL hLe

end CertifiedGaps.SpinOneHeisenberg
