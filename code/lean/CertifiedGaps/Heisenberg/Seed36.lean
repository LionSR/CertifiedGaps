/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Certificate
import CertifiedGaps.Heisenberg.Main
import CertifiedGaps.Heisenberg.SectorTransfer
import CertifiedGaps.Heisenberg.SzBridge
import CertifiedGaps.Heisenberg.TraceTable
import CertifiedGaps.Heisenberg.TrialBound

/-!
# The seed at `(36, 49/4)` of the spin-one chain

This file proves the two starting defects of `Main.lean`,

  `1 - P_s(36, 49/4) ≤ 0.047916`  and  `1 - P_τ(72, 49/4) ≤ 0.181521`,

as an instance of the generic `TableCertificate.purities`. The rational data are `seed36`,
checked by `decide`: the `D₂` scheme (`SectorScheme.d2`), the trace table of `TraceTable.lean` at
the lengths `n = 4, …, 12` within `5 × 10⁻⁵`, two degree-12 filters (the invariant one uses odd
lengths, whose moments carry the sign of the dominant eigenvalue), the caps `1.00136` and
`0.87125`, and the dominance threshold `ℓ = 0.999`. The table is certified to `3 × 10⁻⁵`; the
filters accept `5 × 10⁻⁵`.

The inputs, each a statement about the Lean Heisenberg ring:

* `SectorInput`: the transfer spectrum in the `D₂` sectors (proved: `sectorInput_seed`);
* `TraceTable`: the shifted traces `e^{-βκn} Z_n(β)` and `e^{-βκn} Z_n(β, P)` at `β = 49/4` lie
  within `3 × 10⁻⁵` of the table (Appendix A of the certificate note; open);
* a trial state with `500000 Tr (T⁷⁰ O) < -700741 Tr T⁷²` (eq. (c-trial)), which gives
  `e^{-βκ·72} Z₇₂(β) > 1` (`one_lt_shiftedZ_of_mps`).
-/

open Finset Real Matrix

namespace CertifiedGaps.SpinOneHeisenberg

/-- The energy shift `κ = 1.401482` per bond. -/
def kappa : ℚ := 700741 / 500000

/-- The shifted partition function `e^{-βκL} Z_L(β)`. -/
noncomputable def shiftedZ (L : ℕ) (β : ℝ) : ℝ := shiftedTrace β kappa L (partitionFunction L β)

/-- The shifted twisted partition function `e^{-βκL} Z_L(β, P)`. -/
noncomputable def shiftedZP (L : ℕ) (β : ℝ) : ℝ :=
  shiftedTrace β kappa L (twistedPartitionFunction L β)

/-- **Input 1: the sector representation** of the untwisted and twisted ring traces at inverse
temperature `β`, in the `D₂` sectors. -/
def SectorInput (β : ℝ) (ns : Finset ℕ) : Prop :=
  D2Approx (fun n ↦ partitionFunction n β) (fun n ↦ twistedPartitionFunction n β) ns

/-- **Input 2: the trace table.** The shifted traces at `β = 49/4` lie within `3 × 10⁻⁵` of the
table centers. -/
def TraceTable : Prop :=
  ∀ n ∈ tableLengths, |shiftedZ n (49 / 4) - traceCenter n| ≤ 3 / 100000 ∧
    |shiftedZP n (49 / 4) - twistedCenter n| ≤ 3 / 100000

/-- The rational data and parameters of the seed at `(36, 49/4)`. -/
def seed36 : TableCertificate where
  scheme := SectorScheme.d2
  table := TraceBounds.centered ![traceCenter, twistedCenter] tableLengths (1 / 20000)
  filter := ![some (2, [4 / 591, 51 / 3466, -649 / 2330, -267 / 1582, 1]),
    some (2, [27 / 2258, 0, -489 / 1238, 0, 1])]
  seed :=
    { q := 12, a := 3, mult := ![1, 3]
      cap := ![some (10013641597 / 10000000000), some (4356234043 / 5000000000)]
      mass := ![2131481 / 2000000, 556621 / 2000000]
      ell := 999 / 1000
      spatialDefect := 47916 / 10 ^ 6 }
  beta := 49 / 4
  kappa := kappa
  thermalDefect := 181521 / 10 ^ 6

/-! The conditions of `seed36.Valid`, each checked by its own `decide`; one `decide` on the whole
conjunction exceeds the kernel's recursion depth. -/

theorem seed36_scheme_valid : seed36.scheme.Valid := by decide +kernel
theorem seed36_seed_valid : seed36.seed.Valid := by decide +kernel
theorem seed36_mult : ∀ s, (seed36.seed.mult s : ℤ) = seed36.scheme.coeff 0 s := by decide +kernel
theorem seed36_filters : ∀ s, seed36.FilterCheck s := by decide +kernel
theorem seed36_masses :
    ∀ s, OptLE (seed36.table.hiS seed36.scheme s seed36.seed.q) (seed36.seed.mass s) := by
  decide +kernel
theorem seed36_thermalCheck : (1 - seed36.thermalDefect) * seed36.seed.thermalBound ^ 2 ≤ 1 := by
  decide +kernel

theorem seed36_valid : seed36.Valid :=
  ⟨seed36_scheme_valid, seed36_seed_valid, seed36_mult, seed36_filters, seed36_masses,
    by decide, by decide, by decide +kernel, by decide +kernel, seed36_thermalCheck⟩

/-- The ring lengths that the seed uses: the table and the lengths `36` and `72`. -/
def seedLengths : Finset ℕ := seed36.lengths

theorem partitionFunction_pos (L : ℕ) (β : ℝ) : 0 < partitionFunction L β := by
  rw [partitionFunction_eq_sum]
  exact sum_pos (fun _ _ ↦ exp_pos _) univ_nonempty

/-- The tabulated traces: untwisted and twisted. -/
noncomputable def seedTraces : Fin 2 → ℕ → ℝ → ℝ :=
  ![fun n β ↦ partitionFunction n β, fun n β ↦ twistedPartitionFunction n β]

/-- The seed purities in the certificate's own terms. -/
theorem seed36_purities (hsec : SectorInput (49 / 4) seedLengths) (htab : TraceTable)
    (hmps : 1 < shiftedZ 72 (49 / 4)) (hmps2 : 1 < shiftedZ 72 (49 / 2)) :
    1 - (seed36.seed.spatialDefect : ℝ) ≤ seedTraces 0 (2 * seed36.n₀) (seed36.beta : ℝ) /
        seedTraces 0 seed36.n₀ (seed36.beta : ℝ) ^ 2 ∧
      1 - (seed36.thermalDefect : ℝ) ≤ seedTraces 0 (2 * seed36.n₀) (2 * (seed36.beta : ℝ)) /
        seedTraces 0 (2 * seed36.n₀) (seed36.beta : ℝ) ^ 2 := by
  have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
  have hn : seed36.n₀ = 36 := rfl
  have h2b : (2 : ℝ) * (49 / 4) = 49 / 2 := by norm_num
  refine TableCertificate.purities seed36_valid ?_ ?_ ?_ ?_ ?_
  · rw [hb]
    convert hsec.toSectorApprox using 1
    all_goals first | rfl | (apply heq_of_eq; funext g n; fin_cases g <;> rfl)
  · rw [hb]
    refine TraceBounds.holdsStrict_centered (εt := 3 / 100000) (by norm_num) fun n hn g ↦ ?_
    have e : ((3 / 100000 : ℚ) : ℝ) = 3 / 100000 := by norm_num
    rw [e]
    fin_cases g
    · exact (htab n hn).1
    · exact (htab n hn).2
  · rw [hb]; exact partitionFunction_pos _ _
  · rw [hb, hn, show 2 * 36 = 72 from rfl]; exact hmps
  · rw [hb, hn, h2b, show 2 * 36 = 72 from rfl]; exact hmps2

/-- **The seed at `(36, 49/4)`.** The three inputs give the starting defects of `Main.lean`, by
the generic `TableCertificate.purities`. -/
theorem seed36_of_inputs (hsec : SectorInput (49 / 4) seedLengths) (htab : TraceTable)
    (hmps : 1 < shiftedZ 72 (49 / 4)) (hmps2 : 1 < shiftedZ 72 (49 / 2)) :
    1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4) ∧
      1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4) := by
  obtain ⟨h1, h2⟩ := seed36_purities hsec htab hmps hmps2
  have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
  have hn : seed36.n₀ = 36 := rfl
  rw [hb, hn] at h1 h2
  rw [show (2 : ℝ) * (49 / 4) = 49 / 2 by norm_num] at h2
  have e1 : ((seed36.seed.spatialDefect : ℚ) : ℝ) = 47916 / 10 ^ 6 := by
    rw [show seed36.seed.spatialDefect = 47916 / 10 ^ 6 from rfl]; push_cast; ring
  have e2 : ((seed36.thermalDefect : ℚ) : ℝ) = 181521 / 10 ^ 6 := by
    rw [show seed36.thermalDefect = 181521 / 10 ^ 6 from rfl]; push_cast; ring
  rw [e1] at h1
  rw [e2] at h2
  refine ⟨by simpa [spatialPurity, seedTraces] using h1, ?_⟩
  rw [thermalPurity, show (2 : ℝ) * (49 / 4) = 49 / 2 by norm_num]
  simpa [seedTraces] using h2

/-! ### The trial-state input -/

/-- **The trial-state input from an integer trace inequality.** If a matrix-product state with
integer tensors has `500000 Tr (T⁷⁰ O) < -700741 Tr T⁷²` (eq. (c-trial)), then
`e^{-βκ·72} Z₇₂(β) > 1` for every `β > 0` (eq. (c-z72)). -/
theorem one_lt_shiftedZ_of_mps {D : Type*} [Fintype D] [DecidableEq D] (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72))
    {β : ℝ} (hβ : 0 < β) : 1 < shiftedZ 72 β := by
  have h := exp_lt_partitionFunctionZ_of_mps_700741 A (N := 72) (by norm_num) hU hβ
  rw [← partitionFunction_eq_partitionFunctionZ (by norm_num)] at h
  have hk : (kappa : ℝ) = 700741 / 500000 := by norm_num [kappa]
  rw [shiftedZ]
  exact one_lt_shiftedTrace_of_exp_lt (by rwa [hk])

/-- **The seed at `(36, 49/4)` from the sector representation, the trace table and a trial
state.** -/
theorem seed36_of_mps (hsec : SectorInput (49 / 4) seedLengths) (htab : TraceTable)
    {D : Type*} [Fintype D] [DecidableEq D] (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72)) :
    1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4) ∧
      1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4) :=
  seed36_of_inputs hsec htab (one_lt_shiftedZ_of_mps A hU (by norm_num))
    (one_lt_shiftedZ_of_mps A hU (by norm_num))

/-! ### The sector input -/

/-- **Input 1 holds**: the single-site transfer matrix gives the sector representation at every
length of the seed (`exists_sector_approx`). -/
theorem sectorInput_seed : SectorInput (49 / 4) seedLengths := fun _ hδ ↦
  exists_sector_approx (by norm_num) seedLengths
    (fun _ hn ↦ TableCertificate.two_le_of_mem_lengths seed36_valid hn) hδ

/-- **The seed at `(36, 49/4)` from the trace table and a trial state.** The only remaining
inputs are the trace table at `n = 4, …, 12` and one integer inequality about the trial
state. -/
theorem seed36_of_table_mps (htab : TraceTable) {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72)) :
    1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4) ∧
      1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4) :=
  seed36_of_mps sectorInput_seed htab A hU

/-- **The gap for every even `L ≥ 36` from the trace table and a trial state.** Every even
`L ≥ 36` has a unique ground state and `Δ_L > log(1.05882)/(49/4) ≈ 0.00467`, assuming only the
trace table at `n = 4, …, 12` and one integer inequality about a matrix-product trial state. -/
theorem gap_thirtysix_of_table_mps (htab : TraceTable) {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72))
    {L : ℕ} (hL : 36 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log (105882 / 100000) / (49 / 4)) :=
  gapGT_thirtysix_of_seeds (seed36_of_table_mps htab A hU).1 (seed36_of_table_mps htab A hU).2
    hL hLe

end CertifiedGaps.SpinOneHeisenberg
