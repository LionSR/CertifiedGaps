/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.LiftedSeed
import CertifiedGaps.Heisenberg.Hamiltonian
import CertifiedGaps.Heisenberg.Gap

/-!
# The seed at `(60, 105/4)` of the spin-one chain

This file proves the two seeds of `Main.lean` at length sixty,

  `P_s(60, 105/4) ≥ 7/8`  and  `P_τ(120, 105/4) ≥ 7/8`

(Proposition 5.4 of OpenAI's *The periodic spin-one Haldane gap*), as an instance of the generic
`LiftedSeed.purities`. The rational data are `seed60`, checked by `decide`: the three-twist table
at `b₀ = 21/2` and `n = 6, 8, 10`, the filters `x⁶(-2 + 10x²)²` (sector `I`) and
`x⁶(-25 + 100x²)²` (sector `N`) with caps `U = 1.00245`, `V = 0.88643`, the masses
`m_I = 1.042655`, `m_O = 0.040041`, `m_N = 0.401625` at order `10`, the lifting `b = (5/2) b₀`,
`D = 1.45`, `E = 1.202`, `τ = 1.225`, `τ₀ = 1.264`, the thresholds `0.88 → 0.99 → 0.998` and
`R_q = 0.0684`.

The inputs, each a statement about the Lean Heisenberg ring and a third family `Zc` of traces
(the ring twisted by the `2π/3` rotation `C`):

* the sector representation of the traces `Z_n(β)`, `Z_n(β, P)`, `Zc n β` in the sectors
  `I, O, N` of `SectorScheme.ion`, at every `β > 0` and every finite set of lengths `≥ 2`;
* `Zc 30` is a Boltzmann sum (the Gibbs trace of a Hermitian twisted Hamiltonian);
* `TraceTable60`: the shifted traces `e^{-b₀κn} Z_n(b₀, g)` lie within `10⁻⁵` of the paper's
  centers (their Lemma "Finite thermal inputs");
* a trial state: one eigenvalue of `H₁₂₀` lies below `-120κ`, `κ = 700741/500000`.

The certificate accepts table errors up to `1.01 × 10⁻⁵`, so the hypothesis `≤ 10⁻⁵` holds
strictly.
-/

open Finset Real Matrix

namespace CertifiedGaps.SpinOneHeisenberg

/-! ### The data -/

/-- The table centers of the shifted traces at `b₀ = 21/2` in units `10⁻⁶`: rows `n = 6, 8, 10`,
columns the twists `1, P, C` (the lemma "Finite thermal inputs" of OpenAI's paper). -/
def table60Rows : List (List ℕ) :=
  [[8945303, 346734, 939987], [3741281, 567167, 984108], [2327547, 721090, 1002613]]

/-- The table center of the shifted trace with twist `g` at the length `n ∈ {6, 8, 10}`. -/
def center60 (g : Fin 3) (n : ℕ) : ℚ :=
  ((table60Rows.getD ((n - 6) / 2) []).getD g 0 : ℚ) / 10 ^ 6

/-- The rational data of the seed at `(60, 105/4)`. -/
def seed60 : LiftedSeed where
  table := TraceBounds.centered center60 [6, 8, 10] (101 / 10 ^ 7)
  filter := ![some (3, [-2, 0, 10]), none, some (3, [-25, 0, 100])]
  cap := ![some (100245 / 100000), none, some (88643 / 100000)]
  q := 10
  mass := ![1042655 / 10 ^ 6, 40041 / 10 ^ 6, 401625 / 10 ^ 6]
  r := 3
  j := 2
  l := 5
  D := 1450 / 1000
  E := 1202 / 1000
  tau := 1225 / 1000
  tau0 := 1264 / 1000
  t₀ := 88 / 100
  steps := [99 / 100, 998 / 1000]
  Rq := 684 / 10000
  spatialDefect := 1 / 8
  thermalDefect := 1 / 8
  beta := 21 / 2
  kappa := 700741 / 500000

/-! The conditions of `seed60.Valid`, each checked by its own `decide`. -/

theorem seed60_filters : ∀ s : Fin (seed60.low.scheme.m + 1), seed60.low.FilterCheck s := by
  decide +kernel
theorem seed60_masses :
    ∀ s : Fin 3, OptLE (seed60.table.hiS SectorScheme.ion s seed60.q) (seed60.mass s) := by
  decide +kernel
theorem seed60_shapeValid : seed60.ShapeValid := by decide +kernel
theorem seed60_liftValid : seed60.LiftValid := by decide +kernel
theorem seed60_spatialValid : seed60.SpatialValid := by decide +kernel
theorem seed60_thermalValid : seed60.ThermalValid := by decide +kernel

theorem seed60_valid : seed60.Valid :=
  ⟨seed60_filters, seed60_masses, seed60_shapeValid, seed60_liftValid, seed60_spatialValid,
    seed60_thermalValid⟩

/-! ### The inputs -/

/-- The three tabulated traces: untwisted, twisted by `P`, and the family `Zc` (twisted by
`C`). -/
noncomputable def traces60 (Zc : ℕ → ℝ → ℝ) : Fin 3 → ℕ → ℝ → ℝ :=
  ![fun n β ↦ partitionFunction n β, fun n β ↦ twistedPartitionFunction n β, Zc]

/-- **The trace table at `b₀ = 21/2`.** The shifted traces `e^{-b₀κn} Z_n(b₀, g)`, with
`κ = 700741/500000`, lie within `10⁻⁵` of `center60 g n` for `n = 6, 8, 10` and the twists
`g = 1, P, C`. -/
def TraceTable60 (Zc : ℕ → ℝ → ℝ) : Prop :=
  ∀ n ∈ [6, 8, 10], ∀ g : Fin 3,
    |shiftedTrace (21 / 2) (700741 / 500000) n (traces60 Zc g n (21 / 2)) - center60 g n| ≤
      1 / 10 ^ 5

/-- The twisted partition function is a Boltzmann sum. -/
theorem boltzmannRep_twistedPartitionFunction (L : ℕ) :
    BoltzmannRep (twistedPartitionFunction L) :=
  boltzmannRep_gibbsTrace (H := twistedRingHamiltonian SpinOne.exchange twistedExchange L)
    (twistedRingHamiltonian_isHermitian SpinOne.exchange_isHermitian twistedExchange_isHermitian L)

/-- The partition function is a Boltzmann sum. -/
theorem boltzmannRep_partitionFunction (L : ℕ) : BoltzmannRep (partitionFunction L) :=
  boltzmannRep_gibbsTrace (hamiltonian_isHermitian L)

/-- The untwisted entry of `traces60` is the partition function. -/
theorem traces60_zero (Zc : ℕ → ℝ → ℝ) : traces60 Zc 0 = fun n β ↦ partitionFunction n β := rfl

/-- The lifted length `k = q r = 30`. -/
theorem seed60_k : seed60.k = 30 := rfl

/-- The table temperature `b₀ = 21/2`. -/
theorem seed60_beta : ((seed60.beta : ℚ) : ℝ) = 21 / 2 := by norm_num [seed60]

/-- The lifted temperature `b = (5/2) b₀ = 105/4`. -/
theorem seed60_liftBeta : ((seed60.liftBeta : ℚ) : ℝ) = 105 / 4 := by
  norm_num [seed60, LiftedSeed.liftBeta]

/-- The sector representation at one `β`, for the traces `traces60 Zc`. -/
theorem sectorApprox_traces60 {Zc : ℕ → ℝ → ℝ} {β : ℝ} {ns : Finset ℕ}
    (h : SectorApprox SectorScheme.ion
      ![fun n ↦ partitionFunction n β, fun n ↦ twistedPartitionFunction n β, fun n ↦ Zc n β] ns) :
    SectorApprox SectorScheme.ion (fun g n ↦ traces60 Zc g n β) ns := by
  convert h using 1
  funext g n
  fin_cases g <;> rfl

/-- Every length of the representation at `b₀` is at least `2`. -/
theorem seed60_lengths₀ : ∀ n ∈ seed60.lengths₀, 2 ≤ n := by
  intro n hn
  simp only [LiftedSeed.lengths₀, seed60, TraceBounds.centered, LiftedSeed.k, mem_union,
    List.mem_toFinset, List.mem_cons, List.not_mem_nil, mem_insert, mem_singleton, or_false] at hn
  omega

/-- Every length of the representation at `b` is at least `2`. -/
theorem seed60_lengths₁ : ∀ n ∈ seed60.lengths₁, 2 ≤ n := by
  intro n hn
  simp only [LiftedSeed.lengths₁, seed60, LiftedSeed.k, mem_insert, mem_singleton] at hn
  omega

/-- The energy shift `κ = 700741/500000`. -/
theorem seed60_kappa : ((seed60.kappa : ℚ) : ℝ) = 700741 / 500000 := by norm_num [seed60]

/-- The table hypothesis in the form of the certificate. -/
theorem seed60_holds {Zc : ℕ → ℝ → ℝ} (htab : TraceTable60 Zc) :
    seed60.table.HoldsStrict
      (fun g n ↦ shiftedTrace seed60.beta seed60.kappa n (traces60 Zc g n seed60.beta)) := by
  rw [seed60_beta, seed60_kappa]
  refine TraceBounds.holdsStrict_centered (εt := 1 / 10 ^ 5) (by norm_num) fun n hn g ↦ ?_
  have e : (((1 / 10 ^ 5 : ℚ)) : ℝ) = 1 / 10 ^ 5 := by norm_num
  rw [e]
  exact htab n hn g

/-- The trial-state hypothesis in the form of the certificate. -/
theorem seed60_trial (htrial : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    ∃ m, spectrum (4 * seed60.k) m < -((seed60.kappa : ℝ) * ((4 * seed60.k : ℕ) : ℝ)) := by
  rw [show 4 * seed60.k = 120 from rfl, seed60_kappa, Nat.cast_ofNat]
  exact htrial

/-- The partition function at length `60` is positive. -/
theorem partitionFunction_pos_sixty : 0 < partitionFunction 60 (105 / 4) := by
  rw [partitionFunction_eq_gibbsTrace]
  exact gibbsTrace_pos (hamiltonian_isHermitian 60) _

/-- **The seed at `(60, 105/4)` from its inputs** (OpenAI's Proposition 5.4). The inputs are the
sector representation of the three twisted traces, a Boltzmann representation of the
`C`-twisted trace at length `30`, the trace table at `b₀ = 21/2`, and one eigenvalue of `H₁₂₀`
below `-120κ`. -/
theorem seed60_of_inputs (Zc : ℕ → ℝ → ℝ)
    (hion : ∀ β : ℝ, 0 < β → ∀ ns : Finset ℕ, (∀ n ∈ ns, 2 ≤ n) →
      SectorApprox SectorScheme.ion
        ![fun n ↦ partitionFunction n β, fun n ↦ twistedPartitionFunction n β, fun n ↦ Zc n β] ns)
    (hZc : BoltzmannRep (Zc 30)) (htab : TraceTable60 Zc)
    (htrial : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    7 / 8 ≤ spatialPurity 60 (105 / 4) ∧ 7 / 8 ≤ thermalPurity 120 (105 / 4) := by
  have hrep : ∀ g, BoltzmannRep (traces60 Zc g seed60.k) := fun g ↦ by
    rw [seed60_k]
    fin_cases g
    · exact boltzmannRep_partitionFunction 30
    · exact boltzmannRep_twistedPartitionFunction 30
    · exact hZc
  have hpos : 0 < traces60 Zc 0 (2 * seed60.k) seed60.liftBeta := by
    rw [show 2 * seed60.k = 60 from rfl, seed60_liftBeta]
    exact partitionFunction_pos_sixty
  obtain ⟨h1, h2⟩ := LiftedSeed.purities seed60_valid (Z := traces60 Zc)
    (sectorApprox_traces60 (hion _ (by rw [seed60_beta]; norm_num) _ seed60_lengths₀))
    (sectorApprox_traces60 (hion _ (by rw [seed60_liftBeta]; norm_num) _ seed60_lengths₁))
    (seed60_holds htab) hrep (spectrum (4 * seed60.k))
    (fun β ↦ partitionFunction_eq_sum (4 * seed60.k) β) (seed60_trial htrial) hpos
  have e1 : ((seed60.spatialDefect : ℚ) : ℝ) = 1 / 8 := by norm_num [seed60]
  have e2 : ((seed60.thermalDefect : ℚ) : ℝ) = 1 / 8 := by norm_num [seed60]
  rw [seed60_liftBeta, show 4 * seed60.k = 120 from rfl, show 2 * seed60.k = 60 from rfl, e1]
    at h1
  rw [seed60_liftBeta, show 4 * seed60.k = 120 from rfl, e2] at h2
  simp only [traces60_zero] at h1 h2
  have e : (7 / 8 : ℝ) = 1 - 1 / 8 := by norm_num
  -- Rewrite `2 · (105/4)` in both places: the numerals `2 : ℝ` of the two statements may carry
  -- different instance terms, which the kernel would unfold.
  have e2b : (2 : ℝ) * (105 / 4) = 105 / 2 := by norm_num
  rw [e2b] at h2
  rw [e, spatialPurity, thermalPurity, show 2 * 60 = 120 from rfl, e2b]
  exact ⟨h1, h2⟩

end CertifiedGaps.SpinOneHeisenberg
