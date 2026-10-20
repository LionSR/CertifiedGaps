/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryInit
import CertifiedGaps.Heisenberg.Reseed120

/-!
# The boundary seeds from the trace table

`Heisenberg/BoundaryInit.lean` derives the boundary inputs `seed_boundary` and
`seed_open_thermal` from the data H1 (`ThermalData`), H2 (`TrialData`) and two ring purities.
This file discharges whatever of these the seed-36 trace table already gives:

* the upper end `λ₀ ≤ U = 1.00139` of the window in H1: the leading shifted modulus satisfies
  `λ₀¹²⁰ ≤ Z̃₁₂₀(b₀) ≤ traceBound 10 ≤ 1.17773 ≤ U¹²⁰` (`shiftedZ_le_traceBound`);
* the ring purity `P_s(72, 49/2) ≥ 1 - .0605003` (the chain row from the seeds at `(36, 49/4)`,
  `spatialPurity_seventytwo`), and `P_s(120, 49) ≥ .948` (`ringSpatial120_of_table`), from the
  table and the trial states of the periodic seed.

What remains of H1 is `ThermalCore`: the witness with `λ₀ ≥ 9/10`, `ρ₀ ≤ .874` and
`c₀ P_*(λ₀)² ≤ .003745`, and `Y₆(b₀) ≤ .027`. The ring traces alone do not bound the residual
ratio `ρ₀`: the paper bounds it sector by sector. `Heisenberg/BoundaryLeading.lean` does so on the
Trotter approximants and proves `ThermalCore` from this table and the boundary trace table
(`thermalCore_of_tables`).

## Main results

* `ThermalData.of_table`: H1 from `ThermalCore` and the trace table.
* `boundaryDefect_le_of_table`, `oppositeFieldThermalPurity_ge_of_table`: the statements of
  `seed_boundary` and `seed_open_thermal` from `ThermalCore`, `TrialData`, the trace table and
  the periodic trial states.
-/

open Matrix Real

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryInit

/-- **The part of H1 the trace table does not give**: a witness of the shifted boundary moment
inequalities at `b₀` with `λ₀ ≥ 9/10`, `ρ₀ ≤ .874` and `c₀ P_*(λ₀)² ≤ .003745`, and the shifted
six-site trace `Y₆(b₀) ≤ .027`. -/
structure ThermalCore : Prop where
  witness : ∃ μ c ρ s : ℝ, OpenMomentWitness (shiftedRingZ (49 / 4)) (shiftedOpenZ (49 / 4))
    μ c ρ s ∧ 9 / 10 ≤ μ ∧ ρ ≤ 874 / 1000 ∧
    c * (μ ^ 3 - 7 / 50 * μ ^ 2 - 9 / 50 * μ) ^ 2 ≤ 3745 / 10 ^ 6
  sixth : shiftFactor (49 / 4) ^ 6 * oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) ≤
    27 / 1000

/-- The shifted ring traces of `BoundaryInit` are those of `Seed36`. -/
theorem shiftedRingZ_eq_shiftedZ (L : ℕ) (β : ℝ) : shiftedRingZ β L = shiftedZ L β := by
  rw [shiftedRingZ, shiftedZ, shiftedTrace_eq_pow, shiftFactor, kappa]
  norm_num

/-- **H1 from its core and the trace table.** The upper window end `λ₀ ≤ 1.00139` follows from
`λ₀¹²⁰ ≤ Z̃₁₂₀(b₀) ≤ 1.17773 ≤ 1.00139¹²⁰`. -/
theorem ThermalData.of_table (htab : TraceTable) (h : ThermalCore) : ThermalData := by
  obtain ⟨μ, c, ρ, s, w, hμl, hρ, hC⟩ := h.witness
  refine ⟨⟨μ, c, ρ, s, w, hμl, ?_, hρ, hC⟩, h.sixth⟩
  have h1 : μ ^ (2 * 60) ≤ shiftedRingZ (49 / 4) (2 * 60) := w.pow_le_Z 60 (by norm_num)
  have h2 := shiftedZ_le_traceBound htab (r := 10) (by norm_num) le_rfl
  have h3 : ((seed36.seed.traceBound 10 : ℚ) : ℝ) ≤ ((117773 / 100000 : ℚ) : ℝ) := by
    exact_mod_cast seed36_traceBound_ten
  rw [shiftedRingZ_eq_shiftedZ] at h1
  have h4 : μ ^ 120 ≤ (100139 / 100000 : ℝ) ^ 120 := by
    have : (117773 / 100000 : ℝ) ≤ (100139 / 100000) ^ 120 := by norm_num
    push_cast at h3
    exact h1.trans (h2.trans (h3.trans this))
  exact (pow_le_pow_iff_left₀ (by linarith) (by norm_num) (by norm_num)).1 h4

/-- **The boundary seed from the trace table** (the statement of `seed_boundary`): from the core
H1 data, the trial bounds H2, the trace table, the matrix-product trial state of the periodic
seed at `72` sites and an eigenvalue of the `120`-site ring below `-120κ`. -/
theorem boundaryDefect_le_of_table (h1 : ThermalCore) (h2 : TrialData) (htab : TraceTable)
    {D : Type*} [Fintype D] [DecidableEq D] (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72))
    (htrial120 : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    boundaryDefect (fun L ↦ partitionFunction L 49)
      (fun n ↦ oppositeFieldPartitionFunction (3 / 5) n 49) 142 ≤ 14 / 100 := by
  have h36 := seed36_of_table_mps htab A hU
  exact boundaryDefect_le_of_data (ThermalData.of_table htab h1) h2
    (spatialPurity_seventytwo h36.1 h36.2) (ringSpatial120_of_table htab A hU htrial120)

/-- **The open-chain thermal seed from the trace table** (the statement of
`seed_open_thermal`): from the core H1 data, the trial bounds H2 and the trace table. -/
theorem oppositeFieldThermalPurity_ge_of_table (h1 : ThermalCore) (h2 : TrialData)
    (htab : TraceTable) : 79 / 100 ≤ oppositeFieldThermalPurity (3 / 5) 240 49 :=
  oppositeFieldThermalPurity_ge_of_data (ThermalData.of_table htab h1) h2

end BoundaryInit

end CertifiedGaps.SpinOneHeisenberg
