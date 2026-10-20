/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Reseed120
import SeedKernel.Heisenberg.Seed60

/-!
# The ring inputs of the boundary-field gap from the trace table alone

The kernel-checked trial state (`trial_check_72`, `trial_check_120`: OpenAI's integer
matrix-product state of bond dimension 26 on 72 and 120 sites) discharges both trial-state inputs
of `ringPurities120_of_table`. What remains is the trace table at `n = 4, …, 12` (`TraceTable`,
shared with the seed at `(36, 49/4)`):

* `ringPurities120_of_traceTable`: `1 - P_s(120, 49) ≤ 0.047` and `1 - P_τ(240, 49) ≤ 0.0179`;
* `seed_ring_spatial_of_table`, `seed_ring_thermal_of_table`: the inputs `ring_spatial`,
  `ring_thermal` of `BoundaryFieldInputs` (`P_s(120, 49), P_τ(240, 49) ≥ 0.948`).
-/

namespace CertifiedGaps.SpinOneHeisenberg

/-- **The ring seed at `(120, 49)` from the trace table.** -/
theorem ringPurities120_of_traceTable (htab : TraceTable) :
    1 - 47 / 1000 ≤ spatialPurity 120 49 ∧ 1 - 179 / 10000 ≤ thermalPurity 240 49 :=
  ringPurities120_of_table htab trialA trial_check_72.2 exists_spectrum_lt_onetwenty

/-- **The ring input `P_s(120, 49) ≥ 0.948` from the trace table.** -/
theorem seed_ring_spatial_of_table (htab : TraceTable) : 948 / 1000 ≤ spatialPurity 120 49 :=
  ringSpatial120_of_table htab trialA trial_check_72.2 exists_spectrum_lt_onetwenty

/-- **The ring input `P_τ(240, 49) ≥ 0.948` from the trace table.** -/
theorem seed_ring_thermal_of_table (htab : TraceTable) : 948 / 1000 ≤ thermalPurity 240 49 :=
  ringThermal120_of_table htab trialA trial_check_72.2 exists_spectrum_lt_onetwenty

end CertifiedGaps.SpinOneHeisenberg
