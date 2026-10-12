/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Gap
import CertifiedGaps.Models.SpinOneLargeD

/-!
# The spin-one Heisenberg gap, unconditionally (with open inputs)

This file states the target theorem with no hypotheses: for every even `L ≥ 2304` the periodic
spin-one Heisenberg ring has a unique ground state and `Δ_L > (log 20)/784`.

The transfer-operator input `spatialMomentBounds` (the moment inequalities at every `β > 0`,
for even lengths `≥ 2`) is proved, with finite Trotter transfer matrices. The proof is complete
except for the two numerical seeds, each closed by `sorry` and each a precise statement about
the traces `Z_L(β) = Re Tr e^{-β H_L}` of the Lean Hamiltonian:

* `seed_spatial`: `P_s(2304, 784) ≥ 0.99`.
* `seed_thermal`: `P_τ(4608, 784) ≥ 0.99`.

The two seeds are Proposition 5.5 of OpenAI's *The periodic spin-one Haldane gap* (2026).
`#print axioms spin_one_heisenberg_gap` lists `sorryAx` until both seeds are proved.
-/

open Real ThermalPurityGap

namespace SpinOneHeisenberg

/-- **Input 1 (proved).** The spatial moment inequalities for the Heisenberg traces at every
`β > 0`, from Trotterized finite transfer matrices and the product formula
(`CertifiedGaps.RPRing.spatialMomentBounds`, the Heisenberg ring being the reflection-positive
ring `SpinOneLargeD.model 0`). -/
theorem spatialMomentBounds (β : ℝ) (hβ : 0 < β) :
    SpatialMomentBounds (fun L ↦ partitionFunction L β) 2 := by
  have h := (CertifiedGaps.SpinOneLargeD.model 0).spatialMomentBounds β hβ
  simpa only [CertifiedGaps.SpinOneLargeD.partitionFunction_zero] using h

/-- **Open input 2.** The spatial seed `P_s(2304, 784) ≥ 0.99`. -/
theorem seed_spatial : 99 / 100 ≤ spatialPurity 2304 784 := by
  sorry

/-- **Open input 3.** The thermal seed `P_τ(4608, 784) ≥ 0.99`. -/
theorem seed_thermal : 99 / 100 ≤ thermalPurity 4608 784 := by
  sorry

/-- The inputs of the thermal-MPO route at `(n₀, β₀, ε₀) = (2304, 784, 1/100)`. -/
theorem thermalMPOInputs : ThermalMPOInputs 2304 784 (1 / 100) where
  moments β hβ := (spatialMomentBounds β hβ).mono (by norm_num)
  seed_spatial := by linarith [SpinOneHeisenberg.seed_spatial]
  seed_thermal := by norm_num; linarith [SpinOneHeisenberg.seed_thermal]

/-- **Theorem 1** of the thermal-MPO note. For every even `L ≥ 2304` the periodic spin-one
antiferromagnetic Heisenberg ring has a unique ground state, and every other level lies more
than `(log 20)/784 ≈ 0.003821` above it. -/
theorem spin_one_heisenberg_gap {L : ℕ} (hL : 2304 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log 20 / 784) :=
  gapGT_log_twenty thermalMPOInputs hL hLe

end SpinOneHeisenberg
