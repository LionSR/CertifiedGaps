/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.PurityCriterion
import CertifiedGaps.Heisenberg.Hamiltonian

/-!
# A uniform gap for the periodic spin-one Heisenberg ring from thermal-MPO inputs

This file states Theorem 1 of the thermal-MPO note for the Lean Hamiltonian
`SpinOneHeisenberg.hamiltonian L = ∑ⱼ Sⱼ · Sⱼ₊₁` on `(ℂ³)^{⊗L}` and proves it from explicit
inputs, without `sorry`.

## The inputs

`ThermalMPOInputs n₀ β₀ ε₀` collects everything the doubling argument needs about the ring. Each
field is a statement about `partitionFunction L β = Re Tr e^{-β H_L}`, a trace of a matrix
exponential that Lean defines:
* `moments`: the spatial moment inequalities `SpatialMomentBounds` at every `β > 0`. These are
  closed inequalities between traces at different lengths.
* `seed_spatial`, `seed_thermal`: the two seed purities `P_s(n₀, β₀), P_τ(2n₀, β₀) ≥ 1 - ε₀`.

`Heisenberg.Main` proves `moments` (the Heisenberg ring is the reflection-positive ring
`SpinOneLargeD.model 0`) and states the two seeds, which are finite numerical statements about
rings of 2304 and 4608 sites (Proposition 5.5 of OpenAI's *The periodic spin-one Haldane gap*,
2026).

## Main results

* `gapGT_of_inputs`: for every even `L ≥ n₀`, a nondegenerate ground level and
  `Δ_L > log (1/(5ε₀)) / β₀`.
* `gapGT_log_twenty`: Theorem 1 of the note, `Δ_L > (log 20)/784` for every even `L ≥ 2304`.
-/

open Real

namespace CertifiedGaps.SpinOneHeisenberg

/-- The spatial purity `P_s(n, β) = Z_{2n}(β) / Z_n(β)²` of the Heisenberg ring. -/
noncomputable def spatialPurity (n : ℕ) (β : ℝ) : ℝ :=
  partitionFunction (2 * n) β / partitionFunction n β ^ 2

/-- The thermal purity `P_τ(L, β) = Z_L(2β) / Z_L(β)² = Tr ρ_L(β)²` of the Heisenberg ring. -/
noncomputable def thermalPurity (L : ℕ) (β : ℝ) : ℝ :=
  partitionFunction L (2 * β) / partitionFunction L β ^ 2

/-- The purities of the Heisenberg ring are those of its partition function. -/
theorem spatialPurity_eq_spatialPurityOf (n : ℕ) (β : ℝ) :
    spatialPurity n β = spatialPurityOf partitionFunction n β := rfl

theorem thermalPurity_eq_thermalPurityOf (L : ℕ) (β : ℝ) :
    thermalPurity L β = thermalPurityOf partitionFunction L β := rfl

/-- The inputs of the thermal-MPO route to the gap, as statements about the traces
`Re Tr e^{-β H_L}`: the spatial moment inequalities at every `β > 0` for even lengths `≥ n₀`,
and the two seed purities at `(n₀, β₀)` with defect `ε₀`. -/
structure ThermalMPOInputs (n₀ : ℕ) (β₀ ε₀ : ℝ) : Prop where
  /-- The spatial moment inequalities, supplied by a self-adjoint transfer operator. -/
  moments : ∀ β : ℝ, 0 < β → SpatialMomentBounds (fun L ↦ partitionFunction L β) n₀
  /-- The spatial seed `P_s(n₀, β₀) ≥ 1 - ε₀`. -/
  seed_spatial : 1 - ε₀ ≤ spatialPurity n₀ β₀
  /-- The thermal seed `P_τ(2n₀, β₀) ≥ 1 - ε₀`. -/
  seed_thermal : 1 - ε₀ ≤ thermalPurity (2 * n₀) β₀

/-- The trace partition function is the spectral one of the abstract criterion. -/
theorem partitionFunction_spectrum :
    CertifiedGaps.partitionFunction spectrum = partitionFunction := by
  funext L β
  exact (partitionFunction_eq_sum L β).symm

/-- **Uniform gap from thermal-MPO inputs.** Let `n₀ > 0` be even, `β₀ > 0`,
`0 < ε₀ ≤ 1/100`. Under `ThermalMPOInputs n₀ β₀ ε₀`, for every even `L ≥ n₀` the Heisenberg ring
has a nondegenerate ground level and every other level lies more than `log (1/(5ε₀)) / β₀`
above it. -/
theorem gapGT_of_inputs {n₀ : ℕ} (hn₀ : 0 < n₀) (hn₀e : Even n₀) {β₀ ε₀ : ℝ} (hβ₀ : 0 < β₀)
    (hε₀ : 0 < ε₀) (hε₁ : ε₀ ≤ 1 / 100) (h : ThermalMPOInputs n₀ β₀ ε₀) {L : ℕ} (hL : n₀ ≤ L)
    (hLe : Even L) : GapGT (spectrum L) (log (1 / (5 * ε₀)) / β₀) := by
  have hZ : HasSpatialMomentBounds spectrum n₀ := by
    intro β hβ
    rw [partitionFunction_spectrum]
    exact h.moments β hβ
  refine gapGT_of_seeds spectrum hn₀ hn₀e hβ₀ hZ hε₀ hε₁ ?_ ?_ hL hLe
  · simpa [CertifiedGaps.spatialPurity, spatialPurity, partitionFunction_spectrum]
      using h.seed_spatial
  · simpa [CertifiedGaps.thermalPurity, thermalPurity, partitionFunction_spectrum]
      using h.seed_thermal

/-- **Theorem 1** of the thermal-MPO note, from its inputs: if the spatial moment inequalities
hold and `P_s(2304, 784), P_τ(4608, 784) ≥ 0.99`, then for every even `L ≥ 2304` the periodic
spin-one Heisenberg ring has a unique ground state and `Δ_L > (log 20)/784 ≈ 0.003821`. -/
theorem gapGT_log_twenty (h : ThermalMPOInputs 2304 784 (1 / 100)) {L : ℕ} (hL : 2304 ≤ L)
    (hLe : Even L) : GapGT (spectrum L) (log 20 / 784) := by
  have := gapGT_of_inputs (by norm_num) (by decide) (by norm_num) (by norm_num) le_rfl h hL hLe
  norm_num at this
  exact this

end CertifiedGaps.SpinOneHeisenberg
