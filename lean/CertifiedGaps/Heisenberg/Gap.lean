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
`SpinOneHeisenberg.hamiltonian L = ∑ⱼ Sⱼ · Sⱼ₊₁` on `(ℂ³)^{⊗L}`. It proves the theorem from
explicit inputs, so that each input can be discharged separately later.

## The inputs

`ThermalMPOInputs n₀ β₀ ε₀` collects everything the doubling argument needs about the ring. Each
field is a statement about `partitionFunction L β = Re Tr e^{-β H_L}`, a trace of a matrix
exponential that Lean defines:
* `moments`: the spatial moment inequalities `SpatialMomentBounds` at every `β > 0`. These are
  closed inequalities between traces at different lengths.
* `seed_spatial`, `seed_thermal`: the two seed purities `P_s(n₀, β₀), P_τ(2n₀, β₀) ≥ 1 - ε₀`.

## Discharging the inputs

* `moments` follows from the moment form `Z_L(β) = ∑ᵢ μᵢ^L`
  (`ThermalMPOInputs.of_hasMomentForm`). The continuous-time thermal MPO of the note gives that
  form with `μᵢ = |λᵢ(β)|`, where the `λᵢ(β)` are the eigenvalues of the self-adjoint
  Hilbert–Schmidt transfer operator `X_β`.
* Infinite-dimensional operators can also be avoided: `moments` follows from finite approximants
  with the moment form that converge at each even length
  (`ThermalMPOInputs.of_tendsto_hasMomentForm`). With the even/odd Trotter splitting into
  `M` slices, `Tr (K_in K_out)^M = Tr (T_o T_e)^{L/2}` with `T_o T_e = T_e^† T_e ⪰ 0`. This gives
  the moment form for each `M`, and the Lie–Trotter formula gives convergence as `M → ∞`.
* The seeds are finite numerical statements about two rings of 2304 and 4608 sites, taken
  from Proposition 5.5 of OpenAI's *The periodic spin-one Haldane gap* (2026).

## Main results

* `gapGT_of_inputs`: for every even `L ≥ n₀`, a nondegenerate ground level and
  `Δ_L > log (1/(5ε₀)) / β₀`.
* `gapGT_log_twenty`: Theorem 1 of the note, `Δ_L > (log 20)/784` for every even `L ≥ 2304`.
-/

open Filter Real ThermalPurityGap

namespace SpinOneHeisenberg

/-- The spatial purity `P_s(n, β) = Z_{2n}(β) / Z_n(β)²` of the Heisenberg ring. -/
noncomputable def spatialPurity (n : ℕ) (β : ℝ) : ℝ :=
  partitionFunction (2 * n) β / partitionFunction n β ^ 2

/-- The thermal purity `P_τ(L, β) = Z_L(2β) / Z_L(β)² = Tr ρ_L(β)²` of the Heisenberg ring. -/
noncomputable def thermalPurity (L : ℕ) (β : ℝ) : ℝ :=
  partitionFunction L (2 * β) / partitionFunction L β ^ 2

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

/-- The trace partition function is the spectral one of `ThermalPurityGap`. -/
theorem thermalPurityGap_partitionFunction :
    ThermalPurityGap.partitionFunction spectrum = partitionFunction := by
  funext L β
  exact (partitionFunction_eq_sum L β).symm

/-- The moment form at every `β > 0` gives the `moments` input. -/
theorem ThermalMPOInputs.of_hasMomentForm {n₀ : ℕ} (hn₀ : 0 < n₀) {β₀ ε₀ : ℝ}
    (hM : ∀ β : ℝ, 0 < β → HasMomentForm (fun L ↦ partitionFunction L β) n₀)
    (hs : 1 - ε₀ ≤ spatialPurity n₀ β₀) (ht : 1 - ε₀ ≤ thermalPurity (2 * n₀) β₀) :
    ThermalMPOInputs n₀ β₀ ε₀ :=
  ⟨fun β hβ ↦ (hM β hβ).spatialMomentBounds hn₀, hs, ht⟩

/-- Finite approximants: if for each `β > 0` there are functions `Zs M` with the moment form
for large `M`, converging to the traces at every even length `L ≥ n₀`, then the `moments` input
holds. This is the interface for Trotterized transfer matrices. -/
theorem ThermalMPOInputs.of_tendsto_hasMomentForm {n₀ : ℕ} (hn₀ : 0 < n₀) {β₀ ε₀ : ℝ}
    (Zs : ℝ → ℕ → ℕ → ℝ)
    (hM : ∀ β : ℝ, 0 < β → ∀ᶠ M in atTop, HasMomentForm (Zs β M) n₀)
    (hlim : ∀ β : ℝ, 0 < β → ∀ L, n₀ ≤ L → Even L →
      Tendsto (fun M ↦ Zs β M L) atTop (nhds (partitionFunction L β)))
    (hs : 1 - ε₀ ≤ spatialPurity n₀ β₀) (ht : 1 - ε₀ ≤ thermalPurity (2 * n₀) β₀) :
    ThermalMPOInputs n₀ β₀ ε₀ :=
  ⟨fun β hβ ↦ SpatialMomentBounds.of_tendsto_hasMomentForm hn₀ (hM β hβ) (hlim β hβ), hs, ht⟩

/-- **Uniform gap from thermal-MPO inputs.** Let `n₀ > 0` be even, `β₀ > 0`,
`0 < ε₀ ≤ 1/100`. Under `ThermalMPOInputs n₀ β₀ ε₀`, for every even `L ≥ n₀` the Heisenberg ring
has a nondegenerate ground level and every other level lies more than `log (1/(5ε₀)) / β₀`
above it. -/
theorem gapGT_of_inputs {n₀ : ℕ} (hn₀ : 0 < n₀) (hn₀e : Even n₀) {β₀ ε₀ : ℝ} (hβ₀ : 0 < β₀)
    (hε₀ : 0 < ε₀) (hε₁ : ε₀ ≤ 1 / 100) (h : ThermalMPOInputs n₀ β₀ ε₀) {L : ℕ} (hL : n₀ ≤ L)
    (hLe : Even L) : GapGT (spectrum L) (log (1 / (5 * ε₀)) / β₀) := by
  have hZ : HasSpatialMomentBounds spectrum n₀ := by
    intro β hβ
    rw [thermalPurityGap_partitionFunction]
    exact h.moments β hβ
  refine gapGT_of_seeds spectrum hn₀ hn₀e hβ₀ hZ hε₀ hε₁ ?_ ?_ hL hLe
  · simpa [ThermalPurityGap.spatialPurity, spatialPurity, thermalPurityGap_partitionFunction]
      using h.seed_spatial
  · simpa [ThermalPurityGap.thermalPurity, thermalPurity, thermalPurityGap_partitionFunction]
      using h.seed_thermal

/-- **Theorem 1** of the thermal-MPO note, from its inputs: if the spatial moment inequalities
hold and `P_s(2304, 784), P_τ(4608, 784) ≥ 0.99`, then for every even `L ≥ 2304` the periodic
spin-one Heisenberg ring has a unique ground state and `Δ_L > (log 20)/784 ≈ 0.003821`. -/
theorem gapGT_log_twenty (h : ThermalMPOInputs 2304 784 (1 / 100)) {L : ℕ} (hL : 2304 ≤ L)
    (hLe : Even L) : GapGT (spectrum L) (log 20 / 784) := by
  have := gapGT_of_inputs (by norm_num) (by decide) (by norm_num) (by norm_num) le_rfl h hL hLe
  norm_num at this
  exact this

end SpinOneHeisenberg
