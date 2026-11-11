/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.LowSpectrum
import SeedKernel.ExpBounds
import CertifiedGaps.Open.Gap

/-!
# Seed inequalities from low spectra, for any reflection-positive ring

The spacetime bootstrap needs a few finite inequalities at a seed scale `(N₀, β₀)`: the ring
purities `P_s` and `P_τ` and, for open chains, a boundary trace inequality and the open-chain
thermal purity. Each follows from low spectra (`LowSpectrum`, an eigenvalue enclosure of radius
`r` around a Ritz value and lower bounds on the rest) of a few small rings and chains. This file
states these implications for any `RPRing` and any boundary, with the remaining numerical
condition written as a decidable proposition about rationals (the exponentials replaced by
`expUpQ`), so that one `decide +kernel` checks it at any rational inverse temperature:

* `spatialPurity_ge_of_check`: `P_s(L, β) ≥ 1 - u` (`SpatialCheck`);
* `thermalPurity_ge_of_check`: `P_τ(L, β) ≥ 1 - η` for the ring (`ThermalCheck`);
* `openThermalPurity_ge_of_check`: the same for an open chain;
* `openTrace_le_of_check`: `G_{n₁} (d(u) Z_L)^m ≤ E G_{n₂}` (`TraceCheck`), the hypothesis of
  `RPRing.openGapGT_of_traces`.

A model then supplies the low spectra (from kernel-checked integer residuals and counting by
comparison) and one decidable check per seed; nothing else in a seed proof depends on the model.
-/

open Real CertifiedGaps

namespace SeedKernel

/-! ### Rational checks -/

/-- A rational upper bound on the excited-state sum `excited N β a₀ a₁ ρ`. -/
def excitedQ (k N : ℕ) (β a₀ a₁ ρ : ℚ) : ℚ :=
  expUpQ k (-(β * (a₁ - a₀))) + N * expUpQ k (-(β * (ρ - a₀)))

theorem excited_le_excitedQ {k N : ℕ} (hk : 0 < k) {β a₀ a₁ ρ : ℚ} (hβ : 0 ≤ β) (h₁ : a₀ ≤ a₁)
    (hρ : a₀ ≤ ρ) : excited N β a₀ a₁ ρ ≤ (excitedQ k N β a₀ a₁ ρ : ℝ) := by
  have hkq : (0 : ℚ) < k := by exact_mod_cast hk
  have e1 := exp_le_expUpQ hk (y := -(β * (a₁ - a₀))) (by nlinarith)
  have e2 := exp_le_expUpQ hk (y := -(β * (ρ - a₀))) (by nlinarith)
  push_cast at e1 e2
  unfold excited excitedQ
  push_cast
  gcongr

/-- The rational condition for `P_τ(L, β) ≥ 1 - η` from a low spectrum
`(θ - r, θ + r, a₁, ρ)` on `N` states. -/
abbrev ThermalCheck (k N : ℕ) (β θ r a₁ ρ η : ℚ) : Prop :=
  0 < k ∧ 0 ≤ β ∧ θ - r ≤ a₁ ∧ θ - r ≤ ρ ∧
    1 - η ≤ (1 - 4 * β * r) / (1 + excitedQ k N β (θ - r) a₁ ρ) ^ 2

/-- The rational condition for `P_s(L, β) ≥ 1 - u` from a low spectrum
`(θ₁ - r, θ₁ + r, a₁, ρ)` on `N` states at length `L` and an eigenvalue within `r` of `θ₂` at
length `2L`. -/
abbrev SpatialCheck (k N : ℕ) (β θ₁ θ₂ r a₁ ρ u : ℚ) : Prop :=
  0 < k ∧ 0 ≤ β ∧ θ₁ - r ≤ a₁ ∧ θ₁ - r ≤ ρ ∧
    1 - u ≤ (1 - β * (θ₂ + r - 2 * (θ₁ - r))) / (1 + excitedQ k N β (θ₁ - r) a₁ ρ) ^ 2

/-- The rational condition for `Z₁ (τ Z₂)^m ≤ E Z₃` from low spectra `(θᵢ - r, θᵢ + r, aᵢ, ρᵢ)`
on `Nᵢ` states (`i = 1, 2`) and an eigenvalue within `r` of `θ₃`. -/
abbrev TraceCheck (k N₁ N₂ m : ℕ) (β θ₁ a₁ ρ₁ θ₂ a₂ ρ₂ θ₃ r τ E : ℚ) : Prop :=
  0 < k ∧ 0 ≤ β ∧ 0 ≤ τ ∧ θ₁ - r ≤ a₁ ∧ θ₁ - r ≤ ρ₁ ∧ θ₂ - r ≤ a₂ ∧ θ₂ - r ≤ ρ₂ ∧
    β * (θ₃ + r - (θ₁ - r) - m * (θ₂ - r)) < k ∧
    expUpQ k (β * (θ₃ + r - (θ₁ - r) - m * (θ₂ - r))) * (1 + excitedQ k N₁ β (θ₁ - r) a₁ ρ₁) *
      ((1 + excitedQ k N₂ β (θ₂ - r) a₂ ρ₂) * τ) ^ m ≤ E

/-! ### The seeds of a reflection-positive ring -/

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] [Nonempty d] (R : RPRing d ι)

omit [DecidableEq d] [Nonempty d] in
theorem card_configs (L : ℕ) : Fintype.card (Fin L → d) = Fintype.card d ^ L := by simp

/-- **The ring thermal seed.** -/
theorem thermalPurity_ge_of_check {L k : ℕ} {β θ r a₁ ρ η : ℚ}
    (h : LowSpectrum (R.spectrum L) (θ - r) (θ + r) a₁ ρ)
    (hc : ThermalCheck k (Fintype.card d ^ L) β θ r a₁ ρ η) :
    1 - (η : ℝ) ≤ R.thermalPurity L β := by
  obtain ⟨hk, hβ, h₁, hρ, hle⟩ := hc
  have hβ' : (0 : ℝ) ≤ β := by exact_mod_cast hβ
  have ht := excited_le_excitedQ (N := Fintype.card d ^ L) hk hβ h₁ hρ
  push_cast at ht
  rw [← card_configs] at ht
  have h' : LowSpectrum (R.spectrum L) ((θ : ℝ) - r) (θ + r) a₁ ρ := by exact_mod_cast h
  have key := thermal_ratio_ge hβ' h' ht
  rw [RPRing.thermalPurity, RPRing.partitionFunction_eq_sum, RPRing.partitionFunction_eq_sum]
  refine le_trans ?_ key
  have := (Rat.cast_le (K := ℝ)).2 hle
  push_cast [card_configs] at this ⊢
  simpa only [mul_assoc] using this

/-- **The open-chain thermal seed**, for any boundary. -/
theorem openThermalPurity_ge_of_check (b : R.Boundary) {n k : ℕ} {β θ r a₁ ρ η : ℚ}
    (h : LowSpectrum (R.openSpectrum b n) (θ - r) (θ + r) a₁ ρ)
    (hc : ThermalCheck k (Fintype.card d ^ n) β θ r a₁ ρ η) :
    1 - (η : ℝ) ≤ R.openThermalPurity b n β := by
  obtain ⟨hk, hβ, h₁, hρ, hle⟩ := hc
  have hβ' : (0 : ℝ) ≤ β := by exact_mod_cast hβ
  have ht := excited_le_excitedQ (N := Fintype.card d ^ n) hk hβ h₁ hρ
  push_cast at ht
  rw [← card_configs] at ht
  have h' : LowSpectrum (R.openSpectrum b n) ((θ : ℝ) - r) (θ + r) a₁ ρ := by exact_mod_cast h
  have key := thermal_ratio_ge hβ' h' ht
  rw [RPRing.openThermalPurity, RPRing.openPartitionFunction_eq_sum,
    RPRing.openPartitionFunction_eq_sum]
  refine le_trans ?_ key
  have := (Rat.cast_le (K := ℝ)).2 hle
  push_cast [card_configs] at this ⊢
  simpa only [mul_assoc] using this

/-- **The ring spatial seed.** -/
theorem spatialPurity_ge_of_check {L k : ℕ} {β θ₁ θ₂ r a₁ ρ u : ℚ}
    (h₁ : LowSpectrum (R.spectrum L) (θ₁ - r) (θ₁ + r) a₁ ρ)
    (h₂ : ∃ i, |R.spectrum (2 * L) i - θ₂| ≤ r)
    (hc : SpatialCheck k (Fintype.card d ^ L) β θ₁ θ₂ r a₁ ρ u) :
    1 - (u : ℝ) ≤ R.spatialPurity L β := by
  obtain ⟨hk, hβ, h1, hρ, hle⟩ := hc
  have hβ' : (0 : ℝ) ≤ β := by exact_mod_cast hβ
  have ht := excited_le_excitedQ (N := Fintype.card d ^ L) hk hβ h1 hρ
  push_cast at ht
  rw [← card_configs] at ht
  have h₁' : LowSpectrum (R.spectrum L) ((θ₁ : ℝ) - r) (θ₁ + r) a₁ ρ := by exact_mod_cast h₁
  have key := spatial_ratio_ge hβ' h₁' h₂ ht
  rw [RPRing.spatialPurity, RPRing.partitionFunction_eq_sum, RPRing.partitionFunction_eq_sum]
  refine le_trans ?_ key
  have := (Rat.cast_le (K := ℝ)).2 hle
  push_cast [card_configs] at this ⊢
  exact this

omit [Nonempty d] in
/-- **The boundary trace seed**, for any boundary: from low spectra of the open chain of `n₁`
sites and the ring of `L` sites and an eigenvalue of the open chain of `n₂` sites,
`G_{n₁}(β) (d(u) Z_L(β))^m ≤ E G_{n₂}(β)`. -/
theorem openTrace_le_of_check (b : R.Boundary) {n₁ L n₂ k m : ℕ}
    {β θ₁ a₁ ρ₁ θ₂ a₂ ρ₂ θ₃ r u E : ℚ}
    (h₁ : LowSpectrum (R.openSpectrum b n₁) (θ₁ - r) (θ₁ + r) a₁ ρ₁)
    (h₂ : LowSpectrum (R.spectrum L) (θ₂ - r) (θ₂ + r) a₂ ρ₂)
    (h₃ : ∃ i, |R.openSpectrum b n₂ i - θ₃| ≤ r)
    (hc : TraceCheck k (Fintype.card d ^ n₁) (Fintype.card d ^ L) m β θ₁ a₁ ρ₁ θ₂ a₂ ρ₂ θ₃ r
      (OpenCertificate.tailQ u) E) :
    R.openPartitionFunction b n₁ β * (tailBound u * R.partitionFunction L β) ^ m ≤
      E * R.openPartitionFunction b n₂ β := by
  obtain ⟨hk, hβ, hτ, ha₁, hρ₁, ha₂, hρ₂, hy, hle⟩ := hc
  have hβ' : (0 : ℝ) ≤ β := by exact_mod_cast hβ
  have hs₁ : excited (Fintype.card (Fin n₁ → d)) β (θ₁ - r) a₁ ρ₁ ≤
      (excitedQ k (Fintype.card d ^ n₁) β (θ₁ - r) a₁ ρ₁ : ℝ) := by
    rw [card_configs]; exact_mod_cast excited_le_excitedQ hk hβ ha₁ hρ₁
  have hs₂ : excited (Fintype.card (Fin L → d)) β (θ₂ - r) a₂ ρ₂ ≤
      (excitedQ k (Fintype.card d ^ L) β (θ₂ - r) a₂ ρ₂ : ℝ) := by
    rw [card_configs]; exact_mod_cast excited_le_excitedQ hk hβ ha₂ hρ₂
  have h₁' : LowSpectrum (R.openSpectrum b n₁) ((θ₁ : ℝ) - r) (θ₁ + r) a₁ ρ₁ := by
    exact_mod_cast h₁
  have h₂' : LowSpectrum (R.spectrum L) ((θ₂ : ℝ) - r) (θ₂ + r) a₂ ρ₂ := by exact_mod_cast h₂
  have hexp := exp_le_expUpQ hk hy
  have hE : exp ((β : ℝ) * (θ₃ + r - (θ₁ - r) - m * (θ₂ - r))) *
      (1 + (excitedQ k (Fintype.card d ^ n₁) β (θ₁ - r) a₁ ρ₁ : ℝ)) *
      ((1 + (excitedQ k (Fintype.card d ^ L) β (θ₂ - r) a₂ ρ₂ : ℝ)) *
        (OpenCertificate.tailQ u : ℝ)) ^ m ≤ E := by
    have := (Rat.cast_le (K := ℝ)).2 hle
    push_cast at this hexp
    refine le_trans ?_ this
    have h0 : (0 : ℝ) ≤ OpenCertificate.tailQ u := by exact_mod_cast hτ
    have hx1 : (0 : ℝ) ≤ 1 + (excitedQ k (Fintype.card d ^ n₁) β (θ₁ - r) a₁ ρ₁ : ℝ) := by
      linarith [(excited_nonneg _ _ _ _ _).trans hs₁]
    have hx2 : (0 : ℝ) ≤ 1 + (excitedQ k (Fintype.card d ^ L) β (θ₂ - r) a₂ ρ₂ : ℝ) := by
      linarith [(excited_nonneg _ _ _ _ _).trans hs₂]
    gcongr
  have key := trace_le m hβ' h₁' h₂' h₃ (by exact_mod_cast hτ) hs₁ hs₂ hE
  rw [RPRing.openPartitionFunction_eq_sum, RPRing.openPartitionFunction_eq_sum,
    RPRing.partitionFunction_eq_sum, OpenCertificate.cast_tailQ]
  exact key

end SeedKernel
