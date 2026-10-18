/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Gap
import CertifiedGaps.Models.SpinOneLargeD

/-!
# The spin-one Heisenberg gap, unconditionally (with open inputs)

This file states both bounds of Theorem 1.1 of OpenAI's *The periodic spin-one Haldane gap*
(2026) with no hypotheses, each from its own initialization:

* `spin_one_heisenberg_gap_sixty`: for every even `L ≥ 60`, a unique ground state and
  `Δ_L > (4/105) log (80/79)`, from the seed at `(n₀, β₀) = (60, 105/4)` (their eq. (1.3));
* `spin_one_heisenberg_gap`: for every even `L ≥ 2304`, `Δ_L > (log 20)/784`, from the seed at
  `(36, 49/4)` and six exact doubling steps to `(2304, 784)` (their eq. (1.4); Theorem 1 of the
  thermal-MPO note).

The transfer-operator input `spatialMomentBounds` (the moment inequalities at every `β > 0`,
for even lengths `≥ 2`) is proved, with finite Trotter transfer matrices. The proofs are
complete except for four numerical seeds, each closed by `sorry` and each a precise statement
about the traces `Z_L(β) = Re Tr e^{-β H_L}` of the Lean Hamiltonian:

* `seed60_spatial`, `seed60_thermal`: `P_s(60, 105/4), P_τ(120, 105/4) ≥ 7/8`
  (OpenAI's Proposition 5.4);
* `seed36_spatial`, `seed36_thermal`: `1 - P_s(36, 49/4) ≤ 0.047916` and
  `1 - P_τ(72, 49/4) ≤ 0.181521` (the starting defects in the proof of OpenAI's
  Proposition 5.5).

The seeds `seed_spatial`, `seed_thermal` at `(2304, 784)` of the note's Theorem 1 are proved from
the second pair by the kernel-checked chain `chain36`.
-/

open Real

namespace CertifiedGaps.SpinOneHeisenberg

/-- **Input 1 (proved).** The spatial moment inequalities for the Heisenberg traces at every
`β > 0`, from Trotterized finite transfer matrices and the product formula
(`CertifiedGaps.RPRing.spatialMomentBounds`, the Heisenberg ring being the reflection-positive
ring `SpinOneLargeD.model 0`). -/
theorem spatialMomentBounds (β : ℝ) (hβ : 0 < β) :
    SpatialMomentBounds (fun L ↦ partitionFunction L β) 2 := by
  have h := (SpinOneLargeD.model 0).spatialMomentBounds β hβ
  simpa only [SpinOneLargeD.partitionFunction_zero] using h

/-- The moment inequalities in the form of the abstract criterion, for `spectrum`. -/
theorem hasSpatialMomentBounds : HasSpatialMomentBounds spectrum 2 := fun β hβ ↦ by
  rw [partitionFunction_spectrum]; exact spatialMomentBounds β hβ

/-- The abstract spatial purity of `spectrum` is `spatialPurity`. -/
theorem spatialPurity_spectrum (n : ℕ) (β : ℝ) :
    CertifiedGaps.spatialPurity spectrum n β = spatialPurity n β := by
  simp [CertifiedGaps.spatialPurity, spatialPurity, partitionFunction_spectrum]

/-- The abstract thermal purity of `spectrum` is `thermalPurity`. -/
theorem thermalPurity_spectrum (L : ℕ) (β : ℝ) :
    CertifiedGaps.thermalPurity spectrum L β = thermalPurity L β := by
  simp [CertifiedGaps.thermalPurity, thermalPurity, partitionFunction_spectrum]

/-! ### The bound for every even `L ≥ 60` -/

/-- **The gap for `L ≥ 60` from its seed** (OpenAI's Theorem 1.1, eq. (1.3)). If
`P_s(60, 105/4), P_τ(120, 105/4) ≥ 7/8`, then every even `L ≥ 60` has a unique ground state and
`Δ_L > (4/105) log (80/79) ≈ 0.000479`. This is `gapGT_of_coefficient` with `u₀ = 1/8` and
`C = 79/10 ≥ G(1/8) = 913952/117649`. -/
theorem gapGT_sixty_of_seeds (hs : 7 / 8 ≤ spatialPurity 60 (105 / 4))
    (ht : 7 / 8 ≤ thermalPurity 120 (105 / 4)) {L : ℕ} (hL : 60 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (4 / 105 * log (80 / 79)) := by
  have h := gapGT_of_coefficient spectrum (n₀ := 60) (by norm_num) (by decide)
    (β₀ := 105 / 4) (by norm_num) (fun β hβ ↦ (hasSpatialMomentBounds β hβ).mono (by norm_num))
    (u₀ := 1 / 8) (C := 79 / 10) (by norm_num) (by norm_num)
    (by unfold doublingCoefficient; norm_num) (by norm_num) (by norm_num)
    (by rw [spatialPurity_spectrum]; linarith)
    (by rw [show 2 * 60 = 120 from rfl, thermalPurity_spectrum]; linarith) hL hLe
  have e : log (1 / (79 / 10 * (1 / 8))) / (105 / 4) = 4 / 105 * log (80 / 79) := by
    norm_num; ring
  rwa [e] at h

/-- **Open input 2.** The spatial seed `P_s(60, 105/4) ≥ 7/8` (OpenAI's Proposition 5.4). -/
theorem seed60_spatial : 7 / 8 ≤ spatialPurity 60 (105 / 4) := by
  sorry

/-- **Open input 3.** The thermal seed `P_τ(120, 105/4) ≥ 7/8` (OpenAI's Proposition 5.4). -/
theorem seed60_thermal : 7 / 8 ≤ thermalPurity 120 (105 / 4) := by
  sorry

/-- **OpenAI's Theorem 1.1, eq. (1.3).** For every even `L ≥ 60` the periodic spin-one
Heisenberg ring has a unique ground state, and every other level lies more than
`(4/105) log (80/79)` above it. -/
theorem spin_one_heisenberg_gap_sixty {L : ℕ} (hL : 60 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (4 / 105 * log (80 / 79)) :=
  gapGT_sixty_of_seeds seed60_spatial seed60_thermal hL hLe

/-! ### The bound for every even `L ≥ 2304` -/

/-- The six rounded doubling steps of OpenAI's Proposition 5.5, from the defects at
`(36, 49/4)` to `(2304, 784)`, each rounded up to a multiple of `10⁻⁷`. Only the step
conditions (`ChainCertificate.StepsValid`) are used, so `r` and `ε` are unused placeholders. -/
def chain36 : ChainCertificate where
  s := [47916e-6, 605003e-7, 686630e-7, 713693e-7, 633146e-7, 376144e-7, 84709e-7]
  t := [181521e-6, 1733838e-7, 1632461e-7, 1447318e-7, 1055663e-7, 446518e-7, 55149e-7]
  r := 1
  ε := 1 / 100

/-- The six steps are valid, checked by the kernel. -/
theorem chain36_stepsValid : chain36.StepsValid := by decide +kernel

/-- **The seeds at `(2304, 784)` from the seeds at `(36, 49/4)`** (OpenAI's Proposition 5.5):
`1 - P_s(36, 49/4) ≤ 0.047916` and `1 - P_τ(72, 49/4) ≤ 0.181521` give
`P_s(2304, 784), P_τ(4608, 784) ≥ 0.99`. -/
theorem seeds2304_of_seeds36 (hs : 1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4))
    (ht : 1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4)) :
    99 / 100 ≤ spatialPurity 2304 784 ∧ 99 / 100 ≤ thermalPurity 4608 784 := by
  have h := purities_ge_chainCertificate spectrum (n₀ := 36) (by decide) (β₀ := 49 / 4)
    (by norm_num) (fun β hβ ↦ (hasSpatialMomentBounds β hβ).mono (by norm_num)) chain36
    chain36_stepsValid
    (by rw [spatialPurity_spectrum]; simp [chain36, ChainCertificate.sAt]; norm_num at hs ⊢
        linarith)
    (by rw [show 2 * 36 = 72 from rfl, thermalPurity_spectrum]
        simp [chain36, ChainCertificate.tAt]; norm_num at ht ⊢; linarith) 6 (by decide)
  rw [spatialPurity_spectrum, thermalPurity_spectrum] at h
  norm_num [chain36, ChainCertificate.sAt, ChainCertificate.tAt] at h
  exact ⟨by linarith [h.1], by linarith [h.2]⟩

/-- **Open input 4.** The spatial defect `1 - P_s(36, 49/4) ≤ 0.047916` (the starting defect
`p₀` in the proof of OpenAI's Proposition 5.5). -/
theorem seed36_spatial : 1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4) := by
  sorry

/-- **Open input 5.** The thermal defect `1 - P_τ(72, 49/4) ≤ 0.181521` (the starting defect
`q₀` in the proof of OpenAI's Proposition 5.5). -/
theorem seed36_thermal : 1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4) := by
  sorry

/-- The spatial seed `P_s(2304, 784) ≥ 0.99` of the note's Theorem 1, from `seed36_*`. -/
theorem seed_spatial : 99 / 100 ≤ spatialPurity 2304 784 :=
  (seeds2304_of_seeds36 seed36_spatial seed36_thermal).1

/-- The thermal seed `P_τ(4608, 784) ≥ 0.99` of the note's Theorem 1, from `seed36_*`. -/
theorem seed_thermal : 99 / 100 ≤ thermalPurity 4608 784 :=
  (seeds2304_of_seeds36 seed36_spatial seed36_thermal).2

/-- The inputs of the thermal-MPO route at `(n₀, β₀, ε₀) = (2304, 784, 1/100)`. -/
theorem thermalMPOInputs : ThermalMPOInputs 2304 784 (1 / 100) where
  moments β hβ := (spatialMomentBounds β hβ).mono (by norm_num)
  seed_spatial := by linarith [SpinOneHeisenberg.seed_spatial]
  seed_thermal := by norm_num; linarith [SpinOneHeisenberg.seed_thermal]

/-- **Theorem 1** of the thermal-MPO note (OpenAI's Theorem 1.1, eq. (1.4)). For every even
`L ≥ 2304` the periodic spin-one antiferromagnetic Heisenberg ring has a unique ground state,
and every other level lies more than `(log 20)/784 ≈ 0.003821` above it. -/
theorem spin_one_heisenberg_gap {L : ℕ} (hL : 2304 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log 20 / 784) :=
  gapGT_log_twenty thermalMPOInputs hL hLe

end CertifiedGaps.SpinOneHeisenberg
