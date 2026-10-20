/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Reseed
import CertifiedGaps.Heisenberg.Seed36
import CertifiedGaps.Heisenberg.Seed60

/-!
# The ring seed at `(120, 49)` of the spin-one chain

This file proves the two ring inputs of the boundary-field gap (`OpenChain.lean`),

  `P_s(120, 49) ≥ 0.953`  and  `P_τ(240, 49) ≥ 0.9821`,

which exceed the required `0.948`, from the seed at `(36, 49/4)`: the trace table `TraceTable`
and the trial states at `72` and `120` sites. It is the step "Reseeding at `(120, b₂)`" of OpenAI's
*A boundary-field gap for the spin-one Heisenberg chain* (Proposition "Concentration at the
doubling scales"), with `b₀ = 49/4`, `b₁ = 49/2`, `b₂ = 49`, as an instance of the generic
`Seed/Reseed.lean`. All traces are shifted by `κ = 1.401482` per bond.

1. **Trace bounds at `b₀`** (`shiftedZ_le_traceBound`). The seed certificate `seed36` caps
   `|λ| < U = 1.0013641597` in the invariant sector and `|λ| < V = 0.87125` in the other three,
   with twelfth-moment masses `m_J = 1.0657405` and `m_N = 0.2783105`. The capped-mass bounds give
   `Z₃₆(b₀) ≤ 1.07340` and `Z₁₂₀(b₀) ≤ 1.17773`.
2. **The spatial seed `P_s(60, b₁) ≥ 1 - 0.0974`** (`spatial60`). The chain row
   `1 - P_s(72, b₁) ≤ 0.0605003` (`chain36`) and the trial state `Z₇₂(b₁) > 1` give one transfer
   eigenvalue with `λ₁⁷² > 31/32`. Physical squaring gives `Z₃₆(b₁) ≤ Z₃₆(b₀)² ≤ 1.15218`, so the
   relative tail at length `36` is at most `T = 0.1707` and every other eigenvalue has
   `|λ|/λ₁ ≤ ρ₁ = 0.9521`; the tail at length `60` is at most `ρ₁²⁴ T < 0.0526`.
3. **The thermal seed `P_τ(120, b₁) ≥ 1 - 0.0604`** (`thermal120`). One Boltzmann factor at length
   `120` exceeds `1` (the trial state) and `Z₁₂₀(b₀) ≤ 1.17773`, so the other factors carry at
   most `R = 0.0316 > 0.17773²` of its weight at `b₁`.
4. **One coupled update** (`chain120`) from `(60, b₁)`: `1 - P_s(120, 49) ≤ 0.047` and
   `1 - P_τ(240, 49) ≤ 0.0179`.

## Differences from the paper

The paper's numbers come from its own filters (`U = 1.00139`, `V = 0.872`, `m_J = 1.0657205`,
`m_N = 0.2783155`); this file uses those of `seed36`. The paper bounds `Z_k(b₀) ≤ U^k (1 + D(k))`
with the residual moments `D(k)` divided by `ℓ^k`, `ℓ = 0.999`; here the capped-mass bound of
every sector bounds `Z_k(b₀)` directly, which is sharper. The paper's tail length `32` needs the
powers `8/3` of the twelfth moments; here the tail length is `36`, a multiple of `12`. The
resulting numbers, Lean (paper):

* the chain row `1 - P_s(72, b₁) ≤ 0.0605003` (`0.0604996`), and `λ₁⁷² > 31/32 = 0.96875` (`0.968`);
* the tail bound `0.1707` at length `36` (`0.2092` at length `32`), `ρ₁ ≤ 0.9521` (`< 0.9525`);
* `1 - P_s(60, b₁) ≤ 0.0974` (`0.0998`), `1 - P_τ(120, b₁) ≤ 0.0604` (`0.0628`);
* `1 - P_s(120, 49) ≤ 0.047` (`0.0502`), `1 - P_τ(240, 49) ≤ 0.0179` (`0.0198`).

## Main results

* `ringPurities120_of_table`: both purities from the table and the two trial states.
* `ringSpatial120_of_table`, `ringThermal120_of_table`: the inputs `ring_spatial`, `ring_thermal`
  of `BoundaryFieldInputs`.
-/

open Finset Real Matrix

namespace CertifiedGaps.SpinOneHeisenberg

/-! ### The rational data -/

/-- The spatial reseeding step at `b₁ = 49/2`: tail length `36`, chain length `72 = 2 · 36`, target
length `60`. -/
def reseed120 : ReseedData where
  k := 36
  a := 2
  n := 60
  p := 605003 / 10 ^ 7
  p' := 6052 / 100000
  c := 1 / 32
  lk := 98425 / 100000
  B := 115218 / 100000
  T := 1707 / 10000
  rho := 9521 / 10000
  spatialDefect := 974 / 10000

/-- One coupled update from `(60, 49/2)` to `(120, 49)`. Only the step conditions are used, so
`r` and `ε` are unused placeholders. -/
def chain120 : ChainCertificate where
  s := [974e-4, 470e-4]
  t := [604e-4, 179e-4]
  r := 1
  ε := 1 / 100

theorem reseed120_valid : reseed120.Valid := by decide +kernel

theorem chain120_stepsValid : chain120.StepsValid := by decide +kernel

/-- Every sector of `seed36`, the invariant one included, admits its capped-mass bound. -/
theorem seed36_traceValid : seed36.seed.TraceValid := by decide +kernel

/-- `Z₃₆(b₀)² ≤ 1.15218`: the bound `B` of `reseed120` after physical squaring. -/
theorem seed36_traceBound_three : seed36.seed.traceBound 3 ^ 2 ≤ 115218 / 100000 := by
  decide +kernel

/-- `Z₁₂₀(b₀) ≤ 1.17773`. -/
theorem seed36_traceBound_ten : seed36.seed.traceBound 10 ≤ 117773 / 100000 := by
  decide +kernel

/-! ### Trace bounds at `b₀ = 49/4` -/

/-- The trace table in the form of the certificate. -/
theorem seed36_holdsStrict (htab : TraceTable) :
    seed36.table.HoldsStrict
      (fun g n ↦ shiftedTrace seed36.beta seed36.kappa n (seedTraces g n seed36.beta)) := by
  have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
  rw [hb]
  refine TraceBounds.holdsStrict_centered (εt := 3 / 100000) (by norm_num) fun n hn g ↦ ?_
  have e : ((3 / 100000 : ℚ) : ℝ) = 3 / 100000 := by norm_num
  rw [e]
  fin_cases g
  · exact (htab n hn).1
  · exact (htab n hn).2

/-- The sector representation at `b₀`, in the form of the certificate. -/
theorem seed36_sectorApprox {ns : Finset ℕ} (hsec : SectorInput (49 / 4) ns) :
    SectorApprox seed36.scheme (fun g n ↦ seedTraces g n seed36.beta) ns := by
  have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
  rw [hb]
  convert hsec.toSectorApprox using 1
  all_goals first | rfl | (apply heq_of_eq; funext g n; fin_cases g <;> rfl)

/-- **The trace bounds at `b₀`.** From the trace table, `e^{-b₀κ·12r} Z_{12r}(b₀) ≤ traceBound r`
for `r = 1, …, 10`. -/
theorem shiftedZ_le_traceBound (htab : TraceTable) {r : ℕ} (hr : 1 ≤ r) (hr10 : r ≤ 10) :
    shiftedZ (12 * r) (49 / 4) ≤ seed36.seed.traceBound r := by
  set ns : Finset ℕ := seedLengths ∪ (Finset.Icc 1 10).image (12 * ·)
  have hns : ∀ n ∈ ns, 2 ≤ n := by
    intro n hn
    rcases mem_union.1 hn with h | h
    · exact TableCertificate.two_le_of_mem_lengths seed36_valid h
    · obtain ⟨j, hj, rfl⟩ := mem_image.1 h
      have := (mem_Icc.1 hj).1
      omega
  have hsec : SectorInput (49 / 4) ns := fun _ hδ ↦ exists_sector_approx (by norm_num) ns hns hδ
  have h := TableCertificate.shiftedTrace_le seed36_valid seed36_traceValid
    (seed36_sectorApprox hsec)
    (fun n hn ↦ mem_union_left _ (mem_union_left _ (List.mem_toFinset.2 hn)))
    (seed36_holdsStrict htab) hr
    (mem_union_right _ (mem_image.2 ⟨r, mem_Icc.2 ⟨hr, hr10⟩, rfl⟩))
  have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
  rw [hb] at h
  exact h

/-! ### Purities and shifted traces -/

/-- The spatial purity in terms of the shifted traces. -/
theorem le_spatialPurity_iff (p : ℝ) (n : ℕ) (β : ℝ) :
    1 - p ≤ spatialPurity n β ↔ (1 - p) * shiftedZ n β ^ 2 ≤ shiftedZ (2 * n) β := by
  have hZ := partitionFunction_pos n β
  have hs : 0 < exp (-(β * kappa)) ^ (2 * n) := by positivity
  rw [spatialPurity, le_div_iff₀ (by positivity), shiftedZ, shiftedZ, shiftedTrace_eq_pow,
    shiftedTrace_eq_pow, mul_pow, ← pow_mul, mul_comm n 2, mul_left_comm]
  exact (mul_le_mul_iff_of_pos_left hs).symm

/-- The shifted traces at every `β ≥ 0` are approximated by power sums of transfer eigenvalues. -/
theorem powerSumApprox_shiftedZ {β : ℝ} (hβ : 0 ≤ β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    PowerSumApprox (fun n ↦ shiftedZ n β) ns := by
  have h : PowerSumApprox (fun n ↦ partitionFunction n β) ns := fun δ hδ ↦ by
    obtain ⟨D, _, x, hx⟩ := (SpinOneLargeD.model 0).exists_transfer_approx rfl hβ ns hns hδ
    exact ⟨D, inferInstance, x, fun n hn ↦ by
      simpa only [SpinOneLargeD.partitionFunction_zero] using hx n hn⟩
  have := h.scale (exp_pos (-(β * kappa)))
  simpa only [shiftedZ, shiftedTrace_eq_pow] using this

/-! ### The seeds at `b₁ = 49/2` -/

/-- The chain row `1 - P_s(72, 49/2) ≤ 0.0605003` (`chain36`, one doubling from `(36, 49/4)`). -/
theorem spatialPurity_seventytwo (hs : 1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4))
    (ht : 1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4)) :
    1 - 605003 / 10 ^ 7 ≤ spatialPurity 72 (49 / 2) := by
  have h := purities_ge_chainCertificate spectrum (n₀ := 36) (by decide) (β₀ := 49 / 4)
    (by norm_num) (fun β hβ ↦ (hasSpatialMomentBounds β hβ).mono (by norm_num)) chain36
    chain36_stepsValid
    (by rw [spatialPurity_spectrum]; simp [chain36, ChainCertificate.sAt]; norm_num at hs ⊢
        linarith)
    (by rw [show 2 * 36 = 72 from rfl, thermalPurity_spectrum]
        simp [chain36, ChainCertificate.tAt]; norm_num at ht ⊢; linarith) 1 (by decide)
  rw [spatialPurity_spectrum] at h
  norm_num [chain36, ChainCertificate.sAt, ChainCertificate.tAt] at h
  norm_num
  linarith [h.1]

/-- **The spatial seed `1 - P_s(60, 49/2) ≤ 0.0974`** from the trace table, the chain row at
`(72, 49/2)` and the trial state at `72` sites. -/
theorem spatial60 (htab : TraceTable) (h72 : 1 - 605003 / 10 ^ 7 ≤ spatialPurity 72 (49 / 2))
    (htrial : 1 < shiftedZ 72 (49 / 2)) : 1 - 974 / 10000 ≤ spatialPurity 60 (49 / 2) := by
  have hlen : ∀ n ∈ reseed120.lengths, 2 ≤ n := by
    simp only [ReseedData.lengths, reseed120, ReseedData.m, mem_insert, mem_singleton]
    omega
  -- `Z₃₆(b₁) ≤ Z₃₆(b₀)² ≤ B`.
  have hk : shiftedZ 36 (49 / 2) ≤ 115218 / 100000 := by
    have h1 := (boltzmannRep_partitionFunction 36).shiftedTrace_two_mul_le_sq (49 / 4) kappa 36
    rw [show (2 : ℝ) * (49 / 4) = 49 / 2 by norm_num] at h1
    have h2 := shiftedZ_le_traceBound htab (r := 3) (by norm_num) (by norm_num)
    have h3 : ((seed36.seed.traceBound 3 ^ 2 : ℚ) : ℝ) ≤ ((115218 / 100000 : ℚ) : ℝ) := by
      exact_mod_cast seed36_traceBound_three
    have h0 : 0 ≤ shiftedZ 36 (49 / 4) := by
      rw [shiftedZ, shiftedTrace]; exact mul_nonneg (exp_pos _).le (partitionFunction_pos _ _).le
    push_cast at h3
    calc shiftedZ 36 (49 / 2) ≤ shiftedZ 36 (49 / 4) ^ 2 := h1
      _ ≤ (seed36.seed.traceBound 3 : ℝ) ^ 2 := pow_le_pow_left₀ h0 h2 2
      _ ≤ _ := h3
  have hpur := (le_spatialPurity_iff _ 72 _).1 h72
  have hn : 0 < shiftedZ 60 (49 / 2) := by
    rw [shiftedZ, shiftedTrace]; exact mul_pos (exp_pos _) (partitionFunction_pos _ _)
  have hB : ((reseed120.B : ℚ) : ℝ) = 115218 / 100000 := by norm_num [reseed120]
  have h := ReseedData.spatial reseed120_valid (Z := fun n ↦ shiftedZ n (49 / 2))
    (powerSumApprox_shiftedZ (by norm_num) _ hlen) htrial
    (by simpa [reseed120, ReseedData.m] using hpur)
    (show shiftedZ 36 (49 / 2) ≤ ((reseed120.B : ℚ) : ℝ) by rw [hB]; exact hk) hn
  refine (le_spatialPurity_iff _ 60 _).2 ?_
  simpa [reseed120] using h

/-- **The thermal seed `1 - P_τ(120, 49/2) ≤ 0.0604`** from the trace table and an eigenvalue of
`H₁₂₀` below `-120κ` (the trial state at `120` sites). -/
theorem thermal120 (htab : TraceTable)
    (htrial : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    1 - 604 / 10000 ≤ thermalPurity 120 (49 / 2) := by
  have hκ : (kappa : ℝ) = 700741 / 500000 := by norm_num [kappa]
  have hq : shiftedTrace (49 / 4) kappa 120 (∑ i, exp (-(49 / 4 * spectrum 120 i))) ≤
      117773 / 100000 := by
    rw [← partitionFunction_eq_sum]
    have h1 := shiftedZ_le_traceBound htab (r := 10) (by norm_num) le_rfl
    have h2 : ((seed36.seed.traceBound 10 : ℚ) : ℝ) ≤ ((117773 / 100000 : ℚ) : ℝ) := by
      exact_mod_cast seed36_traceBound_ten
    push_cast at h2
    exact h1.trans h2
  obtain ⟨m, hm⟩ := htrial
  have h := thermalPurity_ge_of_dominant (spectrum 120) (β := 49 / 4) (κ := kappa)
    (R := 316 / 10000) (p := 604 / 10000) (by norm_num)
    ⟨m, by rw [hκ]; push_cast; linarith⟩ hq (by norm_num) (by norm_num) (by norm_num)
    (by norm_num)
  rw [show (2 : ℝ) * (49 / 4) = 49 / 2 by norm_num, ← partitionFunction_eq_sum,
    ← partitionFunction_eq_sum] at h
  exact h

/-! ### The seeds at `(120, 49)` -/

/-- **One coupled update** from `(60, 49/2)` to `(120, 49)`. -/
theorem ringPurities120_of_seeds (hs : 1 - 974 / 10000 ≤ spatialPurity 60 (49 / 2))
    (ht : 1 - 604 / 10000 ≤ thermalPurity 120 (49 / 2)) :
    1 - 47 / 1000 ≤ spatialPurity 120 49 ∧ 1 - 179 / 10000 ≤ thermalPurity 240 49 := by
  have h := purities_ge_chainCertificate spectrum (n₀ := 60) (by decide) (β₀ := 49 / 2)
    (by norm_num) (fun β hβ ↦ (hasSpatialMomentBounds β hβ).mono (by norm_num)) chain120
    chain120_stepsValid
    (by rw [spatialPurity_spectrum]; simp [chain120, ChainCertificate.sAt]; norm_num at hs ⊢
        linarith)
    (by rw [show 2 * 60 = 120 from rfl, thermalPurity_spectrum]
        simp [chain120, ChainCertificate.tAt]; norm_num at ht ⊢; linarith) 1 (by decide)
  rw [spatialPurity_spectrum, thermalPurity_spectrum] at h
  norm_num [chain120, ChainCertificate.sAt, ChainCertificate.tAt] at h
  exact ⟨by linarith [h.1], by linarith [h.2]⟩

/-- **The ring seed at `(120, 49)` from the trace table and two trial states**: a matrix-product
state at `72` sites with `500000 Tr (T⁷⁰ O) < -700741 Tr T⁷²` (as in `seed36_of_table_mps`) and
an eigenvalue of `H₁₂₀` below `-120κ`. -/
theorem ringPurities120_of_table (htab : TraceTable) {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72))
    (htrial120 : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    1 - 47 / 1000 ≤ spatialPurity 120 49 ∧ 1 - 179 / 10000 ≤ thermalPurity 240 49 :=
  have h36 := seed36_of_table_mps htab A hU
  ringPurities120_of_seeds
    (spatial60 htab (spatialPurity_seventytwo h36.1 h36.2)
      (one_lt_shiftedZ_of_mps A hU (by norm_num)))
    (thermal120 htab htrial120)

/-- **The ring input `P_s(120, 49) ≥ 0.948`** of the boundary-field gap (`seed_ring_spatial`). -/
theorem ringSpatial120_of_table (htab : TraceTable) {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72))
    (htrial120 : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    948 / 1000 ≤ spatialPurity 120 49 := by
  linarith [(ringPurities120_of_table htab A hU htrial120).1]

/-- **The ring input `P_τ(240, 49) ≥ 0.948`** of the boundary-field gap (`seed_ring_thermal`). -/
theorem ringThermal120_of_table (htab : TraceTable) {D : Type*} [Fintype D] [DecidableEq D]
    (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72))
    (htrial120 : ∃ m, spectrum 120 m < -(700741 / 500000 * 120)) :
    948 / 1000 ≤ thermalPurity 240 49 := by
  linarith [(ringPurities120_of_table htab A hU htrial120).2]

end CertifiedGaps.SpinOneHeisenberg
