/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.OpenChain
import CertifiedGaps.Open.Criterion
import CertifiedGaps.Open.Initialization
import CertifiedGaps.Trace.Poisson

/-!
# The boundary seeds of the spin-one chain from short-chain data

This file derives the two boundary inputs of `Heisenberg/OpenChain.lean`,

* `seed_boundary`: the boundary defect `e₁₄₂ ≤ 0.14` of the opposite-field chain at `β = 49`,
* `seed_open_thermal`: the opposite-field thermal purity `P_τ(240, 49) ≥ 0.79`,

from finite data at `b₀ = 49/4`, following the subsection "Starting residual and
physical-purity bounds" of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain*
(2026): Lemma in:seedbounds, then the two boundary steps (in:78) and (in:142), and the purity
bound `1 - A(240, 2) < .21`.

## The data

All shifted quantities use the paper's shift `a = 700741/500000` per site: a shifted trace at
inverse temperature `β` and length `n` is `e^{-βan}` times the trace. The leading modulus `λ₀`
below is that of the shifted traces (`λ₀ = e^{-b₀a} μ₀`); the leading weight `c₀` needs no shift.

* `ThermalData` (the paper's eq. (in:thermal) with the periodic window, hypothesis H1): a
  witness `(λ₀, c₀, ρ₀, s₀)` of the boundary moment inequalities for the shifted ring and
  open-chain traces at `b₀` with `9/10 ≤ λ₀ ≤ U = 1.00139`, `ρ₀ ≤ .874` and
  `c₀ P_*(λ₀)² ≤ .003745`, `P_*(x) = x³ - .14x² - .18x`; and the shifted six-site trace
  `Y₆(b₀) ≤ .027`.
* `TrialData` (eq. (in:trials), hypothesis H2): at `m = 34, 78, 142, 240` the opposite-field
  chain has an eigenvalue `E` with `E + am < α = .4042`.
* Two ring purities (periodic inputs): `P_s(72, 49/2) ≥ 1 - .0605003`, the first row of the
  chain `chain36` from the seeds at `(36, 49/4)`, for `ρ₁`; and `P_s(120, 49) ≥ .948`, the
  input `seed_ring_spatial`, for `ρ₂`.

## Main results

* `boundaryDefect_le_of_data`: the statement of `seed_boundary` from H1, H2 and the two ring
  purities.
* `oppositeFieldThermalPurity_ge_of_data`: the statement of `seed_open_thermal` from H1 and H2.

`Heisenberg/BoundaryInitTable.lean` discharges the window end `λ₀ ≤ U` and both ring purities
from the seed-36 trace table and the periodic trial states.

## Differences from the paper

* The residual bound uses `λ₀^{m-6} ≤ U^{m-6}` instead of the paper's `λ₀^m / ℓ⁶ ≤ U^m / ℓ⁶`:
  `ε(m) = .027 e^{b₀α} (.874 U)^{m-6}`, slightly smaller than eq. (in:epsR) (`ε(34) < .0915`
  instead of `.0928`). The lower window end `ℓ` then enters only through the monotonicity of
  `x^m / P_*(x)²`, which holds on `[9/10, ∞)`. So `ℓ = 9/10` instead of `ℓ = .999`; this
  resolves the issue that the witness trace inequalities alone give only `λ₀ ≥ .9986`.
* The leading coefficient `e^{b₀α} · .003745 / P_*(U)²` is used through the rational bound
  `e^{b₀α} ≤ 141.38` instead of being relaxed to `e^{.132}`, so
  `R(m) = 141.38 (.003745 U^m / P_*(U)² + .027 (.874 U)^{m-6}) - 1`. The purity lower bounds
  are `≥ .6266, .7632, .8050` instead of the paper's `.6180, .7533, .7923`.
* `ρ₁` comes from the ring purity `P_s(72, 49/2) ≥ 1 - .0605003` (`ρ₁⁷² ≤ d(.0605003)`, so
  `ρ₁⁴⁴ ≤ .125`) instead of the paper's `ρ₁ < .9525` (`ρ₁⁴⁴ < .1182`), which needs the
  reseeding argument; `ρ₂⁶⁴ ≤ .1491` comes from `d(.052)` directly instead of from `ρ₂ < .971`.
  The outputs are `e₇₈(b₁) ≤ .1335` (paper `.133`) and `e₁₄₂(b₂) ≤ .1306 < .14`.
* The thermal inequalities are stated with `≤` instead of `<`.
-/

open Real

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryInit

open PolynomialFilter

/-! ### Traces, shifts and the data -/

/-- The cubic `P_*(x) = x³ - .14x² - .18x = x (x - 1/2)(x + 9/25)` of the leading-weight bound,
as a coefficient list (constant term first). -/
def pStar : List ℚ := [0, -9 / 50, -7 / 50, 1]

theorem eval_pStar (x : ℝ) : eval pStar x = x ^ 3 - 7 / 50 * x ^ 2 - 9 / 50 * x := by
  simp only [pStar, eval_cons, eval_nil]; push_cast; ring

/-- The open-chain spectra of the Heisenberg ring with end term `-(3/5) Sᶻ`: the equal-field
chain at odd lengths and the opposite-field chain at even lengths. -/
noncomputable abbrev openSpec : (n : ℕ) → (Fin n → Fin 3) → ℝ :=
  (SpinOneLargeD.model 0).openSpectrum (boundary (3 / 5))

/-- The open-chain traces `⟨v, Xⁿ v⟩` at inverse temperature `β`. -/
noncomputable def openZ (β : ℝ) (n : ℕ) : ℝ :=
  (SpinOneLargeD.model 0).openPartitionFunction (boundary (3 / 5)) n β

/-- The shift factor `e^{-βa}` per site, `a = 700741/500000`. -/
noncomputable def shiftFactor (β : ℝ) : ℝ := exp (-(β * (700741 / 500000)))

/-- The shifted ring traces `e^{-βaL} Z_L(β)`. -/
noncomputable def shiftedRingZ (β : ℝ) (L : ℕ) : ℝ := shiftFactor β ^ L * partitionFunction L β

/-- The shifted open-chain traces `e^{-βan} ⟨v, Xⁿ v⟩`. -/
noncomputable def shiftedOpenZ (β : ℝ) (n : ℕ) : ℝ := shiftFactor β ^ n * openZ β n

/-- The spectral part of H1 for one witness `(λ₀, c₀, ρ₀, s₀)` of the shifted boundary moment
inequalities at `b₀`: `9/10 ≤ λ₀ ≤ 1.00139`, `ρ₀ ≤ .874` and `c₀ P_*(λ₀)² ≤ .003745`. -/
def WitnessData (μ c ρ s : ℝ) : Prop :=
  OpenMomentWitness (shiftedRingZ (49 / 4)) (shiftedOpenZ (49 / 4)) μ c ρ s ∧
    9 / 10 ≤ μ ∧ μ ≤ 100139 / 100000 ∧ ρ ≤ 874 / 1000 ∧
    c * (μ ^ 3 - 7 / 50 * μ ^ 2 - 9 / 50 * μ) ^ 2 ≤ 3745 / 10 ^ 6

/-- **H1: the thermal and spectral data at `b₀ = 49/4`** (the paper's eq. (in:thermal) and the
window `ℓ < λ₀ ≤ U`, `ρ₀ < .874` of Proposition per:concentration), for the shifted traces. -/
structure ThermalData : Prop where
  /-- A witness of the shifted boundary moment inequalities with `9/10 ≤ λ₀ ≤ 1.00139`,
  `ρ₀ ≤ .874` and `c₀ P_*(λ₀)² ≤ .003745`. -/
  witness : ∃ μ c ρ s : ℝ, WitnessData μ c ρ s
  /-- The shifted six-site opposite-field trace `Y₆(b₀) ≤ .027`. -/
  sixth : shiftFactor (49 / 4) ^ 6 * oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) ≤
    27 / 1000

/-- The variational bound at length `m`: the opposite-field chain has an eigenvalue `E` with
`E + am < α`, `a = 700741/500000`, `α = 2021/5000`. -/
def TrialBelow (m : ℕ) : Prop :=
  ∃ i, oppositeFieldSpectrum (3 / 5) m i + 700741 / 500000 * m < 2021 / 5000

/-- **H2: the four trial bounds** `E₀(H_m^{opp}) + am < α` (the paper's eq. (in:trials)). -/
structure TrialData : Prop where
  trial34 : TrialBelow 34
  trial78 : TrialBelow 78
  trial142 : TrialBelow 142
  trial240 : TrialBelow 240

/-! ### Basic facts about the traces -/

theorem partitionFunction_pos' (L : ℕ) (β : ℝ) : 0 < partitionFunction L β := by
  rw [← partitionFunction_spectrum]; exact CertifiedGaps.partitionFunction_pos spectrum L β

theorem partitionFunction_openSpec (n : ℕ) (β : ℝ) :
    CertifiedGaps.partitionFunction openSpec n β = openZ β n := by
  rw [((SpinOneLargeD.model 0).partitionFunction_spectra (boundary (3 / 5))).2]; rfl

theorem openZ_pos (n : ℕ) (β : ℝ) : 0 < openZ β n := by
  rw [← partitionFunction_openSpec]; exact CertifiedGaps.partitionFunction_pos openSpec n β

/-- The ring and open-chain traces satisfy the boundary moment inequalities at every `β > 0`. -/
theorem openMomentBounds {β : ℝ} (hβ : 0 < β) :
    OpenMomentBounds (fun L ↦ partitionFunction L β) (openZ β) := by
  have h := (SpinOneLargeD.model 0).openMomentBounds rfl (boundary (3 / 5)) hβ
  simp only [SpinOneLargeD.partitionFunction_zero] at h
  exact h

/-- Physical squaring for the ring: `Z_L(2β) ≤ Z_L(β)²`. -/
theorem partitionFunction_two_mul_le (L : ℕ) (β : ℝ) :
    partitionFunction L (2 * β) ≤ partitionFunction L β ^ 2 := by
  have h := CertifiedGaps.thermalPurity_le_one spectrum L β
  rwa [CertifiedGaps.thermalPurity, partitionFunction_spectrum,
    div_le_one (pow_pos (partitionFunction_pos' L β) 2)] at h

/-- Physical squaring for the open chains: `G_n(2β) ≤ G_n(β)²`. -/
theorem openZ_two_mul_le (n : ℕ) (β : ℝ) : openZ (2 * β) n ≤ openZ β n ^ 2 := by
  have h := CertifiedGaps.thermalPurity_le_one openSpec n β
  rwa [CertifiedGaps.thermalPurity, partitionFunction_openSpec, partitionFunction_openSpec,
    div_le_one (pow_pos (openZ_pos n β) 2)] at h

theorem shiftFactor_pos (β : ℝ) : 0 < shiftFactor β := exp_pos _

theorem shiftFactor_pow (β : ℝ) (n : ℕ) :
    shiftFactor β ^ n = exp (-(β * (700741 / 500000) * n)) := by
  rw [shiftFactor, ← exp_nat_mul]; ring_nf

/-- `e^{b₀α} ≤ 141.38`, from the Taylor enclosure with 25 terms. -/
theorem exp_le : exp (49 / 4 * (2021 / 5000)) ≤ 14138 / 100 := by
  refine (CertifiedGaps.Trace.exp_mem_Icc_taylor (x := 49 / 4 * (2021 / 5000)) (by norm_num)
    (N := 24) (by norm_num)).2.trans ?_
  norm_num [Finset.sum_range_succ, Nat.factorial]

/-- The ground weight is at least `1/141.38`: `1 ≤ 141.38 e^{-b₀α}`. -/
theorem one_le_mul_exp : 1 ≤ 14138 / 100 * exp (-(49 / 4 * (2021 / 5000))) := by
  have h := exp_le
  rw [exp_neg, ← div_eq_mul_inv, one_le_div (exp_pos _)]
  exact h

/-- **The ground weight** (from a trial bound): the shifted trace at length `2N` exceeds the
shifted Boltzmann factor `e^{-b₀α}` of one eigenvalue. -/
theorem ground_le {N : ℕ} (hN : N ≠ 0) (ht : TrialBelow (2 * N)) :
    exp (-(49 / 4 * (2021 / 5000))) ≤ shiftedOpenZ (49 / 4) (2 * N) := by
  obtain ⟨i, hi⟩ := ht
  have hspec : openSpec (2 * N) = oppositeFieldSpectrum (3 / 5) (2 * N) :=
    openSpectrum_even _ (even_two_mul N) (by omega)
  have hsum : exp (-(49 / 4 * oppositeFieldSpectrum (3 / 5) (2 * N) i)) ≤ openZ (49 / 4) (2 * N)
      := by
    rw [← partitionFunction_openSpec, CertifiedGaps.partitionFunction, hspec]
    exact Finset.single_le_sum (f := fun j ↦ exp (-(49 / 4 * oppositeFieldSpectrum (3 / 5)
      (2 * N) j))) (fun _ _ ↦ (exp_pos _).le) (Finset.mem_univ i)
  rw [shiftedOpenZ, shiftFactor_pow]
  calc exp (-(49 / 4 * (2021 / 5000)))
      ≤ exp (-(49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ))) *
          exp (-(49 / 4 * oppositeFieldSpectrum (3 / 5) (2 * N) i)) := by
        rw [← exp_add]; exact exp_le_exp.2 (by push_cast at hi ⊢; nlinarith)
    _ ≤ _ := mul_le_mul_of_nonneg_left hsum (exp_pos _).le

/-! ### The starting bounds at `b₀` (Lemma in:seedbounds) -/

/-- The shifted six-site boundary moment is the opposite-field trace. -/
theorem shiftedOpenZ_six_le (h6 : shiftFactor (49 / 4) ^ 6 *
    oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) ≤ 27 / 1000) :
    shiftedOpenZ (49 / 4) (2 * 3) ≤ 27 / 1000 := by
  rw [shiftedOpenZ, openZ, openPartitionFunction_even (3 / 5) 3]; exact h6

section Seed

variable {μ c ρ s : ℝ}

/-- **The boundary defect at `(34, b₀)`**, `e₃₄(0) ≤ .0915`, in the shifted normalization. -/
theorem shifted_defect34 (hw : WitnessData μ c ρ s)
    (h6 : shiftFactor (49 / 4) ^ 6 * oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) ≤
      27 / 1000) (ht : TrialBelow 34) :
    (1 - 915 / 10000) * shiftedOpenZ (49 / 4) 34 ≤ c * μ ^ 34 := by
  obtain ⟨w, hμl, hμU, hρ, -⟩ := hw
  have hμ : 0 < μ := by linarith
  have hg := ground_le (N := 17) (by norm_num) ht
  have h1 := one_le_mul_exp
  have hg0 := exp_pos (-(49 / 4 * (2021 / 5000) : ℝ))
  refine w.one_sub_mul_le_of_ground hμ (K := 3) (N := 17) (by norm_num) (by norm_num) hμU hρ
    (shiftedOpenZ_six_le h6) hg0 hg ?_
  -- `(.874 U)^{28} .027 ≤ .0915 / 141.38 ≤ .0915 e^{-b₀α}`.
  have hnum : ((874 / 1000 : ℝ) * (100139 / 100000)) ^ (2 * 17 - 2 * 3) * (27 / 1000) * 1 ≤
      915 / 10000 * (1 / (14138 / 100)) := by norm_num
  nlinarith

/-- **The physical purity from a ground weight** (Lemma in:seedbounds, second inequality): at a
length `2N ≥ 8` with a trial bound and `141.38 (.003745 U^{2N} / P_*(U)² + .027 (.874 U)^{2N-6})
≤ 1 + R`, the opposite-field purity at `2^J b₀` is at least `(1 + R^{2^J})⁻²`. -/
theorem thermalPurity_ge (hw : WitnessData μ c ρ s)
    (h6 : shiftFactor (49 / 4) ^ 6 * oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) ≤
      27 / 1000) {N : ℕ} (hN : 4 ≤ N) (ht : TrialBelow (2 * N)) {R : ℝ} (hR0 : 0 ≤ R)
    (hR : 14138 / 100 * (3745 / 10 ^ 6 * ((100139 / 100000) ^ (2 * N) /
      eval pStar (100139 / 100000) ^ 2) + ((874 / 1000) * (100139 / 100000)) ^ (2 * N - 2 * 3) *
        (27 / 1000)) ≤ 1 + R) (J : ℕ) :
    1 / (1 + R ^ 2 ^ J) ^ 2 ≤ CertifiedGaps.thermalPurity openSpec (2 * N) (2 ^ J * (49 / 4)) := by
  obtain ⟨w, hμl, hμU, hρ, hC⟩ := hw
  have hμ : 0 < μ := by linarith
  have hratio : RatioCheck pStar 4 (9 / 10) := by decide +kernel
  -- The shifted trace at length `2N` is at most `X`, the leading plus the residual part.
  have hX := w.le_add_of_bounds hμ (K := 3) (N := N) (by norm_num) (by omega) hratio
    (by norm_num) hN (by push_cast; linarith) hμU hρ (shiftedOpenZ_six_le h6)
    (C := 3745 / 10 ^ 6) (by rw [eval_pStar]; exact hC)
  set X := 3745 / 10 ^ 6 * ((100139 / 100000 : ℝ) ^ (2 * N) /
      eval pStar (100139 / 100000) ^ 2) + ((874 / 1000) * (100139 / 100000)) ^ (2 * N - 2 * 3) *
        (27 / 1000) with hXdef
  have hX0 : 0 ≤ X := (mul_nonneg (pow_nonneg (shiftFactor_pos _).le _)
    (openZ_pos _ _).le).trans hX
  set g := exp (-(49 / 4 * (2021 / 5000) : ℝ))
  have h1 := one_le_mul_exp
  -- `σ^{2N} G(2N) ≤ X ≤ g (1 + R)`.
  have hYg : shiftedOpenZ (49 / 4) (2 * N) ≤ g * (1 + R) :=
    calc _ ≤ X := hX
      _ ≤ 14138 / 100 * g * X := le_mul_of_one_le_left hX0 h1
      _ = g * (14138 / 100 * X) := by ring
      _ ≤ g * (1 + R) := mul_le_mul_of_nonneg_left hR (exp_pos _).le
  obtain ⟨i, hi⟩ := ht
  have hspec : openSpec (2 * N) = oppositeFieldSpectrum (3 / 5) (2 * N) :=
    openSpectrum_even _ (even_two_mul N) (by omega)
  refine thermalPurity_ge_of_ground openSpec (by norm_num) (i₀ := i)
    (e₀ := 2021 / 5000 - 700741 / 500000 * ((2 * N : ℕ) : ℝ)) (by rw [hspec]; linarith) hR0
    ?_ J
  -- Unshift: `e^{-b₀(α - a·2N)} = g e^{b₀ a 2N}` and `σ^{2N} e^{b₀ a 2N} = 1`.
  rw [partitionFunction_openSpec]
  rw [shiftedOpenZ, shiftFactor_pow] at hYg
  have he : exp (-(49 / 4 * (2021 / 5000 - 700741 / 500000 * ((2 * N : ℕ) : ℝ)))) =
      g * exp (49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ)) := by
    rw [← exp_add]; ring_nf
  have hσ : exp (-(49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ))) *
      exp (49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ)) = 1 := by
    rw [← exp_add]; simp
  rw [he]
  have hE := exp_pos (49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ))
  calc openZ (49 / 4) (2 * N)
      = exp (49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ)) *
          (exp (-(49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ))) * openZ (49 / 4) (2 * N)) := by
        rw [← mul_assoc, mul_comm (exp _) (exp _), hσ, one_mul]
    _ ≤ exp (49 / 4 * (700741 / 500000) * ((2 * N : ℕ) : ℝ)) * (g * (1 + R)) :=
        mul_le_mul_of_nonneg_left hYg hE.le
    _ = _ := by ring

end Seed

/-! ### The two boundary steps -/

/-- The unshifted witness at `b₀` from H1, with the boundary defect `e₃₄(0) ≤ .0915`. -/
theorem exists_witness_b₀ (h1 : ThermalData) (ht : TrialBelow 34) :
    ∃ μ c ρ s, OpenMomentWitness (fun L ↦ partitionFunction L (49 / 4)) (openZ (49 / 4))
      μ c ρ s ∧ (1 - 915 / 10000) * openZ (49 / 4) 34 ≤ c * μ ^ 34 := by
  obtain ⟨μ, c, ρ, s, hw⟩ := h1.witness
  have hd := shifted_defect34 hw h1.sixth ht
  have w := hw.1
  set σ := shiftFactor (49 / 4)
  have hσ : 0 < σ := shiftFactor_pos _
  have w' := w.scale (inv_pos.2 hσ)
  have hZ : (fun L ↦ σ⁻¹ ^ L * shiftedRingZ (49 / 4) L) = fun L ↦ partitionFunction L (49 / 4) :=
    funext fun L ↦ by
      rw [shiftedRingZ, ← mul_assoc, ← mul_pow, inv_mul_cancel₀ hσ.ne', one_pow, one_mul]
  have hG : (fun n ↦ σ⁻¹ ^ n * shiftedOpenZ (49 / 4) n) = openZ (49 / 4) :=
    funext fun n ↦ by
      rw [shiftedOpenZ, ← mul_assoc, ← mul_pow, inv_mul_cancel₀ hσ.ne', one_pow, one_mul]
  rw [hZ, hG] at w'
  refine ⟨_, c, ρ, s, w', ?_⟩
  -- Divide the shifted inequality by `σ^{34}`.
  rw [shiftedOpenZ] at hd
  have hσ34 : 0 < σ ^ 34 := pow_pos hσ _
  have e : c * μ ^ 34 = σ ^ 34 * (c * (σ⁻¹ * μ) ^ 34) := by
    rw [mul_pow, inv_pow]; field_simp
  rw [e] at hd
  have : σ ^ 34 * ((1 - 915 / 10000) * openZ (49 / 4) 34) ≤ σ ^ 34 * (c * (σ⁻¹ * μ) ^ 34) := by
    linarith
  exact le_of_mul_le_mul_left this hσ34

/-- `ρ₁⁴⁴ ≤ .125` at `b₁ = 49/2` from `P_s(72, 49/2) ≥ 1 - .0605003`. -/
theorem rho₁_pow_le {μ c ρ s : ℝ}
    (w : OpenMomentWitness (fun L ↦ partitionFunction L (49 / 2)) (openZ (49 / 2)) μ c ρ s)
    (hS : 1 - 605003 / 10 ^ 7 ≤ spatialPurity 72 (49 / 2)) : ρ ^ 44 ≤ 1250 / 10000 := by
  have hμ := w.mu_pos (partitionFunction_pos' 4 _)
  have ht := w.tail_le hμ (N := 36) (by norm_num) (partitionFunction_pos' _ _)
    (ε := 605003 / 10 ^ 7) (by norm_num) hS
  exact pow_le_of_pow_le_gcd w.rho_nonneg (by norm_num) (a := 44) (by norm_num) ht
    (by norm_num [tailBound])

/-- `ρ₂⁶⁴ ≤ .1491` at `b₂ = 49` from `P_s(120, 49) ≥ .948`. -/
theorem rho₂_pow_le {μ c ρ s : ℝ}
    (w : OpenMomentWitness (fun L ↦ partitionFunction L 49) (openZ 49) μ c ρ s)
    (hS : 948 / 1000 ≤ spatialPurity 120 49) : ρ ^ 64 ≤ 1491 / 10000 := by
  have hμ := w.mu_pos (partitionFunction_pos' 4 _)
  have ht := w.tail_le hμ (N := 60) (by norm_num) (partitionFunction_pos' _ _)
    (ε := 52 / 1000) (by norm_num) (by norm_num at hS ⊢; exact hS)
  exact pow_le_of_pow_le_gcd w.rho_nonneg (by norm_num) (a := 64) (by norm_num) ht
    (by norm_num [tailBound])

/-- `μ' ≤ μ²` for witnesses at `β` and `2β`, from physical squaring of the ring traces. -/
theorem mu_le_sq' {β μ c ρ s μ' c' ρ' s' : ℝ}
    (w : OpenMomentWitness (fun L ↦ partitionFunction L β) (openZ β) μ c ρ s)
    (w' : OpenMomentWitness (fun L ↦ partitionFunction L (2 * β)) (openZ (2 * β)) μ' c' ρ' s') :
    μ' ≤ μ ^ 2 :=
  w.mu_le_sq w' (w.mu_pos (partitionFunction_pos' 4 _))
    fun N _ ↦ partitionFunction_two_mul_le (2 * N) β

/-- **The boundary defect `e₁₄₂ ≤ .14` at `b₂ = 49`** for the witness at `b₂`, from H1, H2 and
the two ring purities: the steps (in:78) and (in:142) of the paper, each one `move_one`. -/
theorem boundary_step (h1 : ThermalData) (h2 : TrialData)
    (hS72 : 1 - 605003 / 10 ^ 7 ≤ spatialPurity 72 (49 / 2))
    (hS120 : 948 / 1000 ≤ spatialPurity 120 49) {μ₂ c₂ ρ₂ s₂ : ℝ}
    (w₂ : OpenMomentWitness (fun L ↦ partitionFunction L 49) (openZ 49) μ₂ c₂ ρ₂ s₂) :
    (1 - 14 / 100) * openZ 49 142 ≤ c₂ * μ₂ ^ 142 := by
  obtain ⟨μ, c, ρ, s, w⟩ := h1.witness
  obtain ⟨μ₀, c₀, ρ₀, s₀, w₀, he₀⟩ := exists_witness_b₀ h1 h2.trial34
  obtain ⟨μ₁, c₁, ρ₁, s₁, w₁⟩ := openMomentBounds (β := 49 / 2) (by norm_num)
  have hb₁ : (49 / 2 : ℝ) = 2 * (49 / 4) := by norm_num
  have hb₂ : (49 : ℝ) = 2 * (49 / 2) := by norm_num
  -- The purities `A(78, 0) ≥ (1 + .2632)⁻²` and `A(142, 1) ≥ (1 + .3804²)⁻²`.
  have hA78 := thermalPurity_ge w h1.sixth (N := 39)
    (by norm_num) h2.trial78 (R := 2632 / 10000) (by norm_num)
    (by rw [eval_pStar]; norm_num) 0
  have hA142 := thermalPurity_ge w h1.sixth (N := 71)
    (by norm_num) h2.trial142 (R := 3804 / 10000) (by norm_num)
    (by rw [eval_pStar]; norm_num) 1
  rw [CertifiedGaps.thermalPurity, partitionFunction_openSpec, partitionFunction_openSpec,
    le_div_iff₀ (pow_pos (openZ_pos _ _) 2)] at hA78 hA142
  norm_num only [pow_zero, one_mul, pow_one] at hA78 hA142
  -- Step (in:78): `b₀ → b₁`, `(K, N) = (17, 39)`.
  have w₁' : OpenMomentWitness (fun L ↦ partitionFunction L (2 * (49 / 4))) (openZ (2 * (49 / 4)))
      μ₁ c₁ ρ₁ s₁ := by rw [← hb₁]; exact w₁
  have hμ₁ := w₁.mu_pos (partitionFunction_pos' 4 _)
  have step₁ := w₀.move_one w₁' hμ₁ (mu_le_sq' w₀ w₁') (K := 17) (N := 39) (by norm_num)
    (by norm_num) (openZ_two_mul_le 34 _) (t := 1 / (1 + 2632 / 10000) ^ 2) (by norm_num)
    (by rw [← hb₁]; norm_num at hA78 ⊢; linarith) (E := 915 / 10000) (by norm_num) he₀
    (q := 1250 / 10000) (by norm_num) (by norm_num) (rho₁_pow_le w₁ hS72)
  rw [← hb₁] at step₁
  have hF₁ : overlapUpdate (1250 / 10000) (915 / 10000) (1 / (1 + 2632 / 10000) ^ 2) ≤
      (1335 / 10000 : ℝ) := by norm_num [overlapUpdate]
  have he₁ : (1 - 1335 / 10000) * openZ (49 / 2) (2 * 39) ≤ c₁ * μ₁ ^ (2 * 39) :=
    le_trans (mul_le_mul_of_nonneg_right (by linarith) (openZ_pos _ _).le) step₁
  -- Step (in:142): `b₁ → b₂`, `(K, N) = (39, 71)`.
  have w₂' : OpenMomentWitness (fun L ↦ partitionFunction L (2 * (49 / 2))) (openZ (2 * (49 / 2)))
      μ₂ c₂ ρ₂ s₂ := by rw [← hb₂]; exact w₂
  have hμ₂ := w₂.mu_pos (partitionFunction_pos' 4 _)
  have step₂ := w₁.move_one w₂' hμ₂ (mu_le_sq' w₁ w₂') (K := 39) (N := 71) (by norm_num)
    (by norm_num) (openZ_two_mul_le 78 _) (t := 1 / (1 + (3804 / 10000) ^ 2) ^ 2) (by norm_num)
    (by rw [← hb₂]; norm_num at hA142 ⊢; linarith) (E := 1335 / 10000) (by norm_num) he₁
    (q := 1491 / 10000) (by norm_num) (by norm_num) (rho₂_pow_le w₂ hS120)
  rw [← hb₂] at step₂
  have hF₂ : overlapUpdate (1491 / 10000) (1335 / 10000) (1 / (1 + (3804 / 10000) ^ 2) ^ 2) ≤
      (14 / 100 : ℝ) := by norm_num [overlapUpdate]
  exact le_trans (mul_le_mul_of_nonneg_right (by linarith) (openZ_pos _ _).le) step₂

/-! ### The two boundary inputs -/

/-- **The boundary seed from short-chain data** (the statement of `seed_boundary`): under H1,
H2, `P_s(72, 49/2) ≥ 1 - .0605003` and `P_s(120, 49) ≥ .948`, the boundary defect of the
opposite-field chain at `(142, 49)` is at most `.14`. -/
theorem boundaryDefect_le_of_data (h1 : ThermalData) (h2 : TrialData)
    (hS72 : 1 - 605003 / 10 ^ 7 ≤ spatialPurity 72 (49 / 2))
    (hS120 : 948 / 1000 ≤ spatialPurity 120 49) :
    boundaryDefect (fun L ↦ partitionFunction L 49)
      (fun n ↦ oppositeFieldPartitionFunction (3 / 5) n 49) 142 ≤ 14 / 100 := by
  obtain ⟨μ₂, c₂, ρ₂, s₂, w₂⟩ := openMomentBounds (β := 49) (by norm_num)
  have hμ₂ := w₂.mu_pos (partitionFunction_pos' 4 _)
  have hρ₂ := w₂.rho_lt_one hμ₂ (N := 60) (by norm_num) (partitionFunction_pos' _ _)
    (ε := 52 / 1000) (by norm_num) (by norm_num at hS120 ⊢; exact hS120)
  have h := w₂.boundaryDefect_le (partitionFunction_pos' 4 _) hρ₂ (openZ_pos 142 49)
    (boundary_step h1 h2 hS72 hS120 w₂)
  rw [show (142 : ℕ) = 2 * 71 from rfl, boundaryDefect_congr
    (G' := fun n ↦ oppositeFieldPartitionFunction (3 / 5) n 49)
    (fun N ↦ openPartitionFunction_even (3 / 5) N 49)] at h
  exact h

/-- **The open-chain thermal seed from short-chain data** (the statement of
`seed_open_thermal`): under H1 and H2, `A(240, 2) ≥ (1 + .5817⁴)⁻² > .79`. -/
theorem oppositeFieldThermalPurity_ge_of_data (h1 : ThermalData) (h2 : TrialData) :
    79 / 100 ≤ oppositeFieldThermalPurity (3 / 5) 240 49 := by
  obtain ⟨μ, c, ρ, s, w⟩ := h1.witness
  have hA := thermalPurity_ge w h1.sixth (N := 120)
    (by norm_num) h2.trial240 (R := 5817 / 10000) (by norm_num) (by rw [eval_pStar]; norm_num) 2
  rw [CertifiedGaps.thermalPurity, partitionFunction_openSpec, partitionFunction_openSpec] at hA
  rw [oppositeFieldThermalPurity, ← openPartitionFunction_even (3 / 5) 120,
    ← openPartitionFunction_even (3 / 5) 120]
  have e : (2 : ℝ) ^ 2 * (49 / 4) = 49 := by norm_num
  rw [e] at hA
  refine le_trans ?_ hA
  norm_num

end BoundaryInit

end CertifiedGaps.SpinOneHeisenberg
