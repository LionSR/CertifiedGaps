/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Certificate
import CertifiedGaps.Seed.Lift

/-!
# Seeds from a lifted trace table

This file proves the initialization at length sixty of OpenAI's *The periodic spin-one Haldane
gap* (Proposition 5.4) for any ring with the three twists `1, P, C` of `SectorScheme.ion`
(sectors `I, O, N` with multiplicities `1, 2, 3`). The table of the twisted traces is computed at
an inverse temperature `β₀` and the seed is needed at the larger `β = (l/j) β₀`; the argument
has four steps.

1. **Masses at `β₀`.** The table at short lengths bounds the `q`th sector moments (`mass`) and,
   through polynomial filters, caps `|λ|` in the capped sectors. The capped-mass bounds of
   `SectorSeed.sectorBound` then bound the sector moments of order `q r'`, and so the traces at
   length `q r'` (`traceBound`, `twistBound`; `LiftedSeed.trace_le_of_table`).
2. **Temperature lifting.** The traces at length `k = q r` are Boltzmann sums, so
   `Z_k(β, g)^j ≤ Z_k(β₀, g)^l` (`BoltzmannRep.lift`): `Z_k(β) < D` and `Z_k(β, P), Z_k(β, C) < E`.
   By the sector inversion, the `I` and `I + 2O` parts of the `k`th moment are at most
   `(D + 11E)/12 ≤ τ` and `(D + 3E)/4 ≤ τ₀`.
3. **Spatial purity at `β`.** The weights `w = λ^k` have total mass at most `D` and fourth moment
   `Z_{4k}(β) > 1` (a trial state). One weight `w₀ > t₀` exists, and it lies in `I` because
   `2 t₀ ≥ D` (`exists_dominant_ion`). The other weights have squared sum at most
   `R(t) = max_{vertices} x² + y²/2 + (D - t - x - y)²/3` over the polygon of remaining masses
   (`restSum_two_le`), which bootstraps the threshold `t₀ < t₁ < ⋯` (`BootChain`) and gives
   `P_s(2k, β) ≥ (1 + R(t)/t²)⁻²` (`ion_spatial`).
4. **Thermal purity at `β`.** One physical Boltzmann factor at length `4k` exceeds `1` (the trial
   state), and their sum at `β₀` is at most `traceBound (4r) = q'`; if `(q' - 1)^l < R_q^j`, then
   `P_τ(4k, β) ≥ (1 + R_q)⁻²` (`boltzmann_sq_le_of_dominant`).

`LiftedSeed` collects the rational data; `LiftedSeed.Valid` is decidable, so `decide` checks a
concrete certificate. `LiftedSeed.purities` gives both purities from the sector representation at
`β₀` and `β`, the table (strictly), the Boltzmann representations, and one eigenvalue below
`-κ · 4k`.

## Main results

* `polyQ_le_max`: a convex quadratic on the polygon of remaining masses is at most its value at a
  vertex.
* `ion_spatial`: the spatial purity of three-sector weights with multiplicities `1, 2, 3`.
* `LiftedSeed.purities`: the two seed purities from the certificate.
-/

open Finset Real

namespace CertifiedGaps

open PolynomialFilter

/-! ### Convex quadratics on the polygon of remaining masses -/

/-- The squared-weight bound `x² + y²/2 + (D - t - x - y)²/3` for remaining masses `x` (in `I`),
`y` (the `O` pair) and `D - t - x - y` (the `N` triple). -/
def polyQ {K : Type*} [Field K] (D t x y : K) : K := x ^ 2 + y ^ 2 / 2 + (D - t - x - y) ^ 2 / 3

/-- The maximum `R(t)` of `polyQ` over the four vertices `(0, 0)`, `(τ - t, 0)`,
`(τ - t, τ₀ - τ)`, `(0, τ₀ - t)` of the polygon `0 ≤ x ≤ τ - t`, `0 ≤ y`, `x + y ≤ τ₀ - t`. -/
def restBound {K : Type*} [Field K] [LinearOrder K] (D τ τ₀ t : K) : K :=
  max (max (polyQ D t 0 0) (polyQ D t (τ - t) 0))
    (max (polyQ D t (τ - t) (τ₀ - τ)) (polyQ D t 0 (τ₀ - t)))

theorem restBound_cast (D τ τ₀ t : ℚ) :
    ((restBound D τ τ₀ t : ℚ) : ℝ) = restBound (D : ℝ) τ τ₀ t := by
  simp only [restBound, polyQ, Rat.cast_max]; push_cast; rfl

/-- A function below its chord at `s ∈ [l, u]` is at most its larger endpoint value. -/
theorem le_max_of_chord {φ : ℝ → ℝ} {l s u : ℝ} (hl : l ≤ s) (hu : s ≤ u)
    (h : (u - l) * φ s ≤ (u - s) * φ l + (s - l) * φ u) : φ s ≤ max (φ l) (φ u) := by
  rcases hl.eq_or_lt with rfl | hls
  · exact le_max_left _ _
  · have hul : 0 < u - l := by linarith
    have a := le_max_left (φ l) (φ u)
    have b := le_max_right (φ l) (φ u)
    have : (u - l) * φ s ≤ (u - l) * max (φ l) (φ u) := by nlinarith
    exact le_of_mul_le_mul_left this hul

/-- The chord inequality of a quadratic with leading coefficient `α ≥ 0`: the gap is
`α (u - l)(u - s)(s - l)`. -/
theorem chord_of_quadratic {φ : ℝ → ℝ} {α β γ : ℝ} (hφ : ∀ s, φ s = α * s ^ 2 + β * s + γ)
    (hα : 0 ≤ α) {l s u : ℝ} (hl : l ≤ s) (hu : s ≤ u) :
    (u - l) * φ s ≤ (u - s) * φ l + (s - l) * φ u := by
  have e : (u - s) * φ l + (s - l) * φ u - (u - l) * φ s = α * ((u - l) * (u - s) * (s - l)) := by
    rw [hφ, hφ, hφ]; ring
  have : 0 ≤ α * ((u - l) * (u - s) * (s - l)) :=
    mul_nonneg hα (mul_nonneg (mul_nonneg (by linarith) (by linarith)) (by linarith))
  linarith

/-- **The polygon bound.** On `0 ≤ x ≤ τ - t`, `0 ≤ y`, `x + y ≤ τ₀ - t`, the convex quadratic
`polyQ` is at most its largest vertex value `restBound`: it is convex in `y` at fixed `x`, and
convex in `x` along the two edges `y = 0` and `x + y = τ₀ - t`. -/
theorem polyQ_le_max {D τ τ₀ t x y : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ τ - t) (hy0 : 0 ≤ y)
    (hxy : x + y ≤ τ₀ - t) : polyQ D t x y ≤ restBound D τ τ₀ t := by
  have h1 : polyQ D t x y ≤ max (polyQ D t x 0) (polyQ D t x (τ₀ - t - x)) :=
    le_max_of_chord (φ := fun y ↦ polyQ D t x y) hy0 (by linarith)
      (chord_of_quadratic (φ := fun y ↦ polyQ D t x y) (α := 5 / 6) (β := -(2 / 3) * (D - t - x))
        (γ := x ^ 2 + (D - t - x) ^ 2 / 3)
        (fun s ↦ by simp only [polyQ]; ring) (by norm_num) hy0 (by linarith))
  have h2 : polyQ D t x 0 ≤ max (polyQ D t 0 0) (polyQ D t (τ - t) 0) :=
    le_max_of_chord (φ := fun x ↦ polyQ D t x 0) hx0 hx
      (chord_of_quadratic (φ := fun x ↦ polyQ D t x 0) (α := 4 / 3) (β := -(2 / 3) * (D - t))
        (γ := (D - t) ^ 2 / 3)
        (fun s ↦ by simp only [polyQ]; ring) (by norm_num) hx0 hx)
  have h3 : polyQ D t x (τ₀ - t - x) ≤
      max (polyQ D t 0 (τ₀ - t - 0)) (polyQ D t (τ - t) (τ₀ - t - (τ - t))) :=
    le_max_of_chord (φ := fun x ↦ polyQ D t x (τ₀ - t - x)) hx0 hx
      (chord_of_quadratic (φ := fun x ↦ polyQ D t x (τ₀ - t - x)) (α := 3 / 2) (β := -(τ₀ - t))
        (γ := (τ₀ - t) ^ 2 / 2 + (D - τ₀) ^ 2 / 3) (fun s ↦ by simp only [polyQ]; ring)
        (by norm_num) hx0 hx)
  rw [sub_zero, show τ₀ - t - (τ - t) = τ₀ - τ by ring] at h3
  refine h1.trans (max_le (h2.trans (le_max_left _ _)) (h3.trans ?_))
  rw [max_comm]
  exact le_max_right _ _

/-! ### Three sectors with multiplicities `1, 2, 3` -/

section Spatial

variable {ι : Type*} (sec : Fin 3 → Finset ι) (w : ι → ℝ)

/-- The power sum `∑_{i ∈ sec s} wᵢ^r` of one sector. -/
def psum (s : Fin 3) (r : ℕ) : ℝ := ∑ i ∈ sec s, w i ^ r

/-- The trace `I + 2O + 3N` of the `r`th powers, with the multiplicities of `ion`. -/
def ionSum (r : ℕ) : ℝ := psum sec w 0 r + 2 * psum sec w 1 r + 3 * psum sec w 2 r

/-- The untwisted `ion` trace at length `k r` is `ionSum` of the weights `x^k`. -/
theorem ion_trace_zero_mul (x : ι → ℝ) (k r : ℕ) :
    SectorScheme.ion.trace sec x (0 : Fin 3) (k * r) = ionSum sec (fun i ↦ x i ^ k) r := by
  simp only [SectorScheme.ion_trace_zero_fin3, ionSum, psum, pow_mul]

/-- The weights other than `w i₀`: `∑_{sec 0 ∖ i₀} wᵢ^r + 2 P₁(r) + 3 P₂(r)`. -/
def restSum [DecidableEq ι] (i₀ : ι) (r : ℕ) : ℝ :=
  ∑ i ∈ (sec 0).erase i₀, w i ^ r + 2 * psum sec w 1 r + 3 * psum sec w 2 r

variable {sec w}

theorem ionSum_eq [DecidableEq ι] {i₀ : ι} (hi₀ : i₀ ∈ sec 0) (r : ℕ) :
    ionSum sec w r = w i₀ ^ r + restSum sec w i₀ r := by
  simp only [ionSum, restSum, psum, ← add_sum_erase (sec 0) (fun i ↦ w i ^ r) hi₀]
  ring

theorem psum_nonneg (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) (s : Fin 3) (r : ℕ) :
    0 ≤ psum sec w s r :=
  sum_nonneg fun i hi ↦ pow_nonneg (hw s i hi) r

theorem restSum_nonneg [DecidableEq ι] (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) (i₀ : ι) (r : ℕ) :
    0 ≤ restSum sec w i₀ r := by
  have := sum_nonneg fun i (hi : i ∈ (sec 0).erase i₀) ↦ pow_nonneg (hw 0 i (mem_of_mem_erase hi)) r
  have := psum_nonneg hw 1 r
  have := psum_nonneg hw 2 r
  simp only [restSum]
  linarith

/-- `∑ wᵢ⁴ ≤ (∑ wᵢ²)²` for nonnegative weights. -/
theorem sum_four_le {κ : Type*} (s : Finset κ) (w : κ → ℝ) (hw : ∀ i ∈ s, 0 ≤ w i) :
    ∑ i ∈ s, w i ^ 4 ≤ (∑ i ∈ s, w i ^ 2) ^ 2 := by
  have := sum_pow_le_pow_sum s (fun i ↦ w i ^ 2) (fun i hi ↦ pow_nonneg (hw i hi) 2)
    (n := 2) (by norm_num)
  simpa only [← pow_mul] using this

/-- The fourth moment of the remaining weights is at most the square of their second moment. -/
theorem restSum_four_le [DecidableEq ι] (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) (i₀ : ι) :
    restSum sec w i₀ 4 ≤ restSum sec w i₀ 2 ^ 2 := by
  have h0 := sum_four_le ((sec 0).erase i₀) w fun i hi ↦ hw 0 i (mem_of_mem_erase hi)
  have h1 := sum_four_le (sec 1) w (hw 1)
  have h2 := sum_four_le (sec 2) w (hw 2)
  have a0 : 0 ≤ ∑ i ∈ (sec 0).erase i₀, w i ^ 2 :=
    sum_nonneg fun i hi ↦ pow_nonneg (hw 0 i (mem_of_mem_erase hi)) 2
  have a1 := psum_nonneg hw 1 2
  have a2 := psum_nonneg hw 2 2
  simp only [restSum, psum] at *
  nlinarith

/-- **The remaining squared weights.** If `w i₀ ≥ t` in sector `0`, the total mass is at most
`D`, the `I` mass at most `τ` and the `I + 2O` mass at most `τ₀`, then the squared weights other
than `w i₀` sum to at most `restBound D τ τ₀ t`. -/
theorem restSum_two_le [DecidableEq ι] (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) {i₀ : ι}
    (hi₀ : i₀ ∈ sec 0) {D τ τ₀ t : ℝ} (ht : t ≤ w i₀) (h1 : ionSum sec w 1 ≤ D)
    (hI : psum sec w 0 1 ≤ τ) (hJ : psum sec w 0 1 + 2 * psum sec w 1 1 ≤ τ₀) :
    restSum sec w i₀ 2 ≤ restBound D τ τ₀ t := by
  have hw0 : ∀ i ∈ (sec 0).erase i₀, 0 ≤ w i := fun i hi ↦ hw 0 i (mem_of_mem_erase hi)
  set x := ∑ i ∈ (sec 0).erase i₀, w i
  set y := 2 * ∑ i ∈ sec 1, w i
  set z := 3 * ∑ i ∈ sec 2, w i
  have hx0 : 0 ≤ x := sum_nonneg hw0
  have hy0 : 0 ≤ y := by have := sum_nonneg (hw 1); positivity
  have hz0 : 0 ≤ z := by have := sum_nonneg (hw 2); positivity
  have hP0 : psum sec w 0 1 = w i₀ + x := by
    simp only [psum, pow_one, x]; exact (add_sum_erase _ _ hi₀).symm
  have hP1 : 2 * psum sec w 1 1 = y := by simp only [psum, pow_one, y]
  have hP2 : 3 * psum sec w 2 1 = z := by simp only [psum, pow_one, z]
  have hsum : w i₀ + x + y + z ≤ D := by simp only [ionSum] at h1; linarith
  -- Squared weights are at most the squared masses, sector by sector.
  have s0 : ∑ i ∈ (sec 0).erase i₀, w i ^ 2 ≤ x ^ 2 := sum_pow_le_pow_sum _ w hw0 (by norm_num)
  have s1 : psum sec w 1 2 ≤ (y / 2) ^ 2 := by
    have := sum_pow_le_pow_sum _ w (hw 1) (n := 2) (by norm_num)
    simp only [psum, y]; convert this using 2; ring
  have s2 : psum sec w 2 2 ≤ (z / 3) ^ 2 := by
    have := sum_pow_le_pow_sum _ w (hw 2) (n := 2) (by norm_num)
    simp only [psum, z]; convert this using 2; ring
  have hzD : z ^ 2 ≤ (D - t - x - y) ^ 2 := pow_le_pow_left₀ hz0 (by linarith) 2
  have hQ := polyQ_le_max (D := D) (τ := τ) (τ₀ := τ₀) (t := t) hx0 (by linarith) hy0
    (by linarith)
  simp only [polyQ] at hQ
  simp only [restSum]
  nlinarith

/-- **One bootstrap step.** If `w i₀ ≥ t` and `t'⁴ + R(t)² < 1`, then the fourth moment `> 1`
forces `w i₀ > t'`. -/
theorem lt_of_boot (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) {i₀ : ι}
    (hi₀ : i₀ ∈ sec 0) {D τ τ₀ t t' : ℝ} (ht : t ≤ w i₀) (h1 : ionSum sec w 1 ≤ D)
    (hI : psum sec w 0 1 ≤ τ) (hJ : psum sec w 0 1 + 2 * psum sec w 1 1 ≤ τ₀)
    (h4 : 1 < ionSum sec w 4) (hboot : t' ^ 4 + restBound D τ τ₀ t ^ 2 < 1) : t' < w i₀ := by
  classical
  by_contra hc
  push Not at hc
  have hR := restSum_two_le hw hi₀ ht h1 hI hJ
  have hR0 := restSum_nonneg hw i₀ 2
  have h4' := restSum_four_le hw i₀
  have hsq : restSum sec w i₀ 2 ^ 2 ≤ restBound D τ τ₀ t ^ 2 := pow_le_pow_left₀ hR0 hR 2
  have hw4 : w i₀ ^ 4 ≤ t' ^ 4 := pow_le_pow_left₀ (hw 0 i₀ hi₀) hc 4
  rw [ionSum_eq hi₀] at h4
  linarith

/-- **The dominant weight.** If every weight is at most `t₀ ≥ D/2 ≥ 0` with `D t₀³ < 1`, the
fourth moment is below `1`; a weight above `t₀` cannot lie in the repeated sectors `1`, `2`
because `2 t₀ ≥ D`. So some weight of sector `0` exceeds `t₀`. -/
theorem exists_dominant_ion (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) {D t₀ : ℝ} (ht₀ : 0 ≤ t₀)
    (hD2 : D ≤ 2 * t₀) (hdom : D * t₀ ^ 3 < 1) (h1 : ionSum sec w 1 ≤ D)
    (h4 : 1 < ionSum sec w 4) : ∃ i₀ ∈ sec 0, t₀ < w i₀ := by
  by_contra hno
  push Not at hno
  have m0 := psum_nonneg hw 0 1
  have m1 := psum_nonneg hw 1 1
  have m2 := psum_nonneg hw 2 1
  have hall : ∀ s, ∀ i ∈ sec s, w i ≤ t₀ := by
    intro s i hi
    by_contra hc
    push Not at hc
    have hle : w i ≤ psum sec w s 1 := by
      simp only [psum, pow_one]
      exact single_le_sum (fun j hj ↦ hw s j hj) hi
    simp only [ionSum] at h1
    fin_cases s
    · exact absurd (hno i hi) (not_le.2 hc)
    · simp only [Fin.mk_one] at hle; linarith
    · simp only [Fin.reduceFinMk] at hle; linarith
  have h4s (s : Fin 3) : psum sec w s 4 ≤ t₀ ^ 3 * psum sec w s 1 := by
    have := sum_pow_succ_le (sec s) w (hw s) (hall s) 3
    simpa only [psum, pow_one] using this
  have := h4s 0
  have := h4s 1
  have := h4s 2
  have ht3 : 0 ≤ t₀ ^ 3 := pow_nonneg ht₀ 3
  have : ionSum sec w 4 ≤ t₀ ^ 3 * ionSum sec w 1 := by simp only [ionSum]; linarith
  have : t₀ ^ 3 * ionSum sec w 1 ≤ t₀ ^ 3 * D := mul_le_mul_of_nonneg_left h1 ht3
  linarith

/-- The bootstrap thresholds `t₀ < t₁ < ⋯`: each step needs `tₖ₊₁⁴ + R(tₖ)² < 1`. -/
def BootChain (D τ τ₀ : ℚ) : ℚ → List ℚ → Prop
  | _, [] => True
  | t, t' :: l => t' ^ 4 + restBound D τ τ₀ t ^ 2 < 1 ∧ BootChain D τ τ₀ t' l

instance decBootChain (D τ τ₀ : ℚ) : (t : ℚ) → (l : List ℚ) → Decidable (BootChain D τ τ₀ t l)
  | _, [] => isTrue trivial
  | _, t' :: l =>
    @instDecidableAnd _ _ (inferInstanceAs (Decidable (_ < _))) (decBootChain D τ τ₀ t' l)

/-- The last threshold of a chain. -/
def lastT : ℚ → List ℚ → ℚ
  | t, [] => t
  | _, t' :: l => lastT t' l

/-- The bootstrap reaches the last threshold. -/
theorem lt_lastT (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) {i₀ : ι} (hi₀ : i₀ ∈ sec 0)
    {D τ τ₀ : ℚ} (h1 : ionSum sec w 1 ≤ D) (hI : psum sec w 0 1 ≤ τ)
    (hJ : psum sec w 0 1 + 2 * psum sec w 1 1 ≤ τ₀) (h4 : 1 < ionSum sec w 4) :
    ∀ (l : List ℚ) (t : ℚ), BootChain D τ τ₀ t l → (t : ℝ) < w i₀ → (lastT t l : ℝ) < w i₀
  | [], _, _, h => h
  | t' :: l, t, hc, h => by
    have hb : (t' : ℝ) ^ 4 + restBound (D : ℝ) τ τ₀ t ^ 2 < 1 := by
      rw [← restBound_cast]; exact_mod_cast hc.1
    exact lt_lastT hw hi₀ h1 hI hJ h4 l t' hc.2 (lt_of_boot hw hi₀ h.le h1 hI hJ h4 hb)

/-- **The spatial purity of three sectors with multiplicities `1, 2, 3`.** Let `w ≥ 0` be
weights in sectors `I, O, N` with total mass `I + 2O + 3N ≤ D`, `I ≤ τ`, `I + 2O ≤ τ₀`, and
fourth moment `> 1`. If `0 ≤ t₀`, `D ≤ 2 t₀`, `D t₀³ < 1`, the bootstrap chain from `t₀` is valid
and ends at `T > 0`, and `(1 - p)(1 + R(T)/T²)² ≤ 1`, then
`(1 - p) (∑ μ w²)² ≤ ∑ μ w⁴`. -/
theorem ion_spatial (hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i) {D τ τ₀ t₀ p : ℚ} {steps : List ℚ}
    (ht₀ : 0 ≤ t₀) (hD2 : D ≤ 2 * t₀) (hdom : D * t₀ ^ 3 < 1)
    (hchain : BootChain D τ τ₀ t₀ steps) (hlast : 0 < lastT t₀ steps) (hp1 : p ≤ 1)
    (hp : (1 - p) * (1 + restBound D τ τ₀ (lastT t₀ steps) / lastT t₀ steps ^ 2) ^ 2 ≤ 1)
    (h1 : ionSum sec w 1 ≤ D) (hI : psum sec w 0 1 ≤ τ)
    (hJ : psum sec w 0 1 + 2 * psum sec w 1 1 ≤ τ₀) (h4 : 1 < ionSum sec w 4) :
    (1 - p : ℝ) * ionSum sec w 2 ^ 2 ≤ ionSum sec w 4 := by
  classical
  obtain ⟨i₀, hi₀, hdom₀⟩ := exists_dominant_ion hw (D := D) (t₀ := t₀) (by exact_mod_cast ht₀)
    (by exact_mod_cast hD2) (by exact_mod_cast hdom) h1 h4
  set T := lastT t₀ steps
  have hT : (T : ℝ) < w i₀ := lt_lastT hw hi₀ h1 hI hJ h4 steps t₀ hchain hdom₀
  have hT0 : (0 : ℝ) < T := by exact_mod_cast hlast
  have hp' : (1 - (p : ℝ)) * (1 + restBound (D : ℝ) τ τ₀ T / (T : ℝ) ^ 2) ^ 2 ≤ 1 := by
    rw [← restBound_cast]; exact_mod_cast hp
  have hR := restSum_two_le hw hi₀ hT.le h1 hI hJ (D := D)
  set R : ℝ := restBound (D : ℝ) τ τ₀ T
  have hR0 := restSum_nonneg hw i₀ 2
  have hp0 : 0 ≤ 1 - (p : ℝ) := by have : (p : ℝ) ≤ 1 := by exact_mod_cast hp1
                                   linarith
  set a := w i₀ ^ 2
  have hTa : (T : ℝ) ^ 2 ≤ a := pow_le_pow_left₀ hT0.le hT.le 2
  have hT2 : (0 : ℝ) < (T : ℝ) ^ 2 := by positivity
  -- `R ≤ a R / T²`, so `∑ μ w² ≤ a (1 + R/T²)`.
  have hRa : R ≤ a * (R / (T : ℝ) ^ 2) := by
    rw [mul_div_assoc', le_div_iff₀ hT2, mul_comm a]
    exact mul_le_mul_of_nonneg_left hTa (hR0.trans hR)
  have h2 : ionSum sec w 2 ≤ a * (1 + R / (T : ℝ) ^ 2) := by
    rw [ionSum_eq hi₀]; nlinarith
  have h20 : 0 ≤ ionSum sec w 2 := by
    rw [ionSum_eq hi₀]; have := sq_nonneg (w i₀); linarith
  have h4' : a ^ 2 ≤ ionSum sec w 4 := by
    rw [ionSum_eq hi₀, show w i₀ ^ 4 = a ^ 2 by ring]
    linarith [restSum_nonneg hw i₀ 4]
  calc (1 - (p : ℝ)) * ionSum sec w 2 ^ 2 ≤ (1 - p) * (a * (1 + R / (T : ℝ) ^ 2)) ^ 2 := by gcongr
    _ = a ^ 2 * ((1 - p) * (1 + R / (T : ℝ) ^ 2) ^ 2) := by ring
    _ ≤ a ^ 2 * 1 := by gcongr
    _ ≤ ionSum sec w 4 := by linarith

end Spatial


/-! ### Certificates -/

/-- The rational data of a seed from a table at `β₀`, lifted to `β = (l/j) β₀`, for a ring with
the twists `1, P, C` of `SectorScheme.ion`. The seed lengths are `n₀ = 2k` and `2n₀ = 4k` with
`k = q r`. -/
structure LiftedSeed where
  /-- The table of bounds on the shifted twisted traces at `β₀`. -/
  table : TraceBounds 2
  /-- The filter `X^(2a) p²` of each sector, as `(a, p)`, if the sector is capped. -/
  filter : Fin 3 → Option (ℕ × List ℚ)
  /-- The cap on `|λ|` in each sector, if a filter provides one. -/
  cap : Fin 3 → Option ℚ
  /-- The order of the tabulated sector masses. -/
  q : ℕ
  /-- The bound on the `q`th moment of each sector at `β₀`. -/
  mass : Fin 3 → ℚ
  /-- The lifted length is `k = q r`. -/
  r : ℕ
  /-- The lifting `β = (l/j) β₀`: denominator. -/
  j : ℕ
  /-- The lifting `β = (l/j) β₀`: numerator. -/
  l : ℕ
  /-- The bound on `Z_k(β)`. -/
  D : ℚ
  /-- The bound on the twisted traces `Z_k(β, P)`, `Z_k(β, C)`. -/
  E : ℚ
  /-- The bound on the `I` part of the `k`th moment at `β`. -/
  tau : ℚ
  /-- The bound on the `I + 2O` part of the `k`th moment at `β`. -/
  tau0 : ℚ
  /-- The first dominance threshold. -/
  t₀ : ℚ
  /-- The bootstrapped thresholds. -/
  steps : List ℚ
  /-- The bound on the non-dominant Boltzmann mass at length `4k` and inverse temperature `β`. -/
  Rq : ℚ
  /-- The spatial defect `1 - P_s(2k, β)` to certify. -/
  spatialDefect : ℚ
  /-- The thermal defect `1 - P_τ(4k, β)` to certify. -/
  thermalDefect : ℚ
  /-- The inverse temperature `β₀` of the table. -/
  beta : ℚ
  /-- The energy shift `κ` per bond of the table. -/
  kappa : ℚ

namespace LiftedSeed

variable (C : LiftedSeed)

/-- The masses and caps as a `SectorSeed`, to reuse `SectorSeed.sectorBound`. -/
def masses : SectorSeed 2 where
  q := C.q
  a := C.r
  mult := ![1, 2, 3]
  cap := C.cap
  mass := C.mass
  ell := 0
  spatialDefect := 0

/-- The table and filters as a `TableCertificate`, to reuse `TableCertificate.FilterCheck`. -/
def low : TableCertificate where
  scheme := SectorScheme.ion
  table := C.table
  filter := C.filter
  seed := C.masses
  beta := C.beta
  kappa := C.kappa
  thermalDefect := 0

/-- The lifted length `k = q r`. -/
def k : ℕ := C.q * C.r

/-- The capped-mass bound on `Z_{q r'}(β₀) = I + 2O + 3N`. -/
def traceBound (r' : ℕ) : ℚ :=
  C.masses.sectorBound 0 r' + 2 * C.masses.sectorBound 1 r' + 3 * C.masses.sectorBound 2 r'

/-- The capped-mass bound `I + 2O` on `Z_{q r'}(β₀, P)` and `Z_{q r'}(β₀, C)`. -/
def twistBound (r' : ℕ) : ℚ := C.masses.sectorBound 0 r' + 2 * C.masses.sectorBound 1 r'

/-- The lifted inverse temperature `β = (l/j) β₀`. -/
def liftBeta : ℚ := C.l / C.j * C.beta

/-- The last bootstrap threshold. -/
def last : ℚ := lastT C.t₀ C.steps

/-- The lengths of the sector representation at `β₀`. -/
def lengths₀ : Finset ℕ := C.table.lengths.toFinset ∪ {C.k, 4 * C.k}

/-- The lengths of the sector representation at `β`. -/
def lengths₁ : Finset ℕ := {C.k, 2 * C.k, 4 * C.k}

/-- The shape conditions. -/
def ShapeValid : Prop :=
  Even C.q ∧ 1 ≤ C.r ∧ (∀ s, 0 ≤ C.mass s ∧ SectorSeed.OneBar C.q (C.cap s) (C.mass s)) ∧
    (∀ n ∈ C.table.lengths, 2 ≤ n) ∧ 2 ≤ C.k ∧ 0 < C.beta ∧ 1 ≤ C.j ∧ C.j ≤ C.l

/-- The lifting conditions: `m^l < D^j`, `p^l < E^j`, and the sector inversion at `β`. -/
def LiftValid : Prop :=
  0 ≤ C.D ∧ 0 ≤ C.E ∧ C.traceBound C.r ^ C.l < C.D ^ C.j ∧ C.twistBound C.r ^ C.l < C.E ^ C.j ∧
    (C.D + 11 * C.E) / 12 ≤ C.tau ∧ (C.D + 3 * C.E) / 4 ≤ C.tau0

/-- The spatial conditions: dominance, the bootstrap and the spatial defect. -/
def SpatialValid : Prop :=
  0 ≤ C.t₀ ∧ C.D ≤ 2 * C.t₀ ∧ C.D * C.t₀ ^ 3 < 1 ∧ BootChain C.D C.tau C.tau0 C.t₀ C.steps ∧
    0 < C.last ∧ C.spatialDefect ≤ 1 ∧
    (1 - C.spatialDefect) * (1 + restBound C.D C.tau C.tau0 C.last / C.last ^ 2) ^ 2 ≤ 1

/-- The thermal conditions: `(q' - 1)^l < R_q^j` and the thermal defect. -/
def ThermalValid : Prop :=
  0 ≤ C.Rq ∧ (C.traceBound (4 * C.r) - 1) ^ C.l < C.Rq ^ C.j ∧ C.thermalDefect ≤ 1 ∧
    (1 - C.thermalDefect) * (1 + C.Rq) ^ 2 ≤ 1

/-- The rational conditions of a lifted seed. -/
def Valid : Prop :=
  (∀ s : Fin (C.low.scheme.m + 1), C.low.FilterCheck s) ∧
    (∀ s : Fin 3, OptLE (C.table.hiS SectorScheme.ion s C.q) (C.mass s)) ∧
    C.ShapeValid ∧ C.LiftValid ∧ C.SpatialValid ∧ C.ThermalValid

instance : Decidable C.ShapeValid := by unfold ShapeValid; infer_instance
instance : Decidable C.LiftValid := by unfold LiftValid; infer_instance
instance : Decidable C.SpatialValid := by unfold SpatialValid; infer_instance
instance : Decidable C.ThermalValid := by unfold ThermalValid; infer_instance
instance : Decidable (∀ s : Fin (C.low.scheme.m + 1), C.low.FilterCheck s) := inferInstance
instance : Decidable (∀ s : Fin 3, OptLE (C.table.hiS SectorScheme.ion s C.q) (C.mass s)) :=
  inferInstance
instance : Decidable C.Valid := by unfold Valid; infer_instance

variable {C}

/-- **Step 1: the capped-mass bounds at `β₀`.** Sector data that obey the table have traces at
length `q r'` at most `traceBound r'` (untwisted) and `twistBound r'` (twisted by `P` or `C`). -/
theorem trace_le_of_table (hC : C.Valid) {ι : Type*} {sec : Fin 3 → Finset ι} {x : ι → ℝ}
    (h : C.table.Holds (SectorScheme.ion.trace sec x)) {r' : ℕ} (hr : 1 ≤ r') :
    SectorScheme.ion.trace sec x (0 : Fin 3) (C.q * r') ≤ C.traceBound r' ∧
      SectorScheme.ion.trace sec x (1 : Fin 3) (C.q * r') ≤ C.twistBound r' ∧
      SectorScheme.ion.trace sec x (2 : Fin 3) (C.q * r') ≤ C.twistBound r' := by
  obtain ⟨hfilter, hmass, ⟨hq, -, hbar, -⟩, -⟩ := hC
  have hM := TraceBounds.momentBounds SectorScheme.ion_valid h
  -- Caps from the filters.
  have hcap : ∀ s v, C.cap s = some v → ∀ i ∈ sec s, |x i| < v := by
    intro s v hv i hi
    have hf := hfilter s
    unfold TableCertificate.FilterCheck at hf
    have hv' : C.low.seed.cap s = some v := hv
    rw [hv'] at hf
    rcases hfs : C.filter s with _ | ⟨a, qf⟩
    · have hfs' : C.low.filter s = none := hfs
      rw [hfs'] at hf; exact hf.elim
    · have hfs' : C.low.filter s = some (a, qf) := hfs
      rw [hfs'] at hf
      obtain ⟨B, hB, hcapB⟩ := CapCheckOpt.exists hf
      exact abs_lt_of_filter_bounds (sec s) x (eval_filterPoly_nonneg _ _) (hM s) hB hcapB hi
  -- Masses from the table.
  have hmass' : ∀ s, ∑ i ∈ sec s, x i ^ C.q ≤ C.mass s := by
    intro s
    obtain ⟨b, hb, hbm⟩ := (hmass s).exists
    exact ((hM s C.q b).2 hb).trans (by exact_mod_cast hbm)
  have hsec (s : Fin 3) : ∑ i ∈ sec s, x i ^ (C.q * r') ≤ C.masses.sectorBound s r' :=
    SectorSeed.sum_le_sectorBound (S := C.masses) hq (hbar s).1 (hbar s).2 (hcap s) (hmass' s) hr
  have hnn (s : Fin 3) : 0 ≤ ∑ i ∈ sec s, x i ^ (C.q * r') :=
    sum_nonneg fun i _ ↦ (hq.mul_right r').pow_nonneg _
  rw [SectorScheme.ion_trace_zero_fin3, SectorScheme.ion_trace_one_fin3,
    SectorScheme.ion_trace_two_fin3]
  simp only [traceBound, twistBound]
  push_cast
  have := hsec 0; have := hsec 1; have := hsec 2
  have := hnn 0; have := hnn 1; have := hnn 2
  exact ⟨by linarith, by linarith, by linarith⟩

/-- **Step 1 for the traces.** If the traces `Z g n` have the sector representation and obey the
table strictly, then at every length `q r'` of `lengths₀` they obey the capped-mass bounds. -/
theorem le_of_table (hC : C.Valid) {Z : Fin 3 → ℕ → ℝ}
    (happ : SectorApprox SectorScheme.ion Z C.lengths₀) (htab : C.table.HoldsStrict Z)
    {r' : ℕ} (hr : 1 ≤ r') (hmem : C.q * r' ∈ C.lengths₀) :
    Z 0 (C.q * r') ≤ C.traceBound r' ∧ Z 1 (C.q * r') ≤ C.twistBound r' ∧
      Z 2 (C.q * r') ≤ C.twistBound r' := by
  obtain ⟨μ, hμ, hmargin⟩ := C.table.exists_margin C.table.lengths htab
  have hδ : ∀ δ, 0 < δ → δ < μ →
      Z 0 (C.q * r') ≤ C.traceBound r' + 1 * δ ∧ Z 1 (C.q * r') ≤ C.twistBound r' + 1 * δ ∧
        Z 2 (C.q * r') ≤ C.twistBound r' + 1 * δ := by
    intro δ hδ hδμ
    obtain ⟨ι, _, sec, x, hx⟩ := happ δ hδ
    have htable : C.table.Holds (SectorScheme.ion.trace sec x) := by
      intro n hn g b
      have a := abs_le.1 (hx n (mem_union_left _ (List.mem_toFinset.2 hn)) g)
      have hm := hmargin n hn g b
      exact ⟨fun hb ↦ by linarith [hm.1 hb], fun hb ↦ by linarith [hm.2 hb]⟩
    obtain ⟨h0, h1, h2⟩ := trace_le_of_table hC htable hr
    have c0 := abs_le.1 (hx _ hmem (0 : Fin 3))
    have c1 := abs_le.1 (hx _ hmem (1 : Fin 3))
    have c2 := abs_le.1 (hx _ hmem (2 : Fin 3))
    exact ⟨by linarith, by linarith, by linarith⟩
  exact ⟨le_of_forall_small one_pos hμ fun δ h h' ↦ (hδ δ h h').1,
    le_of_forall_small one_pos hμ fun δ h h' ↦ (hδ δ h h').2.1,
    le_of_forall_small one_pos hμ fun δ h h' ↦ (hδ δ h h').2.2⟩

theorem mem_lengths₀_k : C.q * C.r ∈ C.lengths₀ := by simp [lengths₀, k]

theorem mem_lengths₀_four : C.q * (4 * C.r) ∈ C.lengths₀ := by
  simp only [lengths₀, k, mem_union, mem_insert, mem_singleton]
  right; right; ring

/-- **The lifted seed.** Let `Z g n β` be the twisted traces of a ring with the twists `1, P, C`.
Assume the sector representation at `β₀` (on `lengths₀`) and at `β = (l/j) β₀` (on
`lengths₁`), the table of the shifted traces at `β₀` (strictly), Boltzmann representations of the
three traces at length `k`, and a Boltzmann representation of `Z_{4k}` with one energy below
`-κ · 4k` (a trial state). Then `1 - p ≤ P_s(2k, β)` and `1 - q ≤ P_τ(4k, β)`. -/
theorem purities (hC : C.Valid) {Z : Fin 3 → ℕ → ℝ → ℝ}
    (happ₀ : SectorApprox SectorScheme.ion (fun g n ↦ Z g n C.beta) C.lengths₀)
    (happ : SectorApprox SectorScheme.ion (fun g n ↦ Z g n C.liftBeta) C.lengths₁)
    (htab : C.table.HoldsStrict (fun g n ↦ shiftedTrace C.beta C.kappa n (Z g n C.beta)))
    (hrep : ∀ g, BoltzmannRep (Z g C.k)) {ι : Type*} [Fintype ι] (En : ι → ℝ)
    (hE : ∀ β, Z 0 (4 * C.k) β = ∑ i, exp (-(β * En i)))
    (htrial : ∃ i, En i < -((C.kappa : ℝ) * ((4 * C.k : ℕ) : ℝ)))
    (hpos : 0 < Z 0 (2 * C.k) C.liftBeta) :
    1 - (C.spatialDefect : ℝ) ≤ Z 0 (4 * C.k) C.liftBeta / Z 0 (2 * C.k) C.liftBeta ^ 2 ∧
      1 - (C.thermalDefect : ℝ) ≤
        Z 0 (4 * C.k) (2 * C.liftBeta) / Z 0 (4 * C.k) C.liftBeta ^ 2 := by
  have hC' := hC
  obtain ⟨-, -, ⟨hq, hr, -, -, -, hβ, hj, hjl⟩, ⟨hD0, hE0, hDm, hEm, hτ, hτ₀⟩,
    ⟨ht₀, hD2, hdom, hchain, hlast, hp1, hp⟩, ⟨hRq, hqR, hqd1, hqd⟩⟩ := hC'
  set β₀ : ℝ := (C.beta : ℝ) with hβ₀def
  set κ : ℝ := (C.kappa : ℝ)
  set b : ℝ := (C.liftBeta : ℝ)
  have hj0 : (0 : ℝ) < C.j := by exact_mod_cast (show 0 < C.j by omega)
  have hl0 : (0 : ℝ) < C.l := by exact_mod_cast (show 0 < C.l by omega)
  have hb : b = C.l / C.j * β₀ := by simp only [b, β₀, liftBeta]; push_cast; ring
  have hβ₀ : 0 < β₀ := by rw [hβ₀def]; exact_mod_cast hβ
  have hb0 : 0 < b := by rw [hb]; positivity
  obtain ⟨i₀, hi₀⟩ := htrial
  -- Step 1: the capped-mass bounds at `β₀`.
  set s₀ := exp (-(β₀ * κ))
  have hsh₀ (n : ℕ) (z : ℝ) : shiftedTrace β₀ κ n z = s₀ ^ n * z := shiftedTrace_eq_pow β₀ κ n z
  have hs₀ : 0 < s₀ := exp_pos _
  have htab' : C.table.HoldsStrict (fun (g : Fin 3) n ↦ s₀ ^ n * Z g n β₀) := by
    intro n hn g c
    have := htab n hn g c
    simp only [hsh₀] at this
    exact this
  obtain ⟨hm, hm1, hm2⟩ := le_of_table hC (happ₀.scale hs₀) htab' hr mem_lengths₀_k
  have hq4 := (le_of_table hC (happ₀.scale hs₀) htab' (by omega) mem_lengths₀_four).1
  rw [show C.q * (4 * C.r) = 4 * C.k by simp only [k]; ring] at hq4
  -- Step 2: the temperature lifting of the traces at length `k`.
  have hlift (g : Fin 3) {M Dg : ℚ} (hM : s₀ ^ C.k * Z g C.k β₀ ≤ M) (hD : 0 ≤ Dg)
      (hMD : M ^ C.l < Dg ^ C.j) : shiftedTrace b κ C.k (Z g C.k b) < Dg := by
    have := ((hrep g).shift κ C.k).lt_of_lift hj hjl (β := β₀) (m := M) (D := Dg)
      (by rw [hsh₀]; exact hM) (by exact_mod_cast hD) (by exact_mod_cast hMD)
    simpa only [← hb] using this
  -- Step 3: the spatial purity at `β`.
  set s := exp (-(b * κ))
  have hs0 : 0 < s := exp_pos _
  have hsh (n : ℕ) (z : ℝ) : shiftedTrace b κ n z = s ^ n * z := shiftedTrace_eq_pow b κ n z
  set Zs : Fin 3 → ℕ → ℝ := fun g n ↦ s ^ n * Z g n b
  have hZ0 : Zs 0 C.k < C.D := by
    change s ^ C.k * Z 0 C.k b < C.D; rw [← hsh]; exact hlift 0 hm hD0 hDm
  have hZ1 : Zs 1 C.k < C.E := by
    change s ^ C.k * Z 1 C.k b < C.E; rw [← hsh]; exact hlift 1 hm1 hE0 hEm
  have hZ2 : Zs 2 C.k < C.E := by
    change s ^ C.k * Z 2 C.k b < C.E; rw [← hsh]; exact hlift 2 hm2 hE0 hEm
  have hZ4 : 1 < Zs 0 (4 * C.k) := by
    change 1 < s ^ (4 * C.k) * Z 0 (4 * C.k) b
    rw [← hsh, hE]; exact one_lt_shiftedTrace_of_lt En hb0 hi₀
  have hZ2k : 0 < Zs 0 (2 * C.k) := mul_pos (pow_pos hs0 _) hpos
  set μ := min (min (C.D - Zs 0 C.k) (C.E - Zs 1 C.k)) (min (C.E - Zs 2 C.k) (Zs 0 (4 * C.k) - 1))
  have hμ : 0 < μ :=
    lt_min (lt_min (by linarith) (by linarith)) (lt_min (by linarith) (by linarith))
  have hp1' : (C.spatialDefect : ℝ) ≤ 1 := by exact_mod_cast hp1
  have happrox : ∀ δ, 0 < δ → δ < min μ (Zs 0 (2 * C.k)) →
      (1 - (C.spatialDefect : ℝ)) * (Zs 0 (2 * C.k) - δ) ^ 2 ≤ Zs 0 (4 * C.k) + δ ∧
        Zs 0 (4 * C.k) ≤ Zs 0 (4 * C.k) + δ := by
    intro δ hδ hδμ
    have hδ' : δ < μ := lt_of_lt_of_le hδμ (min_le_left _ _)
    have hδ2 : δ < Zs 0 (2 * C.k) := lt_of_lt_of_le hδμ (min_le_right _ _)
    have m1 : δ < C.D - Zs 0 C.k := hδ'.trans_le ((min_le_left _ _).trans (min_le_left _ _))
    have m2 : δ < C.E - Zs 1 C.k := hδ'.trans_le ((min_le_left _ _).trans (min_le_right _ _))
    have m3 : δ < C.E - Zs 2 C.k := hδ'.trans_le ((min_le_right _ _).trans (min_le_left _ _))
    have m4 : δ < Zs 0 (4 * C.k) - 1 :=
      hδ'.trans_le ((min_le_right _ _).trans (min_le_right _ _))
    obtain ⟨ι', _, sec, x, hx⟩ := happ.scale hs0 δ hδ
    have c (g : Fin 3) {n : ℕ} (hn : n ∈ C.lengths₁) := abs_le.1 (hx n hn g)
    have hk1 : C.k ∈ C.lengths₁ := by simp [lengths₁]
    have hk2 : 2 * C.k ∈ C.lengths₁ := by simp [lengths₁]
    have hk4 : 4 * C.k ∈ C.lengths₁ := by simp [lengths₁]
    have c0 := c 0 hk1
    have c1 := c 1 hk1
    have c2 := c 2 hk1
    have c22 := c 0 hk2
    have c4 := c 0 hk4
    -- The weights `w = x^k` and their sector moments.
    set w : ι' → ℝ := fun i ↦ x i ^ C.k
    have hkeven : Even C.k := hq.mul_right C.r
    have hw : ∀ s, ∀ i ∈ sec s, 0 ≤ w i := fun _ i _ ↦ hkeven.pow_nonneg (x i)
    have e1 : ionSum sec w 1 = SectorScheme.ion.trace sec x (0 : Fin 3) C.k := by
      have := ion_trace_zero_mul sec x C.k 1; rw [mul_one] at this; exact this.symm
    have e2 : ionSum sec w 2 = SectorScheme.ion.trace sec x (0 : Fin 3) (2 * C.k) := by
      have := ion_trace_zero_mul sec x C.k 2; rw [mul_comm] at this; exact this.symm
    have e4 : ionSum sec w 4 = SectorScheme.ion.trace sec x (0 : Fin 3) (4 * C.k) := by
      have := ion_trace_zero_mul sec x C.k 4; rw [mul_comm] at this; exact this.symm
    have t0 := SectorScheme.ion_trace_zero_fin3 sec x C.k
    have t1 := SectorScheme.ion_trace_one_fin3 sec x C.k
    have t2 := SectorScheme.ion_trace_two_fin3 sec x C.k
    have p0 : psum sec w 0 1 = ∑ i ∈ sec (0 : Fin 3), x i ^ C.k := by simp only [psum, pow_one, w]
    have p1 : psum sec w 1 1 = ∑ i ∈ sec (1 : Fin 3), x i ^ C.k := by simp only [psum, pow_one, w]
    have hτ' : ((C.D + 11 * C.E) / 12 : ℝ) ≤ C.tau := by exact_mod_cast hτ
    have hτ₀' : ((C.D + 3 * C.E) / 4 : ℝ) ≤ C.tau0 := by exact_mod_cast hτ₀
    have key := ion_spatial hw ht₀ hD2 hdom hchain hlast hp1 hp (by rw [e1]; linarith)
      (by rw [p0]; linarith) (by rw [p0, p1]; linarith) (by rw [e4]; linarith)
    rw [e2, e4] at key
    refine ⟨?_, by linarith⟩
    have hp0 : 0 ≤ 1 - (C.spatialDefect : ℝ) := by linarith
    calc (1 - (C.spatialDefect : ℝ)) * (Zs 0 (2 * C.k) - δ) ^ 2
        ≤ (1 - C.spatialDefect) * SectorScheme.ion.trace sec x (0 : Fin 3) (2 * C.k) ^ 2 :=
          mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by linarith) (by linarith) 2) hp0
      _ ≤ _ := key
      _ ≤ Zs 0 (4 * C.k) + δ := by linarith
  have hS := (seed_limit hp1' hZ2k (lt_min hμ hZ2k) happrox).1
  -- Step 4: the thermal purity from the dominant Boltzmann factor at length `4k`.
  have hE' (β : ℝ) : shiftedTrace β κ (4 * C.k) (Z 0 (4 * C.k) β) =
      ∑ i, exp (-(β * (En i + κ * ((4 * C.k : ℕ) : ℝ)))) := by
    rw [hE, BoltzmannRep.shiftedTrace_eq]
  have hq' : ∑ i, exp (-(β₀ * (En i + κ * ((4 * C.k : ℕ) : ℝ)))) ≤ C.traceBound (4 * C.r) := by
    rw [← hE', hsh₀]; exact hq4
  have h1' : 1 < exp (-(β₀ * (En i₀ + κ * ((4 * C.k : ℕ) : ℝ)))) := by
    rw [← exp_zero]; exact exp_lt_exp.2 (by nlinarith)
  have hT := boltzmann_sq_le_of_dominant (fun i ↦ En i + κ * ((4 * C.k : ℕ) : ℝ)) hj hjl h1'
    (R := (C.Rq : ℝ))
    (by exact_mod_cast hRq) hq' (by exact_mod_cast hqR)
  rw [← hb, ← hE', ← hE', shiftedTrace_sq_le_iff] at hT
  -- Back to the traces.
  have hZb : 0 < Z 0 (4 * C.k) b := by
    rw [hE]; exact sum_pos (fun _ _ ↦ exp_pos _) ⟨i₀, mem_univ _⟩
  have hZ2b : 0 ≤ Z 0 (4 * C.k) (2 * b) := by
    rw [hE]; exact sum_nonneg fun _ _ ↦ (exp_pos _).le
  refine ⟨?_, ?_⟩
  · rw [le_div_iff₀ (by positivity)]
    have e : s ^ (4 * C.k) = (s ^ (2 * C.k)) ^ 2 := by ring
    have hS' : (s ^ (2 * C.k)) ^ 2 * ((1 - C.spatialDefect) * Z 0 (2 * C.k) b ^ 2) ≤
        (s ^ (2 * C.k)) ^ 2 * Z 0 (4 * C.k) b := by
      simp only [Zs] at hS; rw [e] at hS; nlinarith
    exact le_of_mul_le_mul_left hS' (by positivity)
  · rw [le_div_iff₀ (by positivity)]
    have hqd' : (1 - (C.thermalDefect : ℝ)) * (1 + C.Rq) ^ 2 ≤ 1 := by exact_mod_cast hqd
    have hqd0 : 0 ≤ 1 - (C.thermalDefect : ℝ) := by
      have : (C.thermalDefect : ℝ) ≤ 1 := by exact_mod_cast hqd1
      linarith
    calc (1 - (C.thermalDefect : ℝ)) * Z 0 (4 * C.k) b ^ 2
        ≤ (1 - C.thermalDefect) * ((1 + C.Rq) ^ 2 * Z 0 (4 * C.k) (2 * b)) :=
          mul_le_mul_of_nonneg_left hT hqd0
      _ = ((1 - C.thermalDefect) * (1 + C.Rq) ^ 2) * Z 0 (4 * C.k) (2 * b) := by ring
      _ ≤ 1 * Z 0 (4 * C.k) (2 * b) := mul_le_mul_of_nonneg_right hqd' hZ2b
      _ = Z 0 (4 * C.k) (2 * b) := one_mul _

end LiftedSeed

end CertifiedGaps
