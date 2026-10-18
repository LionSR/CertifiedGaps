/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.ChainCheck
import CertifiedGaps.Open.Criterion

/-!
# Rational certificates for the boundary-field gap

An `OpenCertificate` is the finite rational data that `gapGT_of_openChain` needs besides the
four seeds: the seed half-lengths `N₀` (ring) and `K₀` (boundary defect), the ring defect `u₀`
with its doubling coefficient `C`, the base steps `(qᵢ, Eᵢ, ηᵢ)` and the hand-over `x`.
`OpenCertificate.Valid` is the list of rational inequalities, which is decidable, so a concrete
certificate is checked by `decide +kernel`. `gapGT_of_openCertificate` turns a valid certificate
and the four seed inequalities into the gap bound `log(1/(10x)) / (β₀ 2^m)` for every odd chain
of at least `4 N₀ 2^m` sites.

The data are produced by `lean/scripts/open_certificate.jl` from the seeds; that script only
searches, and Lean re-checks every inequality.

## Main results

* `gapGT_of_openCertificate`: the boundary-field gap from a valid certificate and four seeds.
* `boundaryFieldCertificate`: the certificate of the note's eqs. (C1)–(C6).
* `gapGT_of_openSeeds`: the boundary-field gap `Δ_n > (log 10)/392` for every odd `n ≥ 1921`
  from the four seeds of the note's eq. (C1).
-/

open Real

namespace CertifiedGaps

/-- The rational data of an open-chain certificate with `m` base steps (`m` the length of `q`).
The defects `E` and `η` have `m + 1` entries: the seeds at index `0`, then one per step. -/
structure OpenCertificate where
  /-- Half the ring seed length: the ring seeds are `P_s(2N₀, β₀)` and `P_τ(4N₀, β₀)`. -/
  N₀ : ℕ
  /-- Half the boundary seed length: the boundary seed is `e_{2K₀}(β₀)`. -/
  K₀ : ℕ
  /-- The ring seed defect. -/
  u₀ : ℚ
  /-- The ring doubling coefficient, with `G(u₀) ≤ C ≤ 6`. -/
  C : ℚ
  /-- The suppressions `qᵢ ≥ ρᵢ₊₁^{4Nᵢ - 2Kᵢ}` of the base steps. -/
  q : List ℚ
  /-- Upper bounds on the boundary defects `e_{2Kᵢ}(βᵢ)`. -/
  E : List ℚ
  /-- Upper bounds on the open-chain thermal defects `1 - P_τ(4Nᵢ, βᵢ)`. -/
  η : List ℚ
  /-- The hand-over defect `≤ 1/100`: the gap bound is `log(1/(10x)) / (β₀ 2^m)`. -/
  x : ℚ

namespace OpenCertificate

/-- The tail bound `d(ε) = ε / (2 - 3ε)` over `ℚ`. -/
def tailQ (ε : ℚ) : ℚ := ε / (2 - 3 * ε)

/-- The overlap update `F(q, e, t)` over `ℚ`. -/
def overlapQ (q e t : ℚ) : ℚ := q / (1 - q) * (1 / ((1 - e) ^ 2 * t) - 1)

/-- The doubling coefficient `G(u)` over `ℚ`. -/
def coefficientQ (u : ℚ) : ℚ := ((1 - u)⁻¹ + (1 - u)⁻¹ ^ 2 + (1 - u)⁻¹ ^ 3) ^ 2 / 2

variable (c : OpenCertificate)

/-- The number of base steps. -/
def m : ℕ := c.q.length

/-- The suppression of step `i`. -/
def qAt (i : ℕ) : ℚ := c.q.getD i 0

/-- The boundary defect bound at scale `i`. -/
def EAt (i : ℕ) : ℚ := c.E.getD i 0

/-- The open-chain thermal defect bound at scale `i`. -/
def ηAt (i : ℕ) : ℚ := c.η.getD i 0

/-- The ring defect `εᵢ = (C u₀)^{2^i} / C` at scale `i`. -/
def εAt (i : ℕ) : ℚ := (c.C * c.u₀) ^ 2 ^ i / c.C

/-- The open-chain length `4Nᵢ` at scale `i`. -/
def lengthAt (i : ℕ) : ℕ := 2 * (2 * (c.N₀ * 2 ^ i))

/-- The exponent `4Nᵢ - 2Kᵢ` of the suppression at step `i`. -/
def gapAt (i : ℕ) : ℕ := c.lengthAt i - 2 * boundaryHalf c.N₀ c.K₀ i

/-- The conditions of base step `i` (`OpenStep`), over `ℚ`. -/
def StepValid (i : ℕ) : Prop :=
  0 ≤ c.qAt i ∧ c.qAt i < 1 ∧
    tailQ (c.εAt (i + 1)) ^ (c.gapAt i / (c.gapAt i).gcd (c.lengthAt i)) ≤
      c.qAt i ^ (c.lengthAt i / (c.gapAt i).gcd (c.lengthAt i)) ∧
    c.EAt i < 1 ∧ c.ηAt i < 1 ∧ overlapQ (c.qAt i) (c.EAt i) (1 - c.ηAt i) ≤ c.EAt (i + 1) ∧
    c.EAt (i + 1) ≤ 1 ∧
    0 < (1 - c.ηAt i) * ((1 - c.εAt i) * (1 - c.εAt (i + 1))) * (1 - c.EAt (i + 1)) ∧
    ChainCertificate.fQ
      (1 - (1 - c.ηAt i) * ((1 - c.εAt i) * (1 - c.εAt (i + 1))) * (1 - c.EAt (i + 1))) ≤
        c.ηAt (i + 1)

/-- The rational conditions of `gapGT_of_openChain`. -/
def Valid : Prop :=
  1 ≤ c.N₀ ∧ 1 ≤ c.K₀ ∧ c.K₀ ≤ 2 * c.N₀ ∧ 1 ≤ c.m ∧ 0 < c.u₀ ∧ c.u₀ ≤ 1 / 10 ∧
    coefficientQ c.u₀ ≤ c.C ∧ c.C ≤ 6 ∧ c.C * c.u₀ ≤ 1 ∧ 0 < c.x ∧ c.x ≤ 1 / 100 ∧
    c.EAt c.m ≤ c.x ∧ c.ηAt c.m ≤ c.x ∧ c.εAt c.m ≤ c.x ∧ ∀ i < c.m, c.StepValid i

/-- The conditions are decidable, so `decide` checks a concrete certificate. -/
instance : Decidable c.Valid := by unfold Valid StepValid; infer_instance

/-- `tailQ` computes the tail bound. -/
theorem cast_tailQ (ε : ℚ) : tailBound ε = (tailQ ε : ℝ) := by simp [tailBound, tailQ]

/-- `overlapQ` computes the overlap update. -/
theorem cast_overlapQ (q e t : ℚ) : overlapUpdate q e t = (overlapQ q e t : ℝ) := by
  simp [overlapUpdate, overlapQ]

/-- `coefficientQ` computes the doubling coefficient. -/
theorem cast_coefficientQ (u : ℚ) : doublingCoefficient u = (coefficientQ u : ℝ) := by
  simp [doublingCoefficient, coefficientQ]

/-- `εAt` computes the ring defects. -/
theorem cast_εAt (i : ℕ) : coefficientDefect (c.C : ℝ) c.u₀ i = (c.εAt i : ℝ) := by
  simp [coefficientDefect, εAt]

/-- A valid step gives `OpenStep`. -/
theorem StepValid.openStep {c : OpenCertificate} {i : ℕ} (h : c.StepValid i) :
    OpenStep (coefficientDefect (c.C : ℝ) c.u₀ i) (coefficientDefect (c.C : ℝ) c.u₀ (i + 1))
      (c.qAt i) (c.EAt i) (c.ηAt i) (c.EAt (i + 1)) (c.ηAt (i + 1)) (c.gapAt i)
      (c.lengthAt i) := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9⟩ := h
  rw [cast_εAt, cast_εAt]
  have hprod : ((1 - c.ηAt i) * ((1 - c.εAt i) * (1 - c.εAt (i + 1))) * (1 - c.EAt (i + 1)) :
      ℚ) = (1 - (c.ηAt i : ℝ)) * ((1 - c.εAt i) * (1 - c.εAt (i + 1))) * (1 - c.EAt (i + 1)) := by
    push_cast; ring
  refine ⟨by exact_mod_cast h1, by exact_mod_cast h2, ?_, by exact_mod_cast h4,
    by exact_mod_cast h5, ?_, by exact_mod_cast h7, ?_, ?_⟩
  · rw [cast_tailQ]; exact_mod_cast h3
  · rw [show (1 : ℝ) - c.ηAt i = ((1 - c.ηAt i : ℚ) : ℝ) by push_cast; ring, cast_overlapQ]
    exact_mod_cast h6
  · rw [← hprod]; exact_mod_cast h8
  · rw [← hprod, show (1 : ℝ) - ((1 - c.ηAt i) * ((1 - c.εAt i) * (1 - c.εAt (i + 1))) *
        (1 - c.EAt (i + 1)) : ℚ) = ((1 - (1 - c.ηAt i) * ((1 - c.εAt i) * (1 - c.εAt (i + 1))) *
        (1 - c.EAt (i + 1)) : ℚ) : ℝ) by push_cast; ring, ChainCertificate.cast_fQ]
    exact_mod_cast h9

end OpenCertificate

section Main

variable {ιr ιo : ℕ → Type*} [∀ L, Fintype (ιr L)] [∀ L, Nonempty (ιr L)]
  [∀ n, Fintype (ιo n)] [∀ n, Nonempty (ιo n)] (Er : (L : ℕ) → ιr L → ℝ)
  (Eo : (n : ℕ) → ιo n → ℝ)

open OpenCertificate in
/-- **The boundary-field gap from a rational certificate.** Let the ring spectra `Er` satisfy
the spatial moment inequalities and the ring and open-chain traces the boundary moment
inequalities at every `β > 0`. A valid certificate `c` and the four seeds

* `P_s(2N₀, β₀), P_τ(4N₀, β₀) ≥ 1 - u₀` for the ring,
* the boundary defect `e_{2K₀}(β₀) ≤ E₀`, and
* the open-chain thermal purity `P_τ(4N₀, β₀) ≥ 1 - η₀`

give `Δ_n > log(1/(10x)) / (β₀ 2^m)` for every odd open chain of `n ≥ 4 N₀ 2^m` sites. -/
theorem gapGT_of_openCertificate (hring : HasSpatialMomentBounds Er 2)
    (hopen : ∀ β : ℝ, 0 < β → OpenMomentBounds (fun L ↦ partitionFunction Er L β)
      (fun n ↦ partitionFunction Eo n β))
    (c : OpenCertificate) (hc : c.Valid) {β₀ : ℝ} (hβ₀ : 0 < β₀)
    (hs : 1 - (c.u₀ : ℝ) ≤ spatialPurity Er (2 * c.N₀) β₀)
    (ht : 1 - (c.u₀ : ℝ) ≤ thermalPurity Er (2 * (2 * c.N₀)) β₀)
    (he : boundaryDefect (fun L ↦ partitionFunction Er L β₀) (fun n ↦ partitionFunction Eo n β₀)
      (2 * c.K₀) ≤ c.EAt 0)
    (hη : 1 - (c.ηAt 0 : ℝ) ≤ thermalPurity Eo (2 * (2 * c.N₀)) β₀) {n : ℕ}
    (hn : 4 * c.N₀ * 2 ^ c.m ≤ n) (hno : Odd n) :
    GapGT (Eo n) (log (1 / (10 * (c.x : ℝ))) / (β₀ * 2 ^ c.m)) := by
  obtain ⟨hN₀, hK₀, hK₀N, hm, hu₀, hu₁, hG, hC, hCu, hx0, hx, hEm, hηm, hεm, hstep⟩ := hc
  refine gapGT_of_openChain Er Eo hN₀ hK₀ hK₀N hm (fun β hβ ↦ (hring β hβ).mono (by omega))
    hopen hβ₀ (u₀ := c.u₀) (C := c.C) (by exact_mod_cast hu₀)
    (by rw [show (1 / 10 : ℝ) = ((1 / 10 : ℚ) : ℝ) by norm_num]; exact_mod_cast hu₁)
    (by rw [cast_coefficientQ]; exact_mod_cast hG) (by exact_mod_cast hC)
    (by exact_mod_cast hCu) (fun i ↦ c.EAt i) (fun i ↦ c.ηAt i) (fun i ↦ c.qAt i)
    (fun i hi ↦ (hstep i hi).openStep) (by exact_mod_cast hx0)
    (by rw [show (1 / 100 : ℝ) = ((1 / 100 : ℚ) : ℝ) by norm_num]; exact_mod_cast hx)
    (by exact_mod_cast hEm) (by exact_mod_cast hηm) (by rw [cast_εAt]; exact_mod_cast hεm)
    hs ht he hη hn hno

end Main

/-- The certificate of the note's eqs. (C1)–(C6): ring seeds at `n = 120` and `240` with
`u₀ = 0.052` and `C = 6`, the boundary seed at length `142`, three base steps and the hand-over
`x = 1/100` at scale `3`. -/
def boundaryFieldCertificate : OpenCertificate where
  N₀ := 60
  K₀ := 71
  u₀ := 13 / 250
  C := 6
  q := [71 / 500, 141 / 5000, 7 / 2500]
  E := [14 / 100, 12 / 100, 16 / 1000, 1 / 1000]
  η := [21 / 100, 155 / 1000, 26 / 1000, 1 / 1000]
  x := 1 / 100

/-- The certificate of the note is valid. -/
theorem boundaryFieldCertificate_valid : boundaryFieldCertificate.Valid := by decide +kernel

/-- A sharper certificate from the same seeds as `boundaryFieldCertificate`, found by
`lean/scripts/open_certificate.jl 60 71 0.052 0.14 0.21 49 3`: the coefficient `C = 5.583`
in place of `6`, four-digit bounds along the three base steps, and the hand-over at the actual
defect `x = 2.241 · 10⁻⁴` in place of `1/100`. The gap bound becomes
`log(1/(10x)) / 392 = log(10⁶/2241) / 392 ≈ 0.01556`. -/
def sharpBoundaryFieldCertificate : OpenCertificate where
  N₀ := 60
  K₀ := 71
  u₀ := 52e-3
  C := 5583e-3
  q := [1373e-4, 2525e-5, 2126e-6]
  E := [14e-2, 1133e-4, 1241e-5, 9706e-8]
  η := [21e-2, 1399e-4, 1938e-5, 2241e-7]
  x := 2241e-7

/-- The sharper certificate is valid. -/
theorem sharpBoundaryFieldCertificate_valid : sharpBoundaryFieldCertificate.Valid := by
  decide +kernel

section Main

variable {ιr ιo : ℕ → Type*} [∀ L, Fintype (ιr L)] [∀ L, Nonempty (ιr L)]
  [∀ n, Fintype (ιo n)] [∀ n, Nonempty (ιo n)] (Er : (L : ℕ) → ιr L → ℝ)
  (Eo : (n : ℕ) → ιo n → ℝ)

/-- **The boundary-field gap criterion** (Theorem 2 of the note). Let the ring spectra `Er`
satisfy the spatial moment inequalities and the ring and open-chain traces the boundary moment
inequalities at every `β > 0`. If

* `P_s(120, 49) ≥ 0.948` and `P_τ(240, 49) ≥ 0.948` for the ring,
* the boundary defect `e₁₄₂ ≤ 0.14` at `β = 49`, and
* the open-chain thermal purity `P_τ(240, 49) ≥ 0.79`,

then every odd open chain of `n ≥ 1921` sites has a nondegenerate ground level and gap
`Δ_n > (log 10)/392 ≈ 0.005874`. -/
theorem gapGT_of_openSeeds (hring : HasSpatialMomentBounds Er 2)
    (hopen : ∀ β : ℝ, 0 < β → OpenMomentBounds (fun L ↦ partitionFunction Er L β)
      (fun n ↦ partitionFunction Eo n β))
    (hs : 948 / 1000 ≤ spatialPurity Er 120 49) (ht : 948 / 1000 ≤ thermalPurity Er 240 49)
    (he : boundaryDefect (fun L ↦ partitionFunction Er L 49) (fun n ↦ partitionFunction Eo n 49)
      142 ≤ 14 / 100)
    (hη : 79 / 100 ≤ thermalPurity Eo 240 49) {n : ℕ} (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (Eo n) (log 10 / 392) := by
  have h := gapGT_of_openCertificate Er Eo hring hopen boundaryFieldCertificate
    boundaryFieldCertificate_valid (β₀ := 49) (by norm_num)
    (by norm_num [boundaryFieldCertificate]; linarith)
    (by norm_num [boundaryFieldCertificate]; linarith)
    (by norm_num [boundaryFieldCertificate, OpenCertificate.EAt]; linarith)
    (by norm_num [boundaryFieldCertificate, OpenCertificate.ηAt]; linarith)
    (n := n) (by simp [boundaryFieldCertificate, OpenCertificate.m]; omega) hno
  rw [show boundaryFieldCertificate.m = 3 from rfl,
    show (boundaryFieldCertificate.x : ℝ) = 1 / 100 by norm_num [boundaryFieldCertificate]] at h
  convert h using 2 <;> norm_num

end Main

end CertifiedGaps
