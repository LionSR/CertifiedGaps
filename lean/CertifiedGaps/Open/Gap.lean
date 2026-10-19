/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.Bounds
import CertifiedGaps.Open.Certificate

/-!
# The boundary-field gap for open chains of reflection-positive rings

For a reflection-positive ring without internal reflection and boundary weights `b`, the ring
satisfies the spatial moment inequalities (`RPRing.hasSpatialMomentBounds`) and the ring and
open-chain traces satisfy the boundary moment inequalities (`RPRing.openMomentBounds`). The
criterion `gapGT_of_openCertificate` then turns a valid rational certificate and four seed
inequalities into a uniform gap for every long enough chain, odd or even; `gapGT_of_openSeeds` is
the case of the note's eq. (C1).

## Main results

* `RPRing.openGapGT_of_certificate`: the boundary-field gap from a certificate and four seeds.
* `RPRing.openGapGT_of_traces`: the same, with the boundary seed replaced by finite traces.
* `RPRing.freeGapGT_of_traces`: the gap with open boundary conditions (free ends).
* `RPRing.openGapGT_of_seeds`: the boundary-field gap from the four seeds of the note.
-/

open Matrix Real

namespace CertifiedGaps

namespace RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- The open-chain thermal purity `P_τ(n, β) = G_n(2β) / G_n(β)²`. -/
noncomputable def openThermalPurity (b : R.Boundary) (n : ℕ) (β : ℝ) : ℝ :=
  R.openPartitionFunction b n (2 * β) / R.openPartitionFunction b n β ^ 2

/-- The ring and open-chain partition functions are the Boltzmann sums of their spectra. -/
theorem partitionFunction_spectra (b : R.Boundary) :
    CertifiedGaps.partitionFunction (ι := fun L ↦ Fin L → d) R.spectrum =
        (fun L β ↦ R.partitionFunction L β) ∧
      CertifiedGaps.partitionFunction (ι := fun n ↦ Fin n → d) (R.openSpectrum b) =
        fun n β ↦ R.openPartitionFunction b n β := by
  constructor <;> funext L β
  · exact (R.partitionFunction_eq_sum L β).symm
  · exact (R.openPartitionFunction_eq_sum b L β).symm

/-- **The boundary-field gap from a rational certificate.** Let the ring have no internal
reflection and let `b` be boundary weights. A valid certificate `c` and the seeds

* `P_s(2N₀, β₀), P_τ(4N₀, β₀) ≥ 1 - u₀` for the ring,
* the boundary defect `e_{2K₀}(β₀) ≤ E₀`, and
* the open-chain thermal purity `P_τ(4N₀, β₀) ≥ 1 - η₀`

give `Δ_n > log(1/(10x)) / (β₀ 2^m)` for every open chain of `n ≥ 4 N₀ 2^m` sites. -/
theorem openGapGT_of_certificate [Nonempty d] (hrefl : R.refl = 1) (b : R.Boundary)
    (c : OpenCertificate) (hc : c.Valid) {β₀ : ℝ} (hβ₀ : 0 < β₀)
    (hs : 1 - (c.u₀ : ℝ) ≤ R.spatialPurity (2 * c.N₀) β₀)
    (ht : 1 - (c.u₀ : ℝ) ≤ R.thermalPurity (2 * (2 * c.N₀)) β₀)
    (he : boundaryDefect (fun L ↦ R.partitionFunction L β₀)
      (fun n ↦ R.openPartitionFunction b n β₀) (2 * c.K₀) ≤ c.EAt 0)
    (hη : 1 - (c.ηAt 0 : ℝ) ≤ R.openThermalPurity b (2 * (2 * c.N₀)) β₀) {n : ℕ}
    (hn : 4 * c.N₀ * 2 ^ c.m ≤ n) :
    GapGT (R.openSpectrum b n) (log (1 / (10 * (c.x : ℝ))) / (β₀ * 2 ^ c.m)) := by
  obtain ⟨hZ, hG⟩ := R.partitionFunction_spectra b
  refine gapGT_of_openCertificate (ιr := fun L ↦ Fin L → d) (ιo := fun n ↦ Fin n → d) R.spectrum
    (R.openSpectrum b) R.hasSpatialMomentBounds (fun β hβ ↦ ?_) c hc hβ₀ ?_ ?_ ?_ ?_ hn
  · rw [hZ, hG]; exact R.openMomentBounds hrefl b hβ
  · rw [CertifiedGaps.spatialPurity, hZ]; exact hs
  · rw [CertifiedGaps.thermalPurity, hZ]; exact ht
  · rw [hZ, hG]; exact he
  · rw [CertifiedGaps.thermalPurity, hG]; exact hη

/-- **The open-chain gap from a certificate and finite traces.** As `openGapGT_of_certificate`,
with the boundary seed replaced by the trace inequality
`G(2K₀ - 2N₀ r) (d(u₀) Z(2N₀))^r ≤ E₀ G(2K₀)` at `β₀` (`boundaryDefect_le_of_traces`). All four
seeds are then inequalities between partition functions of finite rings and open chains. -/
theorem openGapGT_of_traces [Nonempty d] (hrefl : R.refl = 1) (b : R.Boundary)
    (c : OpenCertificate) (hc : c.Valid) {β₀ : ℝ} (hβ₀ : 0 < β₀)
    (hs : 1 - (c.u₀ : ℝ) ≤ R.spatialPurity (2 * c.N₀) β₀)
    (ht : 1 - (c.u₀ : ℝ) ≤ R.thermalPurity (2 * (2 * c.N₀)) β₀) {r : ℕ}
    (hr : c.N₀ * r + 1 ≤ c.K₀)
    (he : R.openPartitionFunction b (2 * (c.K₀ - c.N₀ * r)) β₀ *
        (tailBound c.u₀ * R.partitionFunction (2 * c.N₀) β₀) ^ r ≤
      c.EAt 0 * R.openPartitionFunction b (2 * c.K₀) β₀)
    (hη : 1 - (c.ηAt 0 : ℝ) ≤ R.openThermalPurity b (2 * (2 * c.N₀)) β₀) {n : ℕ}
    (hn : 4 * c.N₀ * 2 ^ c.m ≤ n) :
    GapGT (R.openSpectrum b n) (log (1 / (10 * (c.x : ℝ))) / (β₀ * 2 ^ c.m)) := by
  obtain ⟨hZ, hG⟩ := R.partitionFunction_spectra b
  refine gapGT_of_openCertificate_traces (ιr := fun L ↦ Fin L → d) (ιo := fun n ↦ Fin n → d)
    R.spectrum (R.openSpectrum b) R.hasSpatialMomentBounds (fun β hβ ↦ ?_) c hc hβ₀ ?_ ?_ hr ?_
    ?_ hn
  · rw [hZ, hG]; exact R.openMomentBounds hrefl b hβ
  · rw [CertifiedGaps.spatialPurity, hZ]; exact hs
  · rw [CertifiedGaps.thermalPurity, hZ]; exact ht
  · rw [hZ, hG]; exact he
  · rw [CertifiedGaps.thermalPurity, hG]; exact hη

/-- **The gap with open boundary conditions.** For a reflection-positive ring without internal
reflection, the open chains with free ends (`freeBoundary`, `openHamiltonian_free`) have a gap
`Δ_n > log(1/(10x)) / (β₀ 2^m)` at every length `n ≥ 4 N₀ 2^m` from a valid certificate and four
seeds on finite rings and free chains. -/
theorem freeGapGT_of_traces [Nonempty d] (hrefl : R.refl = 1) (c : OpenCertificate)
    (hc : c.Valid) {β₀ : ℝ} (hβ₀ : 0 < β₀)
    (hs : 1 - (c.u₀ : ℝ) ≤ R.spatialPurity (2 * c.N₀) β₀)
    (ht : 1 - (c.u₀ : ℝ) ≤ R.thermalPurity (2 * (2 * c.N₀)) β₀) {r : ℕ}
    (hr : c.N₀ * r + 1 ≤ c.K₀)
    (he : R.openPartitionFunction R.freeBoundary (2 * (c.K₀ - c.N₀ * r)) β₀ *
        (tailBound c.u₀ * R.partitionFunction (2 * c.N₀) β₀) ^ r ≤
      c.EAt 0 * R.openPartitionFunction R.freeBoundary (2 * c.K₀) β₀)
    (hη : 1 - (c.ηAt 0 : ℝ) ≤ R.openThermalPurity R.freeBoundary (2 * (2 * c.N₀)) β₀)
    {n : ℕ} (hn : 4 * c.N₀ * 2 ^ c.m ≤ n) :
    GapGT (R.openSpectrum R.freeBoundary n)
      (log (1 / (10 * (c.x : ℝ))) / (β₀ * 2 ^ c.m)) :=
  R.openGapGT_of_traces hrefl _ c hc hβ₀ hs ht hr he hη hn

/-- **The boundary-field gap for an open chain.** Let the ring have no internal reflection and
let `b` be boundary weights. If

* `P_s(120, 49) ≥ 0.948` and `P_τ(240, 49) ≥ 0.948` for the ring,
* the boundary defect `e₁₄₂ ≤ 0.14` at `β = 49`, and
* the open-chain thermal purity `P_τ(240, 49) ≥ 0.79`,

then every open chain of `n ≥ 1920` sites has a nondegenerate ground state and gap
`Δ_n > (log 10)/392`. -/
theorem openGapGT_of_seeds [Nonempty d] (hrefl : R.refl = 1) (b : R.Boundary)
    (hs : 948 / 1000 ≤ R.spatialPurity 120 49) (ht : 948 / 1000 ≤ R.thermalPurity 240 49)
    (he : boundaryDefect (fun L ↦ R.partitionFunction L 49)
      (fun n ↦ R.openPartitionFunction b n 49) 142 ≤ 14 / 100)
    (hη : 79 / 100 ≤ R.openThermalPurity b 240 49) {n : ℕ} (hn : 1920 ≤ n) :
    GapGT (R.openSpectrum b n) (log 10 / 392) := by
  obtain ⟨hZ, hG⟩ := R.partitionFunction_spectra b
  refine gapGT_of_openSeeds (ιr := fun L ↦ Fin L → d) (ιo := fun n ↦ Fin n → d) R.spectrum
    (R.openSpectrum b) R.hasSpatialMomentBounds (fun β hβ ↦ ?_) ?_ ?_ ?_ ?_ hn
  · rw [hZ, hG]; exact R.openMomentBounds hrefl b hβ
  · rw [CertifiedGaps.spatialPurity, hZ]; exact hs
  · rw [CertifiedGaps.thermalPurity, hZ]; exact ht
  · rw [hZ, hG]; exact he
  · rw [CertifiedGaps.thermalPurity, hG]; exact hη

end RPRing

end CertifiedGaps
