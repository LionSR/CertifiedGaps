/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.Bounds
import CertifiedGaps.Open.Criterion

/-!
# The boundary-field gap for open chains of reflection-positive rings

For a reflection-positive ring without internal reflection and boundary weights `b`, the ring
satisfies the spatial moment inequalities (`RPRing.hasSpatialMomentBounds`) and the ring and
open-chain traces satisfy the boundary moment inequalities (`RPRing.openMomentBounds`). The
criterion `gapGT_of_openSeeds` then turns the four seed inequalities of the note's eq. (C1)
into a uniform gap for every odd chain of at least `1921` sites.

## Main results

* `RPRing.openGapGT_of_seeds`: the boundary-field gap from the four seeds.
-/

open Real

namespace CertifiedGaps

namespace RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- The open-chain thermal purity `P_τ(n, β) = G_n(2β) / G_n(β)²`. -/
noncomputable def openThermalPurity (b : R.Boundary) (n : ℕ) (β : ℝ) : ℝ :=
  R.openPartitionFunction b n (2 * β) / R.openPartitionFunction b n β ^ 2

/-- **The boundary-field gap for an open chain.** Let the ring have no internal reflection and
let `b` be boundary weights. If

* `P_s(120, 49) ≥ 0.948` and `P_τ(240, 49) ≥ 0.948` for the ring,
* the boundary defect `e₁₄₂ ≤ 0.14` at `β = 49`, and
* the open-chain thermal purity `P_τ(240, 49) ≥ 0.79`,

then every odd open chain of `n ≥ 1921` sites has a nondegenerate ground state and gap
`Δ_n > (log 10)/392`. -/
theorem openGapGT_of_seeds [Nonempty d] (hrefl : R.refl = 1) (b : R.Boundary)
    (hs : 948 / 1000 ≤ R.spatialPurity 120 49) (ht : 948 / 1000 ≤ R.thermalPurity 240 49)
    (he : boundaryDefect (fun L ↦ R.partitionFunction L 49)
      (fun n ↦ R.openPartitionFunction b n 49) 142 ≤ 14 / 100)
    (hη : 79 / 100 ≤ R.openThermalPurity b 240 49) {n : ℕ} (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (R.openSpectrum b n) (log 10 / 392) := by
  have hZ : CertifiedGaps.partitionFunction (ι := fun L ↦ Fin L → d) R.spectrum =
      fun L β ↦ R.partitionFunction L β := by
    funext L β; exact (R.partitionFunction_eq_sum L β).symm
  have hG : CertifiedGaps.partitionFunction (ι := fun n ↦ Fin n → d) (R.openSpectrum b) =
      fun n β ↦ R.openPartitionFunction b n β := by
    funext n β; exact (R.openPartitionFunction_eq_sum b n β).symm
  refine gapGT_of_openSeeds (ιr := fun L ↦ Fin L → d) (ιo := fun n ↦ Fin n → d) R.spectrum
    (R.openSpectrum b) R.hasSpatialMomentBounds (fun β hβ ↦ ?_) ?_ ?_ ?_ ?_ hn hno
  · rw [hZ, hG]; exact R.openMomentBounds hrefl b hβ
  · rw [CertifiedGaps.spatialPurity, hZ]; exact hs
  · rw [CertifiedGaps.thermalPurity, hZ]; exact ht
  · rw [hZ, hG]; exact he
  · rw [CertifiedGaps.thermalPurity, hG]; exact hη

end RPRing

end CertifiedGaps
