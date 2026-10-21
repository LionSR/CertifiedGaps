/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Gap
import CertifiedGaps.Seed.Certificate
import CertifiedGaps.Transfer.SingleSiteD2

/-!
# Seeds of reflection-positive rings from trace tables

The entry point for a model. `RPRing.seeds_of_certificate` turns a valid `TableCertificate`, the
sector representation of the ring's (twisted) partition functions, the table bounds and the
trial-state bounds into the seed purities, in exactly the form of a model's `seeds` lemma
(`lean/README.md`, "Adding a model"):

  `1 - s₀ ≤ R.spatialPurity n₀ β₀`  and  `1 - t₀ ≤ R.thermalPurity (2 n₀) β₀`.

The sector representation is free for a ring without an internal cell reflection:

* `RPRing.sectorApprox_plain`: the untwisted traces, one sector (`SectorScheme.plain`);
* `RPRing.sectorApprox_d2`: the untwisted and `s_z`-twisted traces in the `D₂` sectors
  (`SectorScheme.d2`), from the model's `SingleSite.D2Symmetry`.

What remains model-specific is numerical: the table bounds and the trial-state bounds.
-/

open Finset Real

namespace CertifiedGaps.RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- The untwisted traces of a ring without an internal cell reflection, as one sector. -/
theorem sectorApprox_plain (h : R.refl = 1) {β : ℝ} (hβ : 0 ≤ β) (ns : Finset ℕ)
    (hns : ∀ n ∈ ns, 2 ≤ n) :
    SectorApprox SectorScheme.plain (fun _ n ↦ R.partitionFunction n β) ns := by
  intro δ hδ
  obtain ⟨D, _, x, hx⟩ := R.exists_transfer_approx h hβ ns hns hδ
  refine ⟨D, inferInstance, fun _ ↦ univ, x, fun n hn g ↦ ?_⟩
  fin_cases g
  show |SectorScheme.plain.trace (fun _ ↦ univ) x 0 n - R.partitionFunction n β| ≤ δ
  rw [SectorScheme.plain_trace]
  exact hx n hn

/-- The untwisted and `s_z`-twisted traces of a ring without an internal cell reflection, in the
`D₂` sectors of its `D2Symmetry`. -/
theorem sectorApprox_d2 (h : R.refl = 1) (S : SingleSite.D2Symmetry R.C R.o) {β : ℝ}
    (hβ : 0 ≤ β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    SectorApprox SectorScheme.d2
      ![fun n ↦ R.partitionFunction n β, fun n ↦ SingleSite.partitionFunction R.C R.o n S.sz β]
      ns :=
  (R.d2Approx h S hβ ns hns).toSectorApprox

/-- **The seeds of a reflection-positive ring from a trace table.** `Z g n β` are the tabulated
twisted partition functions, with `Z 0` the ring's own. -/
theorem seeds_of_certificate [Nonempty d] {C : TableCertificate} (hC : C.Valid)
    {Z : Fin (C.scheme.t + 1) → ℕ → ℝ → ℝ} (hZ : ∀ n β, Z 0 n β = R.partitionFunction n β)
    (happ : SectorApprox C.scheme (fun g n ↦ Z g n C.beta) C.lengths)
    (htab : C.table.HoldsStrict (fun g n ↦ shiftedTrace C.beta C.kappa n (Z g n C.beta)))
    (htrial : 1 < shiftedTrace C.beta C.kappa (2 * C.n₀)
      (R.partitionFunction (2 * C.n₀) C.beta))
    (htrial₂ : 1 < shiftedTrace (2 * C.beta) C.kappa (2 * C.n₀)
      (R.partitionFunction (2 * C.n₀) (2 * C.beta))) :
    1 - (C.seed.spatialDefect : ℝ) ≤ R.spatialPurity C.n₀ C.beta ∧
      1 - (C.thermalDefect : ℝ) ≤ R.thermalPurity (2 * C.n₀) C.beta := by
  have h := C.purities hC happ htab (by rw [hZ]; exact R.partitionFunction_pos _ _)
    (by rw [hZ]; exact htrial) (by rw [hZ]; exact htrial₂)
  simp only [hZ] at h
  exact h

end CertifiedGaps.RPRing
