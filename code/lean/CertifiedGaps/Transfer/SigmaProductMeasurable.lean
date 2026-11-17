import Mathlib.MeasureTheory.MeasurableSpace.Prod

/-!
# Joint measurability on countable disjoint unions

Componentwise joint measurability extends to the product of two countable measurable sigma
types. The proof uses measurable component embeddings and a countable union of their images.
-/

open Set Function

namespace CertifiedGaps.SigmaProductMeasurable

variable {I J : Type*} {X : I → Type*} {Y : J → Type*}
variable [∀ i, MeasurableSpace (X i)] [∀ j, MeasurableSpace (Y j)]

/-- A set in a measurable disjoint union is measurable exactly when all its sections are. -/
theorem measurableSet_sigma_iff {s : Set (Sigma X)} :
    MeasurableSet s ↔ ∀ i, MeasurableSet (Sigma.mk i ⁻¹' s) := by
  exact MeasurableSpace.measurableSet_iInf

/-- Each component inclusion into a measurable disjoint union is measurable. -/
theorem measurable_sigma_mk (i : I) : Measurable (Sigma.mk i : X i → Sigma X) := by
  intro s hs
  exact measurableSet_sigma_iff.mp hs i

/-- Each component inclusion into a measurable disjoint union is a measurable embedding. -/
theorem measurableEmbedding_sigma_mk (i : I) :
    MeasurableEmbedding (Sigma.mk i : X i → Sigma X) where
  injective := sigma_mk_injective
  measurable := measurable_sigma_mk i
  measurableSet_image' := by
    intro s hs
    rw [measurableSet_sigma_iff]
    intro j
    by_cases h : j = i
    · subst j
      simpa only [sigma_mk_injective.preimage_image] using hs
    · have hempty : Sigma.mk j ⁻¹' (Sigma.mk i '' s) = ∅ := by
        ext x
        simp only [mem_preimage, mem_image, mem_empty_iff_false, iff_false]
        rintro ⟨y, _, hy⟩
        exact h (congrArg Sigma.fst hy).symm
      rw [hempty]
      exact MeasurableSet.empty

/-- Componentwise joint measurability gives joint measurability on products of countable
measurable disjoint unions. The codomain can be any measurable space. -/
theorem measurable_of_sigma_prod [Countable I] [Countable J]
    {Z : Type*} [MeasurableSpace Z] {f : Sigma X × Sigma Y → Z}
    (hf : ∀ i j, Measurable (fun p : X i × Y j => f (⟨i, p.1⟩, ⟨j, p.2⟩))) :
    Measurable f := by
  intro s hs
  have hdecomp : f ⁻¹' s = ⋃ i, ⋃ j,
      (Prod.map (Sigma.mk i) (Sigma.mk j)) ''
        ((fun p : X i × Y j => f (⟨i, p.1⟩, ⟨j, p.2⟩)) ⁻¹' s) := by
    ext p
    rcases p with ⟨⟨i, x⟩, ⟨j, y⟩⟩
    constructor
    · intro h
      exact mem_iUnion.mpr ⟨i, mem_iUnion.mpr ⟨j, ⟨(x, y), h, rfl⟩⟩⟩
    · intro h
      rcases mem_iUnion.mp h with ⟨i', h⟩
      rcases mem_iUnion.mp h with ⟨j', q, hq, heq⟩
      change f (Prod.map (Sigma.mk i') (Sigma.mk j') q) ∈ s at hq
      change f (⟨i, x⟩, ⟨j, y⟩) ∈ s
      simpa only [heq] using hq
  rw [hdecomp]
  exact MeasurableSet.iUnion fun i => MeasurableSet.iUnion fun j =>
    ((measurableEmbedding_sigma_mk i).prodMap (measurableEmbedding_sigma_mk j)).measurableSet_image'
      (hf i j hs)

/-- On countable measurable disjoint unions, joint measurability is equivalent to
joint measurability on each pair of components. -/
theorem measurable_sigma_prod_iff [Countable I] [Countable J]
    {Z : Type*} [MeasurableSpace Z] {f : Sigma X × Sigma Y → Z} :
    Measurable f ↔
      ∀ i j, Measurable (fun p : X i × Y j => f (⟨i, p.1⟩, ⟨j, p.2⟩)) := by
  constructor
  · intro hf i j
    exact hf.comp ((measurable_sigma_mk i).prodMap (measurable_sigma_mk j))
  · exact measurable_of_sigma_prod

end CertifiedGaps.SigmaProductMeasurable
