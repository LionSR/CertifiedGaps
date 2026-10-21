import CertifiedGaps.Transfer.PoissonHistory
import CertifiedGaps.Transfer.SigmaProductMeasurable
import Mathlib.MeasureTheory.Measure.Prod

/-!
# Null collisions between histories

Every ordered event-time coordinate has an atomless marginal. Fubini and the countable
history decomposition show that two independently integrated histories almost surely have
distinct merged event times, including all empty-history cases.
-/

open MeasureTheory Set
open scoped ENNReal NNReal

namespace CertifiedGaps.PoissonHistory

open TimeOrdering SigmaProductMeasurable

/-- Prescribing a single ordered event time cuts out a Lebesgue-null hyperplane. -/
theorem orderedTimesMeasure_coordinate_eq_zero (b : ℝ≥0) (n : ℕ)
    (i : Fin n) (a : ℝ) :
    orderedTimesMeasure b n {t | t.val i = a} = 0 := by
  rw [orderedTimesMeasure, comap_subtype_coe_apply (measurableSet_timeSimplex b n)]
  exact measure_mono_null (by rintro x ⟨t, ht, rfl⟩; exact ht)
    (Measure.pi_hyperplane (fun _ : Fin n => (volume : Measure ℝ)) i a)

/-- Every simplex coordinate has an atomless marginal, including the empty simplex. -/
theorem ae_orderedTimes_coordinate_ne (b : ℝ≥0) (n : ℕ) (i : Fin n) (a : ℝ) :
    ∀ᵐ t ∂orderedTimesMeasure b n, t.val i ≠ a :=
  compl_mem_ae_iff.mpr (orderedTimesMeasure_coordinate_eq_zero b n i a)

/-- A fixed finite history almost surely has no event times in common with a second history. -/
theorem ae_history_times_disjoint (b : ℝ≥0) (m : ℕ) (h : History b m) :
    ∀ᵐ k ∂historyMeasure b m, ∀ i j, h.2.val i ≠ k.2.val j := by
  rw [historyMeasure, Measure.ae_sum_iff]
  intro w
  rw [(measurableEmbedding_sigma_mk w).ae_map_iff]
  change ∀ᵐ k ∂orderedTimesMeasure b w.1,
    ∀ i : Fin h.1.1, ∀ j : Fin w.1, h.2.val i ≠ k.val j
  apply ae_all_iff.mpr
  intro i
  apply ae_all_iff.mpr
  intro j
  exact (ae_orderedTimes_coordinate_ne b w.1 j (h.2.val i)).mono fun _ ht => ht.symm

/-- Collision-free pairs form a measurable subset of the product history space. -/
theorem measurableSet_merged_history_injective (b : ℝ≥0) (m : ℕ) :
    MeasurableSet {p : History b m × History b m |
      Function.Injective (Sum.elim p.1.2.val p.2.2.val)} := by
  rw [measurableSet_setOfPred]
  apply measurable_of_sigma_prod
  intro w v
  rw [← measurableSet_setOfPred]
  have hm : MeasurableSet {p : OrderedTimes b w.1 × OrderedTimes b v.1 |
      ∀ i j, p.1.val i ≠ p.2.val j} := by
    simp only [ofPred_forall]
    exact MeasurableSet.iInter fun i => MeasurableSet.iInter fun j =>
      (measurableSet_eq_fun
        ((measurable_pi_apply i).comp (measurable_subtype_coe.comp measurable_fst))
        ((measurable_pi_apply j).comp (measurable_subtype_coe.comp measurable_snd))).compl
  convert hm using 1
  ext p
  simp only [mem_ofPred_eq, Sum.elim_injective]
  exact ⟨fun h => h.2.2, fun h => ⟨p.1.property.2.injective, p.2.property.2.injective, h⟩⟩

/-- Independently integrated finite histories almost surely have distinct merged event times. -/
theorem ae_merged_history_injective (b : ℝ≥0) (m : ℕ) :
    ∀ᵐ p ∂(historyMeasure b m).prod (historyMeasure b m),
      Function.Injective (Sum.elim p.1.2.val p.2.2.val) := by
  apply (Measure.ae_prod_iff_ae_ae (measurableSet_merged_history_injective b m)).mpr
  apply ae_of_all
  intro h
  exact (ae_history_times_disjoint b m h).mono fun k hk =>
    h.2.property.2.injective.sumElim k.2.property.2.injective hk

end CertifiedGaps.PoissonHistory
