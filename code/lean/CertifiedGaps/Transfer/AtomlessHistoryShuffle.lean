import CertifiedGaps.Transfer.AtomlessTimeOrdering

/-!
# Finite-family shuffles for atomless time measures

Chronological interleaving is a product-measure change of variables for any atomless time law.
The event order remains inside the integrand, so these identities apply to noncommuting values.
-/

open MeasureTheory Set

namespace CertifiedGaps.AtomlessTimeOrdering

open TimeOrdering (PreservesTimeOrder TimeShuffle historyEventEquiv splitFamilyHistoryTimes
  splitFamilyHistoryTimes_apply familyHistoryOrder)

/-- An arbitrary measurable time domain with finitely many strict precedence constraints. -/
def orderedTimeDomain (s : Set ℝ) {n : ℕ} (r : Fin n → Fin n → Prop) :
    Set (Fin n → ℝ) := timeCube s n ∩ {t | ∀ i j, r i j → t i < t j}

theorem measurableSet_orderedTimeDomain {s : Set ℝ} (hs : MeasurableSet s) {n : ℕ}
    (r : Fin n → Fin n → Prop) : MeasurableSet (orderedTimeDomain s r) := by
  refine (measurableSet_timeCube hs n).inter ?_
  simp only [ofPred_forall]
  exact MeasurableSet.iInter fun i => MeasurableSet.iInter fun j =>
    MeasurableSet.iInter fun _ => measurableSet_lt (measurable_pi_apply i)
      (measurable_pi_apply j)

theorem mem_orderedTimeDomain_iff {s : Set ℝ} {n : ℕ} {r : Fin n → Fin n → Prop}
    {σ : Equiv.Perm (Fin n)} {t : Fin n → ℝ} (ht : t ∈ timeChamber s n σ) :
    t ∈ orderedTimeDomain s r ↔ PreservesTimeOrder r σ := by
  constructor
  · intro h i j hij
    apply ht.2.lt_iff_lt.mp
    simpa only [Function.comp_apply, Equiv.apply_symm_apply] using h.2 i j hij
  · intro h
    refine ⟨ht.1, fun i j hij => ?_⟩
    simpa only [Function.comp_apply, Equiv.apply_symm_apply] using ht.2 (h i j hij)

theorem timeChamber_subset_orderedTimeDomain {s : Set ℝ} {n : ℕ}
    {r : Fin n → Fin n → Prop} (σ : TimeShuffle r) :
    timeChamber s n σ.val ⊆ orderedTimeDomain s r :=
  fun _ ht => (mem_orderedTimeDomain_iff ht).mpr σ.property

/-- The admissible shuffles cover the constrained domain outside null coordinate collisions. -/
theorem iUnion_timeShuffle_ae_eq (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    (s : Set ℝ) {n : ℕ} (r : Fin n → Fin n → Prop) :
    (⋃ σ : TimeShuffle r, timeChamber s n σ.val) =ᵐ[timeMeasure ν n]
      orderedTimeDomain s r := by
  filter_upwards [ae_injective_times ν n] with t ht
  apply propext
  constructor
  · intro h
    obtain ⟨σ, hσ⟩ := Set.mem_iUnion.mp h
    exact timeChamber_subset_orderedTimeDomain σ hσ
  · intro h
    obtain ⟨σ, hσ⟩ := Set.mem_iUnion.mp (mem_iUnion_timeChamber h.1 ht)
    exact Set.mem_iUnion.mpr ⟨⟨σ, (mem_orderedTimeDomain_iff hσ).mp h⟩, hσ⟩

/-- A constrained time integral is the sum over chronological orders satisfying the constraints. -/
theorem integral_orderedTimeDomain_eq_sum_timeSimplex
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    {s : Set ℝ} (hs : MeasurableSet s) {n : ℕ} {r : Fin n → Fin n → Prop}
    {f : (Fin n → ℝ) → E} (hf : IntegrableOn f (orderedTimeDomain s r) (timeMeasure ν n)) :
    ∫ t in orderedTimeDomain s r, f t ∂timeMeasure ν n = ∑ σ : TimeShuffle r,
      ∫ t in timeSimplex s n, f (t ∘ σ.val.symm) ∂timeMeasure ν n := by
  classical
  rw [← setIntegral_congr_set (iUnion_timeShuffle_ae_eq ν s r)]
  have hdis : Pairwise (fun σ τ : TimeShuffle r =>
      Disjoint (timeChamber s n σ.val) (timeChamber s n τ.val)) := by
    intro σ τ hστ
    exact pairwise_disjoint_timeChamber s n (fun h => hστ (Subtype.ext h))
  rw [integral_iUnion_fintype (fun σ => measurableSet_timeChamber hs n σ.val) hdis
    (fun σ => hf.mono_set (timeChamber_subset_orderedTimeDomain σ))]
  exact Finset.sum_congr rfl fun σ _ => integral_timeChamber ν s n σ.val f

section FiniteFamily

variable {ι : Type*} [Fintype ι] (counts : ι → ℕ)

/-- Flattening independent finite time tuples preserves the actual product measure. -/
theorem measurePreserving_uncurryHistoryTimes (ν : Measure ℝ) [SigmaFinite ν] :
    MeasurePreserving (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ)).symm
      (Measure.pi fun i => timeMeasure ν (counts i))
      (Measure.pi fun _ : Σ i, Fin (counts i) => ν) where
  measurable := (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ)).symm.measurable
  map_eq := by
    refine (Measure.pi_eq fun u _ => ?_).symm
    rw [MeasurableEquiv.map_apply]
    have hpre : (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ)).symm ⁻¹'
        (Set.univ.pi u) = Set.univ.pi (fun i => Set.univ.pi (fun j => u ⟨i, j⟩)) := by
      ext x
      simp [Set.mem_pi, MeasurableEquiv.coe_piCurry_symm, Sigma.forall, Sigma.uncurry]
    rw [hpre, Measure.pi_pi]
    simp only [timeMeasure, Measure.pi_pi]
    rw [Fintype.prod_sigma]

/-- Splitting globally indexed events into their original histories is measure preserving. -/
theorem measurePreserving_splitFamilyHistoryTimes (ν : Measure ℝ) [SigmaFinite ν] :
    MeasurePreserving (splitFamilyHistoryTimes counts) (timeMeasure ν (∑ i, counts i))
      (Measure.pi fun i => timeMeasure ν (counts i)) :=
  (measurePreserving_uncurryHistoryTimes counts ν).symm.comp
    (MeasureTheory.measurePreserving_piCongrLeft (fun _ : Σ i, Fin (counts i) => ν)
      (historyEventEquiv counts).symm)

/-- The product of the chronological simplexes is precisely the family precedence domain. -/
theorem preimage_pi_timeSimplex (s : Set ℝ) :
    splitFamilyHistoryTimes counts ⁻¹'
        (Set.univ.pi fun i => timeSimplex s (counts i)) =
      orderedTimeDomain s (familyHistoryOrder counts) := by
  ext t
  simp only [Set.mem_preimage, Set.mem_pi, Set.mem_univ, forall_true_left]
  constructor
  · intro h
    refine ⟨?_, ?_⟩
    · intro p _
      obtain ⟨⟨i, j⟩, rfl⟩ := (historyEventEquiv counts).surjective p
      simpa only [splitFamilyHistoryTimes_apply] using (h i).1 j (Set.mem_univ _)
    · rintro p q ⟨i, j, k, hjk, rfl, rfl⟩
      simpa only [splitFamilyHistoryTimes_apply] using (h i).2 hjk
  · rintro ⟨hc, ho⟩ i
    refine ⟨?_, ?_⟩
    · intro j _
      simpa only [splitFamilyHistoryTimes_apply] using hc _ (Set.mem_univ _)
    · intro j k hjk
      simpa only [splitFamilyHistoryTimes_apply] using ho _ _ ⟨i, j, k, hjk, rfl, rfl⟩

/-- Every finite product of ordered time domains is the sum of its genuine global chronological
interleavings, for an arbitrary atomless time measure and arbitrary Banach-valued integrands. -/
theorem integral_pi_timeSimplex_eq_sum_shuffles
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    {s : Set ℝ} (hs : MeasurableSet s)
    {f : ((i : ι) → Fin (counts i) → ℝ) → E}
    (hf : IntegrableOn f (Set.univ.pi fun i => timeSimplex s (counts i))
      (Measure.pi fun i => timeMeasure ν (counts i))) :
    ∫ p in Set.univ.pi (fun i => timeSimplex s (counts i)), f p
      ∂(Measure.pi fun i => timeMeasure ν (counts i)) =
      ∑ σ : TimeShuffle (familyHistoryOrder counts),
        ∫ t in timeSimplex s (∑ i, counts i),
          f (fun i j => t (σ.val.symm (historyEventEquiv counts ⟨i, j⟩)))
            ∂timeMeasure ν (∑ i, counts i) := by
  have hmp := measurePreserving_splitFamilyHistoryTimes counts ν
  have hfi := (hmp.integrableOn_comp_preimage
    (splitFamilyHistoryTimes counts).measurableEmbedding).mpr hf
  rw [preimage_pi_timeSimplex] at hfi
  have h := hmp.setIntegral_preimage_emb (splitFamilyHistoryTimes counts).measurableEmbedding f
    (Set.univ.pi fun i => timeSimplex s (counts i))
  rw [preimage_pi_timeSimplex] at h
  rw [← h]
  convert integral_orderedTimeDomain_eq_sum_timeSimplex ν hs hfi using 1
  · rfl
  · apply Finset.sum_congr rfl
    intro σ _
    congr 1
    funext t
    congr 1
    funext i j
    simp only [splitFamilyHistoryTimes_apply, Function.comp_apply]

end FiniteFamily

end CertifiedGaps.AtomlessTimeOrdering
