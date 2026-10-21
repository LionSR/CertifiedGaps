import CertifiedGaps.Transfer.WeightedPoissonHistory
import CertifiedGaps.Transfer.AtomlessHistoryShuffle
import Mathlib.MeasureTheory.Integral.Pi

/-!
# Finite products of labeled Poisson histories

The product history measure is decomposed along its countable measurable partition by
bondwise label words. This gives an absolutely convergent integral expansion over actual
product ordered-simplex measures; chronological shuffle identities then apply to each term.
-/

open MeasureTheory Set Function
open scoped ENNReal NNReal

namespace CertifiedGaps.WeightedPoissonHistory

open AtomlessTimeOrdering

/-- The component with one specified complete event word. -/
def wordFiber {s : TimeDomain} {A : Type*} (w : EventWord A) : Set (History s A) :=
  {h | h.1 = w}

theorem measurableSet_wordFiber {s : TimeDomain} {A : Type*} (w : EventWord A) :
    MeasurableSet (wordFiber (s := s) w) := by
  apply MeasurableSpace.measurableSet_iInf.mpr
  intro v
  change MeasurableSet ((fun _ : OrderedTimes s v.1 => v) ⁻¹' {w})
  by_cases h : v = w <;> simp [h]

/-- Restricting the history measure to a word recovers its actual ordered-simplex measure. -/
theorem restrict_historyMeasure_wordFiber {A : Type*} (ν : Measure ℝ) (s : TimeDomain)
    (q : A → ℝ≥0) (w : EventWord A) :
    (historyMeasure ν s q).restrict (wordFiber w) =
      (wordTimesMeasure ν s q w).map (Sigma.mk w) := by
  classical
  rw [historyMeasure, Measure.restrict_sum _ (measurableSet_wordFiber w)]
  have hterm (v : EventWord A) :
      ((wordTimesMeasure ν s q v).map (Sigma.mk v)).restrict (wordFiber w) =
        if v = w then (wordTimesMeasure ν s q w).map (Sigma.mk w) else 0 := by
    rw [Measure.restrict_map (measurable_historyMk s v) (measurableSet_wordFiber w)]
    by_cases h : v = w
    · subst v
      simp [wordFiber]
    · have he : Sigma.mk v ⁻¹' wordFiber (s := s) w = ∅ := by
        ext t
        simp [wordFiber, h]
      simp [he, h]
  simp_rw [hterm]
  ext s hs
  rw [Measure.sum_apply _ hs, tsum_eq_single w]
  · simp
  · intro v hv
    simp [hv]

section FiniteFamily

variable {ι : Type*} [Fintype ι] (ν : Measure ℝ) [IsFiniteMeasure ν]
  [NullSingletonClass ν] (s : TimeDomain) (A : ι → Type*) [∀ i, Finite (A i)]
  (q : (i : ι) → A i → ℝ≥0)

/-- The actual finite product of the bond-history measures. -/
noncomputable def historyFamilyMeasure : Measure ((i : ι) → History s (A i)) :=
  Measure.pi (fun i => historyMeasure ν s (q i))

/-- Insert one tuple of ordered times into a specified family of words. -/
def historyFamilyMk (w : (i : ι) → EventWord (A i))
    (t : (i : ι) → OrderedTimes s (w i).1) (i : ι) : History s (A i) :=
  ⟨w i, t i⟩

omit [Fintype ι] [∀ i, Finite (A i)] in
theorem measurable_historyFamilyMk (w : (i : ι) → EventWord (A i)) :
    Measurable (historyFamilyMk s A w) := by
  exact Measurable.of_eval fun i =>
    (measurable_historyMk s (w i)).comp (measurable_pi_apply i)

/-- Product-simplex measure for a specified family of label words. -/
noncomputable def familyTimesMeasure (w : (i : ι) → EventWord (A i)) :
    Measure ((i : ι) → OrderedTimes s (w i).1) :=
  Measure.pi (fun i => wordTimesMeasure ν s (q i) (w i))

/-- Each word-vector part of the finite product is a pushforward of the product of the
ordered-time component measures. -/
theorem restrict_historyFamilyMeasure (w : (i : ι) → EventWord (A i)) :
    (historyFamilyMeasure ν s A q).restrict (Set.univ.pi fun i => wordFiber (w i)) =
      (familyTimesMeasure ν s A q w).map (historyFamilyMk s A w) := by
  rw [historyFamilyMeasure, Measure.restrict_pi_pi]
  simp_rw [restrict_historyMeasure_wordFiber]
  exact (Measure.pi_map_pi (fun i =>
    (measurable_historyMk s (w i)).aemeasurable)).symm

/-- A finite product of countable history disjoint unions is genuinely the countable sum
of its product-simplex components. This identity is proved from the measurable partition,
not imposed as a hypothesis. -/
theorem historyFamilyMeasure_eq_sum :
    historyFamilyMeasure ν s A q = Measure.sum fun w : (i : ι) → EventWord (A i) =>
      (familyTimesMeasure ν s A q w).map (historyFamilyMk s A w) := by
  classical
  have hcover : (⋃ w : (i : ι) → EventWord (A i),
      Set.univ.pi fun i => wordFiber (s := s) (w i)) = Set.univ := by
    ext h
    simp only [Set.mem_iUnion, Set.mem_univ, iff_true]
    exact ⟨fun i => (h i).1, fun i _ => rfl⟩
  have hdisj : Pairwise (Disjoint on (fun w : (i : ι) → EventWord (A i) =>
      Set.univ.pi fun i => wordFiber (s := s) (w i))) := by
    intro w v hwv
    apply Set.disjoint_left.mpr
    intro h hw hv
    apply hwv
    funext i
    exact (hw i (Set.mem_univ i)).symm.trans (hv i (Set.mem_univ i))
  have hm (w : (i : ι) → EventWord (A i)) :
      MeasurableSet (Set.univ.pi fun i => wordFiber (s := s) (w i)) :=
    MeasurableSet.univ_pi fun i => measurableSet_wordFiber (w i)
  have h := Measure.restrict_iUnion (μ := historyFamilyMeasure ν s A q) hdisj hm
  rw [hcover, Measure.restrict_univ] at h
  simpa only [restrict_historyFamilyMeasure] using h


/-- The product-history integral expands into product ordered-simplex integrals. -/
theorem hasSum_integral_historyFamilyMeasure
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    HasSum (fun w : (i : ι) → EventWord (A i) =>
      ∫ t, f (historyFamilyMk s A w t) ∂familyTimesMeasure ν s A q w)
      (∫ h, f h ∂historyFamilyMeasure ν s A q) := by
  rw [historyFamilyMeasure_eq_sum] at hf ⊢
  have h := hasSum_integral_measure hf
  convert h using 1
  funext w
  exact (integral_map (measurable_historyFamilyMk s A w).aemeasurable
    (hf.mono_measure (Measure.le_sum _ w)).aestronglyMeasurable).symm

/-- Absolute convergence is proved at the level of the integrals of the norms, which also
justifies every subsequent finite regrouping and chronological subdivision. -/
theorem summable_integral_norm_historyFamilyMeasure
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    Summable (fun w : (i : ι) → EventWord (A i) =>
      ∫ t, ‖f (historyFamilyMk s A w t)‖ ∂familyTimesMeasure ν s A q w) := by
  rw [historyFamilyMeasure_eq_sum] at hf
  convert hf.summable_integral using 1
  funext w
  exact (integral_map (measurable_historyFamilyMk s A w).aemeasurable
    (hf.mono_measure (Measure.le_sum _ w)).norm.aestronglyMeasurable).symm

/-- Separate bondwise complete words into their count vector and finite label-word family. -/
def familyWordsEquiv : ((i : ι) → EventWord (A i)) ≃
    Σ counts : ι → ℕ, (i : ι) → Fin (counts i) → A i where
  toFun w := ⟨fun i => (w i).1, fun i => (w i).2⟩
  invFun v i := ⟨v.1 i, v.2 i⟩
  left_inv w := by rfl
  right_inv v := by cases v; rfl

/-- Count-vector and label-word expansion of the actual finite-product history integral. -/
theorem hasSum_integral_historyFamilyMeasure_by_counts [∀ i, Fintype (A i)] [DecidableEq ι]
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    HasSum (fun counts : ι → ℕ => ∑ a : (i : ι) → Fin (counts i) → A i,
      ∫ t, f (historyFamilyMk s A (fun i => ⟨counts i, a i⟩) t)
        ∂familyTimesMeasure ν s A q (fun i => ⟨counts i, a i⟩))
      (∫ h, f h ∂historyFamilyMeasure ν s A q) := by
  have h := ((familyWordsEquiv A).symm.hasSum_iff).mpr
    (hasSum_integral_historyFamilyMeasure ν s A q hf)
  exact h.sigma fun _ => hasSum_fintype _

/-- Product exponential event-count growth is integrable for any finite bond family. -/
theorem integrable_family_eventCount_pow {R : ι → ℝ} (hR : ∀ i, 0 ≤ R i) :
    Integrable (fun h : (i : ι) → History s (A i) =>
      ∏ i, (R i) ^ eventCount (h i)) (historyFamilyMeasure ν s A q) :=
  Integrable.fintype_prod_dep fun i => integrable_eventCount_pow (ν := ν) (s := s) (w := q i) (hR i)

/-- Exact product exponential moment, obtained from the actual finite product measure. -/
theorem integral_family_eventCount_pow [∀ i, Fintype (A i)] {R : ι → ℝ} (hR : ∀ i, 0 ≤ R i) :
    ∫ h : (i : ι) → History s (A i), ∏ i, (R i) ^ eventCount (h i)
      ∂historyFamilyMeasure ν s A q =
        ∏ i, Real.exp ((totalIntensity (q i) : ℝ) * timeMass ν s * R i) := by
  rw [historyFamilyMeasure]
  calc
    _ = ∏ i, ∫ h : History s (A i), R i ^ eventCount h ∂historyMeasure ν s (q i) :=
      integral_fintype_prod_eq_prod (fun i (h : History s (A i)) => R i ^ eventCount h)
    _ = _ := Finset.prod_congr rfl fun i _ =>
      integral_eventCount_pow (ν := ν) (s := s) (w := q i) (hR i)

/-- An explicit bondwise event-count bound supplies the absolute integrability required by
the product-history series, without assuming convergence of that series. -/
theorem integrable_of_family_eventCount_bound
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E} {C : ℝ} {R : ι → ℝ}
    (hf : AEStronglyMeasurable f (historyFamilyMeasure ν s A q)) (hR : ∀ i, 0 ≤ R i)
    (hbound : ∀ᵐ h ∂historyFamilyMeasure ν s A q,
      ‖f h‖ ≤ C * ∏ i, (R i) ^ eventCount (h i)) :
    Integrable f (historyFamilyMeasure ν s A q) :=
  ((integrable_family_eventCount_pow ν s A q hR).const_mul C).mono' hf hbound



/-- The product intensity of a family of complete words. -/
noncomputable def familyWeight (w : (i : ι) → EventWord (A i)) : ℝ≥0 :=
  ∏ i, wordWeight (q i) (w i).2

omit [∀ i, Finite (A i)] in
/-- The weighted product component is exactly its total word weight times the unweighted
product-simplex measure. Zero-weight words are included without any positivity assumption. -/
theorem familyTimesMeasure_eq_smul (w : (i : ι) → EventWord (A i)) :
    familyTimesMeasure ν s A q w = familyWeight A q w •
      (Measure.pi fun i => orderedTimesMeasure ν s (w i).1) := by
  apply Measure.pi_eq
  intro u _
  simp only [Measure.smul_apply, Measure.pi_pi, wordTimesMeasure,
    familyWeight, ENNReal.smul_def, smul_eq_mul, ENNReal.ofNNReal_finsetProd,
    Finset.prod_mul_distrib]

/-- Forget the proofs that every history's time tuple is chronological. -/
def familyTimesCoe (counts : ι → ℕ) (t : (i : ι) → OrderedTimes s (counts i))
    (i : ι) : Fin (counts i) → ℝ := (t i).val

omit [Fintype ι] [∀ i, Finite (A i)] in
/-- The product subtype inclusion is a measurable embedding. -/
theorem measurableEmbedding_familyTimesCoe [Finite ι] (counts : ι → ℕ) :
    MeasurableEmbedding (familyTimesCoe s counts) := by
  have hs : MeasurableSet {t : (i : ι) → Fin (counts i) → ℝ |
      ∀ i, t i ∈ timeSimplex s.val (counts i)} := by
    simp only [Set.ofPred_forall]
    exact MeasurableSet.iInter fun i =>
      (measurableSet_timeSimplex s.property (counts i)).preimage (measurable_pi_apply i)
  exact (MeasurableEmbedding.subtype_coe hs).comp
    (MeasurableEquiv.subtypePiEquivPi
      (p := fun i t => t ∈ timeSimplex s.val (counts i))).symm.measurableEmbedding

/-- The unweighted product component maps to the actual product time law restricted to its
product of chronological simplexes. -/
theorem measurePreserving_familyTimesCoe (counts : ι → ℕ) :
    MeasurePreserving (familyTimesCoe s counts)
      (Measure.pi fun i => orderedTimesMeasure ν s (counts i))
      ((Measure.pi fun i => timeMeasure ν (counts i)).restrict
        (Set.univ.pi fun i => timeSimplex s.val (counts i))) where
  measurable := (measurableEmbedding_familyTimesCoe s counts).measurable
  map_eq := by
    change (Measure.pi fun i => orderedTimesMeasure ν s (counts i)).map
      (fun t i => (t i).val) = _
    have (i : ι) : SigmaFinite ((orderedTimesMeasure ν s (counts i)).map
        (Subtype.val : OrderedTimes s (counts i) → (Fin (counts i) → ℝ))) := inferInstance
    rw [Measure.pi_map_pi (X := fun i => OrderedTimes s (counts i))
      (Y := fun i => Fin (counts i) → ℝ)
      (μ := fun i => orderedTimesMeasure ν s (counts i))
      (f := fun i => (Subtype.val : OrderedTimes s (counts i) → (Fin (counts i) → ℝ)))
      (fun i => measurable_subtype_coe.aemeasurable)]
    simp_rw [orderedTimesMeasure,
      map_comap_subtype_coe (measurableSet_timeSimplex s.property _)]
    exact (Measure.restrict_pi_pi _ _).symm

/-- Extend a word-family observable to ambient time coordinates by zero outside its simplexes. -/
noncomputable def familyObservable
    {E : Type*} [Zero E] (f : ((i : ι) → History s (A i)) → E)
    (w : (i : ι) → EventWord (A i)) (t : (i : ι) → Fin (w i).1 → ℝ) : E := by
  classical
  exact if ht : ∀ i, t i ∈ timeSimplex s.val (w i).1 then
    f (historyFamilyMk s A w (fun i => ⟨t i, ht i⟩)) else 0

omit [∀ i, Finite (A i)] in
@[simp]
theorem familyObservable_familyTimesCoe
    {E : Type*} [Zero E] (f : ((i : ι) → History s (A i)) → E)
    (w : (i : ι) → EventWord (A i)) (t : (i : ι) → OrderedTimes s (w i).1) :
    familyObservable s A f w (familyTimesCoe s (fun i => (w i).1) t) =
      f (historyFamilyMk s A w t) := by
  simp [familyObservable, familyTimesCoe, fun i => (t i).property]

/-- Integrability descends to every weighted product-simplex component. -/
theorem integrable_historyFamilyMk {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) (w : (i : ι) → EventWord (A i)) :
    Integrable (fun t => f (historyFamilyMk s A w t)) (familyTimesMeasure ν s A q w) := by
  rw [historyFamilyMeasure_eq_sum] at hf
  exact (hf.mono_measure (Measure.le_sum _ w)).comp_measurable
    (measurable_historyFamilyMk s A w)

/-- A component of nonzero weight has an integrable ambient extension. No positive-weight
assumption is imposed on the label system; zero-weight components are handled separately. -/
theorem integrableOn_familyObservable_of_weight_ne_zero
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) (w : (i : ι) → EventWord (A i))
    (hw : familyWeight A q w ≠ 0) :
    IntegrableOn (familyObservable s A f w)
      (Set.univ.pi fun i => timeSimplex s.val (w i).1)
      (Measure.pi fun i => timeMeasure ν (w i).1) := by
  apply ((measurePreserving_familyTimesCoe ν s (fun i => (w i).1)).integrable_comp_emb
    (measurableEmbedding_familyTimesCoe s (fun i => (w i).1))).mp
  have hi := integrable_historyFamilyMk ν s A q hf w
  rw [familyTimesMeasure_eq_smul, ENNReal.smul_def] at hi
  have hb := (integrable_smul_measure (by exact_mod_cast hw) ENNReal.coe_ne_top).mp hi
  simpa only [Function.comp_def, familyObservable_familyTimesCoe] using hb

/-- Each weighted product-simplex integral becomes its exact chronological shuffle sum.
This includes zero-intensity words, without assuming their unweighted observables integrable. -/
theorem integral_familyTimesMeasure_eq_sum_shuffles
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) (w : (i : ι) → EventWord (A i)) :
    ∫ t, f (historyFamilyMk s A w t) ∂familyTimesMeasure ν s A q w =
      ∑ σ : TimeOrdering.TimeShuffle (TimeOrdering.familyHistoryOrder (fun i => (w i).1)),
        (familyWeight A q w : ℝ) •
          ∫ t in timeSimplex s.val (∑ i, (w i).1),
            familyObservable s A f w
              (fun i j => t (σ.val.symm
                (TimeOrdering.historyEventEquiv (fun i => (w i).1) ⟨i, j⟩)))
            ∂timeMeasure ν (∑ i, (w i).1) := by
  classical
  rw [familyTimesMeasure_eq_smul, integral_smul_nnreal_measure, NNReal.smul_def]
  by_cases hw : familyWeight A q w = 0
  · simp [hw]
  have he := (measurePreserving_familyTimesCoe ν s (fun i => (w i).1)).integral_comp
    (measurableEmbedding_familyTimesCoe s (fun i => (w i).1)) (familyObservable s A f w)
  simp only [familyObservable_familyTimesCoe] at he
  rw [he, integral_pi_timeSimplex_eq_sum_shuffles _ ν s.property
    (integrableOn_familyObservable_of_weight_ne_zero ν s A q hf w hw), Finset.smul_sum]

/-- An explicit product exponential majorant bounds the finite-product history L¹ norm. -/
theorem integral_norm_le_of_family_eventCount_bound [∀ i, Fintype (A i)]
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E} {C : ℝ} {R : ι → ℝ}
    (hf : AEStronglyMeasurable f (historyFamilyMeasure ν s A q)) (hR : ∀ i, 0 ≤ R i)
    (hbound : ∀ᵐ h ∂historyFamilyMeasure ν s A q,
      ‖f h‖ ≤ C * ∏ i, (R i) ^ eventCount (h i)) :
    ∫ h, ‖f h‖ ∂historyFamilyMeasure ν s A q ≤
      C * ∏ i, Real.exp ((totalIntensity (q i) : ℝ) * timeMass ν s * R i) := by
  calc
    _ ≤ ∫ h, C * ∏ i, (R i) ^ eventCount (h i) ∂historyFamilyMeasure ν s A q :=
      integral_mono_ae (integrable_of_family_eventCount_bound ν s A q hf hR hbound).norm
        ((integrable_family_eventCount_pow ν s A q hR).const_mul C) hbound
    _ = _ := by rw [integral_const_mul, integral_family_eventCount_pow ν s A q hR]

/-- The finite-product history integral is an absolutely convergent count, label, and
chronological-shuffle expansion with the exact product of event intensities. -/
theorem hasSum_integral_historyFamilyMeasure_by_counts_shuffles
    [∀ i, Fintype (A i)] [DecidableEq ι]
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    HasSum (fun counts : ι → ℕ => ∑ a : (i : ι) → Fin (counts i) → A i,
      ∑ σ : TimeOrdering.TimeShuffle (TimeOrdering.familyHistoryOrder counts),
        (familyWeight A q (fun i => ⟨counts i, a i⟩) : ℝ) •
          ∫ t in timeSimplex s.val (∑ i, counts i),
            familyObservable s A f (fun i => ⟨counts i, a i⟩)
              (fun i j => t (σ.val.symm (TimeOrdering.historyEventEquiv counts ⟨i, j⟩)))
            ∂timeMeasure ν (∑ i, counts i))
      (∫ h, f h ∂historyFamilyMeasure ν s A q) := by
  have h := hasSum_integral_historyFamilyMeasure_by_counts ν s A q hf
  convert h using 1
  funext counts
  apply Finset.sum_congr rfl
  intro a _
  exact (integral_familyTimesMeasure_eq_sum_shuffles ν s A q hf
    (fun i => ⟨counts i, a i⟩)).symm

omit [∀ i, Finite (A i)] in
@[simp]
theorem familyObservable_norm {E : Type*} [NormedAddCommGroup E]
    (f : ((i : ι) → History s (A i)) → E) (w : (i : ι) → EventWord (A i))
    (t : (i : ι) → Fin (w i).1 → ℝ) :
    familyObservable s A (fun h => ‖f h‖) w t = ‖familyObservable s A f w t‖ := by
  classical
  unfold familyObservable
  split_ifs <;> simp

/-- One index chooses the counts, every label, and one admissible chronological interleaving. -/
abbrev FamilyShuffleIndex := Σ counts : ι → ℕ,
  ((i : ι) → Fin (counts i) → A i) ×
    TimeOrdering.TimeShuffle (TimeOrdering.familyHistoryOrder counts)

/-- One fully ungrouped chronological term, carrying its actual product intensity. -/
noncomputable def shuffledIntegral
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (f : ((i : ι) → History s (A i)) → E) (v : FamilyShuffleIndex A) : E :=
  (familyWeight A q (fun i => ⟨v.1 i, v.2.1 i⟩) : ℝ) •
    ∫ t in timeSimplex s.val (∑ i, v.1 i),
      familyObservable s A f (fun i => ⟨v.1 i, v.2.1 i⟩)
        (fun i j => t (v.2.2.val.symm (TimeOrdering.historyEventEquiv v.1 ⟨i, j⟩)))
      ∂timeMeasure ν (∑ i, v.1 i)

/-- Absolute convergence holds for every individual chronological norm integral. -/
theorem summable_shuffledIntegral_norm {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    Summable (shuffledIntegral ν s A q (fun h => ‖f h‖)) := by
  classical
  let (i : ι) : Fintype (A i) := Fintype.ofFinite (A i)
  have hn (v : FamilyShuffleIndex A) :
      0 ≤ shuffledIntegral ν s A q (fun h => ‖f h‖) v := by
    unfold shuffledIntegral
    rw [smul_eq_mul]
    apply mul_nonneg (NNReal.coe_nonneg _)
    exact integral_nonneg fun t => by rw [familyObservable_norm]; exact norm_nonneg _
  apply (summable_sigma_of_nonneg hn).mpr
  refine ⟨fun _ => (hasSum_fintype _).summable, ?_⟩
  have h := (hasSum_integral_historyFamilyMeasure_by_counts_shuffles ν s A q hf.norm).summable
  simpa only [tsum_fintype, Fintype.sum_prod_type, shuffledIntegral] using h

/-- The norms of the individual Banach-valued chronological terms are summable. -/
theorem summable_norm_shuffledIntegral
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    Summable (fun v => ‖shuffledIntegral ν s A q f v‖) := by
  apply Summable.of_nonneg_of_le (fun _ => norm_nonneg _)
    (fun v => ?_) (summable_shuffledIntegral_norm ν s A q hf)
  simp only [shuffledIntegral, norm_smul, Real.norm_eq_abs,
    abs_of_nonneg (NNReal.coe_nonneg _), familyObservable_norm, smul_eq_mul]
  exact mul_le_mul_of_nonneg_left (norm_integral_le_integral_norm _) (NNReal.coe_nonneg _)

/-- Fully ungrouped weighted chronological expansion of the actual finite-product measure. -/
theorem hasSum_shuffledIntegral
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    HasSum (shuffledIntegral ν s A q f) (∫ h, f h ∂historyFamilyMeasure ν s A q) := by
  classical
  let (i : ι) : Fintype (A i) := Fintype.ofFinite (A i)
  have hs := (summable_norm_shuffledIntegral ν s A q hf).of_norm
  apply hs.hasSum_iff.mpr
  rw [hs.tsum_sigma]
  simpa only [tsum_fintype, Fintype.sum_prod_type, shuffledIntegral] using
    (hasSum_integral_historyFamilyMeasure_by_counts_shuffles ν s A q hf).tsum_eq

/-- The sum of all individual chronological norm integrals is exactly the original L¹ norm. -/
theorem tsum_shuffledIntegral_norm {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History s (A i)) → E}
    (hf : Integrable f (historyFamilyMeasure ν s A q)) :
    ∑' v, shuffledIntegral ν s A q (fun h => ‖f h‖) v =
      ∫ h, ‖f h‖ ∂historyFamilyMeasure ν s A q :=
  (hasSum_shuffledIntegral ν s A q hf.norm).tsum_eq

/-- An explicit exponential bound controls the absolutely convergent ungrouped shuffle series. -/
theorem tsum_norm_shuffledIntegral_le_of_family_eventCount_bound [∀ i, Fintype (A i)]
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : ((i : ι) → History s (A i)) → E} {C : ℝ} {R : ι → ℝ}
    (hf : AEStronglyMeasurable f (historyFamilyMeasure ν s A q)) (hR : ∀ i, 0 ≤ R i)
    (hbound : ∀ᵐ h ∂historyFamilyMeasure ν s A q,
      ‖f h‖ ≤ C * ∏ i, (R i) ^ eventCount (h i)) :
    ∑' v, ‖shuffledIntegral ν s A q f v‖ ≤
      C * ∏ i, Real.exp ((totalIntensity (q i) : ℝ) * timeMass ν s * R i) := by
  have hi := integrable_of_family_eventCount_bound ν s A q hf hR hbound
  calc
    _ ≤ ∑' v, shuffledIntegral ν s A q (fun h => ‖f h‖) v := by
      apply Summable.tsum_le_tsum _ (summable_norm_shuffledIntegral ν s A q hi)
        (summable_shuffledIntegral_norm ν s A q hi)
      intro v
      simp only [shuffledIntegral, norm_smul, Real.norm_eq_abs,
        abs_of_nonneg (NNReal.coe_nonneg _), familyObservable_norm, smul_eq_mul]
      exact mul_le_mul_of_nonneg_left (norm_integral_le_integral_norm _) (NNReal.coe_nonneg _)
    _ = ∫ h, ‖f h‖ ∂historyFamilyMeasure ν s A q := tsum_shuffledIntegral_norm ν s A q hi
    _ ≤ _ := integral_norm_le_of_family_eventCount_bound ν s A q hf hR hbound

end FiniteFamily

end CertifiedGaps.WeightedPoissonHistory
