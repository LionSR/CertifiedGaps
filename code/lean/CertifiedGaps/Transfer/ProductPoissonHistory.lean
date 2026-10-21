import CertifiedGaps.Transfer.PoissonHistory
import CertifiedGaps.Transfer.HistoryShuffle
import Mathlib.MeasureTheory.Integral.Pi

/-!
# Finite products of labeled Poisson histories

The product history measure is decomposed along its countable measurable partition by
bondwise label words. This gives an absolutely convergent integral expansion over actual
product ordered-simplex measures; chronological shuffle identities then apply to each term.
-/

open MeasureTheory Set Function
open scoped ENNReal NNReal

namespace CertifiedGaps.PoissonHistory

open TimeOrdering

/-- The component with one specified complete event word. -/
def wordFiber {b : ℝ≥0} {m : ℕ} (w : EventWord m) : Set (History b m) :=
  {h | h.1 = w}

theorem measurableSet_wordFiber {b : ℝ≥0} {m : ℕ} (w : EventWord m) :
    MeasurableSet (wordFiber (b := b) w) := by
  apply MeasurableSpace.measurableSet_iInf.mpr
  intro v
  change MeasurableSet ((fun _ : OrderedTimes b v.1 => v) ⁻¹' {w})
  by_cases h : v = w <;> simp [h]

/-- Restricting the history measure to a word recovers its actual ordered-simplex measure. -/
theorem restrict_historyMeasure_wordFiber (b : ℝ≥0) (m : ℕ) (w : EventWord m) :
    (historyMeasure b m).restrict (wordFiber w) =
      (orderedTimesMeasure b w.1).map (Sigma.mk w) := by
  classical
  rw [historyMeasure, Measure.restrict_sum _ (measurableSet_wordFiber w)]
  have hterm (v : EventWord m) :
      ((orderedTimesMeasure b v.1).map (Sigma.mk v)).restrict (wordFiber w) =
        if v = w then (orderedTimesMeasure b w.1).map (Sigma.mk w) else 0 := by
    rw [Measure.restrict_map (measurable_historyMk b m v) (measurableSet_wordFiber w)]
    by_cases h : v = w
    · subst v
      simp [wordFiber]
    · have he : Sigma.mk v ⁻¹' wordFiber (b := b) w = ∅ := by
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

variable {ι : Type*} [Fintype ι] (b : ℝ≥0) (labels : ι → ℕ)

/-- The actual finite product of the bond-history measures. -/
noncomputable def historyFamilyMeasure : Measure ((i : ι) → History b (labels i)) :=
  Measure.pi (fun i => historyMeasure b (labels i))

/-- Insert one tuple of ordered times into a specified family of words. -/
def historyFamilyMk (w : (i : ι) → EventWord (labels i))
    (t : (i : ι) → OrderedTimes b (w i).1) (i : ι) : History b (labels i) :=
  ⟨w i, t i⟩

omit [Fintype ι] in
theorem measurable_historyFamilyMk (w : (i : ι) → EventWord (labels i)) :
    Measurable (historyFamilyMk b labels w) := by
  exact Measurable.of_eval fun i =>
    (measurable_historyMk b (labels i) (w i)).comp (measurable_pi_apply i)

/-- Product-simplex measure for a specified family of label words. -/
noncomputable def familyTimesMeasure (w : (i : ι) → EventWord (labels i)) :
    Measure ((i : ι) → OrderedTimes b (w i).1) :=
  Measure.pi (fun i => orderedTimesMeasure b (w i).1)

/-- Each word-vector part of the finite product is a pushforward of the product of the
ordered-time component measures. -/
theorem restrict_historyFamilyMeasure (w : (i : ι) → EventWord (labels i)) :
    (historyFamilyMeasure b labels).restrict (Set.univ.pi fun i => wordFiber (w i)) =
      (familyTimesMeasure b labels w).map (historyFamilyMk b labels w) := by
  rw [historyFamilyMeasure, Measure.restrict_pi_pi]
  simp_rw [restrict_historyMeasure_wordFiber]
  exact (Measure.pi_map_pi (fun i =>
    (measurable_historyMk b (labels i) (w i)).aemeasurable)).symm

/-- A finite product of countable history disjoint unions is genuinely the countable sum
of its product-simplex components. This identity is proved from the measurable partition,
not imposed as a hypothesis. -/
theorem historyFamilyMeasure_eq_sum :
    historyFamilyMeasure b labels = Measure.sum fun w : (i : ι) → EventWord (labels i) =>
      (familyTimesMeasure b labels w).map (historyFamilyMk b labels w) := by
  classical
  have hcover : (⋃ w : (i : ι) → EventWord (labels i),
      Set.univ.pi fun i => wordFiber (b := b) (w i)) = Set.univ := by
    ext h
    simp only [Set.mem_iUnion, Set.mem_univ, iff_true]
    exact ⟨fun i => (h i).1, fun i _ => rfl⟩
  have hdisj : Pairwise (Disjoint on (fun w : (i : ι) → EventWord (labels i) =>
      Set.univ.pi fun i => wordFiber (b := b) (w i))) := by
    intro w v hwv
    apply Set.disjoint_left.mpr
    intro h hw hv
    apply hwv
    funext i
    exact (hw i (Set.mem_univ i)).symm.trans (hv i (Set.mem_univ i))
  have hm (w : (i : ι) → EventWord (labels i)) :
      MeasurableSet (Set.univ.pi fun i => wordFiber (b := b) (w i)) :=
    MeasurableSet.univ_pi fun i => measurableSet_wordFiber (w i)
  have h := Measure.restrict_iUnion (μ := historyFamilyMeasure b labels) hdisj hm
  rw [hcover, Measure.restrict_univ] at h
  simpa only [restrict_historyFamilyMeasure] using h


/-- The product-history integral expands into product ordered-simplex integrals. -/
theorem hasSum_integral_historyFamilyMeasure
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    HasSum (fun w : (i : ι) → EventWord (labels i) =>
      ∫ t, f (historyFamilyMk b labels w t) ∂familyTimesMeasure b labels w)
      (∫ h, f h ∂historyFamilyMeasure b labels) := by
  rw [historyFamilyMeasure_eq_sum] at hf ⊢
  have h := hasSum_integral_measure hf
  convert h using 1
  funext w
  exact (integral_map (measurable_historyFamilyMk b labels w).aemeasurable
    (hf.mono_measure (Measure.le_sum _ w)).aestronglyMeasurable).symm

/-- Absolute convergence is proved at the level of the integrals of the norms, which also
justifies every subsequent finite regrouping and chronological subdivision. -/
theorem summable_integral_norm_historyFamilyMeasure
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    Summable (fun w : (i : ι) → EventWord (labels i) =>
      ∫ t, ‖f (historyFamilyMk b labels w t)‖ ∂familyTimesMeasure b labels w) := by
  rw [historyFamilyMeasure_eq_sum] at hf
  convert hf.summable_integral using 1
  funext w
  exact (integral_map (measurable_historyFamilyMk b labels w).aemeasurable
    (hf.mono_measure (Measure.le_sum _ w)).norm.aestronglyMeasurable).symm

/-- Separate bondwise complete words into their count vector and finite label-word family. -/
def familyWordsEquiv : ((i : ι) → EventWord (labels i)) ≃
    Σ counts : ι → ℕ, (i : ι) → Fin (counts i) → Fin (labels i) where
  toFun w := ⟨fun i => (w i).1, fun i => (w i).2⟩
  invFun v i := ⟨v.1 i, v.2 i⟩
  left_inv w := by rfl
  right_inv v := by cases v; rfl

/-- Count-vector and label-word expansion of the actual finite-product history integral. -/
theorem hasSum_integral_historyFamilyMeasure_by_counts [DecidableEq ι]
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    HasSum (fun counts : ι → ℕ => ∑ a : (i : ι) → Fin (counts i) → Fin (labels i),
      ∫ t, f (historyFamilyMk b labels (fun i => ⟨counts i, a i⟩) t)
        ∂familyTimesMeasure b labels (fun i => ⟨counts i, a i⟩))
      (∫ h, f h ∂historyFamilyMeasure b labels) := by
  have h := ((familyWordsEquiv labels).symm.hasSum_iff).mpr
    (hasSum_integral_historyFamilyMeasure b labels hf)
  exact h.sigma fun _ => hasSum_fintype _

/-- Product exponential event-count growth is integrable for any finite bond family. -/
theorem integrable_family_eventCount_pow {R : ι → ℝ} (hR : ∀ i, 0 ≤ R i) :
    Integrable (fun h : (i : ι) → History b (labels i) =>
      ∏ i, (R i) ^ eventCount (h i)) (historyFamilyMeasure b labels) :=
  Integrable.fintype_prod_dep fun i => integrable_eventCount_pow b (labels i) (hR i)

/-- Exact product exponential moment, obtained from the actual finite product measure. -/
theorem integral_family_eventCount_pow {R : ι → ℝ} (hR : ∀ i, 0 ≤ R i) :
    ∫ h : (i : ι) → History b (labels i), ∏ i, (R i) ^ eventCount (h i)
      ∂historyFamilyMeasure b labels = ∏ i, Real.exp (labels i * b * R i) := by
  rw [historyFamilyMeasure]
  calc
    _ = ∏ i, ∫ h : History b (labels i), R i ^ eventCount h ∂historyMeasure b (labels i) :=
      integral_fintype_prod_eq_prod (fun i (h : History b (labels i)) => R i ^ eventCount h)
    _ = _ := Finset.prod_congr rfl fun i _ => integral_eventCount_pow b (labels i) (hR i)

/-- An explicit bondwise event-count bound supplies the absolute integrability required by
the product-history series, without assuming convergence of that series. -/
theorem integrable_of_family_eventCount_bound
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E} {C : ℝ} {R : ι → ℝ}
    (hf : AEStronglyMeasurable f (historyFamilyMeasure b labels)) (hR : ∀ i, 0 ≤ R i)
    (hbound : ∀ᵐ h ∂historyFamilyMeasure b labels,
      ‖f h‖ ≤ C * ∏ i, (R i) ^ eventCount (h i)) :
    Integrable f (historyFamilyMeasure b labels) :=
  ((integrable_family_eventCount_pow b labels hR).const_mul C).mono' hf hbound


/-- Forget the proofs that each bond's time tuple is ordered. -/
def familyTimesCoe (counts : ι → ℕ) (t : (i : ι) → OrderedTimes b (counts i))
    (i : ι) : Fin (counts i) → ℝ := (t i).val

omit [Fintype ι] in
/-- Product subtype inclusion is a measurable embedding into the ambient time coordinates. -/
theorem measurableEmbedding_familyTimesCoe [Finite ι] (counts : ι → ℕ) :
    MeasurableEmbedding (familyTimesCoe b counts) := by
  have hs : MeasurableSet {t : (i : ι) → Fin (counts i) → ℝ |
      ∀ i, t i ∈ timeSimplex b (counts i)} := by
    simp only [Set.ofPred_forall]
    exact MeasurableSet.iInter fun i =>
      (measurableSet_timeSimplex b (counts i)).preimage (measurable_pi_apply i)
  exact (MeasurableEmbedding.subtype_coe hs).comp
    (MeasurableEquiv.subtypePiEquivPi
      (p := fun i t => t ∈ timeSimplex b (counts i))).symm.measurableEmbedding

/-- The product of the simplex subtype measures maps to product Lebesgue measure restricted
to the product of the simplexes. -/
theorem measurePreserving_familyTimesCoe (counts : ι → ℕ) :
    MeasurePreserving (familyTimesCoe b counts)
      (Measure.pi fun i => orderedTimesMeasure b (counts i))
      (volume.restrict (Set.univ.pi fun i => timeSimplex b (counts i))) where
  measurable := (measurableEmbedding_familyTimesCoe b counts).measurable
  map_eq := by
    change (Measure.pi fun i => orderedTimesMeasure b (counts i)).map
      (fun t i => (t i).val) = _
    have (i : ι) : SigmaFinite ((orderedTimesMeasure b (counts i)).map
        (Subtype.val : OrderedTimes b (counts i) → (Fin (counts i) → ℝ))) := inferInstance
    rw [Measure.pi_map_pi (X := fun i => OrderedTimes b (counts i))
      (Y := fun i => Fin (counts i) → ℝ)
      (μ := fun i => orderedTimesMeasure b (counts i))
      (f := fun i => (Subtype.val : OrderedTimes b (counts i) → (Fin (counts i) → ℝ)))
      (fun i => measurable_subtype_coe.aemeasurable)]
    simp_rw [orderedTimesMeasure,
      map_comap_subtype_coe (measurableSet_timeSimplex b _)]
    exact (Measure.restrict_pi_pi _ _).symm

/-- The canonical ambient extension of an observable in one word-vector component.
It vanishes off the product of the ordered simplexes. -/
noncomputable def familyObservable
    {E : Type*} [Zero E] (f : ((i : ι) → History b (labels i)) → E)
    (w : (i : ι) → EventWord (labels i))
    (t : (i : ι) → Fin (w i).1 → ℝ) : E := by
  classical
  exact if ht : ∀ i, t i ∈ timeSimplex b (w i).1 then
    f (historyFamilyMk b labels w (fun i => ⟨t i, ht i⟩)) else 0

@[simp]
theorem familyObservable_familyTimesCoe
    {E : Type*} [Zero E] (f : ((i : ι) → History b (labels i)) → E)
    (w : (i : ι) → EventWord (labels i))
    (t : (i : ι) → OrderedTimes b (w i).1) :
    familyObservable b labels f w (familyTimesCoe b (fun i => (w i).1) t) =
      f (historyFamilyMk b labels w t) := by
  simp [familyObservable, familyTimesCoe, fun i => (t i).property]

/-- Integrability descends to each product-simplex component of the actual history measure. -/
theorem integrable_historyFamilyMk
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels))
    (w : (i : ι) → EventWord (labels i)) :
    Integrable (fun t => f (historyFamilyMk b labels w t)) (familyTimesMeasure b labels w) := by
  rw [historyFamilyMeasure_eq_sum] at hf
  exact (hf.mono_measure (Measure.le_sum _ w)).comp_measurable
    (measurable_historyFamilyMk b labels w)

/-- The canonical ambient extension is integrable on the product of the bond simplexes. -/
theorem integrableOn_familyObservable
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels))
    (w : (i : ι) → EventWord (labels i)) :
    IntegrableOn (familyObservable b labels f w)
      (Set.univ.pi fun i => timeSimplex b (w i).1) := by
  apply ((measurePreserving_familyTimesCoe b (fun i => (w i).1)).integrable_comp_emb
    (measurableEmbedding_familyTimesCoe b (fun i => (w i).1))).mp
  simpa only [Function.comp_def, familyObservable_familyTimesCoe, familyTimesMeasure] using
    integrable_historyFamilyMk b labels hf w

/-- Each actual product-simplex integral becomes precisely the sum over its chronological
shuffles, with no exchange or commutation of the observable's values. -/
theorem integral_familyTimesMeasure_eq_sum_shuffles
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels))
    (w : (i : ι) → EventWord (labels i)) :
    ∫ t, f (historyFamilyMk b labels w t) ∂familyTimesMeasure b labels w =
      ∑ σ : TimeShuffle (familyHistoryOrder (fun i => (w i).1)),
        ∫ t in timeSimplex b (∑ i, (w i).1),
          familyObservable b labels f w
            (fun i j => t (σ.val.symm (historyEventEquiv (fun i => (w i).1) ⟨i, j⟩))) := by
  have he := (measurePreserving_familyTimesCoe b (fun i => (w i).1)).integral_comp
    (measurableEmbedding_familyTimesCoe b (fun i => (w i).1))
    (familyObservable b labels f w)
  simp only [familyObservable_familyTimesCoe] at he
  exact he.trans (integral_pi_timeSimplex_eq_sum_shuffles _
    (integrableOn_familyObservable b labels hf w))

/-- The actual finite-product history integral equals an absolutely convergent count-vector,
label-word and chronological-shuffle series. Only ordinary integrability is assumed; the
measure decomposition and chronological change of variables are proved above. -/
theorem hasSum_integral_historyFamilyMeasure_by_counts_shuffles [DecidableEq ι]
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    HasSum (fun counts : ι → ℕ => ∑ a : (i : ι) → Fin (counts i) → Fin (labels i),
      ∑ σ : TimeShuffle (familyHistoryOrder counts),
        ∫ t in timeSimplex b (∑ i, counts i),
          familyObservable b labels f (fun i => ⟨counts i, a i⟩)
            (fun i j => t (σ.val.symm (historyEventEquiv counts ⟨i, j⟩))))
      (∫ h, f h ∂historyFamilyMeasure b labels) := by
  have h := hasSum_integral_historyFamilyMeasure_by_counts b labels hf
  convert h using 1
  funext counts
  apply Finset.sum_congr rfl
  intro a _
  exact (integral_familyTimesMeasure_eq_sum_shuffles b labels hf
    (fun i => ⟨counts i, a i⟩)).symm


/-- Quantitative absolute convergence from bondwise event-count bounds. -/
theorem integral_norm_le_of_family_eventCount_bound
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E} {C : ℝ} {R : ι → ℝ}
    (hf : AEStronglyMeasurable f (historyFamilyMeasure b labels)) (hR : ∀ i, 0 ≤ R i)
    (hbound : ∀ᵐ h ∂historyFamilyMeasure b labels,
      ‖f h‖ ≤ C * ∏ i, (R i) ^ eventCount (h i)) :
    ∫ h, ‖f h‖ ∂historyFamilyMeasure b labels ≤
      C * ∏ i, Real.exp (labels i * b * R i) := by
  calc
    _ ≤ ∫ h, C * ∏ i, (R i) ^ eventCount (h i) ∂historyFamilyMeasure b labels :=
      integral_mono_ae (integrable_of_family_eventCount_bound b labels hf hR hbound).norm
        ((integrable_family_eventCount_pow b labels hR).const_mul C) hbound
    _ = _ := by rw [integral_const_mul, integral_family_eventCount_pow b labels hR]

@[simp]
theorem familyObservable_norm
    {E : Type*} [NormedAddCommGroup E]
    (f : ((i : ι) → History b (labels i)) → E)
    (w : (i : ι) → EventWord (labels i))
    (t : (i : ι) → Fin (w i).1 → ℝ) :
    familyObservable b labels (fun h => ‖f h‖) w t =
      ‖familyObservable b labels f w t‖ := by
  classical
  unfold familyObservable
  split_ifs <;> simp

/-- One index records a count vector, one label word for every bond, and one admissible
chronological interleaving. Every choice occurs once, with no extra factorial weight. -/
abbrev FamilyShuffleIndex := Σ counts : ι → ℕ,
  ((i : ι) → Fin (counts i) → Fin (labels i)) × TimeShuffle (familyHistoryOrder counts)

/-- One term of the completely ungrouped chronological series. -/
noncomputable def shuffledIntegral
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (f : ((i : ι) → History b (labels i)) → E) (v : FamilyShuffleIndex labels) : E :=
  ∫ t in timeSimplex b (∑ i, v.1 i),
    familyObservable b labels f (fun i => ⟨v.1 i, v.2.1 i⟩)
      (fun i j => t (v.2.2.val.symm (historyEventEquiv v.1 ⟨i, j⟩)))

/-- The integrals of the norms are summable over every count, label and shuffle individually,
not merely after cancellation inside a fixed count vector. -/
theorem summable_shuffledIntegral_norm
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    Summable (shuffledIntegral b labels (fun h => ‖f h‖)) := by
  classical
  have hn (v : FamilyShuffleIndex labels) :
      0 ≤ shuffledIntegral b labels (fun h => ‖f h‖) v := by
    unfold shuffledIntegral
    exact integral_nonneg fun t => by rw [familyObservable_norm]; exact norm_nonneg _
  apply (summable_sigma_of_nonneg hn).mpr
  refine ⟨fun _ => (hasSum_fintype _).summable, ?_⟩
  have h := (hasSum_integral_historyFamilyMeasure_by_counts_shuffles b labels hf.norm).summable
  simpa only [tsum_fintype, Fintype.sum_prod_type, shuffledIntegral] using h

/-- Absolute convergence of the individual Banach-valued chronological terms. -/
theorem summable_norm_shuffledIntegral
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    Summable (fun v => ‖shuffledIntegral b labels f v‖) := by
  apply Summable.of_nonneg_of_le (fun _ => norm_nonneg _)
    (fun v => ?_) (summable_shuffledIntegral_norm b labels hf)
  simp only [shuffledIntegral, familyObservable_norm]
  exact norm_integral_le_integral_norm _

/-- Fully ungrouped chronological expansion of the finite-product history integral. Its
absolute convergence is the preceding theorem, and event-count bounds imply its hypothesis. -/
theorem hasSum_shuffledIntegral
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    HasSum (shuffledIntegral b labels f) (∫ h, f h ∂historyFamilyMeasure b labels) := by
  classical
  have hs := (summable_norm_shuffledIntegral b labels hf).of_norm
  apply hs.hasSum_iff.mpr
  rw [hs.tsum_sigma]
  simpa only [tsum_fintype, Fintype.sum_prod_type, shuffledIntegral] using
    (hasSum_integral_historyFamilyMeasure_by_counts_shuffles b labels hf).tsum_eq


/-- The sum of the individual chronological norm integrals equals the original L¹ norm. -/
theorem tsum_shuffledIntegral_norm
    {E : Type*} [NormedAddCommGroup E]
    {f : ((i : ι) → History b (labels i)) → E}
    (hf : Integrable f (historyFamilyMeasure b labels)) :
    ∑' v, shuffledIntegral b labels (fun h => ‖f h‖) v =
      ∫ h, ‖f h‖ ∂historyFamilyMeasure b labels :=
  (hasSum_shuffledIntegral b labels hf.norm).tsum_eq

/-- Explicit absolute-convergence majorant for the fully ungrouped chronological series. -/
theorem tsum_norm_shuffledIntegral_le_of_family_eventCount_bound
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : ((i : ι) → History b (labels i)) → E} {C : ℝ} {R : ι → ℝ}
    (hf : AEStronglyMeasurable f (historyFamilyMeasure b labels)) (hR : ∀ i, 0 ≤ R i)
    (hbound : ∀ᵐ h ∂historyFamilyMeasure b labels,
      ‖f h‖ ≤ C * ∏ i, (R i) ^ eventCount (h i)) :
    ∑' v, ‖shuffledIntegral b labels f v‖ ≤
      C * ∏ i, Real.exp (labels i * b * R i) := by
  have hi := integrable_of_family_eventCount_bound b labels hf hR hbound
  calc
    _ ≤ ∑' v, shuffledIntegral b labels (fun h => ‖f h‖) v := by
      apply Summable.tsum_le_tsum _ (summable_norm_shuffledIntegral b labels hi)
        (summable_shuffledIntegral_norm b labels hi)
      intro v
      simp only [shuffledIntegral, familyObservable_norm]
      exact norm_integral_le_integral_norm _
    _ = ∫ h, ‖f h‖ ∂historyFamilyMeasure b labels :=
      tsum_shuffledIntegral_norm b labels hi
    _ ≤ _ := integral_norm_le_of_family_eventCount_bound b labels hf hR hbound

end FiniteFamily

end CertifiedGaps.PoissonHistory
