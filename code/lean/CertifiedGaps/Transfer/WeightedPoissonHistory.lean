import CertifiedGaps.Transfer.AtomlessTimeOrdering
import CertifiedGaps.Transfer.PoissonHistory

/-!
# Weighted histories over arbitrary atomless time laws

A finite label type carries arbitrary nonnegative intensities. Each label word receives
its product intensity times the actual ordered-time measure. The exponential total mass,
Poisson count law, and convergence bounds are consequences of this construction.
-/

open MeasureTheory Set
open scoped ENNReal NNReal

namespace CertifiedGaps.WeightedPoissonHistory

open AtomlessTimeOrdering

/-- A measurable domain for event times. -/
abbrev TimeDomain := {s : Set ℝ // MeasurableSet s}

/-- Every real time is allowed. -/
def fullTimeDomain : TimeDomain := ⟨Set.univ, MeasurableSet.univ⟩

/-- The open interval used by the original Lebesgue history construction. -/
def intervalTimeDomain (b : ℝ≥0) : TimeDomain := ⟨Ioo 0 b, measurableSet_Ioo⟩

/-- An event count and a word in an arbitrary label type. -/
abbrev EventWord (A : Type*) := Σ n : ℕ, Fin n → A

/-- A strictly chronological time tuple in the specified measurable domain. -/
abbrev OrderedTimes (s : TimeDomain) (n : ℕ) :=
  {t : Fin n → ℝ // t ∈ timeSimplex s.val n}

/-- A finite labeled chronological history, including the empty history. -/
abbrev History (s : TimeDomain) (A : Type*) := Σ w : EventWord A, OrderedTimes s w.1

/-- The product of the intensities of all events in a word. -/
noncomputable def wordWeight {A : Type*} (w : A → ℝ≥0) {n : ℕ} (a : Fin n → A) : ℝ≥0 :=
  ∏ j, w (a j)

/-- Total label intensity. -/
noncomputable def totalIntensity {A : Type*} [Fintype A] (w : A → ℝ≥0) : ℝ≥0 := ∑ a, w a

/-- Finite nonnegative mass of the time domain. -/
noncomputable def timeMass (ν : Measure ℝ) (s : TimeDomain) : ℝ≥0 := (ν s.val).toNNReal

/-- Product time measure on the chronological simplex subtype. -/
noncomputable def orderedTimesMeasure (ν : Measure ℝ) (s : TimeDomain) (n : ℕ) :
    Measure (OrderedTimes s n) := (timeMeasure ν n).comap Subtype.val

/-- A word's time component includes its actual product intensity. -/
noncomputable def wordTimesMeasure {A : Type*} (ν : Measure ℝ) (s : TimeDomain)
    (w : A → ℝ≥0) (a : EventWord A) : Measure (OrderedTimes s a.1) :=
  wordWeight w a.2 • orderedTimesMeasure ν s a.1

/-- The weighted history measure is a countable disjoint sum of ordered-time components. -/
noncomputable def historyMeasure {A : Type*} (ν : Measure ℝ) (s : TimeDomain)
    (w : A → ℝ≥0) : Measure (History s A) :=
  Measure.sum fun a : EventWord A => (wordTimesMeasure ν s w a).map (Sigma.mk a)

theorem measurable_historyMk {A : Type*} (s : TimeDomain) (a : EventWord A) :
    Measurable (fun t : OrderedTimes s a.1 => (⟨a, t⟩ : History s A)) :=
  measurable_iff_le_map.mpr (iInf_le _ a)

@[simp]
theorem coe_timeMass (ν : Measure ℝ) [IsFiniteMeasure ν] (s : TimeDomain) :
    (timeMass ν s : ℝ≥0∞) = ν s.val := ENNReal.coe_toNNReal (measure_ne_top _ _)

/-- The ordered component has the factorial-normalized power of the time-domain mass. -/
theorem orderedTimesMeasure_univ (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    (s : TimeDomain) (n : ℕ) :
    orderedTimesMeasure ν s n Set.univ = ν s.val ^ n / n.factorial := by
  rw [orderedTimesMeasure, comap_subtype_coe_apply (measurableSet_timeSimplex s.property n)]
  simp only [Set.image_univ, Subtype.range_coe, measure_timeSimplex ν s.property]

instance (ν : Measure ℝ) [IsFiniteMeasure ν] [NullSingletonClass ν]
    (s : TimeDomain) (n : ℕ) : IsFiniteMeasure (orderedTimesMeasure ν s n) := by
  constructor
  rw [orderedTimesMeasure_univ]
  finiteness

instance {A : Type*} (ν : Measure ℝ) [IsFiniteMeasure ν] [NullSingletonClass ν]
    (s : TimeDomain) (w : A → ℝ≥0) (a : EventWord A) :
    IsFiniteMeasure (wordTimesMeasure ν s w a) := by
  unfold wordTimesMeasure
  infer_instance

/-- Summing all words of length `n` gives the `n`th power of the total intensity. -/
theorem sum_wordWeight {A : Type*} [Fintype A] (w : A → ℝ≥0) (n : ℕ) :
    ∑ a : Fin n → A, wordWeight w a = totalIntensity w ^ n :=
  (Fintype.sum_pow w n).symm

/-- Word weights and the time law separate exactly, including zero-intensity labels. -/
theorem wordTimesMeasure_univ {A : Type*} (ν : Measure ℝ)
    [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0)
    (a : EventWord A) :
    wordTimesMeasure ν s w a Set.univ =
      (wordWeight w a.2 : ℝ≥0∞) * (ν s.val ^ a.1 / a.1.factorial) := by
  simp only [wordTimesMeasure, Measure.smul_apply,
    orderedTimesMeasure_univ, ENNReal.smul_def, smul_eq_mul]

/-- Count expansion of the constructed weighted history measure. -/
theorem historyMeasure_univ_series {A : Type*} [Fintype A] (ν : Measure ℝ)
    [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) :
    historyMeasure ν s w Set.univ =
      ∑' n : ℕ, ((totalIntensity w : ℝ≥0∞) * ν s.val) ^ n / n.factorial := by
  classical
  simp only [historyMeasure, Measure.sum_apply _ MeasurableSet.univ,
    Measure.map_apply (measurable_historyMk s _) MeasurableSet.univ, Set.preimage_univ,
    wordTimesMeasure_univ]
  rw [ENNReal.tsum_sigma (fun n (a : Fin n → A) =>
    (wordWeight w a : ℝ≥0∞) * (ν s.val ^ n / n.factorial))]
  congr 1 with n
  simp only [tsum_fintype, ← Finset.sum_mul, ← ENNReal.ofNNReal_finsetSum,
    sum_wordWeight, ENNReal.coe_pow, mul_pow, mul_div_assoc]

/-- The full history mass is the exponential of total intensity times time mass. -/
theorem historyMeasure_univ {A : Type*} [Fintype A] (ν : Measure ℝ)
    [IsFiniteMeasure ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) :
    historyMeasure ν s w Set.univ =
      ENNReal.ofReal (Real.exp ((totalIntensity w : ℝ) * timeMass ν s)) := by
  rw [historyMeasure_univ_series, ← coe_timeMass ν s]
  have hs := NormedSpace.expSeries_div_hasSum_exp ((totalIntensity w : ℝ) * timeMass ν s)
  rw [← Real.exp_eq_exp_ℝ] at hs
  rw [← hs.tsum_eq, ENNReal.ofReal_tsum_of_nonneg (fun n => by positivity) hs.summable]
  congr 1 with n
  simp [ENNReal.ofReal_div_of_pos (Nat.cast_pos.mpr (Nat.factorial_pos n)),
    ENNReal.ofReal_pow (by positivity : 0 ≤ (totalIntensity w : ℝ) * timeMass ν s),
    ENNReal.ofReal_mul]

instance {A : Type*} [Finite A] (ν : Measure ℝ) [IsFiniteMeasure ν]
    [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) :
    IsFiniteMeasure (historyMeasure ν s w) := by
  let := Fintype.ofFinite A
  constructor
  rw [historyMeasure_univ]
  exact ENNReal.ofReal_lt_top

/-- Exponential normalization turns weighted histories into a probability law. -/
noncomputable def poissonHistoryMeasure {A : Type*} [Fintype A] (ν : Measure ℝ)
    (s : TimeDomain) (w : A → ℝ≥0) : Measure (History s A) :=
  ENNReal.ofReal (Real.exp (-((totalIntensity w : ℝ) * timeMass ν s))) • historyMeasure ν s w

instance {A : Type*} [Fintype A] (ν : Measure ℝ) [IsFiniteMeasure ν]
    [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) :
    IsProbabilityMeasure (poissonHistoryMeasure ν s w) := by
  constructor
  simp only [poissonHistoryMeasure, Measure.smul_apply, smul_eq_mul, historyMeasure_univ]
  rw [← ENNReal.ofReal_mul (by positivity), ← Real.exp_add]
  simp

/-- The history's event count. -/
def eventCount {s : TimeDomain} {A : Type*} (h : History s A) : ℕ := h.1.1

theorem measurable_eventCount (s : TimeDomain) (A : Type*) :
    Measurable (eventCount : History s A → ℕ) := by
  intro u hu
  apply MeasurableSpace.measurableSet_iInf.mpr
  intro a
  change MeasurableSet ((fun _ : OrderedTimes s a.1 => a.1) ⁻¹' u)
  exact measurable_const hu

/-- The exact unnormalized mass of degree `n`. -/
theorem historyMeasure_eventCount {A : Type*} [Fintype A] (ν : Measure ℝ)
    [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) (n : ℕ) :
    historyMeasure ν s w {h | eventCount h = n} =
      ((totalIntensity w : ℝ≥0∞) * ν s.val) ^ n / n.factorial := by
  classical
  have hm : MeasurableSet {h : History s A | eventCount h = n} :=
    (measurable_eventCount s A) (measurableSet_singleton n)
  simp only [historyMeasure, Measure.sum_apply _ hm,
    Measure.map_apply (measurable_historyMk s _) hm]
  have hp (a : EventWord A) :
      (Sigma.mk a ⁻¹' {h : History s A | eventCount h = n}) =
        if a.1 = n then Set.univ else ∅ := by
    ext t
    simp [eventCount]
  simp only [hp, apply_ite, wordTimesMeasure_univ, measure_empty]
  rw [ENNReal.tsum_sigma (fun k (a : Fin k → A) =>
    if k = n then (wordWeight w a : ℝ≥0∞) * (ν s.val ^ k / k.factorial) else 0)]
  simp only [tsum_fintype]
  simp_rw [Finset.sum_ite_irrel, Finset.sum_const_zero]
  rw [tsum_ite_eq]
  simp only [← Finset.sum_mul, ← ENNReal.ofNNReal_finsetSum, sum_wordWeight,
    ENNReal.coe_pow, mul_pow, mul_div_assoc]

/-- The event count has the Poisson distribution with the derived total rate. -/
theorem map_eventCount_poissonHistoryMeasure {A : Type*} [Fintype A] (ν : Measure ℝ)
    [IsFiniteMeasure ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) :
    (poissonHistoryMeasure ν s w).map eventCount =
      ProbabilityTheory.poissonMeasure (totalIntensity w * timeMass ν s) := by
  apply Measure.ext_of_singleton
  intro n
  rw [Measure.map_apply (measurable_eventCount s A) (measurableSet_singleton n)]
  change poissonHistoryMeasure ν s w {h | eventCount h = n} = _
  rw [poissonHistoryMeasure, Measure.smul_apply, historyMeasure_eventCount,
    ProbabilityTheory.poissonMeasure_singleton, ← coe_timeMass ν s]
  simp only [NNReal.coe_mul, smul_eq_mul]
  rw [ENNReal.ofReal_div_of_pos (Nat.cast_pos.mpr (Nat.factorial_pos n)),
    ENNReal.ofReal_mul (Real.exp_nonneg _),
    ENNReal.ofReal_pow (by positivity : 0 ≤ (totalIntensity w : ℝ) * timeMass ν s),
    ENNReal.ofReal_mul (NNReal.coe_nonneg _)]
  simp [mul_div_assoc]

section Integral

variable {A : Type*} [Fintype A] {ν : Measure ℝ} [IsFiniteMeasure ν]
  [NullSingletonClass ν] {s : TimeDomain} {w : A → ℝ≥0}

omit [Fintype A] [IsFiniteMeasure ν] [NullSingletonClass ν] in
/-- Integrability gives an absolutely convergent expansion over complete weighted words. -/
theorem hasSum_integral_historyMeasure
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : History s A → E} (hf : Integrable f (historyMeasure ν s w)) :
    HasSum (fun a : EventWord A => ∫ t, f ⟨a, t⟩ ∂wordTimesMeasure ν s w a)
      (∫ h, f h ∂historyMeasure ν s w) := by
  have h := hasSum_integral_measure hf
  convert h using 1
  · funext a
    exact (integral_map (measurable_historyMk s a).aemeasurable
      (hf.mono_measure (Measure.le_sum _ a)).aestronglyMeasurable).symm
  · rfl

omit [IsFiniteMeasure ν] [NullSingletonClass ν] in
/-- Regroup the weighted history expansion by event count. -/
theorem hasSum_integral_historyMeasure_by_count
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {f : History s A → E} (hf : Integrable f (historyMeasure ν s w)) :
    HasSum (fun n : ℕ => ∑ a : Fin n → A,
      ∫ t, f ⟨⟨n, a⟩, t⟩ ∂wordTimesMeasure ν s w ⟨n, a⟩)
      (∫ h, f h ∂historyMeasure ν s w) := by
  exact (hasSum_integral_historyMeasure hf).sigma (fun _ => hasSum_fintype _)

omit [Fintype A] in
/-- Exponential event-count growth is integrable, proved from the derived Poisson law. -/
theorem integrable_eventCount_pow [Finite A] {R : ℝ} (hR : 0 ≤ R) :
    Integrable (fun h : History s A => R ^ eventCount h) (historyMeasure ν s w) := by
  let := Fintype.ofFinite A
  have hp : Integrable (fun n : ℕ => R ^ n)
      (ProbabilityTheory.poissonMeasure (totalIntensity w * timeMass ν s)) := by
    rw [ProbabilityTheory.integrable_poissonMeasure_iff]
    have hs := (NormedSpace.expSeries_div_hasSum_exp
      ((totalIntensity w : ℝ) * timeMass ν s * R)).summable.mul_left
        (Real.exp (-((totalIntensity w : ℝ) * timeMass ν s)))
    convert hs using 1
    funext n
    simp only [NNReal.coe_mul, norm_pow, Real.norm_eq_abs, abs_of_nonneg hR, mul_pow]
    ring
  have hmap : Integrable (fun n : ℕ => R ^ n)
      ((poissonHistoryMeasure ν s w).map eventCount) := by
    rwa [map_eventCount_poissonHistoryMeasure]
  have hi := (integrable_map_measure (measurable_of_countable _).aestronglyMeasurable
    (measurable_eventCount s A).aemeasurable).mp hmap
  exact (integrable_smul_measure (by simp [Real.exp_pos]) (by finiteness)).mp hi

/-- Exact exponential moment of the constructed unnormalized history measure. -/
theorem integral_eventCount_pow {R : ℝ} (hR : 0 ≤ R) :
    ∫ h : History s A, R ^ eventCount h ∂historyMeasure ν s w =
      Real.exp ((totalIntensity w : ℝ) * timeMass ν s * R) := by
  have hs := hasSum_integral_historyMeasure_by_count
    (integrable_eventCount_pow (ν := ν) (s := s) (w := w) hR)
  have he := NormedSpace.expSeries_div_hasSum_exp
    ((totalIntensity w : ℝ) * timeMass ν s * R)
  rw [← Real.exp_eq_exp_ℝ] at he
  apply HasSum.unique _ he
  convert hs using 1
  funext n
  simp only [eventCount, integral_const, measureReal_def, wordTimesMeasure, Measure.smul_apply,
    orderedTimesMeasure_univ, ENNReal.smul_def,
    ENNReal.toReal_mul, ENNReal.toReal_div, ENNReal.toReal_pow, ENNReal.coe_toReal,
    ENNReal.toReal_natCast, smul_eq_mul]
  rw [← Finset.sum_mul]
  simp_rw [mul_assoc]
  rw [← Finset.sum_mul, ← NNReal.coe_sum, sum_wordWeight]
  simp only [NNReal.coe_pow, timeMass, ENNReal.toReal]
  ring

omit [Fintype A] in
/-- Exponential domination supplies absolute convergence rather than presupposing it. -/
theorem integrable_of_eventCount_bound [Finite A] {E : Type*} [NormedAddCommGroup E]
    {f : History s A → E} {C R : ℝ}
    (hf : AEStronglyMeasurable f (historyMeasure ν s w)) (hR : 0 ≤ R)
    (hbound : ∀ᵐ h ∂historyMeasure ν s w, ‖f h‖ ≤ C * R ^ eventCount h) :
    Integrable f (historyMeasure ν s w) :=
  ((integrable_eventCount_pow hR).const_mul C).mono' hf hbound

/-- The exponential majorant quantitatively bounds the observable's L¹ norm. -/
theorem integral_norm_le_of_eventCount_bound {E : Type*} [NormedAddCommGroup E]
    {f : History s A → E} {C R : ℝ}
    (hf : AEStronglyMeasurable f (historyMeasure ν s w)) (hR : 0 ≤ R)
    (hbound : ∀ᵐ h ∂historyMeasure ν s w, ‖f h‖ ≤ C * R ^ eventCount h) :
    ∫ h, ‖f h‖ ∂historyMeasure ν s w ≤
      C * Real.exp ((totalIntensity w : ℝ) * timeMass ν s * R) := by
  calc
    _ ≤ ∫ h, C * R ^ eventCount h ∂historyMeasure ν s w :=
      integral_mono_ae (integrable_of_eventCount_bound hf hR hbound).norm
        ((integrable_eventCount_pow hR).const_mul C) hbound
    _ = _ := by rw [integral_const_mul, integral_eventCount_pow hR]

end Integral

/-- Restricting the time law to the allowed domain leaves the actual ordered component unchanged. -/
theorem orderedTimesMeasure_restrict (ν : Measure ℝ) [SigmaFinite ν]
    (s : TimeDomain) (n : ℕ) :
    orderedTimesMeasure (ν.restrict s.val) s n = orderedTimesMeasure ν s n := by
  ext u _
  simp only [orderedTimesMeasure,
    comap_subtype_coe_apply (measurableSet_timeSimplex s.property n)]
  change (Measure.pi fun _ : Fin n => ν.restrict s.val) (Subtype.val '' u) = _
  rw [← Measure.restrict_pi_pi (fun _ : Fin n => ν) (fun _ => s.val)]
  rw [Measure.restrict_apply' (MeasurableSet.univ_pi fun _ : Fin n => s.property)]
  congr 1
  apply Set.inter_eq_left.mpr
  rintro _ ⟨t, _, rfl⟩
  exact t.property.1

/-- The original ordered-time Lebesgue measure is recovered by actual measure equality. -/
theorem orderedTimesMeasure_interval (b : ℝ≥0) (n : ℕ) :
    orderedTimesMeasure (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b) n =
      PoissonHistory.orderedTimesMeasure b n := by
  rw [show Ioo 0 (b : ℝ) = (intervalTimeDomain b).val from rfl,
    orderedTimesMeasure_restrict]
  rfl

/-- Unit label intensities and restricted Lebesgue time recover the original history measure. -/
theorem historyMeasure_interval_unit (b : ℝ≥0) (m : ℕ) :
    historyMeasure (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b)
      (fun _ : Fin m => (1 : ℝ≥0)) = PoissonHistory.historyMeasure b m := by
  unfold historyMeasure PoissonHistory.historyMeasure
  congr 1
  funext a
  unfold wordTimesMeasure
  rw [show wordWeight (fun _ : Fin m => (1 : ℝ≥0)) a.2 = 1 by simp [wordWeight]]
  rw [one_smul, orderedTimesMeasure_interval]
  rfl

/-- Domain restriction leaves the entire weighted history measure unchanged. -/
theorem historyMeasure_restrict {A : Type*} (ν : Measure ℝ) [SigmaFinite ν]
    (s : TimeDomain) (w : A → ℝ≥0) :
    historyMeasure (ν.restrict s.val) s w = historyMeasure ν s w := by
  simp only [historyMeasure, wordTimesMeasure, orderedTimesMeasure_restrict]

/-- The interval specialization has exactly its original elapsed-time parameter. -/
theorem timeMass_interval (b : ℝ≥0) :
    timeMass (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b) = b := by
  simp [timeMass, intervalTimeDomain, Measure.restrict_apply, Real.volume_Ioo]

/-- The normalized original Poisson history law is also recovered by actual measure equality. -/
theorem poissonHistoryMeasure_interval_unit (b : ℝ≥0) (m : ℕ) :
    poissonHistoryMeasure (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b)
      (fun _ : Fin m => (1 : ℝ≥0)) = PoissonHistory.poissonHistoryMeasure b m := by
  rw [poissonHistoryMeasure, historyMeasure_interval_unit, timeMass_interval]
  simp only [totalIntensity, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul, mul_one, NNReal.coe_natCast, PoissonHistory.poissonHistoryMeasure]
  rfl

end CertifiedGaps.WeightedPoissonHistory
