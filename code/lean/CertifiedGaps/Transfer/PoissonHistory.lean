import CertifiedGaps.Transfer.TimeOrdering
import Mathlib.Probability.Distributions.Poisson.Basic

/-!
# Finite labeled event histories

The history measure is a countable disjoint union of labeled increasing-time simplexes.
It is constructed from Lebesgue measure; its exponential mass and Poisson event-count law
are theorems, not inputs. There are no spin matrices or partition-function assumptions here.
-/

open MeasureTheory Set
open scoped ENNReal NNReal

namespace CertifiedGaps.PoissonHistory

open TimeOrdering

/-- An event count together with the word of event labels. -/
abbrev EventWord (m : ℕ) := Σ n : ℕ, Fin n → Fin m

/-- Increasing event times in the open interval `(0,b)`. -/
abbrev OrderedTimes (b : ℝ≥0) (n : ℕ) := {t : Fin n → ℝ // t ∈ timeSimplex b n}

/-- A finite labeled chronological history, including the empty history. -/
abbrev History (b : ℝ≥0) (m : ℕ) := Σ w : EventWord m, OrderedTimes b w.1

/-- Lebesgue measure on one ordered-time simplex. -/
noncomputable def orderedTimesMeasure (b : ℝ≥0) (n : ℕ) : Measure (OrderedTimes b n) :=
  volume.comap Subtype.val

/-- The history measure gives each label word its ordered-simplex Lebesgue measure. -/
noncomputable def historyMeasure (b : ℝ≥0) (m : ℕ) : Measure (History b m) :=
  Measure.sum fun w : EventWord m => (orderedTimesMeasure b w.1).map (Sigma.mk w)

/-- Inclusion of a labeled simplex is measurable for the disjoint-union sigma algebra. -/
theorem measurable_historyMk (b : ℝ≥0) (m : ℕ) (w : EventWord m) :
    Measurable (fun t : OrderedTimes b w.1 => (⟨w, t⟩ : History b m)) := by
  exact measurable_iff_le_map.mpr (iInf_le _ w)

/-- The measure of an entire ordered-time component. -/
theorem orderedTimesMeasure_univ (b : ℝ≥0) (n : ℕ) :
    orderedTimesMeasure b n Set.univ = (b : ℝ≥0∞) ^ n / n.factorial := by
  rw [orderedTimesMeasure, comap_subtype_coe_apply (measurableSet_timeSimplex b n)]
  simp only [Set.image_univ, Subtype.range_coe, volume_timeSimplex, ENNReal.ofReal_coe_nnreal]

instance (b : ℝ≥0) (n : ℕ) : IsFiniteMeasure (orderedTimesMeasure b n) := by
  constructor
  rw [orderedTimesMeasure_univ]
  finiteness

/-- Expanding the history mass by event count and label words. -/
theorem historyMeasure_univ_series (b : ℝ≥0) (m : ℕ) :
    historyMeasure b m Set.univ =
      ∑' n : ℕ, ((m : ℝ≥0∞) * b) ^ n / n.factorial := by
  simp only [historyMeasure, Measure.sum_apply _ MeasurableSet.univ,
    Measure.map_apply (measurable_historyMk b m _) MeasurableSet.univ, Set.preimage_univ,
    orderedTimesMeasure_univ]
  rw [ENNReal.tsum_sigma (fun n (_ : Fin n → Fin m) => (b : ℝ≥0∞) ^ n / n.factorial)]
  congr 1 with n
  simp [tsum_fintype, nsmul_eq_mul, mul_pow, mul_div_assoc]

/-- Summing the factorial weights gives the exponential mass of the history space. -/
theorem historyMeasure_univ (b : ℝ≥0) (m : ℕ) :
    historyMeasure b m Set.univ = ENNReal.ofReal (Real.exp (m * b)) := by
  rw [historyMeasure_univ_series]
  have hs := NormedSpace.expSeries_div_hasSum_exp ((m : ℝ) * b)
  rw [← Real.exp_eq_exp_ℝ] at hs
  rw [← hs.tsum_eq, ENNReal.ofReal_tsum_of_nonneg (fun n => by positivity) hs.summable]
  congr 1 with n
  simp [ENNReal.ofReal_div_of_pos (Nat.cast_pos.mpr (Nat.factorial_pos n)),
    ENNReal.ofReal_pow (by positivity : 0 ≤ (m : ℝ) * b), ENNReal.ofReal_mul,
    ENNReal.ofReal_natCast]

instance (b : ℝ≥0) (m : ℕ) : IsFiniteMeasure (historyMeasure b m) := by
  constructor
  rw [historyMeasure_univ]
  exact ENNReal.ofReal_lt_top

/-- Dividing the history measure by its exponential mass gives a probability measure. -/
noncomputable def poissonHistoryMeasure (b : ℝ≥0) (m : ℕ) : Measure (History b m) :=
  ENNReal.ofReal (Real.exp (-(m * b))) • historyMeasure b m

instance (b : ℝ≥0) (m : ℕ) : IsProbabilityMeasure (poissonHistoryMeasure b m) := by
  constructor
  simp only [poissonHistoryMeasure, Measure.smul_apply, smul_eq_mul, historyMeasure_univ]
  rw [← ENNReal.ofReal_mul (by positivity), ← Real.exp_add]
  simp

/-- The number of events in a history. -/
def eventCount {b : ℝ≥0} {m : ℕ} (h : History b m) : ℕ := h.1.1

/-- Event count is measurable on the disjoint union of ordered simplexes. -/
theorem measurable_eventCount (b : ℝ≥0) (m : ℕ) :
    Measurable (eventCount : History b m → ℕ) := by
  intro s hs
  apply MeasurableSpace.measurableSet_iInf.mpr
  intro w
  change MeasurableSet ((fun _ : OrderedTimes b w.1 => w.1) ⁻¹' s)
  exact measurable_const hs

/-- The unnormalized mass of histories containing exactly `n` events. -/
theorem historyMeasure_eventCount (b : ℝ≥0) (m n : ℕ) :
    historyMeasure b m {h | eventCount h = n} =
      ((m : ℝ≥0∞) * b) ^ n / n.factorial := by
  have hm : MeasurableSet {h : History b m | eventCount h = n} :=
    (measurable_eventCount b m) (measurableSet_singleton n)
  simp only [historyMeasure, Measure.sum_apply _ hm,
    Measure.map_apply (measurable_historyMk b m _) hm]
  have hp (w : EventWord m) :
      (Sigma.mk w ⁻¹' {h : History b m | eventCount h = n}) =
        if w.1 = n then Set.univ else ∅ := by
    ext t
    simp [eventCount]
  simp only [hp, apply_ite, orderedTimesMeasure_univ, measure_empty]
  rw [ENNReal.tsum_sigma (fun k (_ : Fin k → Fin m) =>
    if k = n then (b : ℝ≥0∞) ^ k / k.factorial else 0)]
  simp [tsum_fintype, nsmul_eq_mul, mul_pow, mul_div_assoc]

/-- The normalized history count is exactly Poisson with rate `m * b`. -/
theorem map_eventCount_poissonHistoryMeasure (b : ℝ≥0) (m : ℕ) :
    (poissonHistoryMeasure b m).map eventCount =
      ProbabilityTheory.poissonMeasure (m * b) := by
  apply Measure.ext_of_singleton
  intro n
  rw [Measure.map_apply (measurable_eventCount b m) (measurableSet_singleton n)]
  change poissonHistoryMeasure b m {h | eventCount h = n} = _
  rw [poissonHistoryMeasure, Measure.smul_apply, historyMeasure_eventCount,
    ProbabilityTheory.poissonMeasure_singleton]
  simp only [NNReal.coe_mul, NNReal.coe_natCast, smul_eq_mul]
  rw [ENNReal.ofReal_div_of_pos (Nat.cast_pos.mpr (Nat.factorial_pos n)),
    ENNReal.ofReal_mul (Real.exp_nonneg _),
    ENNReal.ofReal_pow (by positivity : 0 ≤ (m : ℝ) * b),
    ENNReal.ofReal_mul (Nat.cast_nonneg m)]
  simp [mul_div_assoc]

/-- The unique degree-zero component has mass one, even when `b = 0` or there are no labels. -/
theorem historyMeasure_empty (b : ℝ≥0) (m : ℕ) :
    historyMeasure b m {h | eventCount h = 0} = 1 := by
  simp [historyMeasure_eventCount]

section Integral

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- Integrating an absolutely integrable history observable is the sum of its ordered-time
integrals over all finite label words. -/
theorem integral_historyMeasure [CompleteSpace E] {b : ℝ≥0} {m : ℕ} {f : History b m → E}
    (hf : Integrable f (historyMeasure b m)) :
    ∫ h, f h ∂historyMeasure b m = ∑' w : EventWord m,
      ∫ t, f ⟨w, t⟩ ∂orderedTimesMeasure b w.1 := by
  rw [historyMeasure, integral_sum_measure hf]
  congr 1 with w
  exact integral_map (measurable_historyMk b m w).aemeasurable
    (hf.mono_measure (Measure.le_sum _ w)).aestronglyMeasurable

/-- The ordered-time expansion is summable whenever the history observable is integrable. -/
theorem hasSum_integral_historyMeasure [CompleteSpace E]
    {b : ℝ≥0} {m : ℕ} {f : History b m → E}
    (hf : Integrable f (historyMeasure b m)) :
    HasSum (fun w : EventWord m => ∫ t, f ⟨w, t⟩ ∂orderedTimesMeasure b w.1)
      (∫ h, f h ∂historyMeasure b m) := by
  have h := hasSum_integral_measure hf
  convert h using 1
  · funext w
    exact (integral_map (measurable_historyMk b m w).aemeasurable
      (hf.mono_measure (Measure.le_sum _ w)).aestronglyMeasurable).symm
  · rfl

/-- Regroup the absolutely convergent history expansion first by count and then by label word. -/
theorem hasSum_integral_historyMeasure_by_count [CompleteSpace E]
    {b : ℝ≥0} {m : ℕ} {f : History b m → E}
    (hf : Integrable f (historyMeasure b m)) :
    HasSum (fun n : ℕ => ∑ a : Fin n → Fin m,
      ∫ t, f ⟨⟨n, a⟩, t⟩ ∂orderedTimesMeasure b n)
      (∫ h, f h ∂historyMeasure b m) := by
  exact (hasSum_integral_historyMeasure hf).sigma (fun _ => hasSum_fintype _)

end Integral

/-- Exponential growth in event count is integrable. This supplies absolute convergence for
ordered operator expansions whose factors have a common norm bound. -/
theorem integrable_eventCount_pow (b : ℝ≥0) (m : ℕ) {R : ℝ} (hR : 0 ≤ R) :
    Integrable (fun h : History b m => R ^ eventCount h) (historyMeasure b m) := by
  have hp : Integrable (fun n : ℕ => R ^ n) (ProbabilityTheory.poissonMeasure (m * b)) := by
    rw [ProbabilityTheory.integrable_poissonMeasure_iff]
    have hs := (NormedSpace.expSeries_div_hasSum_exp ((m : ℝ) * b * R)).summable.mul_left
      (Real.exp (-((m : ℝ) * b)))
    convert hs using 1
    funext n
    simp only [NNReal.coe_mul, NNReal.coe_natCast, norm_pow, Real.norm_eq_abs, abs_of_nonneg hR,
      mul_pow]
    ring
  have hmap : Integrable (fun n : ℕ => R ^ n)
      ((poissonHistoryMeasure b m).map eventCount) := by
    rwa [map_eventCount_poissonHistoryMeasure]
  have hi := (integrable_map_measure (measurable_of_countable _).aestronglyMeasurable
    (measurable_eventCount b m).aemeasurable).mp hmap
  exact (integrable_smul_measure (by simp [Real.exp_pos]) (by finiteness)).mp hi

/-- The generating function of event count under the unnormalized history measure. -/
theorem integral_eventCount_pow (b : ℝ≥0) (m : ℕ) {R : ℝ} (hR : 0 ≤ R) :
    ∫ h : History b m, R ^ eventCount h ∂historyMeasure b m =
      Real.exp (m * b * R) := by
  have hs := hasSum_integral_historyMeasure_by_count (integrable_eventCount_pow b m hR)
  have he := NormedSpace.expSeries_div_hasSum_exp ((m : ℝ) * b * R)
  rw [← Real.exp_eq_exp_ℝ] at he
  apply HasSum.unique _ he
  convert hs using 1
  funext n
  simp only [eventCount, integral_const, measureReal_def, orderedTimesMeasure_univ,
    ENNReal.toReal_div, ENNReal.toReal_pow, ENNReal.coe_toReal, ENNReal.toReal_natCast,
    smul_eq_mul, Finset.sum_const, Finset.card_univ, Fintype.card_fun, Fintype.card_fin,
    nsmul_eq_mul, Nat.cast_pow]
  ring

/-- An observable dominated by exponential event-count growth is absolutely integrable. -/
theorem integrable_of_eventCount_bound {E : Type*} [NormedAddCommGroup E]
    {b : ℝ≥0} {m : ℕ} {f : History b m → E} {C R : ℝ}
    (hf : AEStronglyMeasurable f (historyMeasure b m)) (hR : 0 ≤ R)
    (hbound : ∀ᵐ h ∂historyMeasure b m, ‖f h‖ ≤ C * R ^ eventCount h) :
    Integrable f (historyMeasure b m) :=
  ((integrable_eventCount_pow b m hR).const_mul C).mono' hf hbound

/-- Quantitative absolute convergence of the history expansion. -/
theorem integral_norm_le_of_eventCount_bound {E : Type*} [NormedAddCommGroup E]
    {b : ℝ≥0} {m : ℕ} {f : History b m → E} {C R : ℝ}
    (hf : AEStronglyMeasurable f (historyMeasure b m)) (hR : 0 ≤ R)
    (hbound : ∀ᵐ h ∂historyMeasure b m, ‖f h‖ ≤ C * R ^ eventCount h) :
    ∫ h, ‖f h‖ ∂historyMeasure b m ≤ C * Real.exp (m * b * R) := by
  calc
    ∫ h, ‖f h‖ ∂historyMeasure b m ≤
        ∫ h, C * R ^ eventCount h ∂historyMeasure b m :=
      integral_mono_ae (integrable_of_eventCount_bound hf hR hbound).norm
        ((integrable_eventCount_pow b m hR).const_mul C) hbound
    _ = C * Real.exp (m * b * R) := by
      rw [integral_const_mul, integral_eventCount_pow b m hR]

end CertifiedGaps.PoissonHistory
