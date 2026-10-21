import CertifiedGaps.Transfer.WeightedProductPoissonHistory
import Mathlib.Logic.Equiv.Fin.Rotate
import Mathlib.MeasureTheory.Function.LocallyIntegrable

/-!
# Independent weighted atomless-history regression audit

These checks cover zero and infinite masses, unequal intensities, exact factorial factors,
empty labels and families, inverse permutations, actual original-measure adapters, and
nonzero fully ungrouped chronological expansions. The final command rejects every
transitive axiom outside the permitted three-axiom kernel base.
-/

open MeasureTheory Set
open scoped ENNReal NNReal
open CertifiedGaps.AtomlessTimeOrdering CertifiedGaps.WeightedPoissonHistory

namespace CertifiedGaps.WeightedHistoryAudit

instance : NullSingletonClass (0 : Measure ℝ) := ⟨by simp⟩

-- Empty chronological tuples always have unit mass, even for zero measure and empty domain.
theorem zero_measure_empty_tuple :
    timeMeasure (0 : Measure ℝ) 0 (timeSimplex ∅ 0) = 1 := by
  rw [measure_timeSimplex _ MeasurableSet.empty]
  simp

theorem zero_measure_positive_tuple (n : ℕ) :
    timeMeasure (0 : Measure ℝ) (n + 1) (timeSimplex Set.univ (n + 1)) = 0 := by
  rw [measure_timeSimplex _ MeasurableSet.univ]
  simp

theorem empty_domain_positive_tuple (ν : Measure ℝ) [SigmaFinite ν]
    [NullSingletonClass ν] (n : ℕ) :
    timeMeasure ν (n + 1) (timeSimplex ∅ (n + 1)) = 0 := by
  rw [measure_timeSimplex _ MeasurableSet.empty]
  simp

-- The simplex law does not silently require finite time mass.
theorem infinite_time_simplex :
    timeMeasure volume 3 (timeSimplex Set.univ 3) = ∞ := by
  rw [measure_timeSimplex _ MeasurableSet.univ]
  simp [ENNReal.top_div]

-- A negative starting point and n = 3 independently expose the factorial normalization.
theorem shifted_domain_factorial :
    timeMeasure volume 3 (timeSimplex (Ioo (-2) 1) 3) = 9 / 2 := by
  rw [measure_timeSimplex _ measurableSet_Ioo]
  norm_num
  apply (ENNReal.toReal_eq_toReal_iff' (by finiteness) (by finiteness)).mp
  norm_num

-- Genuine non-Lebesgue time intensity: three times restricted Lebesgue has total mass six.
noncomputable def timeLaw : Measure ℝ := (3 : ℝ≥0) • volume.restrict (Ioo (0 : ℝ) 2)

instance : IsFiniteMeasure timeLaw := by unfold timeLaw; infer_instance
instance : NullSingletonClass timeLaw := by
  constructor
  intro x
  simp [timeLaw, Measure.smul_apply]

@[simp] theorem timeLaw_univ : timeLaw Set.univ = 6 := by
  norm_num [timeLaw, Measure.smul_apply, ENNReal.smul_def, smul_eq_mul]

@[simp] theorem timeMass_timeLaw : timeMass timeLaw fullTimeDomain = 6 := by
  apply ENNReal.coe_injective
  rw [coe_timeMass]
  exact timeLaw_univ

-- Unequal labels have intensities 1 and 3; their weights must be multiplied, not counted.
def intensities : Fin 2 → ℝ≥0 := ![1, 3]

@[simp] theorem totalIntensity_intensities : totalIntensity intensities = 4 := by
  norm_num [totalIntensity, intensities, Fin.sum_univ_succ]

theorem unequal_word_weight : wordWeight intensities ![0, 1, 1] = 9 := by
  norm_num [wordWeight, intensities, Fin.prod_univ_succ]

theorem unequal_word_sum : (∑ a : Fin 3 → Fin 2, wordWeight intensities a) = 64 := by
  rw [sum_wordWeight, totalIntensity_intensities]
  norm_num

-- The chosen two-letter word has mass 3 * 6^2 / 2! = 54.
theorem unequal_component_mass :
    wordTimesMeasure timeLaw fullTimeDomain intensities ⟨2, ![0, 1]⟩ Set.univ = 54 := by
  rw [wordTimesMeasure_univ]
  norm_num [wordWeight, intensities, Fin.prod_univ_succ, fullTimeDomain]
  apply (ENNReal.toReal_eq_toReal_iff' (by finiteness) (by finiteness)).mp
  norm_num

-- Sum over every word of degree two: (4 * 6)^2 / 2! = 288.
theorem unequal_count_two_mass :
    historyMeasure timeLaw fullTimeDomain intensities {h | eventCount h = 2} = 288 := by
  rw [historyMeasure_eventCount]
  norm_num [fullTimeDomain]
  apply (ENNReal.toReal_eq_toReal_iff' (by finiteness) (by finiteness)).mp
  norm_num

theorem unequal_history_mass :
    historyMeasure timeLaw fullTimeDomain intensities Set.univ =
      ENNReal.ofReal (Real.exp 24) := by
  rw [historyMeasure_univ]
  norm_num

theorem unequal_poisson_count_law :
    (poissonHistoryMeasure timeLaw fullTimeDomain intensities).map eventCount =
      ProbabilityTheory.poissonMeasure 24 := by
  rw [map_eventCount_poissonHistoryMeasure]
  norm_num

-- Zero intensities, zero time measure, empty domain, and empty label type are all distinct cases.
theorem zero_intensities_mass (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    (s : TimeDomain) :
    historyMeasure ν s (fun _ : Fin 2 => (0 : ℝ≥0)) Set.univ = 1 := by
  rw [historyMeasure_univ_series]
  simp only [totalIntensity, Finset.sum_const_zero, ENNReal.coe_zero, zero_mul]
  rw [tsum_eq_single 0]
  · simp
  · intro n hn
    simp [hn]

theorem zero_intensities_positive_count (ν : Measure ℝ) [SigmaFinite ν]
    [NullSingletonClass ν] (s : TimeDomain) (n : ℕ) :
    historyMeasure ν s (fun _ : Fin 2 => (0 : ℝ≥0)) {h | eventCount h = n + 1} = 0 := by
  simp [historyMeasure_eventCount, totalIntensity]

theorem zero_time_measure_mass :
    historyMeasure (0 : Measure ℝ) fullTimeDomain intensities Set.univ = 1 := by
  rw [historyMeasure_univ]
  simp [timeMass]

theorem empty_time_domain_mass :
    historyMeasure volume ⟨∅, MeasurableSet.empty⟩ intensities Set.univ = 1 := by
  rw [historyMeasure_univ_series]
  simp only [measure_empty, mul_zero]
  rw [tsum_eq_single 0]
  · simp
  · intro n hn
    simp [hn]

theorem empty_labels_mass (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    (s : TimeDomain) (w : Fin 0 → ℝ≥0) : historyMeasure ν s w Set.univ = 1 := by
  rw [historyMeasure_univ_series]
  simp only [totalIntensity, Finset.univ_eq_empty, Finset.sum_empty,
    ENNReal.coe_zero, zero_mul]
  rw [tsum_eq_single 0]
  · simp
  · intro n hn
    simp [hn]

theorem empty_labels_positive_count (n : ℕ) :
    historyMeasure timeLaw fullTimeDomain (fun i : Fin 0 => Fin.elim0 i)
      {h | eventCount h = n + 1} = 0 := by
  simp [historyMeasure_eventCount, totalIntensity]

theorem infinite_history_mass :
    historyMeasure volume fullTimeDomain (fun _ : Fin 1 => (1 : ℝ≥0)) Set.univ = ∞ := by
  have h := measure_mono (μ := historyMeasure volume fullTimeDomain
    (fun _ : Fin 1 => (1 : ℝ≥0))) (Set.subset_univ {h | eventCount h = 1})
  rw [historyMeasure_eventCount] at h
  simpa [totalIntensity, fullTimeDomain] using h

-- Growth zero selects the unit empty history, while growth two gives exp(48).
theorem zero_growth_moment :
    ∫ h : History fullTimeDomain (Fin 2), (0 : ℝ) ^ eventCount h
      ∂historyMeasure timeLaw fullTimeDomain intensities = 1 := by
  rw [integral_eventCount_pow (by norm_num)]
  simp

theorem unequal_exponential_moment :
    ∫ h : History fullTimeDomain (Fin 2), (2 : ℝ) ^ eventCount h
      ∂historyMeasure timeLaw fullTimeDomain intensities = Real.exp 48 := by
  rw [integral_eventCount_pow (by norm_num)]
  norm_num

-- Actual measure equality, rather than just equality of total masses.
theorem original_measure_adapter (b : ℝ≥0) (m : ℕ) :
    historyMeasure (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b)
      (fun _ : Fin m => (1 : ℝ≥0)) = CertifiedGaps.PoissonHistory.historyMeasure b m :=
  historyMeasure_interval_unit b m

theorem original_probability_adapter (b : ℝ≥0) (m : ℕ) :
    poissonHistoryMeasure (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b)
      (fun _ : Fin m => (1 : ℝ≥0)) =
        CertifiedGaps.PoissonHistory.poissonHistoryMeasure b m :=
  poissonHistoryMeasure_interval_unit b m

-- The finite-domain reduction really proves the mass formula without global finiteness.
theorem finite_domain_exponential_mass {A : Type*} [Fintype A]
    (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain)
    (hs : ν s.val ≠ ∞) (w : A → ℝ≥0) :
    historyMeasure ν s w Set.univ =
      ENNReal.ofReal (Real.exp ((totalIntensity w : ℝ) * timeMass ν s)) := by
  let : IsFiniteMeasure (ν.restrict s.val) := ⟨by simpa using hs.lt_top⟩
  have he := historyMeasure_univ (ν.restrict s.val) s w
  rw [historyMeasure_restrict] at he
  simpa [timeMass, Measure.restrict_apply s.property] using he

-- The same finite-domain reduction proves the normalized law for the original measure.
theorem finite_domain_poisson_law {A : Type*} [Fintype A]
    (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain)
    (hs : ν s.val ≠ ∞) (w : A → ℝ≥0) :
    (poissonHistoryMeasure ν s w).map eventCount =
      ProbabilityTheory.poissonMeasure (totalIntensity w * timeMass ν s) := by
  let : IsFiniteMeasure (ν.restrict s.val) := ⟨by simpa using hs.lt_top⟩
  have he := map_eventCount_poissonHistoryMeasure (ν.restrict s.val) s w
  simpa [poissonHistoryMeasure, historyMeasure_restrict, timeMass,
    Measure.restrict_apply s.property] using he

-- Three-cycle transport detects an incorrect forward-permutation convention.
theorem nonsymmetric_cycle (ν : Measure ℝ) [SigmaFinite ν] (s : Set ℝ) :
    ∫ t in timeChamber s 3 (finRotate 3), (2 * t 0 + 3 * t 1 + 5 * t 2)
      ∂timeMeasure ν 3 =
      ∫ t in timeSimplex s 3, (2 * t 2 + 3 * t 0 + 5 * t 1) ∂timeMeasure ν 3 := by
  simpa [Function.comp_def, finRotate_symm_apply, show (-1 : Fin 3) = 2 by decide] using
    integral_timeChamber ν s 3 (finRotate 3) (fun t => 2 * t 0 + 3 * t 1 + 5 * t 2)

-- The concrete nonsymmetric integrand is independently integrable for the finite time law.
theorem nonsymmetric_integrable :
    IntegrableOn (fun t : Fin 3 → ℝ => 2 * t 0 + 3 * t 1 + 5 * t 2)
      (timeCube (Ioo 0 2) 3) (timeMeasure timeLaw 3) := by
  have hc : Continuous (fun t : Fin 3 → ℝ => 2 * t 0 + 3 * t 1 + 5 * t 2) := by
    fun_prop
  apply (hc.integrableOn_Icc (a := fun _ => 0) (b := fun _ => 2)
    (μ := timeMeasure timeLaw 3)).mono_set
  intro t ht
  exact ⟨fun i => (ht i (mem_univ i)).1.le, fun i => (ht i (mem_univ i)).2.le⟩

-- Empty families preserve the singleton empty tuple even when each alphabet is empty.
theorem empty_family_mass :
    historyFamilyMeasure timeLaw fullTimeDomain (fun _ : Fin 0 => Fin 0)
      (fun _ i => Fin.elim0 i) Set.univ = 1 := by
  simp [historyFamilyMeasure]

theorem zero_intensity_family_mass :
    historyFamilyMeasure timeLaw fullTimeDomain (fun _ : Fin 3 => Fin 2)
      (fun _ _ => 0) Set.univ = 1 := by
  simp [historyFamilyMeasure, Measure.pi_univ, zero_intensities_mass]

-- Unequal fixed word weights and an empty middle bond give product mass 6 * 1 * 54 = 324.
def sampleWords (i : Fin 3) : EventWord (Fin 2) :=
  ![⟨1, ![0]⟩, ⟨0, ![]⟩, ⟨2, ![0, 1]⟩] i

theorem sample_family_mass :
    familyTimesMeasure timeLaw fullTimeDomain (fun _ : Fin 3 => Fin 2)
      (fun _ => intensities) sampleWords Set.univ = 324 := by
  rw [familyTimesMeasure, Measure.pi_univ]
  simp only [Fin.prod_univ_succ]
  change wordTimesMeasure timeLaw fullTimeDomain intensities ⟨1, ![0]⟩ Set.univ *
    (wordTimesMeasure timeLaw fullTimeDomain intensities ⟨0, ![]⟩ Set.univ *
      (wordTimesMeasure timeLaw fullTimeDomain intensities ⟨2, ![0, 1]⟩ Set.univ * 1)) = 324
  rw [unequal_component_mass]
  norm_num [wordTimesMeasure_univ, wordWeight, intensities, Fin.prod_univ_succ,
    fullTimeDomain]

-- A nonzero fully ungrouped expansion reaches the actual weighted product-history integral.
theorem unequal_hasSum_shuffled :
    HasSum (shuffledIntegral timeLaw fullTimeDomain (fun _ : Fin 2 => Fin 2)
      (fun _ => intensities) (fun _ => (1 : ℝ))) (Real.exp 24 * Real.exp 24) := by
  have hi : Integrable (fun _ : Fin 2 → History fullTimeDomain (Fin 2) => (1 : ℝ))
      (historyFamilyMeasure timeLaw fullTimeDomain (fun _ : Fin 2 => Fin 2)
        (fun _ => intensities)) := by
    simpa using integrable_family_eventCount_pow timeLaw fullTimeDomain
      (fun _ : Fin 2 => Fin 2) (fun _ => intensities) (R := fun _ => 1) (fun _ => zero_le_one)
  have h := hasSum_shuffledIntegral timeLaw fullTimeDomain (fun _ : Fin 2 => Fin 2)
    (fun _ => intensities) hi
  have he := integral_family_eventCount_pow timeLaw fullTimeDomain
    (fun _ : Fin 2 => Fin 2) (fun _ => intensities) (R := fun _ => 1) (fun _ => zero_le_one)
  simp only [one_pow, Finset.prod_const_one] at he
  rw [he] at h
  norm_num [Fin.prod_univ_succ, pow_two] at h ⊢
  exact h

-- The concrete (1,0,2) weighted component has three genuine interleavings and total 324.
theorem sample_component_shuffle_sum :
    (∑ σ : TimeOrdering.TimeShuffle
        (TimeOrdering.familyHistoryOrder (fun i => (sampleWords i).1)),
      (familyWeight (fun _ : Fin 3 => Fin 2) (fun _ => intensities) sampleWords : ℝ) •
        ∫ t in timeSimplex fullTimeDomain.val (∑ i, (sampleWords i).1),
          familyObservable fullTimeDomain (fun _ : Fin 3 => Fin 2) (fun _ => (1 : ℝ))
            sampleWords (fun i j => t (σ.val.symm
              (TimeOrdering.historyEventEquiv (fun i => (sampleWords i).1) ⟨i, j⟩)))
          ∂timeMeasure timeLaw (∑ i, (sampleWords i).1)) = 324 := by
  have hi : Integrable (fun _ : Fin 3 → History fullTimeDomain (Fin 2) => (1 : ℝ))
      (historyFamilyMeasure timeLaw fullTimeDomain (fun _ : Fin 3 => Fin 2)
        (fun _ => intensities)) := by
    simpa using integrable_family_eventCount_pow timeLaw fullTimeDomain
      (fun _ : Fin 3 => Fin 2) (fun _ => intensities) (R := fun _ => 1) (fun _ => zero_le_one)
  rw [← integral_familyTimesMeasure_eq_sum_shuffles timeLaw fullTimeDomain
    (fun _ : Fin 3 => Fin 2) (fun _ => intensities) hi sampleWords]
  simp [integral_const, measureReal_def, sample_family_mass]

theorem three_interleavings :
    Fintype.card (TimeOrdering.TimeShuffle (TimeOrdering.familyHistoryOrder ![1, 0, 2])) = 3 := by
  have hf : IntegrableOn
      (fun _ : ((i : Fin 3) → Fin (![1, 0, 2] i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun i : Fin 3 => timeSimplex Set.univ (![1, 0, 2] i))
      (Measure.pi fun i : Fin 3 => timeMeasure timeLaw (![1, 0, 2] i)) :=
    integrableOn_const (by finiteness)
  have he := integral_pi_timeSimplex_eq_sum_shuffles ![1, 0, 2] timeLaw
    MeasurableSet.univ hf
  norm_num [integral_const, measureReal_def, Measure.pi_pi, measure_timeSimplex,
    Fin.prod_univ_succ, Fin.sum_univ_succ] at he
  exact_mod_cast (by linarith :
    (Fintype.card (TimeOrdering.TimeShuffle
      (TimeOrdering.familyHistoryOrder ![1, 0, 2])) : ℝ) = 3)

theorem empty_shuffle_integral :
    (∑ _σ : TimeOrdering.TimeShuffle
        (TimeOrdering.familyHistoryOrder (fun _ : Fin 0 => 0)),
      ∫ _t in timeSimplex Set.univ (∑ _i : Fin 0, 0), (7 : ℝ)
        ∂timeMeasure timeLaw (∑ _i : Fin 0, 0)) = 7 := by
  have hf : IntegrableOn
      (fun _ : ((i : Fin 0) → Fin ((fun _ : Fin 0 => 0) i) → ℝ) => (7 : ℝ))
      (Set.univ.pi fun _ : Fin 0 => timeSimplex Set.univ 0)
      (Measure.pi fun _ : Fin 0 => timeMeasure timeLaw 0) :=
    integrableOn_const (by finiteness)
  rw [← integral_pi_timeSimplex_eq_sum_shuffles (fun _ : Fin 0 => 0)
    timeLaw MeasurableSet.univ hf]
  simp [measureReal_def]

-- A zero-weight label kills the actual component for every observable, with no measurability
-- or unweighted-integrability hypothesis on that observable.
theorem zero_weight_component (f : History fullTimeDomain (Fin 2) → ℝ) :
    ∫ t, f ⟨⟨1, ![1]⟩, t⟩
      ∂wordTimesMeasure timeLaw fullTimeDomain ![2, 0] ⟨1, ![1]⟩ = 0 := by
  simp [wordTimesMeasure, wordWeight]

-- Unequal bond growth rates give a fully numeric exponential majorant for the individual
-- weighted shuffle terms, not merely for a grouped series.
theorem unequal_shuffle_exponential_bound :
    ∑' v, ‖shuffledIntegral timeLaw fullTimeDomain (fun _ : Fin 2 => Fin 2)
      (fun _ => intensities)
      (fun h => ∏ i : Fin 2, (![2, 3] i : ℝ) ^ eventCount (h i)) v‖ ≤
        Real.exp 48 * Real.exp 72 := by
  have hr (i : Fin 2) : 0 ≤ (![2, 3] i : ℝ) := by fin_cases i <;> norm_num
  have hi := integrable_family_eventCount_pow timeLaw fullTimeDomain
    (fun _ : Fin 2 => Fin 2) (fun _ => intensities) hr
  have hb : ∀ᵐ h ∂historyFamilyMeasure timeLaw fullTimeDomain
      (fun _ : Fin 2 => Fin 2) (fun _ => intensities),
      ‖∏ i : Fin 2, (![2, 3] i : ℝ) ^ eventCount (h i)‖ ≤
        1 * ∏ i : Fin 2, (![2, 3] i : ℝ) ^ eventCount (h i) := by
    apply ae_of_all
    intro h
    rw [Real.norm_eq_abs, abs_of_nonneg
      (Finset.prod_nonneg fun i _ => pow_nonneg (hr i) _), one_mul]
  have he := tsum_norm_shuffledIntegral_le_of_family_eventCount_bound timeLaw fullTimeDomain
    (fun _ : Fin 2 => Fin 2) (fun _ => intensities) hi.aestronglyMeasurable hr hb
  norm_num [Fin.prod_univ_succ] at he ⊢
  exact he

end CertifiedGaps.WeightedHistoryAudit

open Lean in
run_cmd do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut checked : Nat := 0
  for (name, _) in (← getEnv).constants.toList do
    if (`CertifiedGaps.AtomlessTimeOrdering).isPrefixOf name ||
        (`CertifiedGaps.WeightedPoissonHistory).isPrefixOf name ||
        (`CertifiedGaps.TimeOrdering).isPrefixOf name ||
        (`CertifiedGaps.PoissonHistory).isPrefixOf name ||
        (`CertifiedGaps.WeightedHistoryAudit).isPrefixOf name then
      checked := checked + 1
      let axioms ← collectAxioms name
      let bad := axioms.filter (!allowed.contains ·)
      unless bad.isEmpty do
        throwError "{name} depends on disallowed axioms: {bad}"
  logInfo m!"Axiom whitelist checked {checked} declarations."
