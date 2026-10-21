import CertifiedGaps.Transfer.WeightedProductPoissonHistory

/-!
# Weighted atomless history regression checks

Exercise unequal and zero intensities, nonuniform disconnected time support, empty histories,
actual Lebesgue specialization, and nonzero absolutely convergent weighted shuffle expansions.
The final command enforces a transitive axiom whitelist for every new declaration.
-/

open MeasureTheory Set
open scoped ENNReal NNReal
open CertifiedGaps.AtomlessTimeOrdering CertifiedGaps.WeightedPoissonHistory

namespace CertifiedGaps.WeightedPoissonHistoryAudit

/-- A nonuniform event-time law supported on two disjoint intervals. -/
noncomputable def twoWindowTime : Measure ℝ :=
  volume.restrict (Ioo 0 1) + (2 : ℝ≥0) • volume.restrict (Ioo 2 3)

instance : IsFiniteMeasure twoWindowTime := by unfold twoWindowTime; infer_instance
instance : NullSingletonClass twoWindowTime := by
  constructor
  intro x
  simp [twoWindowTime, Measure.add_apply, Measure.smul_apply]

theorem twoWindowTime_univ : twoWindowTime Set.univ = 3 := by
  norm_num [twoWindowTime, Measure.add_apply, Measure.smul_apply, ENNReal.smul_def,
    Measure.restrict_apply, Real.volume_Ioo]

theorem twoWindowTime_mass : timeMass twoWindowTime fullTimeDomain = 3 := by
  simp [timeMass, fullTimeDomain, twoWindowTime_univ, ENNReal.toNNReal_ofNat]

/-- Unequal label intensities, including a zero-intensity label. -/
def testWeights : Fin 3 → ℝ≥0 := ![0, 2, 3]

theorem total_testWeights : totalIntensity testWeights = 5 := by
  norm_num [totalIntensity, testWeights, Fin.sum_univ_succ]

-- Three globally ordered events retain the factorial law for the nonuniform time law.
theorem nonuniform_three_simplex :
    timeMeasure twoWindowTime 3 (timeSimplex Set.univ 3) = 9 / 2 := by
  rw [measure_timeSimplex twoWindowTime MeasurableSet.univ, twoWindowTime_univ]
  norm_num only [Nat.factorial, Nat.cast_ofNat]
  apply (ENNReal.div_eq_div_iff (by norm_num) (by simp) (by norm_num) (by simp)).mpr
  norm_num

-- A genuinely weighted count mass uses the total intensity, not the label cardinality.
theorem weighted_two_event_mass :
    historyMeasure twoWindowTime fullTimeDomain testWeights {h | eventCount h = 2} = 225 / 2 := by
  rw [historyMeasure_eventCount]
  norm_num [fullTimeDomain, total_testWeights, twoWindowTime_univ]

-- The derived normalized count law sees both nonuniform time mass and unequal label rates.
theorem weighted_nonuniform_poisson :
    (poissonHistoryMeasure twoWindowTime fullTimeDomain testWeights).map eventCount =
      ProbabilityTheory.poissonMeasure 15 := by
  rw [map_eventCount_poissonHistoryMeasure, total_testWeights, twoWindowTime_mass]
  norm_num

-- The empty-history mass is one for every finite label set and every atomless time law.
theorem empty_history_mass {A : Type*} [Finite A] (ν : Measure ℝ)
    [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0) :
    historyMeasure ν s w {h | eventCount h = 0} = 1 := by
  let := Fintype.ofFinite A
  simp [historyMeasure_eventCount]

-- No labels and zero weights both retain the empty history only.
theorem no_labels_mass (ν : Measure ℝ) [IsFiniteMeasure ν] [NullSingletonClass ν]
    (s : TimeDomain) (w : Fin 0 → ℝ≥0) : historyMeasure ν s w Set.univ = 1 := by
  simp [historyMeasure_univ, totalIntensity]

theorem zero_weights_mass (ν : Measure ℝ) [IsFiniteMeasure ν] [NullSingletonClass ν]
    (s : TimeDomain) : historyMeasure ν s (fun _ : Fin 3 => (0 : ℝ≥0)) Set.univ = 1 := by
  simp [historyMeasure_univ, totalIntensity]

-- Any word containing a zero-intensity event carries the zero measure.
theorem zero_label_word :
    wordTimesMeasure twoWindowTime fullTimeDomain testWeights
      ⟨1, fun _ => (0 : Fin 3)⟩ = 0 := by
  simp [wordTimesMeasure, wordWeight, testWeights]

-- There are still no positive-count histories in a null time domain.
theorem null_domain_positive_count {A : Type*} [Finite A] (ν : Measure ℝ)
    [SigmaFinite ν] [NullSingletonClass ν] (s : TimeDomain) (w : A → ℝ≥0)
    (hs : ν s.val = 0) (n : ℕ) :
    historyMeasure ν s w {h | eventCount h = n + 1} = 0 := by
  let := Fintype.ofFinite A
  simp [historyMeasure_eventCount, hs]

-- The existing Lebesgue history is obtained by equality, not merely matching total masses.
theorem original_measure_specialization (b : ℝ≥0) (m : ℕ) :
    historyMeasure (volume.restrict (Ioo 0 (b : ℝ))) (intervalTimeDomain b)
      (fun _ : Fin m => (1 : ℝ≥0)) = PoissonHistory.historyMeasure b m :=
  historyMeasure_interval_unit b m

-- A nonzero full weighted shuffle expansion of two independent three-label histories.
theorem hasSum_weighted_one_shuffled :
    HasSum (shuffledIntegral twoWindowTime fullTimeDomain (fun _ : Fin 2 => Fin 3)
      (fun _ => testWeights) (fun _ => (1 : ℝ)))
      (∏ _i : Fin 2, Real.exp 15) := by
  have hi : Integrable (fun _ : (i : Fin 2) → History fullTimeDomain (Fin 3) => (1 : ℝ))
      (historyFamilyMeasure twoWindowTime fullTimeDomain (fun _ : Fin 2 => Fin 3)
        (fun _ => testWeights)) := by
    simpa using integrable_family_eventCount_pow twoWindowTime fullTimeDomain
      (fun _ : Fin 2 => Fin 3) (fun _ => testWeights) (R := fun _ => 1) (fun _ => zero_le_one)
  have h := hasSum_shuffledIntegral twoWindowTime fullTimeDomain
    (fun _ : Fin 2 => Fin 3) (fun _ => testWeights) hi
  have he := integral_family_eventCount_pow twoWindowTime fullTimeDomain
    (fun _ : Fin 2 => Fin 3) (fun _ => testWeights) (R := fun _ => 1) (fun _ => zero_le_one)
  simp only [one_pow, Finset.prod_const_one, mul_one, total_testWeights,
    twoWindowTime_mass, NNReal.coe_ofNat, show (5 : ℝ) * 3 = 15 by norm_num] at he
  rwa [he] at h

end CertifiedGaps.WeightedPoissonHistoryAudit

open Lean in
run_cmd do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut checked : Nat := 0
  for (name, _) in (← getEnv).constants.toList do
    if (`CertifiedGaps.AtomlessTimeOrdering).isPrefixOf name ||
        (`CertifiedGaps.WeightedPoissonHistory).isPrefixOf name ||
        (`CertifiedGaps.WeightedPoissonHistoryAudit).isPrefixOf name then
      checked := checked + 1
      let axioms ← collectAxioms name
      let bad := axioms.filter (!allowed.contains ·)
      unless bad.isEmpty do
        throwError "{name} depends on disallowed axioms: {bad}"
  logInfo m!"Axiom whitelist checked {checked} declarations."
