import CertifiedGaps.Transfer.PoissonHistory
import Mathlib.Logic.Equiv.Fin.Rotate

/-!
Regression checks for chronological ordering and the empty-history normalization.
Run with `lake build CertifiedGaps.Transfer.PoissonHistoryAudit`.
-/

open MeasureTheory CertifiedGaps.TimeOrdering CertifiedGaps.PoissonHistory
open scoped ENNReal NNReal

-- The zero-dimensional simplex is one point even at zero time.
example : volume (timeSimplex 0 0) = 1 := by simp [volume_timeSimplex]

-- A positive number of events at zero time has zero measure.
example (n : ℕ) : volume (timeSimplex 0 (n + 1)) = 0 := by
  simp [volume_timeSimplex]

-- A concrete factorial normalization in three dimensions.
example : volume.real (timeSimplex 2 3) = 4 / 3 := by
  rw [volumeReal_timeSimplex (by norm_num)]
  norm_num

-- Zero labels, or zero elapsed time, leave exactly the unit empty-history mass.
example (b : ℝ≥0) : historyMeasure b 0 Set.univ = 1 := by
  simp [historyMeasure_univ]

example (m : ℕ) : historyMeasure 0 m Set.univ = 1 := by
  simp [historyMeasure_univ]

example (b : ℝ≥0) (m : ℕ) : historyMeasure b m {h | eventCount h = 0} = 1 :=
  historyMeasure_empty b m

example (b : ℝ≥0) (n : ℕ) : historyMeasure b 0 {h | eventCount h = n + 1} = 0 := by
  simp [historyMeasure_eventCount]

example (m n : ℕ) : historyMeasure 0 m {h | eventCount h = n + 1} = 0 := by
  simp [historyMeasure_eventCount]

-- The three spin-component labels give the stated mass and count law.
example (b : ℝ≥0) : historyMeasure b 3 Set.univ = ENNReal.ofReal (Real.exp (3 * b)) :=
  historyMeasure_univ b 3

example (b : ℝ≥0) :
    (poissonHistoryMeasure b 3).map eventCount = ProbabilityTheory.poissonMeasure (3 * b) :=
  map_eventCount_poissonHistoryMeasure b 3

-- Growth rate zero selects only the empty history; growth rate one recovers total mass.
example (b : ℝ≥0) (m : ℕ) :
    ∫ h : History b m, (0 : ℝ) ^ eventCount h ∂historyMeasure b m = 1 := by
  rw [integral_eventCount_pow b m (by norm_num)]
  simp

example (b : ℝ≥0) (m : ℕ) :
    ∫ h : History b m, (1 : ℝ) ^ eventCount h ∂historyMeasure b m = Real.exp (m * b) := by
  simpa using integral_eventCount_pow b m (R := 1) (by norm_num)

-- A nonsymmetric integrand and a genuine three-cycle detect the inverse convention.
example (b : ℝ) :
    ∫ t in timeChamber b 3 (finRotate 3), (2 * t 0 + 3 * t 1 + 5 * t 2) =
      ∫ t in timeSimplex b 3, (2 * t 2 + 3 * t 0 + 5 * t 1) := by
  simpa [Function.comp_def, finRotate_symm_apply, show (-1 : Fin 3) = 2 by decide] using
    integral_timeChamber b 3 (finRotate 3) (fun t => 2 * t 0 + 3 * t 1 + 5 * t 2)

#print axioms volume_timeCollision
#print axioms ae_injective_times
#print axioms volume_timeSimplex
#print axioms integral_timeCube_eq_sum_timeSimplex
#print axioms historyMeasure_univ
#print axioms historyMeasure_eventCount
#print axioms map_eventCount_poissonHistoryMeasure
#print axioms hasSum_integral_historyMeasure_by_count
#print axioms integrable_eventCount_pow
#print axioms integral_eventCount_pow
#print axioms integrable_of_eventCount_bound
#print axioms integral_norm_le_of_eventCount_bound

-- Enforce the axiom policy for every declaration in both prerequisite namespaces, not only
-- the selected headline declarations printed above. Kernel replay checks the proof terms too.
open Lean in
run_cmd do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut checked : Nat := 0
  for (name, _) in (← getEnv).constants.toList do
    if (`CertifiedGaps.TimeOrdering).isPrefixOf name ||
        (`CertifiedGaps.PoissonHistory).isPrefixOf name then
      checked := checked + 1
      let axioms ← collectAxioms name
      let bad := axioms.filter (!allowed.contains ·)
      unless bad.isEmpty do
        throwError "{name} depends on disallowed axioms: {bad}"
  logInfo m!"Axiom whitelist checked {checked} declarations."
