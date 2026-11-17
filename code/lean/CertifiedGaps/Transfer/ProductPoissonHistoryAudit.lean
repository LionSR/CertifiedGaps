/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.ProductPoissonHistory

/-!
# Product-history regrouping regression checks

Check empty/zero cases, actual product factorial weights, label and shuffle multiplicities,
nonzero complete-history expansions, and an enforcing transitive axiom whitelist.
-/

open MeasureTheory Set
open scoped ENNReal NNReal
open CertifiedGaps.PoissonHistory CertifiedGaps.TimeOrdering

namespace CertifiedGaps.ProductPoissonHistoryAudit

#print axioms historyFamilyMeasure_eq_sum
#print axioms hasSum_integral_historyFamilyMeasure
#print axioms summable_integral_norm_historyFamilyMeasure
#print axioms hasSum_integral_historyFamilyMeasure_by_counts
#print axioms integrable_family_eventCount_pow
#print axioms integral_family_eventCount_pow
#print axioms integrable_of_family_eventCount_bound

-- Degenerate empty family has one empty tuple with mass one.
theorem empty_family_univ (b : ℝ≥0) :
    historyFamilyMeasure b (fun _ : Fin 0 => 2) Set.univ = 1 := by
  simp [historyFamilyMeasure]

-- At b = 0 the only mass is the all-empty component, hence total product mass one.
theorem zero_time_univ {ι : Type*} [Fintype ι] (labels : ι → ℕ) :
    historyFamilyMeasure 0 labels Set.univ = 1 := by
  simp [historyFamilyMeasure, Measure.pi_univ, historyMeasure_univ]

-- With no labels every factor has only its empty word, for any time length.
theorem zero_labels_univ {ι : Type*} [Fintype ι] (b : ℝ≥0) :
    historyFamilyMeasure b (fun _ : ι => 0) Set.univ = 1 := by
  simp [historyFamilyMeasure, Measure.pi_univ, historyMeasure_univ]

-- A full count vector has precisely product_i m_i^n_i label choices.
theorem label_family_card {ι : Type*} [Fintype ι] [DecidableEq ι] (labels counts : ι → ℕ) :
    Fintype.card ((i : ι) → Fin (counts i) → Fin (labels i)) =
      ∏ i, labels i ^ counts i := by
  simp

-- Empty histories contribute exactly once, even with empty label alphabets.
theorem empty_label_family_card {ι : Type*} (labels : ι → ℕ) :
    Fintype.card ((i : ι) → Fin 0 → Fin (labels i)) = 1 := by
  simp

-- Every individual word-vector carries the product factorial weight, with no extra n!.
theorem family_times_univ {ι : Type*} [Fintype ι] (b : ℝ≥0) (labels : ι → ℕ)
    (w : (i : ι) → EventWord (labels i)) :
    familyTimesMeasure b labels w Set.univ =
      ∏ i, (b : ℝ≥0∞) ^ (w i).1 / ((w i).1).factorial := by
  simp [familyTimesMeasure, Measure.pi_univ, orderedTimesMeasure_univ]

-- Nontrivial scalar generating function specializes to the same finite product.
theorem integral_one_family {ι : Type*} [Fintype ι] (b : ℝ≥0) (labels : ι → ℕ) :
    ∫ _h : (i : ι) → History b (labels i), (1 : ℝ) ∂historyFamilyMeasure b labels =
      ∏ i, Real.exp (labels i * b) := by
  simpa using integral_family_eventCount_pow b labels (R := fun _ => 1) (fun _ => zero_le_one)


-- This is a genuine nonzero expansion of the full, actual three-bond history measure.
theorem hasSum_one_shuffled (b : ℝ≥0) :
    HasSum (shuffledIntegral b (fun _ : Fin 3 => 1) (fun _ => (1 : ℝ)))
      (∏ _i : Fin 3, Real.exp b) := by
  have hi : Integrable (fun _ : (i : Fin 3) → History b 1 => (1 : ℝ))
      (historyFamilyMeasure b (fun _ => 1)) := by
    simpa using integrable_family_eventCount_pow b (fun _ : Fin 3 => 1)
      (R := fun _ => 1) (fun _ => zero_le_one)
  have h := hasSum_shuffledIntegral b (fun _ : Fin 3 => 1) hi
  have he := integral_family_eventCount_pow b (fun _ : Fin 3 => 1)
    (R := fun _ => 1) (fun _ => zero_le_one)
  simp only [one_pow, Finset.prod_const_one, Nat.cast_one, one_mul, mul_one] at he
  rwa [he] at h

-- A fixed nonempty word vector, with an empty middle bond, has a nonzero shuffled value.
theorem sum_shuffles_nonzero :
    (∑ σ : TimeShuffle (familyHistoryOrder ![1, 0, 2]),
      ∫ t in timeSimplex 1 (∑ i : Fin 3, ![1, 0, 2] i),
        familyObservable 1 (fun _ : Fin 3 => 1) (fun _ => (1 : ℝ))
          (fun i => ⟨![1, 0, 2] i, fun _ => 0⟩)
          (fun i j => t (σ.val.symm (historyEventEquiv ![1, 0, 2] ⟨i, j⟩)))) = 1 / 2 := by
  have hi : Integrable (fun _ : (i : Fin 3) → History 1 1 => (1 : ℝ))
      (historyFamilyMeasure 1 (fun _ => 1)) := by
    simpa using integrable_family_eventCount_pow 1 (fun _ : Fin 3 => 1)
      (R := fun _ => 1) (fun _ => zero_le_one)
  have h := integral_familyTimesMeasure_eq_sum_shuffles 1 (fun _ : Fin 3 => 1) hi
    (fun i => ⟨![1, 0, 2] i, fun _ => 0⟩)
  simp only [NNReal.coe_one] at h
  rw [← h]
  norm_num [familyTimesMeasure, measureReal_def, Measure.pi_univ,
    orderedTimesMeasure_univ, Fin.prod_univ_succ]


#print axioms measurePreserving_familyTimesCoe
#print axioms hasSum_integral_historyFamilyMeasure_by_counts_shuffles
#print axioms summable_shuffledIntegral_norm
#print axioms summable_norm_shuffledIntegral
#print axioms hasSum_shuffledIntegral
#print axioms tsum_shuffledIntegral_norm
#print axioms tsum_norm_shuffledIntegral_le_of_family_eventCount_bound

-- Exact chamber-count regression: two single events have two chronological shuffles.
theorem two_single_event_shuffles :
    Fintype.card (TimeShuffle (familyHistoryOrder ![1, 1])) = 2 := by
  have hf : IntegrableOn
      (fun _ : ((i : Fin 2) → Fin (![1, 1] i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun i : Fin 2 => timeSimplex 1 (![1, 1] i)) :=
    integrableOn_const (by norm_num [volume_pi_pi, volume_timeSimplex, Fin.prod_univ_succ])
  have he := integral_pi_timeSimplex_eq_sum_shuffles ![1, 1] hf
  norm_num [integral_const, measureReal_def, volume_pi_pi, volume_timeSimplex,
    Fin.prod_univ_succ, Fin.sum_univ_succ] at he
  exact_mod_cast (by linarith : (Fintype.card (TimeShuffle (familyHistoryOrder ![1, 1])) : ℝ) = 2)

-- A zero-event middle bond creates no duplicate choices; lengths 1,0,2 give exactly three.
theorem zero_middle_bond_shuffles :
    Fintype.card (TimeShuffle (familyHistoryOrder ![1, 0, 2])) = 3 := by
  have hf : IntegrableOn
      (fun _ : ((i : Fin 3) → Fin (![1, 0, 2] i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun i : Fin 3 => timeSimplex 1 (![1, 0, 2] i)) :=
    integrableOn_const (by norm_num [volume_pi_pi, volume_timeSimplex, Fin.prod_univ_succ])
  have he := integral_pi_timeSimplex_eq_sum_shuffles ![1, 0, 2] hf
  norm_num [integral_const, measureReal_def, volume_pi_pi, volume_timeSimplex,
    Fin.prod_univ_succ, Fin.sum_univ_succ] at he
  exact_mod_cast (by linarith :
    (Fintype.card (TimeShuffle (familyHistoryOrder ![1, 0, 2])) : ℝ) = 3)

-- All zero counts admit the unique empty permutation, even for an empty family.
theorem all_zero_counts_shuffles {ι : Type*} [Fintype ι] :
    Fintype.card (TimeShuffle (familyHistoryOrder (fun _ : ι => 0))) = 1 := by
  have hf : IntegrableOn
      (fun _ : ((i : ι) → Fin ((fun _ : ι => 0) i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun _ : ι => timeSimplex 0 0) :=
    integrableOn_const (by simp [volume_pi_pi, volume_timeSimplex])
  have he := integral_pi_timeSimplex_eq_sum_shuffles (fun _ : ι => 0) hf
  simp [integral_const, measureReal_def, volume_pi_pi, volume_timeSimplex] at he
  exact_mod_cast he.symm

end CertifiedGaps.ProductPoissonHistoryAudit

open Lean in
run_cmd do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut checked : Nat := 0
  for (name, _) in (← getEnv).constants.toList do
    if (`CertifiedGaps.PoissonHistory).isPrefixOf name ||
        (`CertifiedGaps.TimeOrdering).isPrefixOf name ||
        (`CertifiedGaps.ProductPoissonHistoryAudit).isPrefixOf name then
      checked := checked + 1
      let axioms ← collectAxioms name
      let bad := axioms.filter (!allowed.contains ·)
      unless bad.isEmpty do
        throwError "{name} depends on disallowed axioms: {bad}"
  logInfo m!"Axiom whitelist checked {checked} declarations."
