/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.HistoryShuffle
import Mathlib.Logic.Equiv.Fin.Rotate
import Mathlib.MeasureTheory.Function.LocallyIntegrable

/-!
# Finite-history shuffle regression checks

Concrete tests cover inverse permutations, forbidden within-history reversals, empty histories,
actual product measures, nonzero integrals, and a nonsymmetric integrand with independently
proved integrability. The final command audits every declaration in the shuffle namespaces.
-/

open MeasureTheory Set CertifiedGaps.TimeOrdering
open scoped ENNReal

namespace CertifiedGaps.HistoryShuffleAudit

-- A genuine three-cycle preserves the two-event second history; its inverse does not.
example : PreservesTimeOrder (twoHistoryOrder 1 2) (finRotate 3) := by
  unfold PreservesTimeOrder twoHistoryOrder
  decide +kernel
example : ¬ PreservesTimeOrder (twoHistoryOrder 1 2) (finRotate 3).symm := by
  unfold PreservesTimeOrder twoHistoryOrder
  decide +kernel

-- Explicit nonsymmetric inverse-permutation test for the (1,2)-shuffle.
example (t : Fin 3 → ℝ) :
    let p := splitHistoryTimes 1 2 (t ∘ (finRotate 3).symm)
    2 * p.1 0 + 3 * p.2 0 + 5 * p.2 1 = 2 * t 2 + 3 * t 0 + 5 * t 1 := by
  simp [Function.comp_def, finRotate_symm_apply]
  rfl

-- The target volume is exactly the product of Lebesgue measures, by definitional equality.
example (n m : ℕ) :
    (volume : Measure ((Fin n → ℝ) × (Fin m → ℝ))) =
      (volume : Measure (Fin n → ℝ)).prod (volume : Measure (Fin m → ℝ)) := rfl

-- Splitting preserves that explicitly displayed product measure.
example (n m : ℕ) :
    MeasurePreserving (splitHistoryTimes n m)
      (Measure.pi (fun _ : Fin (n + m) => (volume : Measure ℝ)))
      ((Measure.pi (fun _ : Fin n => (volume : Measure ℝ))).prod
        (Measure.pi (fun _ : Fin m => (volume : Measure ℝ)))) :=
  measurePreserving_splitHistoryTimes n m

-- No positivity assumption on b is hidden in zero-event normalization.
example (b : ℝ) : volume (timeSimplex b 0 ×ˢ timeSimplex b 0) = 1 := by
  simp [Measure.volume_eq_prod, Measure.prod_prod, volume_timeSimplex]

example (b : ℝ) :
    ∫ _ in timeSimplex b 0 ×ˢ timeSimplex b 0, (7 : ℝ) = 7 := by
  simp [Measure.volume_eq_prod, Measure.prod_prod, volume_timeSimplex,
    measureReal_def]

-- The zero-left-history case is covered by the headline theorem with no extra hypotheses.
example {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    (b : ℝ) (m : ℕ) (f : ((Fin 0 → ℝ) × (Fin m → ℝ)) → E)
    (hf : IntegrableOn f (timeSimplex b 0 ×ˢ timeSimplex b m)) :
    ∫ p in timeSimplex b 0 ×ˢ timeSimplex b m, f p =
      ∑ σ : TimeShuffle (twoHistoryOrder 0 m),
        ∫ t in timeSimplex b (0 + m),
          f (t ∘ σ.val.symm ∘ Fin.castAdd m, t ∘ σ.val.symm ∘ Fin.natAdd 0) :=
  integral_prod_timeSimplex_eq_sum_shuffles hf

#print axioms integral_prod_timeSimplex_eq_sum_shuffles
#print axioms integral_integral_timeSimplex_eq_sum_shuffles
#print axioms iUnion_timeShuffle_ae_eq
#print axioms measurePreserving_splitHistoryTimes


end CertifiedGaps.HistoryShuffleAudit

open MeasureTheory Set CertifiedGaps.TimeOrdering
open scoped ENNReal

namespace CertifiedGaps.HistoryShuffleAudit

-- The constant-integrand identity is nonzero: neither side can pass by the zero convention.
example :
    (∑ _σ : TimeShuffle (twoHistoryOrder 1 2),
      ∫ _ in timeSimplex 1 (1 + 2), (1 : ℝ)) = 1 / 2 := by
  have hf : IntegrableOn (fun _ : ((Fin 1 → ℝ) × (Fin 2 → ℝ)) => (1 : ℝ))
      (timeSimplex 1 1 ×ˢ timeSimplex 1 2) :=
    integrableOn_const (by simp [Measure.volume_eq_prod, Measure.prod_prod, volume_timeSimplex])
  have he := integral_prod_timeSimplex_eq_sum_shuffles hf
  rw [← he]
  norm_num [Measure.volume_eq_prod, Measure.prod_prod, volume_timeSimplex,
    measureReal_def]

/-- Independently discharge the integrability premise for a concrete nonsymmetric integrand. -/
theorem nonsymmetric_integrable (b : ℝ) :
    IntegrableOn
      (fun p : ((Fin 1 → ℝ) × (Fin 2 → ℝ)) => 2 * p.1 0 + 3 * p.2 0 + 5 * p.2 1)
      (timeSimplex b 1 ×ˢ timeSimplex b 2) := by
  have hc : Continuous
      (fun p : ((Fin 1 → ℝ) × (Fin 2 → ℝ)) => 2 * p.1 0 + 3 * p.2 0 + 5 * p.2 1) :=
    by fun_prop
  have hi := hc.integrableOn_Icc (a := (fun _ => 0, fun _ => 0))
    (b := (fun _ => b, fun _ => b)) (μ := volume)
  apply hi.mono_set
  intro p hp
  exact ⟨⟨fun i => (hp.1.1 i (mem_univ _)).1.le,
      fun i => (hp.2.1 i (mem_univ _)).1.le⟩,
    ⟨fun i => (hp.1.1 i (mem_univ _)).2.le,
      fun i => (hp.2.1 i (mem_univ _)).2.le⟩⟩

example (b : ℝ) :
    (∫ x in timeSimplex b 1, ∫ y in timeSimplex b 2, 2 * x 0 + 3 * y 0 + 5 * y 1) =
      ∑ σ : TimeShuffle (twoHistoryOrder 1 2),
        ∫ t in timeSimplex b 3,
          2 * t (σ.val.symm 0) + 3 * t (σ.val.symm 1) + 5 * t (σ.val.symm 2) := by
  simpa [Function.comp_def] using
    integral_integral_timeSimplex_eq_sum_shuffles (nonsymmetric_integrable b)

end CertifiedGaps.HistoryShuffleAudit

open MeasureTheory Set CertifiedGaps.TimeOrdering
open scoped ENNReal

namespace CertifiedGaps.HistoryShuffleAudit

-- Explicitly verify the nested product measure used by the finite-family coordinate change.
example {ι : Type*} [Fintype ι] (counts : ι → ℕ) :
    MeasurePreserving (splitFamilyHistoryTimes counts)
      (Measure.pi (fun _ : Fin (∑ i, counts i) => (volume : Measure ℝ)))
      (Measure.pi (fun i => Measure.pi (fun _ : Fin (counts i) => (volume : Measure ℝ)))) :=
  measurePreserving_splitFamilyHistoryTimes counts

-- For any enumeration, original event (i,j) receives the inverse shuffle coordinate.
example {ι : Type*} [Fintype ι] (counts : ι → ℕ)
    (σ : Equiv.Perm (Fin (∑ i, counts i))) (t : Fin (∑ i, counts i) → ℝ) :
    splitFamilyHistoryTimes counts (t ∘ σ.symm) =
      fun i j => t (σ.symm (historyEventEquiv counts ⟨i, j⟩)) := by
  ext i j
  simp

-- All-zero counts have unit measure for an arbitrary finite index type, even at b < 0.
example {ι : Type*} [Fintype ι] (b : ℝ) :
    volume (Set.univ.pi (fun _ : ι => timeSimplex b 0)) = 1 := by
  simp [volume_pi_pi, volume_timeSimplex]

-- In particular an empty family produces a nonzero unit shuffle integral.
example (b : ℝ) :
    (∑ _σ : TimeShuffle (familyHistoryOrder (fun _ : Fin 0 => 0)),
      ∫ _ in timeSimplex b (∑ _ : Fin 0, 0), (1 : ℝ)) = 1 := by
  have hf : IntegrableOn
      (fun _ : ((i : Fin 0) → Fin ((fun _ : Fin 0 => 0) i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun _ : Fin 0 => timeSimplex b 0) :=
    integrableOn_const (by simp [volume_pi_pi, volume_timeSimplex])
  have he := integral_pi_timeSimplex_eq_sum_shuffles (fun _ : Fin 0 => 0) hf
  rw [← he]
  simp [volume_pi_pi, volume_timeSimplex, measureReal_def]

-- A zero-event bond among nonempty bonds contributes one, not zero, to product mass.
example :
    volume.real (Set.univ.pi (fun i : Fin 3 => timeSimplex 1 (![1, 0, 2] i))) = 1 / 2 := by
  norm_num [measureReal_def, volume_pi_pi, volume_timeSimplex, Fin.prod_univ_succ]

#print axioms measurePreserving_uncurryHistoryTimes
#print axioms integral_pi_timeSimplex_eq_sum_shuffles

end CertifiedGaps.HistoryShuffleAudit


open MeasureTheory Set CertifiedGaps.TimeOrdering
open scoped ENNReal

namespace CertifiedGaps.HistoryShuffleAudit

example {ι : Type*} [Fintype ι] (b : ℝ) :
    (∑ _σ : TimeShuffle (familyHistoryOrder (fun _ : ι => 0)),
      ∫ _ in timeSimplex b (∑ _ : ι, 0), (1 : ℝ)) = 1 := by
  have hf : IntegrableOn
      (fun _ : ((i : ι) → Fin ((fun _ : ι => 0) i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun _ : ι => timeSimplex b 0) :=
    integrableOn_const (by simp [volume_pi_pi, volume_timeSimplex])
  have he := integral_pi_timeSimplex_eq_sum_shuffles (fun _ : ι => 0) hf
  rw [← he]
  simp [volume_pi_pi, volume_timeSimplex, measureReal_def]

example :
    (∑ _σ : TimeShuffle (familyHistoryOrder ![1, 0, 2]),
      ∫ _ in timeSimplex 1 (∑ i : Fin 3, ![1, 0, 2] i), (1 : ℝ)) = 1 / 2 := by
  have hf : IntegrableOn
      (fun _ : ((i : Fin 3) → Fin (![1, 0, 2] i) → ℝ) => (1 : ℝ))
      (Set.univ.pi fun i : Fin 3 => timeSimplex 1 (![1, 0, 2] i)) :=
    integrableOn_const (by norm_num [volume_pi_pi, volume_timeSimplex, Fin.prod_univ_succ])
  have he := integral_pi_timeSimplex_eq_sum_shuffles ![1, 0, 2] hf
  rw [← he]
  norm_num [volume_pi_pi, volume_timeSimplex, measureReal_def, Fin.prod_univ_succ]

end CertifiedGaps.HistoryShuffleAudit

open Lean in
run_cmd do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut checked : Nat := 0
  for (name, _) in (← getEnv).constants.toList do
    if (`CertifiedGaps.TimeOrdering).isPrefixOf name ||
        (`CertifiedGaps.HistoryShuffleAudit).isPrefixOf name then
      checked := checked + 1
      let axioms ← collectAxioms name
      let bad := axioms.filter (!allowed.contains ·)
      unless bad.isEmpty do
        throwError "{name} depends on disallowed axioms: {bad}"
  logInfo m!"Axiom whitelist checked {checked} declarations."
