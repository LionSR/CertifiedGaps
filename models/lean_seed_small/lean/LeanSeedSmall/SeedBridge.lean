import CertifiedGaps.Models.IsingLongitudinal
import LeanSeedSmall.Arithmetic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

/-!
The exact original seed goals, with the still-unproved partition enclosures made explicit.
No seed theorem from the imported model is used. Axiom output alone does NOT make a theorem
unconditional: the enclosure hypotheses below are the remaining mathematical obligations.
-/
namespace LeanSeedSmall
open Real CertifiedGaps

theorem ratio_of_log_bounds {x y eps : ℝ} (hx : 0 < x) (hy : 0 < y)
    (h : -eps ≤ Real.log y - 2 * Real.log x) : 1 - eps ≤ y / x^2 := by
  have hlog : Real.log (y / x^2) = Real.log y - 2 * Real.log x := by
    rw [Real.log_div (ne_of_gt hy) (pow_ne_zero 2 (ne_of_gt hx)), Real.log_pow]
    norm_num
  have hexp := Real.add_one_le_exp (Real.log (y / x^2))
  rw [Real.exp_log (div_pos hy (sq_pos_of_pos hx)), hlog] at hexp
  linarith

noncomputable abbrev selectedModel := IsingLongitudinal.model (1/4) (1/20)

theorem partition_pos (n : ℕ) (b : ℝ) : 0 < selectedModel.partitionFunction n b := by
  rw [RPRing.partitionFunction_eq_sum]
  exact Finset.sum_pos (fun _ _ => Real.exp_pos _) Finset.univ_nonempty

/-- These are precisely the unproved source-to-Lean obligations for the log-envelope route. -/
structure LogEnclosures : Prop where
  small_upper : Real.log (selectedModel.partitionFunction 6 (97/4)) ≤ logSmallUpper
  large_lower : (logLargeLower : ℝ) ≤ Real.log (selectedModel.partitionFunction 12 (97/4))
  large_upper : Real.log (selectedModel.partitionFunction 12 (97/4)) ≤ logLargeUpper
  double_lower : (logDoubleLower : ℝ) ≤ Real.log (selectedModel.partitionFunction 12 (97/2))

/-- Exact spatial seed at the first grid point; conditional on the displayed enclosures. -/
theorem seed_spatial_from_enclosures (h : LogEnclosures) :
    1 - (9064e-9 : ℝ) ≤ (IsingLongitudinal.model (1/4) (1/20)).spatialPurity 6 (97/4) := by
  change 1 - (9064e-9 : ℝ) ≤
    selectedModel.partitionFunction 12 (97/4) / selectedModel.partitionFunction 6 (97/4)^2
  apply ratio_of_log_bounds (partition_pos 6 (97/4)) (partition_pos 12 (97/4))
  have ha : -(9064 / 10^9 : ℝ) ≤ (logLargeLower : ℝ) - 2*(logSmallUpper : ℝ) := by
    have hq := (Rat.cast_le (K := ℝ)).2 spatial_arithmetic
    push_cast at hq
    exact hq
  have h1 := h.small_upper
  have h2 := h.large_lower
  norm_num at ha ⊢
  linarith

/-- Exact thermal seed at the first grid point; conditional on the displayed enclosures. -/
theorem seed_thermal_from_enclosures (h : LogEnclosures) :
    1 - (3141e-8 : ℝ) ≤ (IsingLongitudinal.model (1/4) (1/20)).thermalPurity 12 (97/4) := by
  have hβ : (2 : ℝ) * (97/4) = 97/2 := by norm_num
  simp only [RPRing.thermalPurity, hβ]
  apply ratio_of_log_bounds (partition_pos 12 (97/4)) (partition_pos 12 (97/2))
  have ha : -(3141 / 10^8 : ℝ) ≤ (logDoubleLower : ℝ) - 2*(logLargeUpper : ℝ) := by
    have hq := (Rat.cast_le (K := ℝ)).2 thermal_arithmetic
    push_cast at hq
    exact hq
  have h1 := h.large_upper
  have h2 := h.double_lower
  norm_num at ha ⊢
  linarith

#print axioms spatial_arithmetic
#print axioms thermal_arithmetic
#print axioms seed_spatial_from_enclosures
#print axioms seed_thermal_from_enclosures

/-- The literal first instance of the original thirty-point `seeds` statement. -/
abbrev firstPoint : IsingLongitudinal.Point :=
  IsingLongitudinal.points[0]'(by decide +kernel)

theorem first_point_seed_instance (h : LogEnclosures) :
    1 - (firstPoint.cert.sAt 0 : ℝ) ≤
      (IsingLongitudinal.model firstPoint.g firstPoint.h).spatialPurity 6 firstPoint.β₀ ∧
    1 - (firstPoint.cert.tAt 0 : ℝ) ≤
      (IsingLongitudinal.model firstPoint.g firstPoint.h).thermalPurity 12 firstPoint.β₀ := by
  constructor
  · convert seed_spatial_from_enclosures h using 1 <;>
      norm_num [firstPoint, IsingLongitudinal.points, ChainCertificate.sAt]
  · convert seed_thermal_from_enclosures h using 1 <;>
      norm_num [firstPoint, IsingLongitudinal.points, ChainCertificate.tAt]

#print axioms first_point_seed_instance

end LeanSeedSmall
