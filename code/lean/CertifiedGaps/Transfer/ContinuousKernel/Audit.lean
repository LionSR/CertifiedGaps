import CertifiedGaps.Transfer.ContinuousKernel.Compact

/-!
# Audit of the bounded-kernel operator construction

Every principal conclusion below is derived from the stated kernel and measure hypotheses.
The output must contain only `propext`, `Classical.choice`, and `Quot.sound`.
-/

#print axioms CertifiedGaps.ContinuousKernel.BoundedKernel.operator
#print axioms CertifiedGaps.ContinuousKernel.BoundedKernel.norm_operator_le
#print axioms CertifiedGaps.ContinuousKernel.BoundedKernel.isSelfAdjoint_operator
#print axioms CertifiedGaps.ContinuousKernel.BoundedKernel.hasHilbertSchmidtBound_operator
#print axioms CertifiedGaps.ContinuousKernel.BoundedKernel.isHilbertSchmidt_operator
#print axioms CertifiedGaps.ContinuousKernel.BoundedKernel.isCompactOperator_operator
#print axioms CertifiedGaps.ContinuousKernel.hasSum_rankOne_of_summable_norm_sq
#print axioms CertifiedGaps.ContinuousKernel.IsHilbertSchmidt.isCompactOperator

open MeasureTheory

namespace CertifiedGaps.ContinuousKernel

/- The zero measure is allowed: no hidden positivity or nonzero-mass hypothesis is required. -/
example {α : Type*} [MeasurableSpace α] (K : BoundedKernel α ℂ) :
    K.operator (0 : Measure α) = 0 := by
  apply norm_eq_zero.mp
  exact le_antisymm (by simpa using K.norm_operator_le (0 : Measure α)) (norm_nonneg _)

/- On a one-point probability space, the constant kernel 3 is genuinely the scalar operator 3. -/
private def threeKernel : BoundedKernel Unit ℂ where
  toFun _ _ := 3
  stronglyMeasurable := stronglyMeasurable_const
  bound := 3
  norm_le _ _ := by norm_num

example (f : Lp ℂ 2 (Measure.dirac ())) :
    threeKernel.operator (Measure.dirac ()) f = (3 : ℂ) • f := by
  apply Lp.ext
  filter_upwards [threeKernel.coeFn_apply (Measure.dirac ()) f,
    Lp.coeFn_smul (3 : ℂ) f] with x hx hs
  rw [hs]
  change threeKernel.apply (Measure.dirac ()) f x = (3 : ℂ) * f x
  rw [hx]
  cases x
  simp [BoundedKernel.applyFun, threeKernel]

/- Imaginary off-diagonal entries require conjugation, not plain kernel symmetry. -/
private def imaginaryKernel : BoundedKernel (Fin 2) ℂ where
  toFun x y := if x = y then 0 else if x = 0 then Complex.I else -Complex.I
  stronglyMeasurable := (measurable_of_finite _).stronglyMeasurable
  bound := 1
  norm_le x y := by fin_cases x <;> fin_cases y <;> norm_num

example : IsSelfAdjoint (imaginaryKernel.operator Measure.count) := by
  apply imaginaryKernel.isSelfAdjoint_operator
  intro x y
  fin_cases x <;> fin_cases y <;> simp [imaginaryKernel]

end CertifiedGaps.ContinuousKernel
