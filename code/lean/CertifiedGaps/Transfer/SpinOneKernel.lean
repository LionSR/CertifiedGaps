import CertifiedGaps.Transfer.MergedHistory

/-!
# The continuous spin-one history operator

This is the concrete kernel and the operator part of the spatial-transfer construction in
`transfer.tex`, lines 35–65, pinned at `adc7f1241b42e322a6451854ab7e4b4c146bf78a` in
`openai/math`. The kernel is the trace of the merged chronological word in the real Cartesian
spin-one matrices. Its bound, joint measurability, and symmetry are proved, not assumed.

We construct `K` on complex `L²` and `X a b = exp(-a*b) • K b`. Both are Hilbert–Schmidt,
compact, and self-adjoint. No operator trace, partition-function, moment, rotation-sector, or
twisted-trace identity is asserted here.
-/

noncomputable section

open MeasureTheory
open scoped NNReal ENNReal Matrix.Norms.Operator

namespace CertifiedGaps.SpinOneKernel

open PoissonHistory ContinuousKernel

/-- The complex Hilbert space of finite three-label histories with the unnormalized measure. -/
abbrev Space (b : ℝ≥0) := Lp ℂ 2 (historyMeasure b 3)

/-- The actual spin-one matrix-trace kernel, complexified from its real values. -/
def kernel (b : ℝ≥0) : BoundedKernel (History b 3) ℂ :=
  MergedHistory.boundedKernel SpinOne.G norm_spinOne_G_le_one

/-- The kernel has the paper's explicit merged-chronological trace formula. -/
theorem kernel_apply (b : ℝ≥0) (h k : History b 3) :
    kernel b h k = (MergedHistory.mergedTrace SpinOne.G h k : ℂ) := rfl

/-- The actual spin-one kernel is Hermitian. -/
theorem isHermitian_kernel (b : ℝ≥0) : (kernel b).IsHermitian :=
  MergedHistory.isHermitian_boundedKernel SpinOne.G norm_spinOne_G_le_one

/-- The continuous history integral operator `K_b`. -/
def K (b : ℝ≥0) : Space b →L[ℂ] Space b := (kernel b).operator (historyMeasure b 3)

/-- Its value is the integral against the explicit chronological kernel, almost everywhere. -/
theorem K_apply (b : ℝ≥0) (f : Space b) :
    K b f =ᵐ[historyMeasure b 3] fun h =>
      ∫ k, (MergedHistory.mergedTrace SpinOne.G h k : ℂ) * f k ∂historyMeasure b 3 :=
  (kernel b).coeFn_apply (historyMeasure b 3) f

/-- Real-valued symmetry of the chronological kernel proves self-adjointness of `K_b`. -/
theorem isSelfAdjoint_K (b : ℝ≥0) : IsSelfAdjoint (K b) :=
  (kernel b).isSelfAdjoint_operator (historyMeasure b 3) (isHermitian_kernel b)

/-- The unnormalized history measure has real mass `exp(3*b)`. -/
theorem historyMeasure_real_univ (b : ℝ≥0) :
    (historyMeasure b 3).real Set.univ = Real.exp (3 * b) := by
  rw [measureReal_def, historyMeasure_univ, ENNReal.toReal_ofReal (by positivity)]
  norm_num

/-- The explicit squared Hilbert–Schmidt bound in the spatial-transfer construction. -/
theorem hasHilbertSchmidtBound_K (b : ℝ≥0) :
    HasHilbertSchmidtBound (K b) (9 * Real.exp (6 * b)) := by
  have h := (kernel b).hasHilbertSchmidtBound_operator (historyMeasure b 3)
  have he : Real.exp (3 * (b : ℝ)) ^ 2 = Real.exp (6 * b) := by
    rw [← Real.exp_nat_mul]
    congr 1
    ring_nf
  change HasHilbertSchmidtBound (K b)
    ((3 : ℝ) ^ 2 * ((historyMeasure b 3).real Set.univ) ^ 2) at h
  simpa only [historyMeasure_real_univ, he, show (3 : ℝ) ^ 2 = 9 by norm_num] using h

/-- The actual continuous history operator is Hilbert–Schmidt. -/
theorem isHilbertSchmidt_K (b : ℝ≥0) : IsHilbertSchmidt (K b) :=
  (hasHilbertSchmidtBound_K b).isHilbertSchmidt

/-- The actual continuous history operator is compact. -/
theorem isCompactOperator_K (b : ℝ≥0) : IsCompactOperator (K b) :=
  (isHilbertSchmidt_K b).isCompactOperator

/-- The shifted operator `X_b = exp(-a*b) K_b`, with real shift parameter `a`. -/
def X (a : ℝ) (b : ℝ≥0) : Space b →L[ℂ] Space b := (Real.exp (-a * b) : ℂ) • K b

/-- The normalized operator is integration against the scalar-normalized chronological kernel. -/
theorem X_apply (a : ℝ) (b : ℝ≥0) (f : Space b) :
    X a b f =ᵐ[historyMeasure b 3] fun h =>
      ∫ k, (Real.exp (-a * b) : ℂ) * (MergedHistory.mergedTrace SpinOne.G h k : ℂ) *
        f k ∂historyMeasure b 3 := by
  filter_upwards [Lp.coeFn_smul (Real.exp (-a * b) : ℂ) (K b f), K_apply b f] with h hs hk
  change ((Real.exp (-a * b) : ℂ) • K b f) h = _
  rw [hs, Pi.smul_apply, smul_eq_mul, hk, ← integral_const_mul]
  simp only [mul_assoc]

/-- The empty–empty value is exactly three, including at zero time. -/
theorem kernel_empty (b : ℝ≥0) :
    kernel b (MergedHistory.emptyHistory b 3) (MergedHistory.emptyHistory b 3) = 3 := by
  rw [kernel_apply, MergedHistory.mergedTrace_empty]
  norm_num

/-- A real scalar shift preserves self-adjointness. -/
theorem isSelfAdjoint_X (a : ℝ) (b : ℝ≥0) : IsSelfAdjoint (X a b) :=
  (show IsSelfAdjoint (Real.exp (-a * b) : ℂ) from by
    exact isSelfAdjoint_iff.mpr (Complex.conj_ofReal _)).smul
    (isSelfAdjoint_K b)

/-- Scalar normalization preserves compactness. -/
theorem isCompactOperator_X (a : ℝ) (b : ℝ≥0) : IsCompactOperator (X a b) := by
  exact (isCompactOperator_K b).smul (Real.exp (-a * b) : ℂ)

/-- The normalized operator has the correspondingly rescaled squared Hilbert–Schmidt bound. -/
theorem hasHilbertSchmidtBound_X (a : ℝ) (b : ℝ≥0) :
    HasHilbertSchmidtBound (X a b) (9 * Real.exp ((6 - 2 * a) * b)) := by
  intro n v hv
  have h := hasHilbertSchmidtBound_K b n v hv
  have he : Real.exp (-a * (b : ℝ)) ^ 2 * (9 * Real.exp (6 * b)) =
      9 * Real.exp ((6 - 2 * a) * b) := by
    rw [← Real.exp_nat_mul]
    calc
      _ = 9 * (Real.exp (2 * (-a * b)) * Real.exp (6 * b)) := by ring_nf
      _ = _ := by
        rw [← Real.exp_add, show 2 * (-a * (b : ℝ)) + 6 * b =
          (6 - 2 * a) * b by ring_nf]
  calc
    ∑ i, ‖X a b (v i)‖ ^ 2 =
        Real.exp (-a * b) ^ 2 * ∑ i, ‖K b (v i)‖ ^ 2 := by
      simp only [X, smul_apply, norm_smul, Complex.norm_real, Real.norm_eq_abs,
        abs_of_pos (Real.exp_pos _), mul_pow, Finset.mul_sum]
    _ ≤ Real.exp (-a * b) ^ 2 * (9 * Real.exp (6 * b)) :=
      mul_le_mul_of_nonneg_left h (sq_nonneg _)
    _ = _ := he

/-- The normalized continuous spin-one operator is Hilbert–Schmidt. -/
theorem isHilbertSchmidt_X (a : ℝ) (b : ℝ≥0) : IsHilbertSchmidt (X a b) :=
  (hasHilbertSchmidtBound_X a b).isHilbertSchmidt

end CertifiedGaps.SpinOneKernel
