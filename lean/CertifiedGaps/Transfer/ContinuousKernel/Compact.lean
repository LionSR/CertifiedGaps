import CertifiedGaps.Transfer.ContinuousKernel.Basic
import CertifiedGaps.Transfer.ContinuousKernel.HilbertSchmidt

/-!
# Hilbert–Schmidt bounds and compactness for bounded measurable kernels

The squared norms of the images of every orthonormal family have a uniform finite-sum bound.
This is proved by applying Bessel's inequality to each section of the kernel and integrating.
-/

noncomputable section
open MeasureTheory Filter
open scoped ENNReal NNReal ComplexConjugate

namespace CertifiedGaps.ContinuousKernel.BoundedKernel

variable {α 𝕜 ι : Type*} [MeasurableSpace α] [RCLike 𝕜]
variable (K : BoundedKernel α 𝕜) (μ : Measure α) [IsFiniteMeasure μ]

/-- The conjugated row of the kernel as an `L²` vector. -/
def row (x : α) : Lp 𝕜 2 μ :=
  (MemLp.of_bound
    (StronglyMeasurable.star _ K.stronglyMeasurable.of_uncurry_left).aestronglyMeasurable
    (K.bound : ℝ) (Eventually.of_forall fun y => by
      change ‖star (K x y)‖ ≤ (K.bound : ℝ)
      simpa only [norm_star] using K.norm_le x y)).toLp (fun y => conj (K x y))

lemma coeFn_row (x : α) : K.row μ x =ᵐ[μ] fun y => conj (K x y) :=
  MemLp.coeFn_toLp _

lemma applyFun_eq_inner_row (f : Lp 𝕜 2 μ) (x : α) :
    K.applyFun μ f x = inner 𝕜 (K.row μ x) f := by
  rw [applyFun, L2.inner_def]
  apply integral_congr_ae
  filter_upwards [K.coeFn_row μ x] with y hy
  change K x y * f y = inner 𝕜 (K.row μ x y) (f y)
  rw [RCLike.inner_apply', hy]
  simp only [starRingEnd_apply, star_star]

lemma norm_row_sq_le (x : α) :
    ‖K.row μ x‖ ^ 2 ≤ (K.bound : ℝ) ^ 2 * μ.real Set.univ := by
  have h : ‖K.row μ x‖ ≤ (K.bound : ℝ) * ‖Lp.const 2 μ (1 : ℝ)‖ := by
    apply Lp.norm_le_mul_norm_of_ae_le_mul
    filter_upwards [K.coeFn_row μ x, Lp.coeFn_const 2 μ (1 : ℝ)] with y hy h1
    rw [hy, h1]
    simpa only [RCLike.norm_conj, Function.const_apply, norm_one, mul_one] using K.norm_le x y
  calc
    ‖K.row μ x‖ ^ 2 ≤ ((K.bound : ℝ) * ‖Lp.const 2 μ (1 : ℝ)‖) ^ 2 := by
      exact pow_le_pow_left₀ (norm_nonneg _) h 2
    _ = (K.bound : ℝ) ^ 2 * μ.real Set.univ := by rw [mul_pow, norm_const_one_sq μ]

omit [IsFiniteMeasure μ] in
lemma norm_sq_eq_integral (f : Lp 𝕜 2 μ) : ‖f‖ ^ 2 = ∫ x, ‖f x‖ ^ 2 ∂μ := by
  rw [InnerProductSpace.norm_sq_eq_re_inner (𝕜 := 𝕜), L2.inner_def,
    ← integral_re (L2.integrable_inner f f)]
  simp only [inner_self_eq_norm_sq]

lemma norm_operator_sq_eq_integral (f : Lp 𝕜 2 μ) :
    ‖K.operator μ f‖ ^ 2 = ∫ x, ‖K.applyFun μ f x‖ ^ 2 ∂μ := by
  rw [norm_sq_eq_integral μ]
  exact integral_congr_ae ((K.coeFn_apply μ f).mono fun x hx => by
    change ‖K.apply μ f x‖ ^ 2 = ‖K.applyFun μ f x‖ ^ 2
    rw [hx])

lemma integrable_norm_applyFun_sq (f : Lp 𝕜 2 μ) :
    Integrable (fun x => ‖K.applyFun μ f x‖ ^ 2) μ :=
  (K.memLp_applyFun μ f).integrable_norm_pow (by norm_num)

/-- Bessel's inequality integrated over the rows gives a dimension-independent HS bound. -/
lemma sum_norm_operator_sq_le (v : ι → Lp 𝕜 2 μ) (hv : Orthonormal 𝕜 v) (s : Finset ι) :
    ∑ i ∈ s, ‖K.operator μ (v i)‖ ^ 2 ≤
      (K.bound : ℝ) ^ 2 * (μ.real Set.univ) ^ 2 := by
  classical
  calc
    ∑ i ∈ s, ‖K.operator μ (v i)‖ ^ 2 =
        ∫ x, ∑ i ∈ s, ‖K.applyFun μ (v i) x‖ ^ 2 ∂μ := by
      simp_rw [K.norm_operator_sq_eq_integral μ]
      exact (integral_finsetSum s (fun i _ => K.integrable_norm_applyFun_sq μ (v i))).symm
    _ ≤ ∫ _ : α, (K.bound : ℝ) ^ 2 * μ.real Set.univ ∂μ := by
      apply integral_mono
      · exact integrable_finsetSum _ (fun i _ => K.integrable_norm_applyFun_sq μ (v i))
      · exact integrable_const _
      · intro x
        simp_rw [K.applyFun_eq_inner_row μ, norm_inner_symm]
        exact (hv.sum_inner_products_le (K.row μ x)).trans (K.norm_row_sq_le μ x)
    _ = _ := by rw [integral_const]; simp only [smul_eq_mul]; ring

/-- Every orthonormal family has square-summable images under the kernel operator. -/
lemma summable_norm_operator_sq (v : ι → Lp 𝕜 2 μ) (hv : Orthonormal 𝕜 v) :
    Summable (fun i => ‖K.operator μ (v i)‖ ^ 2) :=
  summable_of_sum_le (fun _ => sq_nonneg _) (K.sum_norm_operator_sq_le μ v hv)

/-- The square sum is bounded by the squared uniform kernel bound times squared total mass. -/
lemma tsum_norm_operator_sq_le (v : ι → Lp 𝕜 2 μ) (hv : Orthonormal 𝕜 v) :
    ∑' i, ‖K.operator μ (v i)‖ ^ 2 ≤ (K.bound : ℝ) ^ 2 * (μ.real Set.univ) ^ 2 :=
  tsum_le_of_sum_le' (mul_nonneg (sq_nonneg _) (sq_nonneg _)) (K.sum_norm_operator_sq_le μ v hv)

/-- The explicit bound on the squared Hilbert–Schmidt norm of the kernel operator. -/
lemma hasHilbertSchmidtBound_operator :
    HasHilbertSchmidtBound (K.operator μ) ((K.bound : ℝ) ^ 2 * (μ.real Set.univ) ^ 2) := by
  intro n v hv
  exact K.sum_norm_operator_sq_le μ v hv Finset.univ

/-- The integral operator of a bounded measurable kernel is Hilbert–Schmidt. -/
lemma isHilbertSchmidt_operator : IsHilbertSchmidt (K.operator μ) :=
  (K.hasHilbertSchmidtBound_operator μ).isHilbertSchmidt

/-- A bounded measurable kernel on a finite measure space defines a compact operator. -/
lemma isCompactOperator_operator : IsCompactOperator (K.operator μ) :=
  (K.isHilbertSchmidt_operator μ).isCompactOperator

end CertifiedGaps.ContinuousKernel.BoundedKernel
