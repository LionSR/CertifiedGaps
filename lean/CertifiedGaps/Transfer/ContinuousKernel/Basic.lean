import Mathlib.MeasureTheory.Function.L2Space
import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.Analysis.InnerProductSpace.Adjoint

/-!
# Bounded measurable kernels on finite measure spaces

This file constructs the integral operator on `L²` from a bounded jointly measurable kernel.
The construction uses actual integration, and is independent of the transfer matrix and its
moments. Hermitian symmetry of the kernel implies self-adjointness of the resulting operator.
-/

noncomputable section

open MeasureTheory Filter
open scoped ENNReal NNReal ComplexConjugate

namespace CertifiedGaps.ContinuousKernel

variable {α 𝕜 : Type*} [MeasurableSpace α] [RCLike 𝕜]

/-- A jointly strongly measurable scalar kernel with an explicit uniform bound. -/
structure BoundedKernel (α 𝕜 : Type*) [MeasurableSpace α] [RCLike 𝕜] where
  /-- The value of the kernel. -/
  toFun : α → α → 𝕜
  /-- Joint measurability, rather than only separate measurability. -/
  stronglyMeasurable : StronglyMeasurable (Function.uncurry toFun)
  /-- An upper bound on the absolute value. -/
  bound : ℝ≥0
  /-- The bound holds pointwise. -/
  norm_le : ∀ x y, ‖toFun x y‖ ≤ bound

namespace BoundedKernel

instance : CoeFun (BoundedKernel α 𝕜) (fun _ => α → α → 𝕜) := ⟨toFun⟩

variable (K : BoundedKernel α 𝕜) (μ : Measure α) [IsFiniteMeasure μ]

/-- The integral formula, before passing to equivalence classes in `L²`. -/
def applyFun (f : α → 𝕜) (x : α) : 𝕜 := ∫ y, K x y * f y ∂μ

omit [IsFiniteMeasure μ] in
lemma integrable_section {f : α → 𝕜} (hf : Integrable f μ) (x : α) :
    Integrable (fun y => K x y * f y) μ :=
  hf.bdd_mul K.stronglyMeasurable.of_uncurry_left.aestronglyMeasurable
    (Eventually.of_forall (K.norm_le x))

lemma aestronglyMeasurable_applyFun {f : α → 𝕜} (hf : AEStronglyMeasurable f μ) :
    AEStronglyMeasurable (K.applyFun μ f) μ :=
  (K.stronglyMeasurable.aestronglyMeasurable.mul hf.comp_snd).integral_prod_right'

omit [IsFiniteMeasure μ] in
lemma norm_applyFun_le {f : α → 𝕜} (hf : Integrable f μ) (x : α) :
    ‖K.applyFun μ f x‖ ≤ K.bound * ∫ y, ‖f y‖ ∂μ := by
  unfold applyFun
  calc
    ‖∫ y, K x y * f y ∂μ‖ ≤ ∫ y, (K.bound : ℝ) * ‖f y‖ ∂μ :=
      norm_integral_le_of_norm_le (hf.norm.const_mul _) (Eventually.of_forall fun y => by
        rw [norm_mul]
        exact mul_le_mul_of_nonneg_right (K.norm_le x y) (norm_nonneg _))
    _ = _ := integral_const_mul _ _

lemma memLp_applyFun (f : Lp 𝕜 2 μ) : MemLp (K.applyFun μ f) 2 μ :=
  MemLp.of_bound (K.aestronglyMeasurable_applyFun μ (Lp.aestronglyMeasurable f))
    (K.bound * ∫ y, ‖f y‖ ∂μ)
    (Eventually.of_forall (K.norm_applyFun_le μ ((Lp.memLp f).integrable (by norm_num))))

/-- The integral operator as a map of `L²` equivalence classes. -/
def apply (f : Lp 𝕜 2 μ) : Lp 𝕜 2 μ := (K.memLp_applyFun μ f).toLp (K.applyFun μ f)

lemma coeFn_apply (f : Lp 𝕜 2 μ) : K.apply μ f =ᵐ[μ] K.applyFun μ f :=
  MemLp.coeFn_toLp _

omit [IsFiniteMeasure μ] in
lemma applyFun_congr {f g : α → 𝕜} (hfg : f =ᵐ[μ] g) :
    K.applyFun μ f = K.applyFun μ g := by
  funext x
  exact integral_congr_ae (hfg.mono fun y h => by simp only [h])

lemma apply_add (f g : Lp 𝕜 2 μ) : K.apply μ (f + g) = K.apply μ f + K.apply μ g := by
  apply Lp.ext
  filter_upwards [K.coeFn_apply μ (f + g), K.coeFn_apply μ f, K.coeFn_apply μ g,
    Lp.coeFn_add (K.apply μ f) (K.apply μ g)] with x hfg hf hg ha
  rw [ha]
  simp only [Pi.add_apply]
  rw [hf, hg, hfg, K.applyFun_congr μ (Lp.coeFn_add f g)]
  simp only [applyFun, Pi.add_apply, mul_add]
  exact integral_add (K.integrable_section μ ((Lp.memLp f).integrable (by norm_num)) x)
    (K.integrable_section μ ((Lp.memLp g).integrable (by norm_num)) x)

lemma apply_smul (c : 𝕜) (f : Lp 𝕜 2 μ) : K.apply μ (c • f) = c • K.apply μ f := by
  apply Lp.ext
  filter_upwards [K.coeFn_apply μ (c • f), K.coeFn_apply μ f,
    Lp.coeFn_smul c (K.apply μ f)] with x hcf hf hs
  rw [hs]
  simp only [Pi.smul_apply]
  rw [hf, hcf, K.applyFun_congr μ (Lp.coeFn_smul c f)]
  simp only [applyFun, Pi.smul_apply, smul_eq_mul]
  simp_rw [← mul_assoc, mul_comm (K x _ ) c, mul_assoc]
  exact integral_const_mul _ _

/-- Integration against a bounded kernel is linear on `L²`. -/
def linearMap : Lp 𝕜 2 μ →ₗ[𝕜] Lp 𝕜 2 μ where
  toFun := K.apply μ
  map_add' := K.apply_add μ
  map_smul' := K.apply_smul μ

/-- The squared `L²` norm of the constant one is the total mass. -/
lemma norm_const_one_sq : ‖Lp.const 2 μ (1 : ℝ)‖ ^ 2 = μ.real Set.univ := by
  rw [← real_inner_self_eq_norm_sq, L2.inner_def]
  calc
    (∫ x, inner ℝ (Lp.const 2 μ (1 : ℝ) x) (Lp.const 2 μ (1 : ℝ) x) ∂μ) =
        ∫ _ : α, (1 : ℝ) ∂μ := integral_congr_ae ((Lp.coeFn_const 2 μ (1 : ℝ)).mono
          fun x hx => by
            change inner ℝ (Lp.const 2 μ (1 : ℝ) x) (Lp.const 2 μ (1 : ℝ) x) = 1
            rw [hx]
            norm_num)
    _ = _ := by simp

/-- Cauchy–Schwarz bounds the `L¹` norm by the `L²` norm on a finite measure space. -/
lemma integral_norm_le (f : Lp 𝕜 2 μ) :
    (∫ x, ‖f x‖ ∂μ) ≤ ‖Lp.const 2 μ (1 : ℝ)‖ * ‖f‖ := by
  let g := (Lp.memLp f).norm.toLp (fun x => ‖f x‖)
  have hg : ‖g‖ = ‖f‖ := by
    rw [Lp.norm_toLp, eLpNorm_norm _ (Lp.aestronglyMeasurable f)]
    rfl
  have hinner : inner ℝ (Lp.const 2 μ (1 : ℝ)) g = ∫ x, ‖f x‖ ∂μ := by
    rw [L2.inner_def]
    apply integral_congr_ae
    filter_upwards [Lp.coeFn_const 2 μ (1 : ℝ), (Lp.memLp f).norm.coeFn_toLp] with x h1 h2
    change inner ℝ (Lp.const 2 μ (1 : ℝ) x) (g x) = ‖f x‖
    rw [Real.inner_apply, h1]
    change 1 * g x = ‖f x‖
    rw [one_mul]
    exact h2
  calc
    (∫ x, ‖f x‖ ∂μ) ≤ ‖∫ x, ‖f x‖ ∂μ‖ := le_abs_self _
    _ = ‖inner ℝ (Lp.const 2 μ (1 : ℝ)) g‖ := congrArg norm hinner.symm
    _ ≤ ‖Lp.const 2 μ (1 : ℝ)‖ * ‖g‖ := norm_inner_le_norm _ _
    _ = _ := by rw [hg]

lemma norm_apply_le (f : Lp 𝕜 2 μ) :
    ‖K.apply μ f‖ ≤ (K.bound : ℝ) * μ.real Set.univ * ‖f‖ := by
  have hf := (Lp.memLp f).integrable (by norm_num : (1 : ℝ≥0∞) ≤ 2)
  have hI : 0 ≤ ∫ y, ‖f y‖ ∂μ := integral_nonneg fun _ => norm_nonneg _
  have hb : ‖K.apply μ f‖ ≤
      ((K.bound : ℝ) * ∫ y, ‖f y‖ ∂μ) * ‖Lp.const 2 μ (1 : ℝ)‖ := by
    apply Lp.norm_le_mul_norm_of_ae_le_mul
    filter_upwards [K.coeFn_apply μ f, Lp.coeFn_const 2 μ (1 : ℝ)] with x hx h1
    rw [hx, h1]
    simpa only [Function.const_apply, norm_one, mul_one] using K.norm_applyFun_le μ hf x
  calc
    ‖K.apply μ f‖ ≤ ((K.bound : ℝ) * ∫ y, ‖f y‖ ∂μ) *
        ‖Lp.const 2 μ (1 : ℝ)‖ := hb
    _ ≤ (K.bound * (‖Lp.const 2 μ (1 : ℝ)‖ * ‖f‖)) *
        ‖Lp.const 2 μ (1 : ℝ)‖ := by gcongr; exact integral_norm_le μ f
    _ = (K.bound : ℝ) * ‖Lp.const 2 μ (1 : ℝ)‖ ^ 2 * ‖f‖ := by ring
    _ = _ := by rw [norm_const_one_sq μ]

/-- The bounded `L²` operator associated with a measurable bounded kernel. -/
def operator : Lp 𝕜 2 μ →L[𝕜] Lp 𝕜 2 μ :=
  (K.linearMap μ).mkContinuous ((K.bound : ℝ) * μ.real Set.univ) (K.norm_apply_le μ)

lemma operator_apply (f : Lp 𝕜 2 μ) : K.operator μ f = K.apply μ f := rfl

/-- The uniform kernel bound gives an explicit operator-norm bound. -/
lemma norm_operator_le : ‖K.operator μ‖ ≤ (K.bound : ℝ) * μ.real Set.univ :=
  LinearMap.mkContinuous_norm_le _ (mul_nonneg K.bound.coe_nonneg (measureReal_nonneg)) _

/-- Hermitian symmetry is the kernel condition needed for self-adjointness. -/
def IsHermitian : Prop := ∀ x y, conj (K x y) = K y x

lemma integrable_inner_integrand (f g : Lp 𝕜 2 μ) :
    Integrable (fun p : α × α => conj (K p.1 p.2 * f p.2) * g p.1) (μ.prod μ) := by
  have hf := (Lp.memLp f).star.integrable (by norm_num : (1 : ℝ≥0∞) ≤ 2)
  have hg := (Lp.memLp g).integrable (by norm_num : (1 : ℝ≥0∞) ≤ 2)
  have h := (hg.mul_prod hf).bdd_mul (c := (K.bound : ℝ))
    (StronglyMeasurable.star _ K.stronglyMeasurable).aestronglyMeasurable
    (Eventually.of_forall fun p : α × α => by
      change ‖star (K p.1 p.2)‖ ≤ (K.bound : ℝ)
      simpa only [norm_star] using K.norm_le p.1 p.2)
  convert h using 1
  ext p
  change conj (K p.1 p.2 * f p.2) * g p.1 =
    conj (K p.1 p.2) * (g p.1 * conj (f p.2))
  simp only [map_mul]
  ring

/-- The integral formula for the scalar pairing with the kernel operator. -/
lemma inner_operator_left (f g : Lp 𝕜 2 μ) :
    inner 𝕜 (K.operator μ f) g =
      ∫ x, ∫ y, conj (K x y * f y) * g x ∂μ ∂μ := by
  rw [L2.inner_def]
  apply integral_congr_ae
  filter_upwards [K.coeFn_apply μ f] with x hx
  change inner 𝕜 (K.apply μ f x) (g x) = _
  rw [hx, RCLike.inner_apply', applyFun, ← integral_conj, ← integral_mul_const]

/-- Hermitian kernels give symmetric operators; the proof is Fubini with absolute integrability. -/
lemma isSymmetric_operator (hK : K.IsHermitian) : (K.operator μ).toLinearMap.IsSymmetric := by
  intro f g
  change inner 𝕜 (K.operator μ f) g = inner 𝕜 f (K.operator μ g)
  rw [K.inner_operator_left μ, integral_integral_swap (K.integrable_inner_integrand μ f g),
    L2.inner_def]
  apply integral_congr_ae
  filter_upwards [K.coeFn_apply μ g] with y hy
  change (∫ x, conj (K x y * f y) * g x ∂μ) = inner 𝕜 (f y) (K.apply μ g y)
  rw [hy, RCLike.inner_apply', applyFun, ← integral_const_mul]
  apply integral_congr_ae
  exact Eventually.of_forall fun x => by
    simp only [map_mul, hK x y]
    ring

/-- The integral operator of a Hermitian bounded measurable kernel is self-adjoint. -/
lemma isSelfAdjoint_operator (hK : K.IsHermitian) : IsSelfAdjoint (K.operator μ) :=
  ContinuousLinearMap.isSelfAdjoint_iff_isSymmetric.mpr (K.isSymmetric_operator μ hK)

end BoundedKernel
end CertifiedGaps.ContinuousKernel
