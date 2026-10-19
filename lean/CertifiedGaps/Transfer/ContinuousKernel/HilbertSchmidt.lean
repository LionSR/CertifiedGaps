import Mathlib.Analysis.InnerProductSpace.l2Space
import Mathlib.Analysis.Normed.Operator.Compact.Basic
import Mathlib.Analysis.Normed.Group.InfiniteSum
import Mathlib.Topology.Algebra.InfiniteSum.ContinuousEval

/-!
# Square-summable basis images give compact operators

An operator whose images of a Hilbert basis have summable squared norms is compact. The proof
constructs the operator-norm convergent series of rank-one truncations, using Cauchy–Schwarz and
Bessel's inequality. The basis index type is arbitrary; no separability assumption is needed.

The predicates `HasHilbertSchmidtBound` and `IsHilbertSchmidt` record the finite orthonormal-family
criterion. The constant in `HasHilbertSchmidtBound T C` bounds the squared Hilbert–Schmidt norm;
we do not construct a Schatten norm or assert that any given bound is optimal.
-/

noncomputable section

open Filter
open scoped Topology

namespace CertifiedGaps.ContinuousKernel

variable {𝕜 ι E F : Type*} [RCLike 𝕜]
  [NormedAddCommGroup E] [InnerProductSpace 𝕜 E]
  [NormedAddCommGroup F]

section HilbertCodomain

variable [InnerProductSpace 𝕜 F]

/-- `C` bounds the squared Hilbert–Schmidt norm of `T`: the squared image norms of every
finite orthonormal family sum to at most `C`. -/
def HasHilbertSchmidtBound (T : E →L[𝕜] F) (C : ℝ) : Prop :=
  ∀ (n : ℕ) (v : Fin n → E), Orthonormal 𝕜 v → ∑ i, ‖T (v i)‖ ^ 2 ≤ C

/-- An operator is Hilbert–Schmidt if its squared image norms have a uniform finite bound over
all finite orthonormal families. -/
def IsHilbertSchmidt (T : E →L[𝕜] F) : Prop := ∃ C, HasHilbertSchmidtBound T C

/-- An explicit squared Hilbert–Schmidt bound implies the Hilbert–Schmidt property. -/
lemma HasHilbertSchmidtBound.isHilbertSchmidt {T : E →L[𝕜] F} {C : ℝ}
    (hT : HasHilbertSchmidtBound T C) : IsHilbertSchmidt T := ⟨C, hT⟩

end HilbertCodomain

section NormedCodomain

variable [NormedSpace 𝕜 F]

/-- Bessel's inequality bounds a finite sum of rank-one operators by the square sum of its
output vectors, when its input vectors are orthonormal. -/
lemma norm_sum_rankOne_le (v : ι → E) (hv : Orthonormal 𝕜 v) (w : ι → F)
    (s : Finset ι) :
    ‖∑ i ∈ s, InnerProductSpace.rankOne 𝕜 (w i) (v i)‖ ≤
      Real.sqrt (∑ i ∈ s, ‖w i‖ ^ 2) := by
  classical
  apply ContinuousLinearMap.opNorm_le_bound _ (Real.sqrt_nonneg _)
  intro x
  calc
    ‖(∑ i ∈ s, InnerProductSpace.rankOne 𝕜 (w i) (v i)) x‖ =
        ‖∑ i ∈ s, inner 𝕜 (v i) x • w i‖ := by simp
    _ ≤ ∑ i ∈ s, ‖inner 𝕜 (v i) x‖ * ‖w i‖ := by
      simpa only [norm_smul] using norm_sum_le s (fun i => inner 𝕜 (v i) x • w i)
    _ ≤ Real.sqrt (∑ i ∈ s, ‖inner 𝕜 (v i) x‖ ^ 2) *
        Real.sqrt (∑ i ∈ s, ‖w i‖ ^ 2) :=
      Real.sum_mul_le_sqrt_mul_sqrt s _ _
    _ ≤ ‖x‖ * Real.sqrt (∑ i ∈ s, ‖w i‖ ^ 2) := by
      apply mul_le_mul_of_nonneg_right _ (Real.sqrt_nonneg _)
      exact (Real.sqrt_le_sqrt (hv.sum_inner_products_le x)).trans_eq
        (Real.sqrt_sq (norm_nonneg _))
    _ = Real.sqrt (∑ i ∈ s, ‖w i‖ ^ 2) * ‖x‖ := mul_comm _ _

/-- Square-summable output vectors and orthonormal input vectors define an operator-norm
summable family of rank-one operators. -/
lemma summable_rankOne [CompleteSpace F] (v : ι → E) (hv : Orthonormal 𝕜 v)
    (w : ι → F) (hw : Summable fun i => ‖w i‖ ^ 2) :
    Summable fun i => InnerProductSpace.rankOne 𝕜 (w i) (v i) := by
  classical
  apply summable_iff_vanishing_norm.mpr
  intro ε hε
  obtain ⟨s, hs⟩ := summable_iff_vanishing_norm.mp hw (ε ^ 2) (sq_pos_of_pos hε)
  refine ⟨s, fun t ht => ?_⟩
  apply (norm_sum_rankOne_le v hv w t).trans_lt
  apply (Real.sqrt_lt' hε).mpr
  exact (le_abs_self _).trans_lt (hs t ht)

/-- A rank-one operator is compact because it factors through the scalar field. -/
lemma isCompactOperator_rankOne (x : F) (y : E) :
    IsCompactOperator (InnerProductSpace.rankOne 𝕜 x y) := by
  rw [InnerProductSpace.rankOne_def']
  exact (isCompactOperator_of_locallyCompactSpace_rng
    (ContinuousLinearMap.toSpanSingleton 𝕜 x)).comp_clm (innerSL 𝕜 y)

/-- The rank-one expansion of an operator with square-summable basis images converges in
operator norm to that operator. -/
lemma hasSum_rankOne_of_summable_norm_sq [CompleteSpace F]
    (b : HilbertBasis ι 𝕜 E) (T : E →L[𝕜] F)
    (hT : Summable fun i => ‖T (b i)‖ ^ 2) :
    HasSum (fun i => InnerProductSpace.rankOne 𝕜 (T (b i)) (b i)) T := by
  let A : ι → E →L[𝕜] F := fun i => InnerProductSpace.rankOne 𝕜 (T (b i)) (b i)
  have hA : Summable A := summable_rankOne b b.orthonormal _ hT
  have heq : ∑' i, A i = T := by
    apply ContinuousLinearMap.ext
    intro x
    have h1 := hasSum_apply hA.hasSum x
    have h2 := (b.hasSum_repr x).mapL T
    apply h1.unique
    simpa only [A, InnerProductSpace.rankOne_apply, map_smul,
      HilbertBasis.repr_apply_apply] using h2
  convert hA.hasSum using 1
  exact heq.symm

/-- An operator with square-summable images of a Hilbert basis is compact. -/
lemma isCompactOperator_of_summable_norm_sq [CompleteSpace F]
    (b : HilbertBasis ι 𝕜 E) (T : E →L[𝕜] F)
    (hT : Summable fun i => ‖T (b i)‖ ^ 2) : IsCompactOperator T := by
  classical
  apply isCompactOperator_of_tendsto (hasSum_rankOne_of_summable_norm_sq b T hT)
  apply Filter.Eventually.of_forall
  intro s
  induction s using Finset.induction_on with
  | empty => simpa using (isCompactOperator_zero : IsCompactOperator (0 : E →L[𝕜] F))
  | @insert i s hi ih =>
    simpa [Finset.sum_insert hi] using (isCompactOperator_rankOne (T (b i)) (b i)).add ih

/-- A uniform bound on the finite sums of squared basis-image norms implies compactness. -/
lemma isCompactOperator_of_sum_norm_sq_le [CompleteSpace F]
    (b : HilbertBasis ι 𝕜 E) (T : E →L[𝕜] F) (C : ℝ)
    (hT : ∀ s : Finset ι, ∑ i ∈ s, ‖T (b i)‖ ^ 2 ≤ C) : IsCompactOperator T :=
  isCompactOperator_of_summable_norm_sq b T (summable_of_sum_le (fun _ => sq_nonneg _) hT)

/-- A bound on the squared image norms of every finite orthonormal family implies compactness.
This criterion does not require selecting a Hilbert basis in its hypotheses. -/
lemma isCompactOperator_of_orthonormal_sum_norm_sq_le [CompleteSpace E] [CompleteSpace F]
    (T : E →L[𝕜] F) (C : ℝ)
    (hT : ∀ (n : ℕ) (v : Fin n → E), Orthonormal 𝕜 v →
      ∑ i, ‖T (v i)‖ ^ 2 ≤ C) : IsCompactOperator T := by
  classical
  obtain ⟨w, b, _⟩ := exists_hilbertBasis 𝕜 E
  apply isCompactOperator_of_sum_norm_sq_le b T C
  intro s
  let e : Fin (Fintype.card s) ≃ s := (Fintype.equivFin s).symm
  have hv : Orthonormal 𝕜 (fun i => b (e i)) :=
    b.orthonormal.comp (fun i => (e i).val) (Subtype.val_injective.comp e.injective)
  have h := hT _ _ hv
  have he : (∑ i : Fin (Fintype.card s), ‖T (b (e i))‖ ^ 2) =
      ∑ i ∈ s, ‖T (b i)‖ ^ 2 := by
    rw [e.sum_comp (fun i : s => ‖T (b i)‖ ^ 2)]
    exact Finset.sum_coe_sort s (fun i => ‖T (b i)‖ ^ 2)
  rwa [he] at h

end NormedCodomain

section HilbertCodomain

variable [InnerProductSpace 𝕜 F]

/-- A Hilbert–Schmidt bound implies compactness when both spaces are complete. -/
lemma HasHilbertSchmidtBound.isCompactOperator [CompleteSpace E] [CompleteSpace F]
    {T : E →L[𝕜] F} {C : ℝ} (hT : HasHilbertSchmidtBound T C) : IsCompactOperator T :=
  isCompactOperator_of_orthonormal_sum_norm_sq_le T C hT

/-- Hilbert–Schmidt operators on complete spaces are compact. -/
lemma IsHilbertSchmidt.isCompactOperator [CompleteSpace E] [CompleteSpace F]
    {T : E →L[𝕜] F} (hT : IsHilbertSchmidt T) : IsCompactOperator T := by
  obtain ⟨C, hC⟩ := hT
  exact hC.isCompactOperator

end HilbertCodomain

end CertifiedGaps.ContinuousKernel
