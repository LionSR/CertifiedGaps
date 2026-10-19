import CertifiedGaps.Local.SpinOne
import Mathlib.Analysis.Matrix.Normed

/-!
# Matrix-word trace bounds

The induced maximum-row-sum norm is submultiplicative. The real Cartesian spin-one factors
are contractions in this norm, giving an absolute trace bound of three for every finite
matrix word in its original order.
-/

noncomputable section
open scoped Matrix.Norms.Operator NNReal

namespace CertifiedGaps.ContinuousKernel

variable {ι Λ : Type*} [Fintype ι] [DecidableEq ι]

/-- An entry is bounded by the induced maximum-row-sum matrix norm. -/
lemma norm_entry_le_opNorm (A : Matrix ι ι ℝ) (i j : ι) : ‖A i j‖ ≤ ‖A‖ := by
  have h : ‖A i j‖₊ ≤ ‖A‖₊ := by
    rw [Matrix.linfty_opNNNorm_def]
    exact (Finset.single_le_sum (fun k _ => (by positivity : (0 : ℝ≥0) ≤ ‖A i k‖₊))
      (Finset.mem_univ j)).trans
      (Finset.le_sup (f := fun i : ι => ∑ j : ι, ‖A i j‖₊) (Finset.mem_univ i))
  exact_mod_cast h

/-- The trace is bounded by the matrix size times a submultiplicative operator norm. -/
lemma norm_trace_le_card_mul_opNorm (A : Matrix ι ι ℝ) :
    ‖Matrix.trace A‖ ≤ Fintype.card ι * ‖A‖ := by
  calc
    ‖Matrix.trace A‖ ≤ ∑ i : ι, ‖A i i‖ := norm_sum_le _ _
    _ ≤ ∑ _ : ι, ‖A‖ := Finset.sum_le_sum fun i _ => norm_entry_le_opNorm A i i
    _ = _ := by simp [nsmul_eq_mul]

/-- Products of contractions remain contractions, in their original list order. -/
lemma norm_prod_le_one [Nonempty ι] (G : Λ → Matrix ι ι ℝ)
    (hG : ∀ a, ‖G a‖ ≤ 1) (l : List Λ) : ‖(l.map G).prod‖ ≤ 1 := by
  induction l with
  | nil => simp
  | cons a l ih =>
    simp only [List.map_cons, List.prod_cons]
    exact (norm_mul_le _ _).trans
      ((mul_le_mul (hG a) ih (norm_nonneg _) zero_le_one).trans_eq (one_mul 1))

/-- Every chronological product of contractions has trace bounded by the local dimension. -/
lemma norm_trace_prod_le [Nonempty ι] (G : Λ → Matrix ι ι ℝ)
    (hG : ∀ a, ‖G a‖ ≤ 1) (l : List Λ) :
    ‖Matrix.trace (l.map G).prod‖ ≤ Fintype.card ι := by
  have hcard : (0 : ℝ) ≤ Fintype.card ι := by positivity
  have hmul := mul_le_mul_of_nonneg_left (norm_prod_le_one G hG l) hcard
  exact (norm_trace_le_card_mul_opNorm _).trans (by simpa using hmul)

/-- The real Cartesian spin-one factors are contractions in the row-sum operator norm. -/
lemma norm_spinOne_G_le_one (a : Fin 3) : ‖SpinOne.G a‖ ≤ 1 := by
  rw [Matrix.linfty_opNorm_def]
  norm_cast
  apply Finset.sup_le
  intro i _
  fin_cases a <;> fin_cases i <;>
    norm_num [SpinOne.G, SpinOneHeisenberg.leviCivita, Fin.sum_univ_succ]

/-- The merged chronological spin-one trace is uniformly bounded by three. -/
lemma norm_trace_spinOne_prod_le_three (l : List (Fin 3)) :
    ‖Matrix.trace (l.map SpinOne.G).prod‖ ≤ 3 := by
  simpa using norm_trace_prod_le SpinOne.G norm_spinOne_G_le_one l

end CertifiedGaps.ContinuousKernel
