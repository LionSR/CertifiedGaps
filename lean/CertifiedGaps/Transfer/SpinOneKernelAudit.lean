import CertifiedGaps.Transfer.SpinOneKernel

/-!
Regression checks for the actual merged-history kernel and normalized complex `L²` operator.
The three-event example distinguishes increasing-time products from their reversals.
-/

noncomputable section

open MeasureTheory CertifiedGaps CertifiedGaps.PoissonHistory CertifiedGaps.MergedHistory
open scoped NNReal ENNReal Matrix.Norms.Operator

namespace CertifiedGaps.SpinOneKernelAudit

/-- A two-event history with another event strictly between its times. -/
def leftHistory : History 4 3 :=
  ⟨⟨2, ![0, 2]⟩, ⟨![1, 3], by
    constructor
    · intro i _
      fin_cases i <;> norm_num
    · intro i j hij
      fin_cases i <;> fin_cases j <;> norm_num at *⟩⟩

/-- The event interleaved between the two events of `leftHistory`. -/
def rightHistory : History 4 3 :=
  ⟨⟨1, ![1]⟩, ⟨![2], by
    constructor
    · intro i _
      fin_cases i
      norm_num
    · intro i j hij
      fin_cases i
      fin_cases j
      norm_num at hij⟩⟩

/-- Enumeration `[left 0, right 0, left 1]` in increasing merged time. -/
def chronologicalEnumeration : Fin 3 ≃ Fin 2 ⊕ Fin 1 :=
  (Equiv.swap 1 2).trans (finSumFinEquiv : Fin 2 ⊕ Fin 1 ≃ Fin (2 + 1)).symm

@[simp] theorem chronologicalEnumeration_zero : chronologicalEnumeration 0 = Sum.inl 0 := by
  decide

@[simp] theorem chronologicalEnumeration_one : chronologicalEnumeration 1 = Sum.inr 0 := by
  decide

@[simp] theorem chronologicalEnumeration_two : chronologicalEnumeration 2 = Sum.inl 1 := by
  decide

-- The real chronological word is genuinely noncommutative: reversal changes its trace.
example : Matrix.trace (SpinOne.G 0 * SpinOne.G 1 * SpinOne.G 2) = 1 := by
  norm_num [SpinOne.G, SpinOneHeisenberg.leviCivita, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_succ]

example : Matrix.trace (SpinOne.G 2 * SpinOne.G 1 * SpinOne.G 0) = -1 := by
  norm_num [SpinOne.G, SpinOneHeisenberg.leviCivita, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_succ]

-- The concrete merged-history kernel uses the original chronological order.
theorem interleaved_trace : mergedTrace SpinOne.G leftHistory rightHistory = 1 := by
  change chamberTrace 3 (Sum.elim (SpinOne.G ∘ ![0, 2]) (SpinOne.G ∘ ![1]))
    (Sum.elim ![1, 3] ![2]) = 1
  rw [chamberTrace_eq _ _ chronologicalEnumeration (by
    intro i j hij
    fin_cases i <;> fin_cases j <;>
      norm_num [Function.comp_def] at *)]
  norm_num [chronologicalProduct, leftHistory, rightHistory,
    Function.comp_def, List.ofFn_succ,
    SpinOne.G, SpinOneHeisenberg.leviCivita, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_succ]

-- Swapping the histories preserves the same noncommuting chronological word.
example : mergedTrace SpinOne.G rightHistory leftHistory = 1 := by
  rw [mergedTrace_swap]
  exact interleaved_trace

-- Equal event times on opposite histories are assigned zero.
example : mergedTrace SpinOne.G rightHistory rightHistory = 0 :=
  mergedTrace_eq_zero_of_collision _ _ _ (i := (0 : Fin 1)) (j := (0 : Fin 1)) rfl

-- Empty histories retain the trace of the identity even at zero time.
example (b : ℝ≥0) :
    SpinOneKernel.kernel b (emptyHistory b 3) (emptyHistory b 3) = 3 :=
  SpinOneKernel.kernel_empty b

example : SpinOneKernel.kernel 0 (emptyHistory 0 3) (emptyHistory 0 3) = 3 :=
  SpinOneKernel.kernel_empty 0

-- The public objects act on complex L² with the actual unnormalized history measure.
example (b : ℝ≥0) :
    IsSelfAdjoint (SpinOneKernel.K b : Lp ℂ 2 (historyMeasure b 3) →L[ℂ]
      Lp ℂ 2 (historyMeasure b 3)) := SpinOneKernel.isSelfAdjoint_K b

example (b : ℝ≥0) : SpinOneKernel.X 0 b = SpinOneKernel.K b := by
  simp [SpinOneKernel.X]

example : ContinuousKernel.HasHilbertSchmidtBound (SpinOneKernel.K 0) 9 := by
  simpa using SpinOneKernel.hasHilbertSchmidtBound_K 0

example (b : ℝ≥0) : ContinuousKernel.HasHilbertSchmidtBound (SpinOneKernel.X 3 b) 9 := by
  convert SpinOneKernel.hasHilbertSchmidtBound_X 3 b using 1
  norm_num

example (a : ℝ) (b : ℝ≥0) : IsCompactOperator (SpinOneKernel.X a b) :=
  SpinOneKernel.isCompactOperator_X a b

end CertifiedGaps.SpinOneKernelAudit

#print axioms CertifiedGaps.PoissonHistory.ae_merged_history_injective
#print axioms CertifiedGaps.MergedHistory.mergedTrace_swap
#print axioms CertifiedGaps.MergedHistory.measurable_mergedTrace
#print axioms CertifiedGaps.MergedHistory.norm_mergedTrace_le
#print axioms CertifiedGaps.SpinOneKernel.K_apply
#print axioms CertifiedGaps.SpinOneKernel.X_apply
#print axioms CertifiedGaps.SpinOneKernel.hasHilbertSchmidtBound_K
#print axioms CertifiedGaps.SpinOneKernel.isSelfAdjoint_X
#print axioms CertifiedGaps.SpinOneKernel.isCompactOperator_X
#print axioms CertifiedGaps.SpinOneKernel.hasHilbertSchmidtBound_X

-- Enforce the transitive axiom whitelist on every declaration in the construction and its
-- analytic dependencies. This rejects `sorryAx`, custom axioms, and native-decision oracles.
open Lean in
run_cmd do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let prefixes : List Name := [`CertifiedGaps.SigmaProductMeasurable,
    `CertifiedGaps.TimeOrdering, `CertifiedGaps.PoissonHistory, `CertifiedGaps.MergedHistory,
    `CertifiedGaps.ContinuousKernel, `CertifiedGaps.SpinOneKernel,
    `CertifiedGaps.SpinOneKernelAudit]
  let mut checked : Nat := 0
  for (name, _) in (← getEnv).constants.toList do
    if prefixes.any (·.isPrefixOf name) then
      checked := checked + 1
      let axioms ← collectAxioms name
      let bad := axioms.filter (!allowed.contains ·)
      unless bad.isEmpty do
        throwError "{name} depends on disallowed axioms: {bad}"
  logInfo m!"Merged-history kernel axiom whitelist checked {checked} declarations."
