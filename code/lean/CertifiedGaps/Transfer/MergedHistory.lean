import CertifiedGaps.Transfer.HistoryCollision
import CertifiedGaps.Transfer.SigmaProductMeasurable
import CertifiedGaps.Transfer.ContinuousKernel.MatrixBounds
import CertifiedGaps.Transfer.ContinuousKernel.Compact

/-!
# Chronologically merged history kernels

A finite chronological chamber determines the ordered matrix word of the union of two
histories. Summing its matrix trace over the disjoint chambers makes the time dependence
explicitly measurable. Colliding times are assigned zero. Swapping histories relabels events,
without reversing or commuting any matrix factors.
-/

noncomputable section

open MeasureTheory
open scoped NNReal Matrix.Norms.Operator

namespace CertifiedGaps.MergedHistory

/-- A finite matrix word, multiplied from earliest to latest event. -/
def chronologicalProduct {E : Type*} {d N : ℕ} (A : E → Matrix (Fin d) (Fin d) ℝ)
    (σ : Fin N ≃ E) : Matrix (Fin d) (Fin d) ℝ :=
  (List.ofFn (A ∘ σ)).prod

/-- Matrix-trace value on the unique chronological chamber; zero on time collisions. -/
def chamberTrace {E : Type*} [Fintype E] {d : ℕ} (N : ℕ)
    (A : E → Matrix (Fin d) (Fin d) ℝ) (t : E → ℝ) : ℝ := by
  classical
  exact ∑ σ : Fin N ≃ E,
    if StrictMono (t ∘ σ) then Matrix.trace (chronologicalProduct A σ) else 0

/-- Relabeling events preserves the chronological product exactly. -/
theorem chronologicalProduct_reindex {E F : Type*} {d N : ℕ}
    (e : E ≃ F) (A : F → Matrix (Fin d) (Fin d) ℝ) (σ : Fin N ≃ E) :
    chronologicalProduct (A ∘ e) σ = chronologicalProduct A (σ.trans e) := rfl

/-- Relabeling events preserves the chamber trace, without reordering the matrix word. -/
theorem chamberTrace_reindex {E F : Type*} [Fintype E] [Fintype F] {d N : ℕ}
    (e : E ≃ F) (A : F → Matrix (Fin d) (Fin d) ℝ) (t : F → ℝ) :
    chamberTrace N (A ∘ e) (t ∘ e) = chamberTrace N A t := by
  classical
  let q : (Fin N ≃ E) ≃ (Fin N ≃ F) := Equiv.equivCongr (Equiv.refl _) e
  unfold chamberTrace
  exact q.sum_comp (fun σ => if StrictMono (t ∘ σ) then
    Matrix.trace (chronologicalProduct A σ) else 0)

/-- Only one enumeration can put distinct event times in increasing order. -/
theorem chronological_unique {E : Type*} {N : ℕ} (t : E → ℝ)
    {σ τ : Fin N ≃ E} (hσ : StrictMono (t ∘ σ)) (hτ : StrictMono (t ∘ τ)) :
    σ = τ := by
  have heq := Tuple.unique_monotone (f := t ∘ σ) (σ := Equiv.refl _)
    (τ := τ.trans σ.symm) (by simpa using hσ.monotone)
    (by simpa [Function.comp_def] using hτ.monotone)
  have hinj : Function.Injective t := hσ.injective.of_comp_right σ.surjective
  apply Equiv.ext
  intro i
  apply hinj
  simpa [Function.comp_def] using congrFun heq i

/-- On any chronological chamber the finite sum is exactly its matrix trace. -/
theorem chamberTrace_eq {E : Type*} [Fintype E] {d N : ℕ}
    (A : E → Matrix (Fin d) (Fin d) ℝ) (t : E → ℝ)
    (σ : Fin N ≃ E) (hσ : StrictMono (t ∘ σ)) :
    chamberTrace N A t = Matrix.trace (chronologicalProduct A σ) := by
  classical
  unfold chamberTrace
  rw [Finset.sum_eq_single σ]
  · simp [hσ]
  · intro τ _ hτ
    simp only [ite_eq_right_iff]
    intro ht
    exact (hτ (chronological_unique t ht hσ)).elim
  · simp

/-- With no strict chronological enumeration, the chamber trace vanishes. -/
theorem chamberTrace_eq_zero {E : Type*} [Fintype E] {d N : ℕ}
    (A : E → Matrix (Fin d) (Fin d) ℝ) (t : E → ℝ)
    (h : ¬ ∃ σ : Fin N ≃ E, StrictMono (t ∘ σ)) : chamberTrace N A t = 0 := by
  classical
  apply Finset.sum_eq_zero
  intro σ _
  exact ite_eq_right (fun hs => h ⟨σ, hs⟩)

/-- Every finite-degree chamber trace is measurable in all event times jointly. -/
theorem measurable_chamberTrace {E : Type*} [Fintype E] {d N : ℕ}
    (A : E → Matrix (Fin d) (Fin d) ℝ) : Measurable (chamberTrace N A) := by
  classical
  apply Finset.measurable_sum
  intro σ _
  apply Measurable.ite _ measurable_const measurable_const
  exact (TimeOrdering.measurableSet_strictMono N).preimage
    (measurable_pi_iff.mpr fun i => measurable_pi_apply (σ i))

/-- A collision-free event family gives exactly its sorted matrix word. -/
theorem chamberTrace_eq_sorted {E : Type*} [Fintype E] {d N : ℕ}
    (A : E → Matrix (Fin d) (Fin d) ℝ) (t : E → ℝ)
    (e : Fin N ≃ E) (ht : Function.Injective t) :
    chamberTrace N A t = Matrix.trace
      (chronologicalProduct A ((Tuple.sort (t ∘ e)).trans e)) := by
  apply chamberTrace_eq
  exact (Tuple.monotone_sort (t ∘ e)).strictMono_of_injective
    (ht.comp (e.injective.comp (Tuple.sort (t ∘ e)).injective))

/-- A collision among distinct events makes every strict chamber inactive. -/
theorem chamberTrace_eq_zero_of_collision {E : Type*} [Fintype E] {d N : ℕ}
    (A : E → Matrix (Fin d) (Fin d) ℝ) (t : E → ℝ)
    (ht : ¬ Function.Injective t) : chamberTrace N A t = 0 := by
  apply chamberTrace_eq_zero
  rintro ⟨σ, hσ⟩
  exact ht (hσ.injective.of_comp_right σ.surjective)

/-- A uniform bound on the finite matrix words bounds the chamber trace without a factorial
loss, since at most one chamber is active. -/
theorem norm_chamberTrace_le {E : Type*} [Fintype E] {d N : ℕ}
    (A : E → Matrix (Fin d) (Fin d) ℝ) (t : E → ℝ) {B : ℝ} (hB : 0 ≤ B)
    (hA : ∀ σ : Fin N ≃ E, ‖Matrix.trace (chronologicalProduct A σ)‖ ≤ B) :
    ‖chamberTrace N A t‖ ≤ B := by
  by_cases h : ∃ σ : Fin N ≃ E, StrictMono (t ∘ σ)
  · obtain ⟨σ, hσ⟩ := h
    rw [chamberTrace_eq A t σ hσ]
    exact hA σ
  · rw [chamberTrace_eq_zero A t h, norm_zero]
    exact hB

open PoissonHistory

/-- The real trace kernel obtained by merging both histories chronologically. Cross-history
collisions are assigned zero. The two histories use the same matrix family. -/
def mergedTrace {b : ℝ≥0} {m d : ℕ} (A : Fin m → Matrix (Fin d) (Fin d) ℝ)
    (h k : History b m) : ℝ :=
  chamberTrace (h.1.1 + k.1.1) (Sum.elim (A ∘ h.1.2) (A ∘ k.1.2))
    (Sum.elim h.2.val k.2.val)

/-- Exchanging the two histories only relabels their union. It never reverses chronological
order and needs no transpose or commutativity assumption on the matrix family. -/
theorem mergedTrace_swap {b : ℝ≥0} {m d : ℕ}
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) (h k : History b m) :
    mergedTrace A h k = mergedTrace A k h := by
  unfold mergedTrace
  rw [Nat.add_comm h.1.1 k.1.1]
  have h := chamberTrace_reindex (N := k.1.1 + h.1.1)
    (Equiv.sumComm (Fin h.1.1) (Fin k.1.1))
    (Sum.elim (A ∘ k.1.2) (A ∘ h.1.2)) (Sum.elim k.2.val h.2.val)
  convert h using 1
  congr 1 <;> funext i <;> cases i <;> rfl

/-- Joint measurability on each pair of ordered-time components. -/
theorem measurable_mergedTrace_component {b : ℝ≥0} {m d : ℕ}
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) (w v : EventWord m) :
    Measurable (fun p : OrderedTimes b w.1 × OrderedTimes b v.1 =>
      mergedTrace A ⟨w, p.1⟩ ⟨v, p.2⟩) := by
  apply (measurable_chamberTrace (N := w.1 + v.1)
    (Sum.elim (A ∘ w.2) (A ∘ v.2))).comp
  apply measurable_pi_iff.mpr
  intro i
  cases i with
  | inl i => exact (measurable_pi_apply i).comp (measurable_subtype_coe.comp measurable_fst)
  | inr i => exact (measurable_pi_apply i).comp (measurable_subtype_coe.comp measurable_snd)

/-- The matrix-trace kernel is jointly measurable on the whole countable disjoint union. -/
theorem measurable_mergedTrace {b : ℝ≥0} {m d : ℕ}
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) :
    Measurable (Function.uncurry (mergedTrace (b := b) A)) :=
  SigmaProductMeasurable.measurable_of_sigma_prod (measurable_mergedTrace_component A)

/-- Contractive factors give a uniform dimension bound on the actual merged trace kernel. -/
theorem norm_mergedTrace_le {b : ℝ≥0} {m d : ℕ} [Nonempty (Fin d)]
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) (hA : ∀ a, ‖A a‖ ≤ 1)
    (h k : History b m) : ‖mergedTrace A h k‖ ≤ d := by
  apply norm_chamberTrace_le _ _ (by positivity)
  intro σ
  have hp := ContinuousKernel.norm_trace_prod_le A hA
    (List.ofFn (Sum.elim h.1.2 k.1.2 ∘ σ))
  have he : Sum.elim (A ∘ h.1.2) (A ∘ k.1.2) = A ∘ Sum.elim h.1.2 k.1.2 := by
    funext i
    cases i <;> rfl
  simpa only [chronologicalProduct, he, List.map_ofFn, Fintype.card_fin,
    Function.comp_assoc] using hp

/-- The complex-valued bounded kernel associated to contractive real matrix factors. Its
analytic properties are proved from the explicit chronological construction. -/
def boundedKernel {b : ℝ≥0} {m d : ℕ} [Nonempty (Fin d)]
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) (hA : ∀ a, ‖A a‖ ≤ 1) :
    ContinuousKernel.BoundedKernel (History b m) ℂ where
  toFun h k := (mergedTrace A h k : ℂ)
  stronglyMeasurable := (Complex.continuous_ofReal.measurable.comp
    (measurable_mergedTrace A)).stronglyMeasurable
  bound := d
  norm_le h k := by
    simpa only [Complex.norm_real, NNReal.coe_natCast] using norm_mergedTrace_le A hA h k

/-- Swap symmetry of the real chronological kernel gives Hermitian symmetry over `ℂ`. -/
theorem isHermitian_boundedKernel {b : ℝ≥0} {m d : ℕ} [Nonempty (Fin d)]
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) (hA : ∀ a, ‖A a‖ ≤ 1) :
    (boundedKernel (b := b) A hA).IsHermitian := by
  intro h k
  change starRingEnd ℂ (mergedTrace A h k : ℂ) = (mergedTrace A k h : ℂ)
  simp [mergedTrace_swap A h k]

/-- The empty history, valid at every nonnegative time and for every label count. -/
def emptyHistory (b : ℝ≥0) (m : ℕ) : History b m :=
  ⟨⟨0, Fin.elim0⟩, ⟨Fin.elim0, by
    constructor
    · intro i
      exact Fin.elim0 i
    · intro i
      exact Fin.elim0 i⟩⟩

/-- Both empty histories contribute the trace of the identity, including at `b = 0`. -/
theorem mergedTrace_empty {b : ℝ≥0} {m d : ℕ}
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) :
    mergedTrace A (emptyHistory b m) (emptyHistory b m) = d := by
  unfold mergedTrace
  rw [chamberTrace_eq _ _ finSumFinEquiv.symm (by intro i; exact Fin.elim0 i)]
  simp [chronologicalProduct, emptyHistory]

/-- A cross-history time collision gives the specified zero value. -/
theorem mergedTrace_eq_zero_of_collision {b : ℝ≥0} {m d : ℕ}
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) (h k : History b m)
    {i : Fin h.1.1} {j : Fin k.1.1} (hij : h.2.val i = k.2.val j) :
    mergedTrace A h k = 0 := by
  apply chamberTrace_eq_zero_of_collision
  intro hinj
  exact Sum.inl_ne_inr (hinj (show Sum.elim h.2.val k.2.val (.inl i) =
    Sum.elim h.2.val k.2.val (.inr j) from hij))

/-- For almost every pair, the chamber construction equals the matrix product of the sorted
union of its events. The initial enumeration concatenates the two component index types. -/
theorem ae_mergedTrace_eq_sorted {b : ℝ≥0} {m d : ℕ}
    (A : Fin m → Matrix (Fin d) (Fin d) ℝ) :
    ∀ᵐ p ∂(historyMeasure b m).prod (historyMeasure b m),
      mergedTrace A p.1 p.2 = Matrix.trace (chronologicalProduct
        (Sum.elim (A ∘ p.1.1.2) (A ∘ p.2.1.2))
        ((Tuple.sort (Sum.elim p.1.2.val p.2.2.val ∘ finSumFinEquiv.symm)).trans
          finSumFinEquiv.symm)) := by
  filter_upwards [ae_merged_history_injective b m] with p hp
  exact chamberTrace_eq_sorted _ _ _ hp

end CertifiedGaps.MergedHistory
