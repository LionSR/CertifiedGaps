/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.PurityCriterion
import Mathlib.Analysis.Matrix.PosDef

/-!
# Moment form from traces of powers of a positive semidefinite matrix

If the even partition functions are traces of powers of a positive semidefinite transfer matrix,
`Z (2 * N) = Re tr (T ^ N)` for `N ≥ 1`, then `Z` has the moment form of
`ThermalPurityGap.HasMomentForm` from `L = 2` on, with weights `μ i = √(λ i)` the square
roots of the eigenvalues of `T` (indexed through `Fintype.equivFin` and padded with zeros).

## Main results

* `CertifiedGaps.trace_pow_eq_sum_eigenvalues`: `tr (T ^ N) = ∑ i, λ i ^ N` for Hermitian `T`.
* `CertifiedGaps.hasMomentForm_of_trace_pow`: the moment form from a positive semidefinite `T`.
-/

open Matrix Unitary
open scoped ComplexOrder

namespace CertifiedGaps

/-- The trace of a power of a Hermitian matrix is the power sum of its eigenvalues. -/
theorem trace_pow_eq_sum_eigenvalues {D : Type*} [Fintype D] [DecidableEq D]
    {T : Matrix D D ℂ} (hT : T.IsHermitian) (N : ℕ) :
    Matrix.trace (T ^ N) = ((∑ i, hT.eigenvalues i ^ N : ℝ) : ℂ) := by
  conv_lhs => rw [hT.spectral_theorem, ← map_pow, diagonal_pow, conjStarAlgAut_apply]
  rw [trace_mul_cycle, coe_star_mul_self, one_mul, trace_diagonal]
  push_cast
  rfl

/-- If `Z (2 * N) = Re tr (T ^ N)` for all `N ≥ 1` with `T` positive semidefinite, then `Z` has
the moment form from `L = 2` on, with weights the square roots of the eigenvalues of `T`. -/
theorem hasMomentForm_of_trace_pow {D : Type*} [Fintype D] [DecidableEq D]
    {T : Matrix D D ℂ} (hT : T.PosSemidef) {Z : ℕ → ℝ}
    (hZ : ∀ N : ℕ, 1 ≤ N → Z (2 * N) = (Matrix.trace (T ^ N)).re) :
    ThermalPurityGap.HasMomentForm Z 2 := by
  set e := Fintype.equivFin D
  let μ : ℕ → ℝ := fun i ↦
    if h : i < Fintype.card D then √(hT.1.eigenvalues (e.symm ⟨i, h⟩)) else 0
  refine ⟨μ, fun i ↦ ?_, fun L hL hLe ↦ ?_⟩
  · simp only [μ]; split_ifs <;> positivity
  obtain ⟨N, rfl⟩ : ∃ N, L = 2 * N := hLe.exists_two_nsmul _ |>.imp fun _ h ↦ h
  have hN : 1 ≤ N := by omega
  have hsum := hasSum_sum_of_ne_finset_zero (f := fun i ↦ μ i ^ (2 * N))
    (L := SummationFilter.unconditional ℕ) (s := Finset.range (Fintype.card D)) fun i hi ↦ by
      simp only [Finset.mem_range, not_lt] at hi
      simp [μ, show ¬ i < Fintype.card D by omega, show 2 * N ≠ 0 by omega]
  convert hsum using 1
  rw [hZ N hN, trace_pow_eq_sum_eigenvalues hT.1, Complex.ofReal_re, Finset.sum_range,
    ← e.symm.sum_comp]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [μ, i.2, dite_true, pow_mul, Real.sq_sqrt (hT.eigenvalues_nonneg _)]

end CertifiedGaps
