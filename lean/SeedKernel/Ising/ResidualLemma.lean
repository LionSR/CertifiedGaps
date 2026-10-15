/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic.Linarith

/-!
# The weighted-residual lemma

The finite scalar step of the residual argument: a small weighted residual forces a nearby
eigenvalue.
-/

namespace SeedKernel.Ising

/-- The scalar spectral-basis step of the residual argument. To use it for a matrix one
must still prove the change of basis and the two norm identities. -/
theorem exists_near_of_weighted_residual {ι : Type*} [Fintype ι]
    (e w : ι → ℝ) (a r : ℝ) (hr : 0 ≤ r)
    (hw : ∀ i, 0 ≤ w i) (hpos : ∃ i, 0 < w i)
    (hres : (∑ i, (e i - a) ^ 2 * w i) ≤ r ^ 2 * ∑ i, w i) :
    ∃ i, |e i-a| ≤ r := by
  classical
  by_contra hn
  push Not at hn
  have hsquare (i : ι) : r^2 < (e i-a)^2 := by
    have hh := hn i
    have ha := abs_nonneg (e i-a)
    nlinarith [sq_abs (e i-a)]
  obtain ⟨i, hi⟩ := hpos
  have hsum := Finset.sum_lt_sum
    (s := Finset.univ)
    (f := fun j => r^2*w j) (g := fun j => (e j-a)^2*w j)
    (fun j _ => mul_le_mul_of_nonneg_right (hsquare j).le (hw j))
    ⟨i, Finset.mem_univ i, mul_lt_mul_of_pos_right (hsquare i) hi⟩
  rw [← Finset.mul_sum] at hsum
  exact (not_lt_of_ge hres) hsum

#print axioms exists_near_of_weighted_residual
end SeedKernel.Ising
