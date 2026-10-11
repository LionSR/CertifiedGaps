/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import ThermalMPOGap.PurityCriterion
import ThermalMPOGap.Transfer.ProductLimit
import ThermalMPOGap.Transfer.MomentForm
import ThermalMPOGap.Transfer.Factorization

/-!
# The spatial moment inequalities for the Heisenberg ring

The transfer-operator input of the thermal-MPO route, proved with finite matrices only. For
fixed `β > 0` and `M` Trotter slices, `trotterTrace M β (2N) = Re Tr (Fᴴ F)^N` for a finite
matrix `F` (`exists_posSemidef_transfer`), so the approximants have the moment form
(`hasMomentForm_of_trace_pow`). As `M → ∞` they converge to `Re Tr e^{-βH_L}`
(`tendsto_list_prod_one_add_smul_pow`). The moment inequalities are closed, so they hold for
the limit (`SpatialMomentBounds.of_tendsto_hasMomentForm`).

## Main result

* `spatialMomentBounds_partitionFunction`: for every `β > 0`, the traces `Re Tr e^{-βH_L}`
  satisfy the moment inequalities above length `2`.
-/

open Filter Matrix ThermalPurityGap ThermalMPOGap
open scoped ComplexOrder

namespace SpinOneHeisenberg

/-- The Trotter approximants have the moment form, for every number of slices. -/
theorem hasMomentForm_trotterTrace (M : ℕ) {β : ℝ} (hβ : 0 < β) :
    HasMomentForm (trotterTrace M β) 2 := by
  obtain ⟨D, _, _, T, hT, htr⟩ := exists_posSemidef_transfer M (β / M) (by positivity)
  exact hasMomentForm_of_trace_pow hT fun N hN ↦ by rw [trotterTrace, htr N hN]

/-- The sum of the bond terms in Trotter order is the Hamiltonian. -/
theorem sum_bondOrder (L : ℕ) (f : Fin L → Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ) :
    ((bondOrder L).map f).sum = ∑ j, f j := by
  have hperm : (bondOrder L).Perm (List.finRange L) := by
    unfold bondOrder
    simpa using List.filter_append_perm (fun j : Fin L ↦ decide (Even j.val)) (List.finRange L)
  rw [(hperm.map f).sum_eq, Fin.sum_univ_def]

/-- The Trotter approximants converge to the partition function. -/
theorem tendsto_trotterTrace (β : ℝ) (L : ℕ) :
    Tendsto (fun M : ℕ ↦ trotterTrace M β L) atTop (nhds (partitionFunction L β)) := by
  have hslice : ∀ τ : ℝ, trotterSlice L τ =
      (((bondOrder L).map fun j ↦ -bond L j).map fun h ↦ 1 + (τ : ℂ) • h).prod := by
    intro τ
    simp only [trotterSlice, List.map_map, Function.comp_def, smul_neg, ← sub_eq_add_neg]
  have hsum : ((bondOrder L).map fun j ↦ -bond L j).sum = -hamiltonian L := by
    rw [sum_bondOrder, Finset.sum_neg_distrib, hamiltonian_eq_sum_bond]
  have h := tendsto_list_prod_one_add_smul_pow ((bondOrder L).map fun j ↦ -bond L j) β
  rw [hsum, smul_neg, ← neg_smul] at h
  have hc : Continuous fun A : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ ↦ (trace A).re :=
    Complex.continuous_re.comp (continuous_id.matrix_trace)
  have := (hc.tendsto _).comp h
  simpa [Function.comp_def, trotterTrace, hslice, partitionFunction] using this

/-- **The transfer-operator input.** For every `β > 0`, the traces `Re Tr e^{-βH_L}` satisfy the
spatial moment inequalities above length `2`. -/
theorem spatialMomentBounds_partitionFunction (β : ℝ) (hβ : 0 < β) :
    SpatialMomentBounds (fun L ↦ partitionFunction L β) 2 :=
  SpatialMomentBounds.of_tendsto_hasMomentForm (Zs := fun M ↦ trotterTrace M β) two_pos
    (Eventually.of_forall fun M ↦ hasMomentForm_trotterTrace M hβ)
    fun L _ _ ↦ tendsto_trotterTrace β L

end SpinOneHeisenberg
