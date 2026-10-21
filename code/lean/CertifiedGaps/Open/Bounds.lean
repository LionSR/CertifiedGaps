/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.MomentBounds
import CertifiedGaps.Open.Transfer
import CertifiedGaps.Ring.Gap

/-!
# Boundary moment inequalities for open chains of reflection-positive rings

For a reflection-positive ring without internal reflection and any boundary weights `b`, the
ring partition functions `Z_L(β)` and the open-chain partition functions `G_n(β)` satisfy the
boundary moment inequalities (`RPRing.openMomentBounds`). For every number of slices `M`, the
Trotter approximants of both are the bulk and boundary moments `Tr X^{2N}` and `⟨v, Xⁿ v⟩` of one
Hermitian matrix `X = F J` (`RPRing.trace_trotterSlice_pow_eq_trace_openTransfer`,
`RPRing.trace_openSlice_pow`), so they satisfy the inequalities
(`OpenMomentBounds.of_hermitian`); the inequalities are closed, so they pass to the limit
`M → ∞` (`OpenMomentBounds.of_tendsto`).

## Main results

* `RPRing.openMomentBounds_trotterTrace`: the inequalities for the Trotter approximants.
* `RPRing.openMomentBounds`: the inequalities for the exact partition functions.
-/

open Filter Matrix

namespace CertifiedGaps

namespace RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- The Trotter approximant `Re Tr (openSlice)^M` of the partition function of the open chain
of `n ≥ 1` sites, with `M` slices of width `β / M`. -/
noncomputable def openTrotterTrace (b : R.Boundary) (M : ℕ) (β : ℝ) (n : ℕ) : ℝ :=
  (trace (R.openSlice b (n - 1) (β / M) ^ M)).re

/-- For every number of slices, the ring and open-chain Trotter approximants satisfy the boundary
moment inequalities: they are the moments of the Hermitian matrix `X = F J`. -/
theorem openMomentBounds_trotterTrace (hrefl : R.refl = 1) (b : R.Boundary) (M : ℕ) {β : ℝ}
    (hβ : 0 < β) : OpenMomentBounds (R.trotterTrace M β) (R.openTrotterTrace b M β) := by
  classical
  have hτ : 0 ≤ β / M := by positivity
  refine OpenMomentBounds.of_hermitian (R.openTransfer_isHermitian hrefl b M (β / M))
    (capVector M (β / M) b.w) (fun N hN ↦ ?_) (fun n hn ↦ ?_)
  · rw [trotterTrace, R.trace_trotterSlice_pow_eq_trace_openTransfer hrefl b hN hτ]
  · obtain ⟨n, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
    rw [openTrotterTrace, Nat.add_sub_cancel, R.trace_openSlice_pow hrefl b n hτ]

/-- **The boundary moment inequalities for an open chain.** For a reflection-positive ring
without internal reflection, any boundary weights `b` and every `β > 0`, the ring partition
functions `Z_L(β)` and the open-chain partition functions `G_n(β)` satisfy the boundary moment
inequalities. -/
theorem openMomentBounds [Nonempty d] (hrefl : R.refl = 1) (b : R.Boundary) {β : ℝ}
    (hβ : 0 < β) :
    OpenMomentBounds (fun L ↦ R.partitionFunction L β) (fun n ↦ R.openPartitionFunction b n β) := by
  have hZ : ∀ L, 0 < R.partitionFunction L β := fun L ↦ by
    rw [partitionFunction_eq_sum]
    exact Finset.sum_pos (fun _ _ ↦ Real.exp_pos _) Finset.univ_nonempty
  refine OpenMomentBounds.of_tendsto (Zs := fun M ↦ R.trotterTrace M β)
    (Gs := fun M ↦ R.openTrotterTrace b M β)
    (Eventually.of_forall fun M ↦ R.openMomentBounds_trotterTrace hrefl b M hβ)
    (fun n _ ↦ R.tendsto_trotterTrace β n) (fun n hn ↦ ?_) (hZ 2) (hZ 4)
  obtain ⟨n, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  simpa [openTrotterTrace] using R.tendsto_openSliceTrace b β n

end RPRing

end CertifiedGaps
