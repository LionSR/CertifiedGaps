/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Gap
import CertifiedGaps.Transfer.SingleSiteLimit

/-!
# Partition functions and purities as sums over transfer eigenvalues

For a reflection-positive ring without an internal cell reflection, the single-site transfer
matrix `X_M = transfer C o M (β / M)` with `M` imaginary-time steps is real symmetric
(`SingleSite.transfer_isHermitian`); let `λ_{M,i}(β)` be its eigenvalues
(`RPRing.transferEigenvalues`). At every ring length `n ≥ 2`, odd lengths included, the
partition function and both purities are limits of power sums of these eigenvalues:

  `Z_n(β) = lim_M ∑ᵢ λ_{M,i}(β)ⁿ`,
  `P_s(n, β) = lim_M ∑ᵢ λ_{M,i}(β)^{2n} / (∑ᵢ λ_{M,i}(β)ⁿ)²`,
  `P_τ(L, β) = lim_M ∑ᵢ λ_{M,i}(2β)^L / (∑ᵢ λ_{M,i}(β)^L)²`.

This is the spectral form of the transfer operator `X_β` of the main note at finite Trotter
step: at each `M` the identity `Tr X_M^n = ∑ᵢ λ_{M,i}ⁿ` is exact
(`SingleSite.trace_stepOp_pow`, `trace_pow_eq_sum_eigenvalues_real`), and `M → ∞` recovers the
ring (`SingleSite.tendsto_trace_stepOp_pow`).

## Main results

* `RPRing.tendsto_sum_transferEigenvalues_pow`: the partition function.
* `RPRing.tendsto_spatialPurity`, `RPRing.tendsto_thermalPurity`: the two purities.
-/

open Filter Matrix Topology

namespace CertifiedGaps.RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] [DecidableEq ι] (R : RPRing d ι)

/-- The eigenvalues `λ_{M,i}(β)` of the single-site transfer matrix with `M` imaginary-time
steps of length `β / M`. -/
noncomputable def transferEigenvalues (M : ℕ) (β : ℝ) : (Fin M → Option ι) → ℝ :=
  (SingleSite.transfer_isHermitian R.C R.o M (β / M)).eigenvalues

/-- **The partition function as a sum over transfer eigenvalues.** At every length `n ≥ 2`,
`Z_n(β) = lim_M ∑ᵢ λ_{M,i}(β)ⁿ`. -/
theorem tendsto_sum_transferEigenvalues_pow (h : R.refl = 1) {β : ℝ} (hβ : 0 ≤ β) {n : ℕ}
    (hn : 2 ≤ n) :
    Tendsto (fun M ↦ ∑ i, R.transferEigenvalues M β i ^ n) atTop
      (𝓝 (R.partitionFunction n β)) := by
  refine Metric.tendsto_nhds.2 fun ε hε ↦ ?_
  have hev := SingleSite.eventually_abs_trace_transfer_sub_le R.C R.o (fun _ ↦ 1) hβ {n}
    (by simpa using hn) (half_pos hε)
  refine hev.mono fun M hM ↦ ?_
  have h1 := hM n (Finset.mem_singleton_self n)
  rw [show SingleSite.labelChar (fun _ ↦ (1 : ℝ)) = fun _ ↦ 1 from
      funext SingleSite.labelChar_one, diagonal_one, Matrix.one_mul,
    trace_pow_eq_sum_eigenvalues_real (SingleSite.transfer_isHermitian R.C R.o M (β / M))] at h1
  rw [Real.dist_eq]
  have h2 : SingleSite.partitionFunction R.C R.o n (fun _ ↦ 1) β = R.partitionFunction n β := by
    rw [partitionFunction, hamiltonian_eq_singleSite R h]; rfl
  rw [h2] at h1
  exact lt_of_le_of_lt h1 (half_lt_self hε)

/-- **The spatial purity as a limit of eigenvalue power sums.** At every length `n ≥ 2`,
`P_s(n, β) = lim_M ∑ᵢ λ_{M,i}^{2n} / (∑ᵢ λ_{M,i}ⁿ)²`. -/
theorem tendsto_spatialPurity [Nonempty d] (h : R.refl = 1) {β : ℝ} (hβ : 0 ≤ β) {n : ℕ}
    (hn : 2 ≤ n) :
    Tendsto (fun M ↦ (∑ i, R.transferEigenvalues M β i ^ (2 * n)) /
        (∑ i, R.transferEigenvalues M β i ^ n) ^ 2) atTop (𝓝 (R.spatialPurity n β)) :=
  (R.tendsto_sum_transferEigenvalues_pow h hβ (by omega)).div
    ((R.tendsto_sum_transferEigenvalues_pow h hβ hn).pow 2)
    (pow_ne_zero 2 (R.partitionFunction_pos n β).ne')

/-- **The thermal purity as a limit of eigenvalue power sums.** At every length `L ≥ 2`,
`P_τ(L, β) = lim_M ∑ᵢ λ_{M,i}(2β)^L / (∑ᵢ λ_{M,i}(β)^L)²`. -/
theorem tendsto_thermalPurity [Nonempty d] (h : R.refl = 1) {β : ℝ} (hβ : 0 ≤ β) {L : ℕ}
    (hL : 2 ≤ L) :
    Tendsto (fun M ↦ (∑ i, R.transferEigenvalues M (2 * β) i ^ L) /
        (∑ i, R.transferEigenvalues M β i ^ L) ^ 2) atTop (𝓝 (R.thermalPurity L β)) :=
  (R.tendsto_sum_transferEigenvalues_pow h (by positivity) hL).div
    ((R.tendsto_sum_transferEigenvalues_pow h hβ hL).pow 2)
    (pow_ne_zero 2 (R.partitionFunction_pos L β).ne')

end CertifiedGaps.RPRing
