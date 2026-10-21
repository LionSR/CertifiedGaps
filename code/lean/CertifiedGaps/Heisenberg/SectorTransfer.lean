/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.SzBridge
import CertifiedGaps.Transfer.SingleSiteD2
import CertifiedGaps.Transfer.SpectralForm
import CertifiedGaps.Heisenberg.Gap
import CertifiedGaps.Models.SpinOneLargeD

/-!
# Sector moment form of the spin-one ring partition functions

The appendix "From partition functions to the transfer spectrum" of the certificate note uses
real numbers `λᵢ` split into the invariant `D₂` sector `J` and one nontrivial sector `N`, with

  `Z_n(β) = ∑_{i ∈ J} λᵢ^n + 3 ∑_{i ∈ N} λᵢ^n`,  `Z_n(β, P) = ∑_{i ∈ J} λᵢ^n - ∑_{i ∈ N} λᵢ^n`

for every ring length `n ≥ 2`, odd lengths included. This file proves the form that the
finite-dimensional argument needs: for any finite set of lengths and any tolerance `δ > 0` there
are such numbers reproducing both partition functions within `δ`.

The argument is the generic one of `Transfer/SingleSiteD2.lean`; this file only supplies the
spin-one data. The bond factors are the real antisymmetric generators `SpinOne.G α = i S^α`,
`(G_α)_{bc} = ε_{αbc}` (`S · S = -∑_α G_α ⊗ G_α`), with no on-site term, and the `D₂` symmetry
(`spinOneD2`) is:

* the `π` rotations about `x` and `y`, which in the Cartesian basis are the diagonal sign
  matrices `SpinOne.rotX`, `SpinOne.rotY` with `R G_α R = r_α G_α` (since `r_b r_c = r_α` whenever
  `ε_{αbc} ≠ 0`); their product is the twist `rotZ`, the `π` rotation `P` about `z`;
* the cyclic relabeling `x → y → z → x`, which preserves `ε`.

The generic ring Hamiltonians with closing-bond signs `1` and `rotZ` are `H_L` and `H_L(P)`
(`hamiltonian_plain`, `hamiltonian_twist`).

## Main results

* `hamiltonian_plain`, `hamiltonian_twist`: the generic ring Hamiltonians are `H_L` and `H_L(P)`.
* `spinOneD2`: the `D₂` symmetry with cyclic relabeling of the spin-one generators.
* `exists_sector_approx`: the approximate moment form for finitely many lengths.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

open SingleSite SpinOne

/-! ### The spin-one instance -/

/-- The plain signed bond of the generators is the exchange `S · S`, since
`G_α ⊗ G_α = -(S^α ⊗ S^α)`. -/
theorem signedBond_G_one : signedBond G (fun _ ↦ 1) = exchange := by
  simp [signedBond, bondC, map_G_kronecker, exchange]

/-- The signed bond of the generators with the signs `rotZ` is the twisted exchange. -/
theorem signedBond_G_rotZ : signedBond G rotZ = twistedExchange := by
  rw [signedBond, twistedExchange, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [bondC, map_G_kronecker, smul_neg, neg_neg, twistSign_eq]

/-- The generic plain ring is the Heisenberg ring. -/
theorem hamiltonian_plain (L : ℕ) :
    SingleSite.hamiltonian G 0 L (fun _ ↦ 1) = hamiltonian L := by
  rw [SingleSite.hamiltonian, Matrix.map_zero _ Complex.ofReal_zero, siteSum_zero, add_zero,
    twistedRingHamiltonian_self, signedBond_G_one]
  rfl

/-- The generic ring with the closing-bond signs `rotZ` is the twisted Heisenberg ring. -/
theorem hamiltonian_twist (L : ℕ) : SingleSite.hamiltonian G 0 L rotZ = twistedHamiltonian L := by
  rw [SingleSite.hamiltonian, Matrix.map_zero _ Complex.ofReal_zero, siteSum_zero, add_zero,
    signedBond_G_one, signedBond_G_rotZ, twistedHamiltonian_eq]

/-! ### The `D₂` symmetry -/

/-- A `π` rotation about a Cartesian axis, given by its sign vector `r` (with `r b r c = r α`
whenever `ε_{αbc} ≠ 0`), is a sign symmetry of the generators with the same signs. -/
noncomputable def axisSymmetry (r : Fin 3 → ℝ) (hr : ∀ α, r α * r α = 1)
    (hε : ∀ α b c, leviCivita α b c ≠ 0 → r b * r c = r α) :
    SignSymmetry G (0 : Matrix (Fin 3) (Fin 3) ℝ) where
  sign := r
  sign_sq := hr
  map := conjAlgHom (signUnitary r hr)
  trace_map := trace_conjAlgHom _
  map_C α := by
    rw [conjAlgHom_signUnitary_apply]
    ext b c
    simp only [G, of_apply, Matrix.smul_apply, smul_eq_mul]
    by_cases h : leviCivita α b c = 0
    · simp [h]
    · rw [← hε α b c h]
  map_o := map_zero _

/-- The cyclic relabeling `x → y → z → x` of the axes preserves the Levi-Civita symbol. -/
noncomputable def cycleSymmetry : LabelSymmetry G (0 : Matrix (Fin 3) (Fin 3) ℝ) where
  perm := finRotate 3
  map := (reindexAlgEquiv ℝ ℝ (finRotate 3)).toAlgHom
  trace_map A := trace_reindex_self _ A
  map_C α := by
    change reindexAlgEquiv ℝ ℝ (finRotate 3) (G α) = G (finRotate 3 α)
    ext b c
    simp only [coe_reindexAlgEquiv, G, reindex_apply, submatrix_apply, of_apply]
    fin_cases α <;> fin_cases b <;> fin_cases c <;> rfl
  map_o := map_zero _

/-- The `D₂` symmetry of the spin-one generators: the `π` rotations about `x` and `y` and the
cyclic relabeling of the axes. -/
noncomputable def spinOneD2 : D2Symmetry G (0 : Matrix (Fin 3) (Fin 3) ℝ) where
  x := axisSymmetry rotX (by intro α; fin_cases α <;> norm_num [rotX])
    (by intro α b c; fin_cases α <;> fin_cases b <;> fin_cases c <;> simp [rotX, leviCivita])
  y := axisSymmetry rotY (by intro α; fin_cases α <;> norm_num [rotY])
    (by intro α b c; fin_cases α <;> fin_cases b <;> fin_cases c <;> simp [rotY, leviCivita])
  cycle := cycleSymmetry
  x_cycle α := by fin_cases α <;> norm_num [axisSymmetry, cycleSymmetry, rotX, rotY]
  y_cycle α := by fin_cases α <;> norm_num [axisSymmetry, cycleSymmetry, rotX, rotY]

/-- The third rotation of `spinOneD2` is the twist. -/
theorem spinOneD2_sz : spinOneD2.sz = rotZ := by
  funext α
  fin_cases α <;> norm_num [D2Symmetry.sz, spinOneD2, axisSymmetry, rotX, rotY, rotZ]

/-! ### The sector moment form -/

/-- **Approximate sector moment form.** For every `β > 0`, every finite set of ring lengths
`n ≥ 2` and every `δ > 0` there are real numbers `x i` and two finite index sets `SJ`, `SN` with
`∑_{SJ} xᵢ^n + 3 ∑_{SN} xᵢ^n` within `δ` of `Z_n(β)` and `∑_{SJ} xᵢ^n - ∑_{SN} xᵢ^n` within `δ`
of `Z_n(β, P)`, for every `n` in the set. -/
theorem exists_sector_approx {β : ℝ} (hβ : 0 < β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n)
    {δ : ℝ} (hδ : 0 < δ) :
    ∃ (ι : Type) (_ : Fintype ι) (SJ SN : Finset ι) (x : ι → ℝ), ∀ n ∈ ns,
      |∑ i ∈ SJ, x i ^ n + 3 * ∑ i ∈ SN, x i ^ n - partitionFunction n β| ≤ δ ∧
      |∑ i ∈ SJ, x i ^ n - ∑ i ∈ SN, x i ^ n - twistedPartitionFunction n β| ≤ δ := by
  obtain ⟨ι, _, SJ, SN, x, hx⟩ := spinOneD2.d2Approx hβ.le ns hns δ hδ
  refine ⟨ι, inferInstance, SJ, SN, x, fun n hn ↦ ?_⟩
  have h := hx n hn
  simp only [spinOneD2_sz, SingleSite.partitionFunction, hamiltonian_plain, hamiltonian_twist]
    at h
  exact h

/-! ### The spectral form of the Heisenberg ring -/

section SpectralForm

open Filter Topology

/-- The Heisenberg partition function as a sum over transfer eigenvalues: at every length
`n ≥ 2`, `Z_n(β) = lim_M ∑ᵢ λ_{M,i}(β)ⁿ`, with `λ_{M,i}` the eigenvalues of the single-site
transfer matrix of the ring `SpinOneLargeD.model 0`. -/
theorem tendsto_sum_transferEigenvalues_pow {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (hn : 2 ≤ n) :
    Tendsto (fun M ↦ ∑ i, (SpinOneLargeD.model 0).transferEigenvalues M β i ^ n) atTop
      (𝓝 (partitionFunction n β)) := by
  simpa [SpinOneLargeD.partitionFunction_zero] using
    (SpinOneLargeD.model 0).tendsto_sum_transferEigenvalues_pow rfl hβ hn

/-- The Heisenberg spatial purity as a limit of transfer-eigenvalue power sums. -/
theorem tendsto_spatialPurity {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (hn : 2 ≤ n) :
    Tendsto (fun M ↦ (∑ i, (SpinOneLargeD.model 0).transferEigenvalues M β i ^ (2 * n)) /
        (∑ i, (SpinOneLargeD.model 0).transferEigenvalues M β i ^ n) ^ 2) atTop
      (𝓝 (spatialPurity n β)) := by
  simpa [RPRing.spatialPurity, spatialPurity, SpinOneLargeD.partitionFunction_zero] using
    (SpinOneLargeD.model 0).tendsto_spatialPurity rfl hβ hn

/-- The Heisenberg thermal purity as a limit of transfer-eigenvalue power sums. -/
theorem tendsto_thermalPurity {β : ℝ} (hβ : 0 ≤ β) {L : ℕ} (hL : 2 ≤ L) :
    Tendsto (fun M ↦ (∑ i, (SpinOneLargeD.model 0).transferEigenvalues M (2 * β) i ^ L) /
        (∑ i, (SpinOneLargeD.model 0).transferEigenvalues M β i ^ L) ^ 2) atTop
      (𝓝 (thermalPurity L β)) := by
  simpa [RPRing.thermalPurity, thermalPurity, SpinOneLargeD.partitionFunction_zero] using
    (SpinOneLargeD.model 0).tendsto_thermalPurity rfl hβ hL

end SpectralForm


end CertifiedGaps.SpinOneHeisenberg
