/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.SectorTransfer
import CertifiedGaps.Local.SpinOneRotation
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Transfer.SingleSiteA4

/-!
# `I`, `O`, `N` sectors of the spin-one ring

The seed of the length-60 initialization uses three twisted traces of the spin-one Heisenberg
ring: untwisted, twisted by the `π` rotation `P` about `z`, and twisted by the rotation `C` by
`2π/3` about `z`. Their sector form is the scheme `SectorScheme.ion`,

  `Z_n = I_n + 2 O_n + 3 N_n`,  `Z_n(P) = I_n + 2 O_n - N_n`,  `Z_n(C) = I_n - O_n`,

with `I`, `O`, `N` power sums of real numbers (`exists_ion_approx`).

**The `C`-twisted ring.** A ring is twisted by a real `3 × 3` matrix `R` through its closing
bond `L - 1 → 0`, `∑_{αβ} R_{αβ} S^α ⊗ S^β` (`rotTwistedHamiltonian`). The twist `C` is the
rotation matrix `cRot` by the angle `θ` with `cos θ = -1/2`, `sin θ = √3/2` about `z`. Since the
second site enters as `∑_β R_{αβ} S^β = g S^α g⁻¹` with `g = Rᵀ`, this is the bond
`∑_α S^α ⊗ g S^α g⁻¹` with the second site rotated by `g = Rᵀ`.

**From `C_z` to the axis cycle.** The generic `A₄` theorem (`Transfer/SingleSiteA4.lean`) twists
the single-site ring of the generators `G_α` by the cyclic relabeling `π : x → y → z → x`
(`cycleSymmetry`), whose closing bond is `-∑_α G_α ⊗ G_{π α} = ∑_α S^α ⊗ S^{π α}`, the ring
twisted by the permutation matrix `permRot π`. That matrix is the rotation by `2π/3` about the
diagonal `-(1, 1, 1)`, conjugate to `cRot` by the proper rotation `cBasis` with rows
`(1, -1, 0)/√2`, `(-1, -1, 2)/√6`, `-(1, 1, 1)/√3`:
`cBasisᵀ cRot cBasis = permRot π` (`cBasis_conj_cRot`). Rotating every site of the ring by
`cBasis` (the spin-one site space is the vector representation, `conj_rotExchange`) turns the
`C`-twisted ring into the `π`-twisted one and leaves the partition function unchanged
(`cTwistedPartitionFunction_eq`).

## Main results

* `rotTwistedHamiltonian`, `cRot`, `cTwistedHamiltonian`, `cTwistedPartitionFunction`.
* `cTwistedPartitionFunction_eq`: the `C`-twisted partition function is that of the single-site
  ring twisted by the axis cycle.
* `spinOneA4`: the `A₄` symmetry of the spin-one generators.
* `exists_ion_approx`: the sector representation of `Z_n`, `Z_n(P)` and `Z_n(C)`.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

open SingleSite SpinOne

/-! ### Rings twisted by a rotation -/

/-- The Heisenberg ring whose closing bond `L - 1 → 0` is twisted by a real matrix `R`:
`∑_{αβ} R_{αβ} S^α_{L-1} S^β_0`. -/
noncomputable def rotTwistedHamiltonian (R : Matrix (Fin 3) (Fin 3) ℝ) (L : ℕ) :
    Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  twistedRingHamiltonian exchange (rotExchange R) L

/-- The rotation by `θ` about `z` with `cos θ = -1/2`, `sin θ = √3/2`, that is `θ = 2π/3`. -/
noncomputable def cRot : Matrix (Fin 3) (Fin 3) ℝ :=
  !![-1 / 2, -(√3 / 2), 0; √3 / 2, -1 / 2, 0; 0, 0, 1]

/-- The Heisenberg ring twisted by the rotation `C` by `2π/3` about `z`. -/
noncomputable def cTwistedHamiltonian (L : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  rotTwistedHamiltonian cRot L

/-- The `C`-twisted partition function `Z_L(β, C) = Re Tr exp (-β H_L(C))`. -/
noncomputable def cTwistedPartitionFunction (L : ℕ) (β : ℝ) : ℝ :=
  gibbsTrace (cTwistedHamiltonian L) β

/-- The `C`-twisted ring is Hermitian. -/
theorem cTwistedHamiltonian_isHermitian (L : ℕ) : (cTwistedHamiltonian L).IsHermitian :=
  twistedRingHamiltonian_isHermitian exchange_isHermitian (rotExchange_isHermitian cRot) L

/-! ### The rotation taking `C_z` to the axis cycle -/

/-- The proper rotation with rows `(1, -1, 0)/√2`, `(-1, -1, 2)/√6`, `-(1, 1, 1)/√3`: it maps
the diagonal `-(1, 1, 1)/√3` to `e_z`. -/
noncomputable def cBasis : Matrix (Fin 3) (Fin 3) ℝ :=
  !![√2 / 2, -(√2 / 2), 0;
     -(√2 * √3 / 6), -(√2 * √3 / 6), √2 * √3 / 3;
     -(√3 / 3), -(√3 / 3), -(√3 / 3)]

/-- `√2² = 2` and `√3² = 3`, the only properties of the roots the checks below use. -/
theorem sqrt_two_sq_and_sqrt_three_sq : √2 ^ 2 = 2 ∧ √3 ^ 2 = 3 :=
  ⟨Real.sq_sqrt (by norm_num), Real.sq_sqrt (by norm_num)⟩

-- One simp set serves all entries; each entry uses only some of the root identities.
set_option linter.unusedSimpArgs false in
/-- `cBasis` is orthogonal. -/
theorem cBasis_transpose_mul_self : cBasisᵀ * cBasis = 1 := by
  obtain ⟨ha, hb⟩ := sqrt_two_sq_and_sqrt_three_sq
  generalize √2 = a at ha ⊢
  generalize √3 = b at hb ⊢
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [cBasis, mul_apply, Fin.sum_univ_three, one_apply] <;> ring_nf <;> simp [ha, hb] <;>
    ring_nf

set_option linter.unusedSimpArgs false in
/-- `cBasis` preserves the cross product: it is a proper rotation. -/
theorem cBasis_preservesCross : PreservesCross cBasis := by
  obtain ⟨ha, hb⟩ := sqrt_two_sq_and_sqrt_three_sq
  generalize √2 = a at ha ⊢
  generalize √3 = b at hb ⊢
  intro α c e
  fin_cases α <;> fin_cases c <;> fin_cases e <;>
    simp [cBasis, Fin.sum_univ_three, leviCivita] <;> ring_nf <;> simp [ha, hb] <;> ring_nf

set_option linter.unusedSimpArgs false in
/-- `cBasis` conjugates the rotation `C_z` to the axis cycle: `cBasisᵀ cRot cBasis = permRot π`
with `π = finRotate 3 : x → y → z → x`. -/
theorem cBasis_conj_cRot : cBasisᵀ * cRot * cBasis = permRot (finRotate 3) := by
  obtain ⟨ha, hb⟩ := sqrt_two_sq_and_sqrt_three_sq
  have hb3 : √3 ^ 3 = 3 * √3 := by rw [pow_succ, hb]
  generalize √2 = a at ha ⊢
  generalize √3 = b at hb hb3 ⊢
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [cBasis, cRot, permRot, mul_apply, Fin.sum_univ_three] <;> ring_nf <;>
    simp [ha, hb, hb3] <;> ring_nf

/-! ### The `C`-twisted ring is the axis-cycle-twisted single-site ring -/

/-- The closing bond of the single-site ring twisted by an axis permutation `π` is the exchange
twisted by `permRot π`. -/
theorem twistBond_G (π : Fin 3 → Fin 3) :
    twistBond G (fun _ ↦ 1) π = rotExchange (permRot π) := by
  rw [twistBond, rotExchange_permRot, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [bondC, bondC, map_G_kronecker_G, Complex.ofReal_one, one_smul, neg_neg]

/-- The single-site ring twisted by an axis permutation `π` is the Heisenberg ring twisted by
`permRot π`. -/
theorem twistHamiltonian_G (L : ℕ) (π : Fin 3 → Fin 3) :
    twistHamiltonian G 0 L (fun _ ↦ 1) π = rotTwistedHamiltonian (permRot π) L := by
  rw [twistHamiltonian, Matrix.map_zero _ Complex.ofReal_zero, siteSum_zero, add_zero,
    signedBond_G_one, twistBond_G, rotTwistedHamiltonian]

/-- **The `C`-twisted partition function** is the partition function of the single-site ring of
the generators twisted by the axis cycle `x → y → z → x`, on every ring of `L ≠ 1` sites. -/
theorem cTwistedPartitionFunction_eq {L : ℕ} (hL : L ≠ 1) (β : ℝ) :
    cTwistedPartitionFunction L β =
      twistPartitionFunction G 0 L (fun _ ↦ 1) (finRotate 3) β := by
  have hV := map_conjTranspose_mul_map cBasis_transpose_mul_self
  have hB : (cBasis.map (↑) ⊗ₖ cBasis.map (↑))ᴴ * exchange *
      (cBasis.map (↑) ⊗ₖ cBasis.map (↑)) = exchange := by
    rw [← rotExchange_one, conj_rotExchange cBasis_preservesCross, Matrix.mul_one,
      cBasis_transpose_mul_self]
  have hBc : (cBasis.map (↑) ⊗ₖ cBasis.map (↑))ᴴ * rotExchange cRot *
      (cBasis.map (↑) ⊗ₖ cBasis.map (↑)) = rotExchange (permRot (finRotate 3)) := by
    rw [conj_rotExchange cBasis_preservesCross, cBasis_conj_cRot]
  rw [twistPartitionFunction, twistHamiltonian_G, cTwistedPartitionFunction,
    cTwistedHamiltonian, rotTwistedHamiltonian, rotTwistedHamiltonian]
  exact (gibbsTrace_twistedRingHamiltonian_conj _ hV hL hB hBc β).symm

/-! ### The `A₄` symmetry and the sector form -/

/-- The `A₄` symmetry of the spin-one generators: the `D₂` symmetry `spinOneD2` (the `π`
rotations about the axes) and the cyclic relabeling of the axes, of order three. -/
noncomputable def spinOneA4 : A4Symmetry G (0 : Matrix (Fin 3) (Fin 3) ℝ) where
  toD2Symmetry := spinOneD2
  cycle_cube α := by fin_cases α <;> rfl

/-- **`I`, `O`, `N` sector form of the spin-one ring.** For every `β > 0` and every finite set
of ring lengths `n ≥ 2` (odd lengths included), the untwisted partition functions `Z_n(β)`, the
`P`-twisted `Z_n(β, P)` and the `C`-twisted `Z_n(β, C)` are, to any accuracy,
`I_n + 2 O_n + 3 N_n`, `I_n + 2 O_n - N_n` and `I_n - O_n` with `I`, `O`, `N` power sums of
finitely many real numbers. -/
theorem exists_ion_approx {β : ℝ} (hβ : 0 < β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    SectorApprox SectorScheme.ion
      ![fun n ↦ partitionFunction n β, fun n ↦ twistedPartitionFunction n β,
        fun n ↦ cTwistedPartitionFunction n β] ns := by
  refine (spinOneA4.sectorApprox_ion hβ.le ns hns).congr fun g n hn ↦ ?_
  fin_cases g
  · change SingleSite.partitionFunction G 0 n (fun _ ↦ 1) β = partitionFunction n β
    rw [SingleSite.partitionFunction, hamiltonian_plain]
    rfl
  · change SingleSite.partitionFunction G 0 n spinOneD2.sz β = twistedPartitionFunction n β
    rw [SingleSite.partitionFunction, spinOneD2_sz, hamiltonian_twist]
    rfl
  · change twistPartitionFunction G 0 n (fun _ ↦ 1) (finRotate 3) β =
      cTwistedPartitionFunction n β
    rw [cTwistedPartitionFunction_eq (by have := hns n hn; omega)]

/-! ### The twist by `C⁻¹` -/

/-- The `π` rotation about `x`, which conjugates `C` to `C⁻¹`. -/
def flipX : Matrix (Fin 3) (Fin 3) ℝ := !![1, 0, 0; 0, -1, 0; 0, 0, -1]

theorem flipX_transpose_mul_self : flipXᵀ * flipX = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [flipX, mul_apply, Fin.sum_univ_three]

theorem flipX_preservesCross : PreservesCross flipX := by
  intro α c e
  fin_cases α <;> fin_cases c <;> fin_cases e <;> simp [flipX, Fin.sum_univ_three, leviCivita]

theorem flipX_conj_cRot : flipXᵀ * cRot * flipX = cRotᵀ := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [flipX, cRot, mul_apply, Fin.sum_univ_three]

/-- **The twist by `C⁻¹` gives the same partition function**: the trace table may use either
orientation of the `2π/3` rotation. -/
theorem gibbsTrace_rotTwisted_cRot_transpose {L : ℕ} (hL : L ≠ 1) (β : ℝ) :
    gibbsTrace (rotTwistedHamiltonian cRotᵀ L) β = cTwistedPartitionFunction L β := by
  have hV := map_conjTranspose_mul_map flipX_transpose_mul_self
  have hB : (flipX.map (↑) ⊗ₖ flipX.map (↑))ᴴ * exchange *
      (flipX.map (↑) ⊗ₖ flipX.map (↑)) = exchange := by
    rw [← rotExchange_one, conj_rotExchange flipX_preservesCross, Matrix.mul_one,
      flipX_transpose_mul_self]
  have hBc : (flipX.map (↑) ⊗ₖ flipX.map (↑))ᴴ * rotExchange cRot *
      (flipX.map (↑) ⊗ₖ flipX.map (↑)) = rotExchange cRotᵀ := by
    rw [conj_rotExchange flipX_preservesCross, flipX_conj_cRot]
  rw [cTwistedPartitionFunction, cTwistedHamiltonian, rotTwistedHamiltonian,
    rotTwistedHamiltonian]
  exact gibbsTrace_twistedRingHamiltonian_conj _ hV hL hB hBc β

end CertifiedGaps.SpinOneHeisenberg
