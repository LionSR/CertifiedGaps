/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Symmetry

/-!
# Phase gauges and bond-changing symmetries of twisted rings

`Ring/Symmetry.lean` treats sign gauges (`p² = 1`) and symmetries that preserve both bond
operators. The ring twisted by the rotation `C` by `2π/3` (`thermal120.tex`) needs two
generalizations:

* **Phase gauges.** The seam of the `C`-twisted ring is moved by the phase `ω^{-m}` on one site,
  which is not an involution. `Intertwines.symm_of_id_of_mul` inverts a gauge with any inverse
  phases, and `intertwines_twisted_rotate_of_inv` is the signed rotation of
  `intertwines_twisted_rotate` for a gauge `p` with inverse `q` (the paper's monomial unitary
  `|x⟩ ↦ e^{-iθ x_{n-1}} |s(x)⟩`).
* **Symmetries that change the bonds.** Reversal and negation of the words take the `C`-twisted
  ring to its complex conjugate. `intertwines_twisted_reverse_swap` and
  `intertwines_twisted_onsite_map` record how the reflection and an on-site relabeling act on
  arbitrary bond operators; the antiunitary symmetry is then the case where the transformed bonds
  are the conjugate ones.
-/

open Matrix

namespace CertifiedGaps

variable {d : Type} [DecidableEq d] {L : ℕ} {R : Type*}

/-- A gauge (`e = id`) with phases `u` that have inverses `v` intertwines in the opposite
direction with `v`. -/
theorem Intertwines.symm_of_id_of_mul {ι : Type*} [CommMonoid R] {u v : ι → R}
    (huv : ∀ σ, u σ * v σ = 1) {H H' : Matrix ι ι R} (h : Intertwines id u H H') :
    Intertwines id v H' H := fun σ τ ↦ by
  have h' := h σ τ
  simp only [id] at h' ⊢
  calc H σ τ * v τ = v σ * (u σ * H σ τ) * v τ := by
        rw [← mul_assoc, mul_comm (v σ), huv, one_mul]
    _ = v σ * (H' σ τ * u τ) * v τ := by rw [h']
    _ = v σ * H' σ τ * (u τ * v τ) := by simp only [mul_assoc]
    _ = v σ * H' σ τ := by rw [huv, mul_one]

/-- **The phase-gauged rotation of a twisted ring.** If the phase `p (σ s)` on one site moves
the seam across that site (`h₁`, `h₂`, as in `intertwines_seam_gauge`) and `q` is the inverse of
`p`, then rotating the ring by one site and applying `q (σ 0)` maps it to itself. For `p² = 1`
this is `intertwines_twisted_rotate`. -/
theorem intertwines_twisted_rotate_of_inv [CommSemiring R] (hL : 2 ≤ L) (p q : d → R)
    (hpq : ∀ a, p a * q a = 1) {B Bc : Matrix (d × d) (d × d) R}
    (h₁ : ∀ a b c e, Bc (a, b) (c, e) * p e = p b * B (a, b) (c, e))
    (h₂ : ∀ a b c e, B (a, b) (c, e) * p c = p a * Bc (a, b) (c, e)) :
    Intertwines (relabel (finRotate L))
      (fun σ ↦ 1 * q (σ ⟨0, by omega⟩)) (twistedRingHamiltonian B Bc L)
      (twistedRingHamiltonian B Bc L) := by
  set last : Fin L := ⟨L - 1, by omega⟩
  set zero : Fin L := ⟨0, by omega⟩
  have hlast : finRotate L last = zero := by
    apply Fin.ext; rw [val_finRotate]; simp [last, zero]; omega
  have hzero : (finRotate L).symm zero = last := by rw [Equiv.symm_apply_eq, hlast]
  rw [twistedRingHamiltonian_eq_seam B Bc last (by simp [last]; omega)]
  -- Gauge the seam from `L - 1 → 0` to `0 → 1` with the inverse phases, then rotate it back.
  have hg := (intertwines_seam_gauge hL p h₁ h₂ zero).symm_of_id_of_mul
    (v := fun σ : Fin L → d ↦ q (σ zero)) (fun σ ↦ hpq _)
  rw [hzero] at hg
  have hr := intertwines_seam_rotate B Bc last
  rw [hlast] at hr
  exact hg.trans hr

/-- **Reflection of a ring with a seam, for arbitrary bonds.** Read through the reflection
`σ ↦ σ ∘ rev`, the ring with bonds `B`, `Bc` and seam `m` is the ring with the swapped bonds
`B₁ = B ∘ swap`, `Bc₁ = Bc ∘ swap` and seam `rev (m + 1)`. -/
theorem intertwines_seam_reverse_swap [CommSemiring R] {B Bc B₁ Bc₁ : Matrix (d × d) (d × d) R}
    (hB : B.submatrix Prod.swap Prod.swap = B₁) (hBc : Bc.submatrix Prod.swap Prod.swap = Bc₁)
    (m : Fin L) :
    Intertwines (relabel Fin.revPerm) (fun _ ↦ 1)
      (seamHamiltonian B₁ Bc₁ L (Fin.rev (finRotate L m))) (seamHamiltonian B Bc L m) :=
  fun σ τ ↦ by
  simp only [seamHamiltonian, Matrix.sum_apply, bondOp_relabel, mul_one, one_mul,
    Fin.revPerm_apply]
  -- The bond `j` becomes the bond `k = rev (j + 1)`, traversed backwards.
  let e : Equiv.Perm (Fin L) := (finRotate L).trans Fin.revPerm
  conv_rhs => rw [← Equiv.sum_comp e]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  have he : e j = Fin.rev (finRotate L j) := rfl
  rw [he, finRotate_rev_finRotate, bondOp_swap (Fin.rev j)]
  have hj : (Fin.rev (finRotate L j) = Fin.rev (finRotate L m)) ↔ j = m := by
    rw [Fin.rev_inj, (finRotate L).injective.eq_iff]
  rw [if_congr hj rfl rfl]
  split_ifs <;> simp only [hB, hBc]

/-- **Reflection of a twisted ring, for arbitrary bonds.** The reflection `j ↦ L - 1 - j` maps
the closing bond to itself and swaps the two sites of every bond: it intertwines the twisted ring
of `B ∘ swap`, `Bc ∘ swap` with the twisted ring of `B`, `Bc`. -/
theorem intertwines_twisted_reverse_swap [CommSemiring R] (hL : 0 < L)
    {B Bc B₁ Bc₁ : Matrix (d × d) (d × d) R} (hB : B.submatrix Prod.swap Prod.swap = B₁)
    (hBc : Bc.submatrix Prod.swap Prod.swap = Bc₁) :
    Intertwines (relabel Fin.revPerm) (fun _ ↦ 1) (twistedRingHamiltonian B₁ Bc₁ L)
      (twistedRingHamiltonian B Bc L) := by
  set last : Fin L := ⟨L - 1, by omega⟩
  have h0 : finRotate L last = ⟨0, hL⟩ := by
    apply Fin.ext; rw [val_finRotate]; simp [last]; omega
  have hlast : Fin.rev (finRotate L last) = last := by
    rw [h0]; apply Fin.ext; simp [last]
  have h := intertwines_seam_reverse_swap hB hBc last
  rwa [hlast, ← twistedRingHamiltonian_eq_seam B₁ Bc₁ last (by simp [last]; omega),
    ← twistedRingHamiltonian_eq_seam B Bc last (by simp [last]; omega)] at h

/-- **An on-site relabeling of a twisted ring, for arbitrary bonds.** Read through
`σ ↦ f ∘ σ`, the twisted ring of `B`, `Bc` is the twisted ring of the relabeled bonds. -/
theorem intertwines_twisted_onsite_map [CommSemiring R] (f : Equiv.Perm d)
    {B Bc B₁ Bc₁ : Matrix (d × d) (d × d) R}
    (hB : B.submatrix (Prod.map f f) (Prod.map f f) = B₁)
    (hBc : Bc.submatrix (Prod.map f f) (Prod.map f f) = Bc₁) :
    Intertwines (onsite f) (fun _ ↦ 1) (twistedRingHamiltonian B₁ Bc₁ L)
      (twistedRingHamiltonian B Bc L) := fun σ τ ↦ by
  simp only [twistedRingHamiltonian, Matrix.sum_apply, bondOp_onsite, mul_one, one_mul]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  split_ifs <;> simp only [hB, hBc]

end CertifiedGaps
