/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Hamiltonian
import CertifiedGaps.Spectral.OperatorBounds

/-!
# Uniform energy bounds for ring Hamiltonians

Rigorous trace enclosures need a priori bounds `a L ≤ H_L ≤ b L` on the spectrum of a ring
Hamiltonian, uniform in the length `L`. This file lifts operator inequalities on one bond
(`Spectral/OperatorBounds.lean`; `a ≤ A` means `∀ x, a * nsq x ≤ qform A x`) to the ring.

## Main results

* Lifting a bond bound (`le_qform_bondOp`, `qform_bondOp_le`): since `bondOp j k` is a
  `*`-algebra embedding for `j ≠ k`, a bound `a ≤ B ≤ b` on a two-site operator holds for
  `bondOp j k B`: `b - B = Cᴴ C` becomes `(bondOp C)ᴴ (bondOp C)`.
* Summing over the `L` bonds gives `a L ≤ H_L ≤ b L` for every twisted or untwisted ring
  (`le_qform_twistedRingHamiltonian`, `qform_twistedRingHamiltonian_le` and the
  `ringHamiltonian` forms).
* Row sums of `bondOp` (`sum_norm_bondOp_mul`): the entries of `bondOp j k B` in the row `σ` are
  those of `B` in the row `(σ j, σ k)`, every other site being frozen. A Schur test
  (`Spectral.abs_qform_le_of_schur`) for an operator supported on a few sites therefore reduces
  to a finite check on those sites, with weights that depend only on their letters.
* Pairs of consecutive bonds (`le_qform_twistedRingHamiltonian_of_pairs`): if every pair
  `B_{j,j+1} + B_{j+1,j+2}` is bounded below by `a`, then `H_L ≥ (a/2) L`, since summing the `L`
  pairs counts every bond twice (`two_mul_qform_twistedRingHamiltonian`).
-/

open Matrix

namespace CertifiedGaps

open Spectral

variable {d : Type} [Fintype d] [DecidableEq d] {L : ℕ}

/-! ### Row sums of an embedded bond -/

/-- **Row sums of an embedded bond.** The weighted absolute row sum of `bondOp j k B` in the row
`σ` runs over the letters `p` of the two bond sites only, every other site being frozen. -/
theorem sum_norm_bondOp_mul {j k : Fin L} (hjk : j ≠ k) (B : Matrix (d × d) (d × d) ℂ)
    (σ : Fin L → d) (q : (Fin L → d) → ℝ) :
    ∑ τ, ‖bondOp j k B σ τ‖ * q τ =
      ∑ p : d × d, ‖B (σ j, σ k) p‖ * q (Function.update (Function.update σ j p.1) k p.2) := by
  simp only [bondOp, of_apply, apply_ite norm, norm_zero, ite_mul, zero_mul]
  rw [sum_ite_agree_off hjk]
  simp only [Function.update_self, Function.update_of_ne hjk]

/-! ### Lifting bond bounds to the ring -/

/-- **Lifting an upper bound.** If a Hermitian two-site operator satisfies `B ≤ b`, then so does
its embedding on two distinct sites: `b - B = Cᴴ C` embeds as `(bondOp C)ᴴ (bondOp C) ≥ 0`. -/
theorem qform_bondOp_le {j k : Fin L} (hjk : j ≠ k) {B : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) {b : ℝ} (h : ∀ y, qform B y ≤ b * nsq y) (x : (Fin L → d) → ℂ) :
    qform (bondOp j k B) x ≤ b * nsq x := by
  have hbH : ((b : ℂ) • (1 : Matrix (d × d) (d × d) ℂ) - B).IsHermitian :=
    (isHermitian_one.smul (Complex.conj_ofReal b)).sub hB
  obtain ⟨C, hC⟩ := exists_eq_conjTranspose_mul_self hbH fun y ↦ by
    have e : qform ((b : ℂ) • 1 - B) y = b * nsq y - qform B y := qform_smul_one_sub b B y
    linarith [h y]
  have hglob : (b : ℂ) • 1 - bondOp j k B = (bondOp j k C)ᴴ * bondOp j k C := by
    rw [bondOp_conjTranspose, bondOp_mul hjk, ← hC, bondOp_sub, bondOp_smul, bondOp_one]
  have hpos : 0 ≤ qform ((b : ℂ) • 1 - bondOp j k B) x := by
    rw [hglob, ← nsq_mulVec]
    exact nsq_nonneg _
  have e : qform ((b : ℂ) • 1 - bondOp j k B) x = b * nsq x - qform (bondOp j k B) x :=
    qform_smul_one_sub b _ x
  linarith

/-- **Lifting a lower bound.** If a Hermitian two-site operator satisfies `a ≤ B`, then so does
its embedding on two distinct sites. -/
theorem le_qform_bondOp {j k : Fin L} (hjk : j ≠ k) {B : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) {a : ℝ} (h : ∀ y, a * nsq y ≤ qform B y) (x : (Fin L → d) → ℂ) :
    a * nsq x ≤ qform (bondOp j k B) x := by
  have hneg : ∀ y, qform (-B) y ≤ -a * nsq y := fun y ↦ by
    rw [qform_neg_left]
    linarith [h y]
  have := qform_bondOp_le hjk hB.neg hneg x
  rw [bondOp_neg, qform_neg_left] at this
  linarith

/-! ### Ring Hamiltonians -/

/-- **Upper bound for a ring.** If both bond operators satisfy `B, Bc ≤ b`, then the twisted
ring of `L ≥ 2` sites satisfies `H ≤ b L`. -/
theorem qform_twistedRingHamiltonian_le (hL : 2 ≤ L) {B Bc : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) (hBc : Bc.IsHermitian) {b : ℝ} (hb : ∀ y, qform B y ≤ b * nsq y)
    (hbc : ∀ y, qform Bc y ≤ b * nsq y) (x : (Fin L → d) → ℂ) :
    qform (twistedRingHamiltonian B Bc L) x ≤ L * b * nsq x := by
  have := qform_sum_le Finset.univ (A := twistedBond B Bc) (b := fun _ ↦ b) (fun j _ y ↦ by
    rw [twistedBond]
    split_ifs
    · exact qform_bondOp_le (ne_finRotate hL j) hBc hbc y
    · exact qform_bondOp_le (ne_finRotate hL j) hB hb y) x
  simpa [twistedRingHamiltonian_eq_sum] using this

/-- **Lower bound for a ring.** If both bond operators satisfy `a ≤ B, Bc`, then the twisted
ring of `L ≥ 2` sites satisfies `a L ≤ H`. -/
theorem le_qform_twistedRingHamiltonian (hL : 2 ≤ L) {B Bc : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) (hBc : Bc.IsHermitian) {a : ℝ} (ha : ∀ y, a * nsq y ≤ qform B y)
    (hac : ∀ y, a * nsq y ≤ qform Bc y) (x : (Fin L → d) → ℂ) :
    L * a * nsq x ≤ qform (twistedRingHamiltonian B Bc L) x := by
  have := le_qform_sum Finset.univ (A := twistedBond B Bc) (a := fun _ ↦ a) (fun j _ y ↦ by
    rw [twistedBond]
    split_ifs
    · exact le_qform_bondOp (ne_finRotate hL j) hBc hac y
    · exact le_qform_bondOp (ne_finRotate hL j) hB ha y) x
  simpa [twistedRingHamiltonian_eq_sum] using this

/-- **Upper bound for a translation-invariant ring**: `B ≤ b` gives `H_L ≤ b L`. -/
theorem qform_ringHamiltonian_le (hL : 2 ≤ L) {B : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) {b : ℝ} (hb : ∀ y, qform B y ≤ b * nsq y) (x : (Fin L → d) → ℂ) :
    qform (ringHamiltonian B L) x ≤ L * b * nsq x := by
  rw [← twistedRingHamiltonian_self]
  exact qform_twistedRingHamiltonian_le hL hB hB hb hb x

/-- **Lower bound for a translation-invariant ring**: `a ≤ B` gives `a L ≤ H_L`. -/
theorem le_qform_ringHamiltonian (hL : 2 ≤ L) {B : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) {a : ℝ} (ha : ∀ y, a * nsq y ≤ qform B y) (x : (Fin L → d) → ℂ) :
    L * a * nsq x ≤ qform (ringHamiltonian B L) x := by
  rw [← twistedRingHamiltonian_self]
  exact le_qform_twistedRingHamiltonian hL hB hB ha ha x

/-- Summing the `L` pairs of consecutive bonds `B_{j,j+1} + B_{j+1,j+2}` counts every bond
twice. -/
theorem two_mul_qform_twistedRingHamiltonian (B Bc : Matrix (d × d) (d × d) ℂ) (L : ℕ)
    (x : (Fin L → d) → ℂ) :
    2 * qform (twistedRingHamiltonian B Bc L) x =
      ∑ j : Fin L, qform (twistedBond B Bc j + twistedBond B Bc (finRotate L j)) x := by
  simp only [qform_add_left, Finset.sum_add_distrib, twistedRingHamiltonian_eq_sum, qform_sum]
  rw [Equiv.sum_comp (finRotate L) (fun j ↦ qform (twistedBond B Bc j) x)]
  ring

/-- **Lower bound from pairs of bonds.** If every pair of consecutive bonds of the twisted ring
satisfies `a ≤ B_{j,j+1} + B_{j+1,j+2}`, then `(a/2) L ≤ H`. -/
theorem le_qform_twistedRingHamiltonian_of_pairs {B Bc : Matrix (d × d) (d × d) ℂ} {a : ℝ}
    (h : ∀ (j : Fin L) (y : (Fin L → d) → ℂ),
      a * nsq y ≤ qform (twistedBond B Bc j + twistedBond B Bc (finRotate L j)) y)
    (x : (Fin L → d) → ℂ) :
    L * a / 2 * nsq x ≤ qform (twistedRingHamiltonian B Bc L) x := by
  have hsum := two_mul_qform_twistedRingHamiltonian B Bc L x
  have hle : ∑ _j : Fin L, a * nsq x ≤
      ∑ j : Fin L, qform (twistedBond B Bc j + twistedBond B Bc (finRotate L j)) x :=
    Finset.sum_le_sum fun j _ ↦ h j x
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hle
  linarith

end CertifiedGaps
