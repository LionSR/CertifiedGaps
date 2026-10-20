/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.Logic.Equiv.Fin.Rotate
import CertifiedGaps.Transfer.PiTensor

/-!
# Ring Hamiltonians and Gibbs traces

Every ring Hamiltonian of the library is built from the same three pieces, and every partition
function is the same functional of its Hamiltonian:

* a two-site operator `B` on the bond `j → j + 1` (`bondOp`), with the cyclic successor
  `finRotate L j = j + 1`;
* a one-site operator `A` on the site `j` (`siteOp`);
* the Gibbs trace `gibbsTrace H β = Re Tr e^{-βH}`.

The ring `∑ⱼ Bⱼ,ⱼ₊₁` may carry a different operator `Bc` on its closing bond `L - 1 → 0`
(`twistedRingHamiltonian B Bc L`). This covers both the translation-invariant ring
(`ringHamiltonian B L`, the case `Bc = B`) and the ring twisted by an on-site symmetry, whose
partition function is the symmetry-resolved trace of the untwisted ring. Writing all of them
through one definition lets a single conjugation lemma (`conj_twistedRingHamiltonian` in
`Ring/Conjugation.lean`) change the local basis of both.

## Definitions

* `bondOp j k B`, `siteOp j A`: two-site and one-site operators on a ring of `L` sites.
* `twistedRingHamiltonian B Bc L`, `ringHamiltonian B L`, `siteSum A L`, and the `j`-th bond
  `twistedBond B Bc j` of the twisted ring.
* `gibbsTrace H β`: the partition function `Re Tr e^{-βH}` of a Hamiltonian `H`.

## Main results

* `bondOp j k` is a linear `*`-algebra embedding for `j ≠ k`: `bondOp_mul`, `bondOp_one`,
  `bondOp_conjTranspose`, `bondOp_add`, `bondOp_sub`, `bondOp_smul`, `bondOp_sum`; its matrix
  entries in a row `σ` run over the letters of the two bond sites only (`sum_ite_agree_off`).
* `gibbsTrace_eq_sum`: `Re Tr e^{-βH} = ∑ₘ e^{-β Eₘ}` for Hermitian `H`.
* `gibbsTrace_pos`: the Gibbs trace of a Hermitian matrix is positive.
* `gibbsTrace_conj`: unitary conjugation leaves the Gibbs trace unchanged; `exp_smul_one`:
  `exp (z • 1) = e^z • 1`.
* `twistedRingHamiltonian_eq_sum`: the twisted ring is the sum of its bonds `twistedBond`.
* `twistedRingHamiltonian_isHermitian`, `ringHamiltonian_isHermitian`, `siteSum_isHermitian`.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace CertifiedGaps

variable {d : Type} [Fintype d] [DecidableEq d]

/-! ### Ring geometry -/

/-- The cyclic successor in coordinates: `finRotate L j = j + 1`, wrapping `L - 1` to `0`. -/
theorem val_finRotate {L : ℕ} (j : Fin L) :
    (finRotate L j : ℕ) = if (j : ℕ) + 1 = L then 0 else (j : ℕ) + 1 := by
  cases L with
  | zero => exact j.elim0
  | succ n =>
    rw [coe_finRotate]
    split_ifs with h₁ h₂ h₂ <;> simp_all [Fin.ext_iff]

/-- On a ring of at least two sites, the two ends of a bond are distinct. -/
theorem ne_finRotate {L : ℕ} (hL : 2 ≤ L) (j : Fin L) : j ≠ finRotate L j := by
  intro h
  have h' := congrArg Fin.val h
  rw [val_finRotate] at h'
  split_ifs at h' <;> omega

/-- On an even ring the two ends of a bond have opposite parity. -/
theorem odd_finRotate_iff {L : ℕ} (hL : Even L) (j : Fin L) :
    Odd (finRotate L j : ℕ) ↔ ¬ Odd (j : ℕ) := by
  rw [val_finRotate]
  obtain ⟨k, hk⟩ := hL
  have := j.isLt
  split_ifs with h <;> simp only [Nat.odd_iff] <;> omega

/-! ### Operators on the ring -/

/-- The operator `B` on two sites `j, k` of a ring of `L` sites, acting as the identity on all
other sites. The entries may lie in any type with a zero (`ℂ` for Hamiltonians, `ℤ` for the
integer matrices of the seed computations). -/
def bondOp {R : Type*} [Zero R] {L : ℕ} (j k : Fin L) (B : Matrix (d × d) (d × d) R) :
    Matrix (Fin L → d) (Fin L → d) R :=
  of fun σ τ ↦ if ∀ i, i ≠ j → i ≠ k → σ i = τ i then B (σ j, σ k) (τ j, τ k) else 0

/-- The operator `A` on the site `j` of a ring of `L` sites, acting as the identity elsewhere. -/
def siteOp {L : ℕ} (j : Fin L) (A : Matrix d d ℂ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  piTensor fun i ↦ if i = j then A else 1

omit [Fintype d] in
/-- The conjugate transpose of an embedded two-site operator is the embedded conjugate
transpose. -/
theorem bondOp_conjTranspose {L : ℕ} (j k : Fin L) (B : Matrix (d × d) (d × d) ℂ) :
    (bondOp j k B)ᴴ = bondOp j k Bᴴ := by
  ext σ τ
  simp only [conjTranspose_apply, bondOp, of_apply]
  split_ifs with h₁ h₂ h₂
  · rfl
  · exact (h₂ fun i hj hk ↦ (h₁ i hj hk).symm).elim
  · exact (h₁ fun i hj hk ↦ (h₂ i hj hk).symm).elim
  · exact star_zero _

omit [Fintype d] in
/-- A Hermitian two-site operator embeds as a Hermitian operator on the ring. -/
theorem bondOp_isHermitian {L : ℕ} (j k : Fin L) {B : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) : (bondOp j k B).IsHermitian :=
  (bondOp_conjTranspose j k B).trans (congrArg _ hB)

omit [Fintype d] [DecidableEq d] in
/-- The conjugate transpose of a tensor product is the tensor product of the conjugate
transposes. -/
theorem conjTranspose_piTensor {ι : Type*} [Fintype ι] (A : ι → Matrix d d ℂ) :
    (piTensor A)ᴴ = piTensor fun i ↦ (A i)ᴴ := by
  ext σ τ
  simp [piTensor_apply, conjTranspose_apply, star_prod]

omit [Fintype d] in
/-- A Hermitian one-site operator embeds as a Hermitian operator on the ring. -/
theorem siteOp_isHermitian {L : ℕ} (j : Fin L) {A : Matrix d d ℂ} (hA : A.IsHermitian) :
    (siteOp j A).IsHermitian := by
  change (siteOp j A)ᴴ = siteOp j A
  rw [siteOp, conjTranspose_piTensor]
  congr 1
  ext1 i
  split_ifs
  · exact hA
  · exact conjTranspose_one

omit [Fintype d] in
/-- The zero operator at one site is zero. -/
theorem siteOp_zero {L : ℕ} (j : Fin L) : siteOp j (0 : Matrix d d ℂ) = 0 :=
  piTensor_eq_zero (i₀ := j) (by simp)

omit [Fintype d] in
/-- `bondOp` is additive in the two-site operator. -/
theorem bondOp_sum {L : ℕ} (j k : Fin L) {X : Type*} (s : Finset X)
    (B : X → Matrix (d × d) (d × d) ℂ) :
    bondOp j k (∑ x ∈ s, B x) = ∑ x ∈ s, bondOp j k (B x) := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.sum_apply]
  split_ifs <;> simp

omit [Fintype d] in
/-- `bondOp` is additive in the two-site operator. -/
theorem bondOp_add {L : ℕ} (j k : Fin L) (B C : Matrix (d × d) (d × d) ℂ) :
    bondOp j k (B + C) = bondOp j k B + bondOp j k C := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.add_apply]
  split_ifs <;> simp

omit [Fintype d] in
/-- `bondOp` is homogeneous in the two-site operator. -/
theorem bondOp_smul {L : ℕ} (j k : Fin L) (c : ℂ) (B : Matrix (d × d) (d × d) ℂ) :
    bondOp j k (c • B) = c • bondOp j k B := by
  ext σ τ
  simp only [bondOp, of_apply, Matrix.smul_apply, smul_eq_mul]
  split_ifs <;> simp

omit [Fintype d] in
/-- `bondOp` commutes with negation. -/
theorem bondOp_neg {L : ℕ} (j k : Fin L) (B : Matrix (d × d) (d × d) ℂ) :
    bondOp j k (-B) = -bondOp j k B := by
  simpa using bondOp_smul j k (-1) B

omit [Fintype d] in
/-- `bondOp` commutes with subtraction. -/
theorem bondOp_sub {L : ℕ} (j k : Fin L) (B C : Matrix (d × d) (d × d) ℂ) :
    bondOp j k (B - C) = bondOp j k B - bondOp j k C := by
  rw [sub_eq_add_neg, bondOp_add, bondOp_neg, ← sub_eq_add_neg]

omit [Fintype d] in
/-- A product operator `A ⊗ B` on two distinct sites is a tensor product over the ring. -/
theorem bondOp_kronecker {L : ℕ} {j k : Fin L} (hjk : j ≠ k) (A B : Matrix d d ℂ) :
    bondOp j k (A ⊗ₖ B) = piTensor fun i ↦ if i = j then A else if i = k then B else 1 := by
  ext σ τ
  have hk : k ∈ Finset.univ.erase j := Finset.mem_erase.2 ⟨hjk.symm, Finset.mem_univ _⟩
  rw [piTensor_apply, ← Finset.mul_prod_erase _ _ (Finset.mem_univ j),
    ← Finset.mul_prod_erase _ _ hk,
    Finset.prod_congr rfl (g := fun i ↦ if σ i = τ i then (1 : ℂ) else 0) fun i hi ↦ ?_,
    Finset.prod_ite_zero, Finset.prod_const_one]
  · simp only [bondOp, of_apply, kroneckerMap_apply, ite_true, Ne.symm hjk, ite_false,
      Finset.mem_erase, Finset.mem_univ, and_true]
    split_ifs with h₁ h₂ h₂ <;> simp_all
  · simp only [Finset.mem_erase, Finset.mem_univ, and_true] at hi
    simp [hi.1, hi.2, one_apply]

/-- Summing over the configurations `τ` that agree with `σ` away from two distinct sites
`j, k` is summing over the pair of letters `(τ j, τ k)`. -/
theorem sum_ite_agree_off {L : ℕ} {M : Type*} [AddCommMonoid M] {j k : Fin L} (hjk : j ≠ k)
    (σ : Fin L → d) (F : (Fin L → d) → M) :
    ∑ τ, (if ∀ i, i ≠ j → i ≠ k → σ i = τ i then F τ else 0) =
      ∑ p : d × d, F (Function.update (Function.update σ j p.1) k p.2) := by
  rw [← Finset.sum_filter]
  refine Finset.sum_bij' (fun τ _ ↦ (τ j, τ k))
    (fun p _ ↦ Function.update (Function.update σ j p.1) k p.2) (fun _ _ ↦ Finset.mem_univ _)
    (fun p _ ↦ ?_) (fun τ hτ ↦ ?_) (fun p _ ↦ ?_) (fun τ hτ ↦ ?_)
  · simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    intro i hj hk
    rw [Function.update_of_ne hk, Function.update_of_ne hj]
  · simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hτ
    funext i
    by_cases hk : i = k
    · subst hk; rw [Function.update_self]
    · by_cases hj : i = j
      · subst hj; rw [Function.update_of_ne hk, Function.update_self]
      · rw [Function.update_of_ne hk, Function.update_of_ne hj, hτ i hj hk]
  · rw [Function.update_of_ne hjk, Function.update_self, Function.update_self]
  · simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hτ
    congr 1
    funext i
    by_cases hk : i = k
    · subst hk; rw [Function.update_self]
    · by_cases hj : i = j
      · subst hj; rw [Function.update_of_ne hk, Function.update_self]
      · rw [Function.update_of_ne hk, Function.update_of_ne hj, hτ i hj hk]

omit [Fintype d] [DecidableEq d] in
/-- Updating the two sites of a bond does not change the other sites. -/
theorem update_update_of_ne {L : ℕ} {j k i : Fin L} (σ : Fin L → d) (p : d × d) (hj : i ≠ j)
    (hk : i ≠ k) : Function.update (Function.update σ j p.1) k p.2 i = σ i := by
  rw [Function.update_of_ne hk, Function.update_of_ne hj]

/-- `bondOp j k` is multiplicative on two distinct sites. -/
theorem bondOp_mul {L : ℕ} {j k : Fin L} (hjk : j ≠ k) (B C : Matrix (d × d) (d × d) ℂ) :
    bondOp j k B * bondOp j k C = bondOp j k (B * C) := by
  ext σ τ
  simp only [mul_apply, bondOp, of_apply, ite_mul, zero_mul]
  rw [sum_ite_agree_off hjk]
  simp only [Function.update_self, Function.update_of_ne hjk]
  by_cases h : ∀ i, i ≠ j → i ≠ k → σ i = τ i
  · rw [ite_eq_left h]
    refine Finset.sum_congr rfl fun p _ ↦ ?_
    rw [ite_eq_left fun i hj hk ↦ (update_update_of_ne σ p hj hk).trans (h i hj hk)]
  · rw [ite_eq_right h]
    refine Finset.sum_eq_zero fun p _ ↦ ?_
    rw [ite_eq_right fun h' ↦ h fun i hj hk ↦ (update_update_of_ne σ p hj hk).symm.trans
      (h' i hj hk), mul_zero]

omit [Fintype d] in
/-- `bondOp j k` is unital. -/
theorem bondOp_one {L : ℕ} (j k : Fin L) :
    bondOp j k (1 : Matrix (d × d) (d × d) ℂ) = 1 := by
  ext σ τ
  simp only [bondOp, of_apply, one_apply, Prod.mk.injEq]
  by_cases hστ : σ = τ
  · subst hστ; simp
  · rw [ite_eq_right hστ]
    split_ifs with h₁ h₂
    · exact (hστ (funext fun i ↦ if hj : i = j then hj ▸ h₂.1 else
        if hk : i = k then hk ▸ h₂.2 else h₁ i hj hk)).elim
    · rfl
    · rfl

/-! ### The trace of a Gibbs operator -/

/-- For a Hermitian matrix `A` and a real scalar `c`, `Tr exp (c A) = ∑ᵢ exp (c λᵢ)`. -/
theorem trace_exp_smul_eq_sum_of_isHermitian {n : Type*} [Fintype n] [DecidableEq n]
    {A : Matrix n n ℂ} (hA : A.IsHermitian) (c : ℝ) :
    trace (exp ((c : ℂ) • A)) = ∑ i, (Real.exp (c * hA.eigenvalues i) : ℂ) := by
  set U : Matrix n n ℂ := (hA.eigenvectorUnitary : Matrix n n ℂ)
  have hU : star U * U = 1 := Unitary.coe_star_mul_self _
  have hUinv : U⁻¹ = star U := inv_eq_left_inv hU
  have hUunit : IsUnit U := (isUnit_iff_isUnit_det U).mpr
    (isUnit_det_of_left_inverse hU)
  have hA' : (c : ℂ) • A = U * diagonal (fun i ↦ ((c * hA.eigenvalues i : ℝ) : ℂ)) * U⁻¹ := by
    conv_lhs => rw [hA.spectral_theorem]
    rw [Unitary.conjStarAlgAut_apply, hUinv, ← smul_mul_assoc, ← mul_smul_comm, ← diagonal_smul]
    congr
    ext i
    simp
  rw [hA', exp_conj _ _ hUunit, exp_diagonal, trace_mul_cycle, hUinv, hU, one_mul,
    trace_diagonal]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [Pi.coe_exp, ← Complex.exp_eq_exp_ℂ, Complex.ofReal_exp]

/-- Equal Hermitian matrices have equal eigenvalue families (the Hermiticity proofs are
irrelevant). -/
theorem IsHermitian.eigenvalues_eq_of_eq {n : Type*} [Fintype n] [DecidableEq n]
    {A B : Matrix n n ℂ} (hA : A.IsHermitian) (hB : B.IsHermitian) (h : A = B) :
    hA.eigenvalues = hB.eigenvalues := by
  subst h; rfl

/-- Unitary conjugation leaves `Tr exp (c H)` unchanged. -/
theorem trace_exp_smul_conj {n : Type*} [Fintype n] [DecidableEq n] (P H : Matrix n n ℂ)
    (hP : Pᴴ * P = 1) (c : ℂ) :
    trace (exp (c • (Pᴴ * H * P))) = trace (exp (c • H)) := by
  have hinv : P⁻¹ = Pᴴ := inv_eq_left_inv hP
  have hunit : IsUnit P := (isUnit_iff_isUnit_det P).mpr (isUnit_det_of_left_inverse hP)
  rw [← hinv, ← Matrix.smul_mul, ← Matrix.mul_smul, exp_conj' _ _ hunit, Matrix.mul_assoc,
    trace_mul_comm, Matrix.mul_assoc, mul_nonsing_inv _ ((isUnit_iff_isUnit_det P).mp hunit),
    Matrix.mul_one]

/-! ### The Gibbs trace -/

section Gibbs

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The exponential of a complex scalar matrix is the scalar matrix of the exponential. -/
theorem exp_smul_one (z : ℂ) : exp (z • (1 : Matrix n n ℂ)) = exp z • (1 : Matrix n n ℂ) := by
  rw [smul_one_eq_diagonal, smul_one_eq_diagonal, exp_diagonal, Pi.exp_def]

/-- The Gibbs trace `Re Tr e^{-βH}`: the partition function of the Hamiltonian `H` at inverse
temperature `β`. The real part is taken so that the value is a real number for every `H`; for
Hermitian `H` the trace is already real (`gibbsTrace_eq_sum`). -/
noncomputable def gibbsTrace (H : Matrix n n ℂ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • H))).re

/-- For Hermitian `H`, `Re Tr e^{-βH} = ∑ₘ e^{-β Eₘ}` over the eigenvalues with
multiplicity. -/
theorem gibbsTrace_eq_sum {H : Matrix n n ℂ} (hH : H.IsHermitian) (β : ℝ) :
    gibbsTrace H β = ∑ i, Real.exp (-(β * hH.eigenvalues i)) := by
  have h := trace_exp_smul_eq_sum_of_isHermitian hH (-β)
  rw [Complex.ofReal_neg] at h
  rw [gibbsTrace, h, ← Complex.ofReal_sum, Complex.ofReal_re]
  simp [neg_mul]

/-- The Gibbs trace of a Hermitian matrix on a nonempty index type is positive. -/
theorem gibbsTrace_pos [Nonempty n] {H : Matrix n n ℂ} (hH : H.IsHermitian) (β : ℝ) :
    0 < gibbsTrace H β := by
  rw [gibbsTrace_eq_sum hH]
  exact Finset.sum_pos (fun _ _ ↦ Real.exp_pos _) Finset.univ_nonempty

/-- Unitary conjugation leaves the Gibbs trace unchanged. -/
theorem gibbsTrace_conj (P H : Matrix n n ℂ) (hP : Pᴴ * P = 1) (β : ℝ) :
    gibbsTrace (Pᴴ * H * P) β = gibbsTrace H β := by
  rw [gibbsTrace, gibbsTrace, trace_exp_smul_conj _ _ hP]

end Gibbs

/-! ### Ring Hamiltonians -/

/-- The ring Hamiltonian `∑ⱼ Bⱼ,ⱼ₊₁` on `L` sites whose closing bond `L - 1 → 0` carries `Bc`
in place of `B`. With `Bc` the bond conjugated by an on-site symmetry on one end, its Gibbs
trace is the twisted partition function `Tr (P e^{-βH_L})` of the untwisted ring. -/
def twistedRingHamiltonian {R : Type*} [AddCommMonoid R] (B Bc : Matrix (d × d) (d × d) R)
    (L : ℕ) : Matrix (Fin L → d) (Fin L → d) R :=
  ∑ j : Fin L, bondOp j (finRotate L j) (if j.val + 1 = L then Bc else B)

/-- The `j`-th bond `B_{j,j+1}` of the twisted ring `twistedRingHamiltonian B Bc L`: `Bc` on the
closing bond `L - 1 → 0` and `B` elsewhere. -/
def twistedBond {R : Type*} [Zero R] {L : ℕ} (B Bc : Matrix (d × d) (d × d) R) (j : Fin L) :
    Matrix (Fin L → d) (Fin L → d) R :=
  bondOp j (finRotate L j) (if j.val + 1 = L then Bc else B)

omit [Fintype d] in
/-- The twisted ring is the sum of its bonds. -/
theorem twistedRingHamiltonian_eq_sum (B Bc : Matrix (d × d) (d × d) ℂ) (L : ℕ) :
    twistedRingHamiltonian B Bc L = ∑ j : Fin L, twistedBond B Bc j := rfl

/-- The translation-invariant ring Hamiltonian `∑ⱼ Bⱼ,ⱼ₊₁` of a two-site operator `B` on `L`
sites, with the cyclic successor `finRotate L j = j + 1`. For `L = 2` both bonds join the sites
`0` and `1`. -/
def ringHamiltonian {R : Type*} [AddCommMonoid R] (B : Matrix (d × d) (d × d) R) (L : ℕ) :
    Matrix (Fin L → d) (Fin L → d) R :=
  ∑ j : Fin L, bondOp j (finRotate L j) B

omit [Fintype d] in
/-- The untwisted ring is the twisted ring with `Bc = B`. -/
theorem twistedRingHamiltonian_self (B : Matrix (d × d) (d × d) ℂ) (L : ℕ) :
    twistedRingHamiltonian B B L = ringHamiltonian B L := by
  simp only [twistedRingHamiltonian, ite_self, ringHamiltonian]

omit [Fintype d] in
/-- A zero-preserving map of entries commutes with embedding a bond: an integer bond embeds to
the cast of its integer embedding. -/
theorem bondOp_map {R S : Type*} [Zero R] [Zero S] {L : ℕ} (j k : Fin L)
    (B : Matrix (d × d) (d × d) R) (f : R → S) (hf : f 0 = 0) :
    (bondOp j k B).map f = bondOp j k (B.map f) := by
  ext σ τ
  simp only [bondOp, map_apply, of_apply]
  split_ifs <;> simp [hf]

omit [Fintype d] in
/-- An additive map of entries commutes with the twisted ring Hamiltonian. -/
theorem twistedRingHamiltonian_map {R S : Type*} [AddCommMonoid R] [AddCommMonoid S]
    (B Bc : Matrix (d × d) (d × d) R) (L : ℕ) (f : R →+ S) :
    (twistedRingHamiltonian B Bc L).map f = twistedRingHamiltonian (B.map f) (Bc.map f) L := by
  rw [twistedRingHamiltonian, ← AddMonoidHom.mapMatrix_apply, map_sum,
    twistedRingHamiltonian]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [AddMonoidHom.mapMatrix_apply, bondOp_map _ _ _ _ f.map_zero]
  split_ifs <;> rfl

/-- The sum `∑ⱼ Aⱼ` of a one-site operator over the `L` sites of a ring. -/
def siteSum (A : Matrix d d ℂ) (L : ℕ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  ∑ j : Fin L, siteOp j A

omit [Fintype d] in
/-- The site sum of the zero operator is zero. -/
theorem siteSum_zero (L : ℕ) : siteSum (0 : Matrix d d ℂ) L = 0 := by
  simp [siteSum, siteOp_zero]

omit [Fintype d] in
/-- The twisted ring of Hermitian bonds is Hermitian. -/
theorem twistedRingHamiltonian_isHermitian {B Bc : Matrix (d × d) (d × d) ℂ}
    (hB : B.IsHermitian) (hBc : Bc.IsHermitian) (L : ℕ) :
    (twistedRingHamiltonian B Bc L).IsHermitian :=
  isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian j (finRotate L j) (by split_ifs <;> assumption)

omit [Fintype d] in
/-- The ring Hamiltonian of a Hermitian bond is Hermitian. -/
theorem ringHamiltonian_isHermitian {B : Matrix (d × d) (d × d) ℂ} (hB : B.IsHermitian)
    (L : ℕ) : (ringHamiltonian B L).IsHermitian :=
  isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian j (finRotate L j) hB

omit [Fintype d] in
/-- The site sum of a Hermitian operator is Hermitian. -/
theorem siteSum_isHermitian {A : Matrix d d ℂ} (hA : A.IsHermitian) (L : ℕ) :
    (siteSum A L).IsHermitian :=
  isSelfAdjoint_sum _ fun j _ ↦ siteOp_isHermitian j hA

end CertifiedGaps
