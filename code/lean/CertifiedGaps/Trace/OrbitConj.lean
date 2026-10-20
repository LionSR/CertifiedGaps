/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Trace.Orbit

/-!
# Orbits of representative columns under antiunitary symmetries

`Trace/Orbit.lean` derives equal column norms of `F = gibbsHalf H b c` from unitary
signed-permutation symmetries of `H`. For the ring twisted by the `2π/3` rotation `C`
(`thermal120.tex`, §"The finite matrices and their columns"), reversal and negation of the words
take `H` to its complex conjugate: composed with complex conjugation they are antiunitary
symmetries. They still preserve the norm of every column of the real function `F` of `H`.

* `IsConjColumnSymmetry H e`: some unit phases `u` make `(e, u)` intertwine `H` with its entrywise
  conjugate `H̄`, i.e. `P H = H̄ P` for the signed permutation `P`.
* `IsConjColumnSymmetry.nsq_gibbsHalf`: for Hermitian `H`, `‖F e_{e r}‖ = ‖F e_r‖`. Since
  `H̄ = Hᵀ`, `Fᵀ = P F P⁻¹`, and the rows and columns of the Hermitian `F` have equal norms.
* `columnNormSymmetries H b c`: the permutations that preserve every column norm of `F`, a
  submonoid; `nsq_gibbsHalf_rep_of_norms` is `nsq_gibbsHalf_rep` for generators of either kind.
* `IsColumnSymmetry.of_map`, `IsConjColumnSymmetry.of_map`: the hypotheses from an intertwining
  identity over any commutative ring with a ring hom to `ℂ` (such as `ℤ[ω]`).
-/

open Matrix NormedSpace
open CertifiedGaps.Spectral (nsq nsq_eq_sum nsq_smul nsq_mulVec_of_unitary)

namespace CertifiedGaps.Trace

variable {n : Type*}

/-- `e` is an **antiunitary column symmetry** of `H`: with some unit phases `u`, the signed
permutation `e_τ ↦ u_τ e_{e τ}` intertwines `H` with its entrywise complex conjugate. -/
def IsConjColumnSymmetry (H : Matrix n n ℂ) (e : Equiv.Perm n) : Prop :=
  ∃ u : n → ℂ, (∀ σ, ‖u σ‖ = 1) ∧ Intertwines e u H (H.map star)

/-- A unitary column symmetry from an intertwining identity over a commutative ring `S` mapped
to `ℂ`. -/
theorem IsColumnSymmetry.of_map {S : Type*} [CommRing S] (f : S →+* ℂ) {A : Matrix n n S}
    {e : Equiv.Perm n} {u : n → S} (hu : ∀ σ, ‖f (u σ)‖ = 1) (h : Intertwines e u A A) :
    IsColumnSymmetry (A.map f) e :=
  ⟨fun σ ↦ f (u σ), hu, h.map f⟩

/-- An antiunitary column symmetry from an intertwining identity over a commutative ring `S`
with a conjugation `c` that `f` carries to complex conjugation. -/
theorem IsConjColumnSymmetry.of_map {S : Type*} [CommRing S] (f : S →+* ℂ) (c : S →+* S)
    (hc : ∀ z, f (c z) = star (f z)) {A : Matrix n n S} {e : Equiv.Perm n} {u : n → S}
    (hu : ∀ σ, ‖f (u σ)‖ = 1) (h : Intertwines e u A (A.map c)) :
    IsConjColumnSymmetry (A.map f) e := by
  refine ⟨fun σ ↦ f (u σ), hu, ?_⟩
  have hmap : (A.map c).map f = (A.map f).map star := by
    ext i j; simp [hc]
  rw [← hmap]
  exact h.map f

variable [Fintype n] [DecidableEq n]

omit [Fintype n] [DecidableEq n] in
/-- The entrywise conjugate of a Hermitian matrix is its transpose. -/
theorem map_star_eq_transpose_of_isHermitian {H : Matrix n n ℂ} (hH : H.IsHermitian) :
    H.map star = Hᵀ := by
  ext i j
  exact hH.apply j i

/-- The rows and the columns of a Hermitian matrix have the same norms:
`‖Fᵀ e_k‖ = ‖F e_k‖`. -/
theorem nsq_transpose_mulVec_single {F : Matrix n n ℂ} (hF : F.IsHermitian) (k : n) :
    nsq (Fᵀ *ᵥ Pi.single k 1) = nsq (F *ᵥ Pi.single k 1) := by
  rw [mulVec_single_one, mulVec_single_one, nsq_eq_sum, nsq_eq_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [col_apply, transpose_apply]
  rw [← hF.apply k i, norm_star]

/-- **Equal column norms along an antiunitary symmetry.** If `P H = H̄ P` for a unitary signed
permutation `P`, then `Fᵀ = P F P⁻¹` for `F = gibbsHalf H b c`, so the column of `F` at `e r` has
the norm of the row (equivalently, the column) at `r`. -/
theorem IsConjColumnSymmetry.nsq_gibbsHalf {H : Matrix n n ℂ} (hH : H.IsHermitian)
    {e : Equiv.Perm n} (h : IsConjColumnSymmetry H e) (b c : ℝ) (r : n) :
    nsq (gibbsHalf H b c *ᵥ Pi.single (e r) 1) = nsq (gibbsHalf H b c *ᵥ Pi.single r 1) := by
  obtain ⟨u, hu, hi⟩ := h
  set P := signedPerm e u with hPdef
  have hP : Pᴴ * P = 1 := signedPerm_conjTranspose_mul_self hu
  have hPinv : P⁻¹ = Pᴴ := inv_eq_left_inv hP
  have hPu : IsUnit P := (isUnit_iff_isUnit_det P).2 (isUnit_det_of_left_inverse hP)
  have hPH : P * H = Hᵀ * P := by
    rw [← map_star_eq_transpose_of_isHermitian hH]; exact signedPerm_mul hi
  set A : Matrix n n ℂ := ((-(b / 2) : ℝ) : ℂ) • (H + (c : ℂ) • 1) with hA
  -- `P A = Aᵀ P`, hence `Aᵀ = P A P⁻¹` and `exp Aᵀ = P (exp A) P⁻¹`.
  have hPA : P * A = Aᵀ * P := by
    rw [hA, transpose_smul, transpose_add, transpose_smul, transpose_one, Matrix.mul_smul,
      Matrix.smul_mul, Matrix.mul_add, Matrix.add_mul, hPH, Matrix.mul_smul, Matrix.smul_mul,
      Matrix.mul_one, Matrix.one_mul]
  have hAt : Aᵀ = P * A * P⁻¹ := by
    rw [hPA, Matrix.mul_assoc, Matrix.mul_nonsing_inv _ ((isUnit_iff_isUnit_det P).1 hPu),
      Matrix.mul_one]
  have hFt : (gibbsHalf H b c)ᵀ = P * gibbsHalf H b c * P⁻¹ := by
    rw [gibbsHalf, ← hA, ← exp_transpose, hAt, exp_conj _ _ hPu]
  have hF : (gibbsHalf H b c).IsHermitian :=
    ((hH.add (isHermitian_one.smul (Complex.conj_ofReal c))).smul
      (Complex.conj_ofReal (-(b / 2)))).exp
  -- `Fᵀ (P e_r) = P (F e_r)`, and `P e_r = u_r e_{e r}`.
  have hcomm : (gibbsHalf H b c)ᵀ * P = P * gibbsHalf H b c := by
    rw [hFt, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ ((isUnit_iff_isUnit_det P).1 hPu),
      Matrix.mul_one]
  have hv := congrArg (fun M ↦ nsq (M *ᵥ Pi.single r 1)) hcomm
  simp only [← mulVec_mulVec, hPdef, signedPerm_mulVec_single, mulVec_smul, nsq_smul, hu,
    one_pow, one_mul] at hv
  rw [← nsq_transpose_mulVec_single hF, hv, nsq_mulVec_of_unitary hP]

/-- The permutations that preserve the norm of every column of `F = gibbsHalf H b c`. -/
def columnNormSymmetries (H : Matrix n n ℂ) (b c : ℝ) : Submonoid (Equiv.Perm n) where
  carrier := {e | ∀ r, nsq (gibbsHalf H b c *ᵥ Pi.single (e r) 1) =
    nsq (gibbsHalf H b c *ᵥ Pi.single r 1)}
  one_mem' _ := rfl
  mul_mem' {e₂ e₁} h₂ h₁ r := by
    simp only [Set.mem_ofPred_eq] at h₁ h₂
    rw [Equiv.Perm.mul_apply, h₂, h₁]

/-- **Representatives from generators of either kind.** If every generator preserves the column
norms of `F` (a unitary or an antiunitary symmetry) and every representative `rep r` is the image
of `r` under a product of generators, then `‖F e_r‖ = ‖F e_{rep r}‖`. -/
theorem nsq_gibbsHalf_rep_of_norms {H : Matrix n n ℂ} {b c : ℝ} (gens : Set (Equiv.Perm n))
    (hgens : ∀ g ∈ gens, ∀ r, nsq (gibbsHalf H b c *ᵥ Pi.single (g r) 1) =
      nsq (gibbsHalf H b c *ᵥ Pi.single r 1))
    (rep : n → n) (hrep : ∀ r, ∃ g ∈ Submonoid.closure gens, g r = rep r) (r : n) :
    nsq (gibbsHalf H b c *ᵥ Pi.single r 1) = nsq (gibbsHalf H b c *ᵥ Pi.single (rep r) 1) := by
  obtain ⟨g, hg, hgr⟩ := hrep r
  have hsym : g ∈ columnNormSymmetries H b c :=
    (Submonoid.closure_le (S := columnNormSymmetries H b c)).2 hgens hg
  rw [← hgr, hsym r]

end CertifiedGaps.Trace
