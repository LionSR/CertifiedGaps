/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Spectral.Counting
import Mathlib.LinearAlgebra.Matrix.Block

/-!
# Inertia from an `L D Lᴴ` factorization

A factorization `T = L D Lᴴ` with `D = diagonal d` real writes the form of `T` as a weighted sum
of squares, `Re ⟪c, T c⟫ = ∑ₗ dₗ ‖(Lᴴ c)ₗ‖²` (`qform_LDL`). Hence:

* `exists_nonneg_subspace_of_LDL`: for **any** (even rectangular) `L`, the form of `T` is `≥ 0` on
  the subspace `{c : (Lᴴ c)ₗ = 0 whenever dₗ < 0}`, of codimension at most the number of
  negative pivots;
* `exists_neg_subspace_of_LDL`: if `L` is square and invertible, the form is `< 0` on a subspace
  of dimension exactly the number of negative pivots (used for lower bounds on counts).

`LDLCert` is the exact rational certificate: entries in `ℚ`, the identity `L D Lᵀ = T` as a
decidable entrywise statement (`LDLCert.Factors`, summed by `List.foldr` so that `decide +kernel`
evaluates it), and the decidable shape condition `LDLCert.UnitLower`. Its theorems
`LDLCert.exists_nonneg_subspace` and `LDLCert.exists_neg_subspace` transport the conclusions to
the cast matrix over any `RCLike` field.
-/

open Matrix Module

namespace CertifiedGaps.Spectral

variable {𝕜 : Type*} [RCLike 𝕜] {k p : Type*} [Fintype k] [Fintype p] [DecidableEq p]

/-- The form of a real diagonal matrix is `∑ₗ dₗ ‖cₗ‖²`. -/
theorem qform_diagonal (d : p → ℝ) (c : p → 𝕜) :
    qform (diagonal fun l ↦ (d l : 𝕜)) c = ∑ l, d l * ‖c l‖ ^ 2 := by
  simp only [qform, mulVec_diagonal, dotProduct, Pi.star_apply, map_sum]
  refine Finset.sum_congr rfl fun l _ ↦ ?_
  rw [mul_left_comm, RCLike.star_def, RCLike.conj_mul, ← RCLike.ofReal_pow,
    ← RCLike.ofReal_mul, RCLike.ofReal_re]

/-- The form of `L D Lᴴ` is the weighted sum of squares `∑ₗ dₗ ‖(Lᴴ c)ₗ‖²`. -/
theorem qform_LDL (L : Matrix k p 𝕜) (d : p → ℝ) (c : k → 𝕜) :
    qform (L * diagonal (fun l ↦ (d l : 𝕜)) * Lᴴ) c = ∑ l, d l * ‖(Lᴴ *ᵥ c) l‖ ^ 2 := by
  rw [← qform_diagonal, qform_mulVec, conjTranspose_conjTranspose]

/-- **Inertia, upper half.** If `T = L D Lᴴ` then the form of `T` is nonnegative on a subspace of
codimension at most `#{l : dₗ < 0}`. No shape condition on `L` is needed. -/
theorem exists_nonneg_subspace_of_LDL (L : Matrix k p 𝕜) (d : p → ℝ) :
    ∃ P : Submodule 𝕜 (k → 𝕜),
      Fintype.card k ≤ finrank 𝕜 P + (Finset.univ.filter fun l ↦ d l < 0).card ∧
      ∀ c ∈ P, 0 ≤ qform (L * diagonal (fun l ↦ (d l : 𝕜)) * Lᴴ) c := by
  classical
  set neg := Finset.univ.filter fun l ↦ d l < 0
  let f : (k → 𝕜) →ₗ[𝕜] (neg → 𝕜) :=
    (LinearMap.pi fun l : neg ↦ LinearMap.proj (R := 𝕜) (φ := fun _ : p ↦ 𝕜) l.1) ∘ₗ
      Lᴴ.mulVecLin
  refine ⟨LinearMap.ker f, ?_, fun c hc ↦ ?_⟩
  · have h := LinearMap.finrank_range_add_finrank_ker f
    have hr : finrank 𝕜 (LinearMap.range f) ≤ neg.card := by
      simpa [Module.finrank_fintype_fun_eq_card] using Submodule.finrank_le (LinearMap.range f)
    rw [Module.finrank_fintype_fun_eq_card] at h
    omega
  · rw [qform_LDL]
    refine Finset.sum_nonneg fun l _ ↦ ?_
    by_cases hl : d l < 0
    · have := congrFun (LinearMap.mem_ker.1 hc) ⟨l, Finset.mem_filter.2 ⟨Finset.mem_univ _, hl⟩⟩
      simp only [f, LinearMap.comp_apply, LinearMap.pi_apply, LinearMap.proj_apply,
        Matrix.mulVecLin_apply, Pi.zero_apply] at this
      simp [this]
    · exact mul_nonneg (not_lt.1 hl) (by positivity)

/-- **Inertia, lower half.** If `T = L D Lᴴ` with `L` square and invertible, then the form of `T`
is negative on every nonzero vector of a subspace of dimension `#{l : dₗ < 0}`. -/
theorem exists_neg_subspace_of_LDL [DecidableEq k] (L : Matrix k k 𝕜) (hL : IsUnit L.det)
    (d : k → ℝ) :
    ∃ P : Submodule 𝕜 (k → 𝕜), finrank 𝕜 P = (Finset.univ.filter fun l ↦ d l < 0).card ∧
      ∀ c ∈ P, c ≠ 0 → qform (L * diagonal (fun l ↦ (d l : 𝕜)) * Lᴴ) c < 0 := by
  set neg := Finset.univ.filter fun l ↦ d l < 0
  have hLH : IsUnit Lᴴ.det := by
    rw [det_conjTranspose]
    exact hL.star
  have hinv : Lᴴ * (Lᴴ)⁻¹ = 1 := mul_nonsing_inv _ hLH
  have hinv' : (Lᴴ)⁻¹ * Lᴴ = 1 := nonsing_inv_mul _ hLH
  let g : (neg → 𝕜) →ₗ[𝕜] (k → 𝕜) := ((Lᴴ)⁻¹).mulVecLin ∘ₗ extendZero neg
  have hg : Function.Injective g := by
    refine fun c c' hcc' ↦ extendZero_injective neg ?_
    have := congrArg (fun x ↦ Lᴴ *ᵥ x) hcc'
    simpa [g, mulVec_mulVec, hinv] using this
  refine ⟨LinearMap.range g, ?_, fun c hc hc0 ↦ ?_⟩
  · rw [LinearMap.finrank_range_of_inj hg, Module.finrank_fintype_fun_eq_card,
      Fintype.card_coe]
  obtain ⟨e, rfl⟩ := hc
  have hLc : Lᴴ *ᵥ g e = extendZero neg e := by
    simp [g, mulVec_mulVec, hinv]
  obtain ⟨j, hj⟩ : ∃ j, extendZero neg e j ≠ 0 := by
    by_contra h
    push Not at h
    apply hc0
    have h0 : extendZero neg e = 0 := funext h
    have := congrArg (fun x ↦ (Lᴴ)⁻¹ *ᵥ x) hLc
    simpa [mulVec_mulVec, hinv', h0] using this
  have hjneg : j ∈ neg := by
    by_contra h
    exact hj (extendZero_apply_of_not_mem neg e h)
  rw [qform_LDL, hLc]
  refine Finset.sum_neg' (s := Finset.univ) ?_ ⟨j, Finset.mem_univ _, ?_⟩
  · -- Only the indices with `dₗ < 0` carry a nonzero entry.
    intro l _
    by_cases hl : l ∈ neg
    · exact mul_nonpos_of_nonpos_of_nonneg (Finset.mem_filter.1 hl).2.le (by positivity)
    · simp [extendZero_apply_of_not_mem neg e hl]
  · exact mul_neg_of_neg_of_pos (Finset.mem_filter.1 hjneg).2 (by positivity)

/-- A unit lower triangular matrix has determinant `1`, so it is invertible. -/
theorem isUnit_det_of_unitLower {R : Type*} [CommRing R] {m : Type*} [Fintype m] [DecidableEq m]
    [LinearOrder m] (L : Matrix m m R) (hdiag : ∀ i, L i i = 1)
    (hlow : ∀ i j, i < j → L i j = 0) : IsUnit L.det := by
  rw [det_of_isLowerTriangular L fun i j hij ↦ hlow i j hij]
  simp [hdiag]

/-! ### Exact rational certificates -/

/-- An exact `L D Lᵀ` certificate of size `m` over `ℚ`: the factor `L` and the pivots `d`. -/
structure LDLCert (m : ℕ) where
  /-- The factor (unit lower triangular when `UnitLower` holds). -/
  L : Matrix (Fin m) (Fin m) ℚ
  /-- The pivots, the diagonal of `D`. -/
  d : Fin m → ℚ

namespace LDLCert

variable {m : ℕ} (C : LDLCert m)

/-- The `(i, j)` entry of `L D Lᵀ`, summed with `List.foldr` so the kernel can evaluate it. -/
def entry (i j : Fin m) : ℚ :=
  (List.finRange m).foldr (fun l a ↦ C.L i l * C.d l * C.L j l + a) 0

/-- `L D Lᵀ = T`, entrywise. -/
def Factors (T : Matrix (Fin m) (Fin m) ℚ) : Prop := ∀ i j, C.entry i j = T i j

instance (T : Matrix (Fin m) (Fin m) ℚ) : Decidable (C.Factors T) :=
  @Nat.decidableForallFin m _ fun i ↦ @Nat.decidableForallFin m _ fun j ↦
    instDecidableEqRat (C.entry i j) (T i j)

/-- `L` is unit lower triangular. -/
def UnitLower : Prop := ∀ i j : Fin m, (i = j → C.L i j = 1) ∧ (i < j → C.L i j = 0)

instance : Decidable C.UnitLower :=
  @Nat.decidableForallFin m _ fun i ↦ @Nat.decidableForallFin m _ fun j ↦
    @instDecidableAnd _ _ (@instDecidableForall _ _ (instDecidableEqFin m i j)
      (instDecidableEqRat (C.L i j) 1))
      (@instDecidableForall _ _ (Fin.decLt i j) (instDecidableEqRat (C.L i j) 0))

/-- The number of negative pivots. -/
def negCount : ℕ := (Finset.univ.filter fun l ↦ C.d l < 0).card

theorem entry_eq_sum (i j : Fin m) : C.entry i j = ∑ l, C.L i l * C.d l * C.L j l := by
  rw [entry, Fin.sum_univ_def, List.sum, List.foldr_map]

/-- The certificate identity, cast to an `RCLike` field: `T = L D Lᴴ`. -/
theorem map_eq_LDL {T : Matrix (Fin m) (Fin m) ℚ} (h : C.Factors T) :
    T.map (fun q : ℚ ↦ (q : 𝕜)) =
      C.L.map (fun q : ℚ ↦ (q : 𝕜)) * diagonal (fun l ↦ ((C.d l : ℝ) : 𝕜)) *
        (C.L.map fun q : ℚ ↦ (q : 𝕜))ᴴ := by
  ext i j
  rw [map_apply, ← h i j, entry_eq_sum, mul_apply]
  simp only [mul_diagonal, map_apply, conjTranspose_apply, star_ratCast, RCLike.ofReal_ratCast]
  push_cast
  rfl

theorem negCount_eq : C.negCount = (Finset.univ.filter fun l ↦ (C.d l : ℝ) < 0).card := by
  simp [negCount]

/-- A valid certificate makes the form of `T` nonnegative on a subspace of codimension at most
`negCount`. -/
theorem exists_nonneg_subspace {T : Matrix (Fin m) (Fin m) ℚ} (h : C.Factors T) :
    ∃ P : Submodule 𝕜 (Fin m → 𝕜), m ≤ finrank 𝕜 P + C.negCount ∧
      ∀ c ∈ P, 0 ≤ qform (T.map fun q : ℚ ↦ (q : 𝕜)) c := by
  obtain ⟨P, hP, hP'⟩ := exists_nonneg_subspace_of_LDL (C.L.map fun q : ℚ ↦ (q : 𝕜))
    (fun l ↦ (C.d l : ℝ))
  refine ⟨P, ?_, fun c hc ↦ ?_⟩
  · rw [negCount_eq]; simpa using hP
  · rw [C.map_eq_LDL h]; exact hP' c hc

/-- A valid unit lower triangular certificate makes the form of `T` negative on a subspace of
dimension `negCount`. -/
theorem exists_neg_subspace {T : Matrix (Fin m) (Fin m) ℚ} (h : C.Factors T)
    (hL : C.UnitLower) :
    ∃ P : Submodule 𝕜 (Fin m → 𝕜), finrank 𝕜 P = C.negCount ∧
      ∀ c ∈ P, c ≠ 0 → qform (T.map fun q : ℚ ↦ (q : 𝕜)) c < 0 := by
  have hdet : IsUnit (C.L.map fun q : ℚ ↦ (q : 𝕜)).det :=
    isUnit_det_of_unitLower _ (fun i ↦ by simp [(hL i i).1 rfl])
      (fun i j hij ↦ by simp [(hL i j).2 hij])
  obtain ⟨P, hP, hP'⟩ := exists_neg_subspace_of_LDL _ hdet (fun l ↦ (C.d l : ℝ))
  refine ⟨P, ?_, fun c hc hc0 ↦ ?_⟩
  · rw [negCount_eq, hP]
  · rw [C.map_eq_LDL h]; exact hP' c hc hc0

end LDLCert

/-- **Rayleigh–Ritz with a rational inertia certificate.** For any witness `Z : n × K`, if a
rational matrix `T` with a unit lower triangular `L D Lᵀ` certificate has form above that of
`Zᴴ A Z - t Zᴴ Z`, then at least `C.negCount` eigenvalues of `A` lie strictly below `t`. -/
theorem le_eigCount_of_compression_LDL {n : Type*} [Fintype n] [DecidableEq n] {A : Matrix n n 𝕜}
    (hA : A.IsHermitian) {K : ℕ} (t : ℝ) (Z : Matrix n (Fin K) 𝕜) (C : LDLCert K)
    (T : Matrix (Fin K) (Fin K) ℚ) (hC : C.Factors T) (hL : C.UnitLower)
    (hT : ∀ c, qform (Zᴴ * A * Z - (t : 𝕜) • (Zᴴ * Z)) c ≤ qform (T.map fun q : ℚ ↦ (q : 𝕜)) c) :
    C.negCount ≤ eigCount hA t := by
  obtain ⟨P, hPdim, hP⟩ := C.exists_neg_subspace (𝕜 := 𝕜) hC hL
  rw [← hPdim]
  refine le_eigCount_of_compression hA Z P fun c hc hc0 ↦ ?_
  have h := (hT c).trans_lt (hP c hc hc0)
  rw [qform_sub_left, qform_smul_left] at h
  linarith

end CertifiedGaps.Spectral
