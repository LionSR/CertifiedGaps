/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Spectral.Basic
import Mathlib.LinearAlgebra.Matrix.PosDef

/-!
# Counting eigenvalues by subspaces

All results are dimension counts against the eigenspaces of a Hermitian matrix `A` (over any
`RCLike` field). Write `N(t) = eigCount hA t = #{i : λᵢ(A) < t}`.

* `eigCount_le_of_forall_le_qform` (**count by subspace**): if `t ‖x‖² ≤ Re ⟪x, A x⟫` on a subspace
  `V` with `card n ≤ dim V + K` (codimension at most `K`), then `N(t) ≤ K`;
  `eigCount_le_of_ker` is the special case `V = ker f` for a linear map `f` to `ι → 𝕜`.
* `le_eigCount_of_forall_qform_lt` (the converse direction): if `Re ⟪x, A x⟫ < t ‖x‖²` for all
  nonzero `x ∈ V`, then `dim V ≤ N(t)`.
* `le_eigCount_of_compression` (**Rayleigh–Ritz with a Gram matrix**): for any `Z : n × k`, if
  `Re ⟪c, Zᴴ A Z c⟫ < t Re ⟪c, Zᴴ Z c⟫` on nonzero `c ∈ P`, then `dim P ≤ N(t)`. No orthonormality
  of the columns of `Z` is needed. `eigCount_compression_le` is Cauchy interlacing in counting form
  (`Zᴴ Z = 1`).
* `eigCount_le_of_compression` (**count through a compression**): if `range Z + S` is the whole
  space and the form of `A - t` is `≥ 0` on every `Z c + u` with `c ∈ P`, `u ∈ S`, then
  `N(t) + dim P ≤ k`.
* `eigCount_le_eigCount_of_qform_le`, `eigCount_le_eigCount_of_posSemidef` (**Weyl monotonicity**
  in counting form): `B ≤ A` in the Loewner order implies `N_A(t) ≤ N_B(t)` for all `t`, i.e. the
  sorted eigenvalues satisfy `λ_k(B) ≤ λ_k(A)`.
  `eigenvalues₀_le_of_posSemidef` is the same statement for the sorted `eigenvalues₀`.
* `exists_eigenvalue_near` (**residual bound**): some eigenvalue satisfies
  `(λᵢ - θ)² ‖x‖² ≤ ‖A x - θ x‖²`.
-/

open Matrix Module
open scoped ComplexOrder

namespace CertifiedGaps.Spectral

variable {𝕜 : Type*} [RCLike 𝕜] {n : Type*} [Fintype n]

/-- Two subspaces whose dimensions add up to more than the ambient dimension meet nontrivially. -/
theorem exists_ne_zero_mem_inf (E V : Submodule 𝕜 (n → 𝕜))
    (h : Fintype.card n < finrank 𝕜 E + finrank 𝕜 V) : ∃ x ∈ E, x ∈ V ∧ x ≠ 0 := by
  have hsum := Submodule.finrank_sup_add_finrank_inf_eq E V
  have hle : finrank 𝕜 ↥(E ⊔ V) ≤ Fintype.card n := by
    simpa [Module.finrank_fintype_fun_eq_card] using Submodule.finrank_le (E ⊔ V)
  have hne : E ⊓ V ≠ ⊥ := by
    intro hbot
    rw [hbot, finrank_bot] at hsum
    omega
  obtain ⟨x, hx, hx0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hne
  exact ⟨x, (Submodule.mem_inf.1 hx).1, (Submodule.mem_inf.1 hx).2, hx0⟩

/-! ### Spans of eigenvectors -/

variable [DecidableEq n]

/-- Extension by zero from a finite set of indices. -/
def extendZero (s : Finset n) : (s → 𝕜) →ₗ[𝕜] (n → 𝕜) where
  toFun c i := if h : i ∈ s then c ⟨i, h⟩ else 0
  map_add' c d := by ext i; by_cases h : i ∈ s <;> simp [h]
  map_smul' a c := by ext i; by_cases h : i ∈ s <;> simp [h]

omit [Fintype n] in
theorem extendZero_apply_of_mem (s : Finset n) (c : s → 𝕜) {i : n} (h : i ∈ s) :
    extendZero s c i = c ⟨i, h⟩ := by simp [extendZero, h]

omit [Fintype n] in
theorem extendZero_apply_of_not_mem (s : Finset n) (c : s → 𝕜) {i : n} (h : i ∉ s) :
    extendZero s c i = 0 := by simp [extendZero, h]

omit [Fintype n] in
theorem extendZero_injective (s : Finset n) : Function.Injective (extendZero (𝕜 := 𝕜) s) := by
  intro c d hcd
  ext ⟨i, hi⟩
  have := congrFun hcd i
  rwa [extendZero_apply_of_mem s c hi, extendZero_apply_of_mem s d hi] at this

variable {A : Matrix n n 𝕜} (hA : A.IsHermitian)

/-- The span of the eigenvectors `hA.eigenvectorBasis i`, `i ∈ s`. -/
noncomputable def eigSpan (s : Finset n) : Submodule 𝕜 (n → 𝕜) :=
  LinearMap.range ((eigU hA).mulVecLin ∘ₗ extendZero s)

theorem finrank_eigSpan (s : Finset n) : finrank 𝕜 (eigSpan hA s) = s.card := by
  have hinj : Function.Injective ((eigU hA).mulVecLin ∘ₗ extendZero (𝕜 := 𝕜) s) := by
    refine fun c d hcd ↦ extendZero_injective s ?_
    have := congrArg (fun x ↦ (eigU hA)ᴴ *ᵥ x) hcd
    simpa [mulVec_mulVec, eigU_star_mul] using this
  rw [eigSpan, LinearMap.finrank_range_of_inj hinj, Module.finrank_fintype_fun_eq_card,
    Fintype.card_coe]

/-- On the span of eigenvectors with eigenvalues `≥ t`, the form is `≥ t ‖x‖²`. -/
theorem le_qform_of_mem_eigSpan {s : Finset n} {t : ℝ} (hs : ∀ i ∈ s, t ≤ hA.eigenvalues i)
    {x : n → 𝕜} (hx : x ∈ eigSpan hA s) : t * nsq x ≤ qform A x := by
  obtain ⟨c, rfl⟩ := hx
  simp only [LinearMap.comp_apply, Matrix.mulVecLin_apply]
  rw [qform_eigU, nsq_eigU, nsq_eq_sum, Finset.mul_sum]
  refine Finset.sum_le_sum fun i _ ↦ ?_
  by_cases hi : i ∈ s
  · exact mul_le_mul_of_nonneg_right (hs i hi) (by positivity)
  · simp [extendZero_apply_of_not_mem s c hi]

/-- On the span of eigenvectors with eigenvalues `< t`, the form is `< t ‖x‖²` for `x ≠ 0`. -/
theorem qform_lt_of_mem_eigSpan {s : Finset n} {t : ℝ} (hs : ∀ i ∈ s, hA.eigenvalues i < t)
    {x : n → 𝕜} (hx : x ∈ eigSpan hA s) (hx0 : x ≠ 0) : qform A x < t * nsq x := by
  obtain ⟨c, rfl⟩ := hx
  simp only [LinearMap.comp_apply, Matrix.mulVecLin_apply] at hx0 ⊢
  set y := extendZero s c
  obtain ⟨j, hj⟩ : ∃ j, y j ≠ 0 := by
    by_contra h
    push Not at h
    exact hx0 (by rw [show y = 0 from funext h, mulVec_zero])
  have hjs : j ∈ s := by
    by_contra h
    exact hj (extendZero_apply_of_not_mem s c h)
  rw [qform_eigU, nsq_eigU, nsq_eq_sum, Finset.mul_sum]
  refine Finset.sum_lt_sum (fun i _ ↦ ?_) ⟨j, Finset.mem_univ _, ?_⟩
  · by_cases hi : i ∈ s
    · exact mul_le_mul_of_nonneg_right (hs i hi).le (by positivity)
    · simp [y, extendZero_apply_of_not_mem s c hi]
  · exact mul_lt_mul_of_pos_right (hs j hjs) (by positivity)

/-! ### Counting -/

/-- **Count by subspace.** If the form of `A - t` is nonnegative on a subspace `V` of codimension
at most `K`, then at most `K` eigenvalues lie strictly below `t`. -/
theorem eigCount_le_of_forall_le_qform {t : ℝ} (V : Submodule 𝕜 (n → 𝕜)) {K : ℕ}
    (hdim : Fintype.card n ≤ finrank 𝕜 V + K) (hV : ∀ x ∈ V, t * nsq x ≤ qform A x) :
    eigCount hA t ≤ K := by
  set s := Finset.univ.filter fun i ↦ hA.eigenvalues i < t
  by_contra hlt
  have hcard : Fintype.card n < finrank 𝕜 (eigSpan hA s) + finrank 𝕜 V := by
    rw [finrank_eigSpan]
    have : K < s.card := not_le.1 hlt
    omega
  obtain ⟨x, hxE, hxV, hx0⟩ := exists_ne_zero_mem_inf _ _ hcard
  exact (hV x hxV).not_gt
    (qform_lt_of_mem_eigSpan hA (fun i hi ↦ (Finset.mem_filter.1 hi).2) hxE hx0)

/-- **Count by a kernel.** If the form of `A - t` is nonnegative on the kernel of a linear map
to `ι → 𝕜`, then at most `card ι` eigenvalues lie strictly below `t`. -/
theorem eigCount_le_of_ker {ι : Type*} [Fintype ι] {t : ℝ} (f : (n → 𝕜) →ₗ[𝕜] (ι → 𝕜))
    (hV : ∀ x, f x = 0 → t * nsq x ≤ qform A x) : eigCount hA t ≤ Fintype.card ι := by
  refine eigCount_le_of_forall_le_qform hA (LinearMap.ker f) ?_ fun x hx ↦ hV x hx
  have h := LinearMap.finrank_range_add_finrank_ker f
  have hr : finrank 𝕜 (LinearMap.range f) ≤ Fintype.card ι := by
    simpa [Module.finrank_fintype_fun_eq_card] using Submodule.finrank_le (LinearMap.range f)
  rw [Module.finrank_fintype_fun_eq_card] at h
  omega

/-- **Count from below.** If the form of `A - t` is negative on every nonzero vector of `V`, then
at least `dim V` eigenvalues lie strictly below `t`. -/
theorem le_eigCount_of_forall_qform_lt {t : ℝ} (V : Submodule 𝕜 (n → 𝕜))
    (hV : ∀ x ∈ V, x ≠ 0 → qform A x < t * nsq x) : finrank 𝕜 V ≤ eigCount hA t := by
  set s := Finset.univ.filter fun i ↦ ¬hA.eigenvalues i < t
  have hs : s.card + eigCount hA t = Fintype.card n := by
    rw [add_comm]
    exact Finset.card_filter_add_card_filter_not _
  by_contra hlt
  have hcard : Fintype.card n < finrank 𝕜 (eigSpan hA s) + finrank 𝕜 V := by
    rw [finrank_eigSpan]
    omega
  obtain ⟨x, hxE, hxV, hx0⟩ := exists_ne_zero_mem_inf _ _ hcard
  exact (hV x hxV hx0).not_ge
    (le_qform_of_mem_eigSpan hA (fun i hi ↦ not_lt.1 (Finset.mem_filter.1 hi).2) hxE)

/-- **Rayleigh–Ritz with a Gram matrix.** For any matrix `Z` (columns need not be orthonormal),
if `Re ⟪c, Zᴴ A Z c⟫ < t Re ⟪c, Zᴴ Z c⟫` for every nonzero `c ∈ P`, then at least `dim P`
eigenvalues of `A` lie strictly below `t`. -/
theorem le_eigCount_of_compression {k : Type*} [Fintype k] {t : ℝ} (Z : Matrix n k 𝕜)
    (P : Submodule 𝕜 (k → 𝕜))
    (hP : ∀ c ∈ P, c ≠ 0 → qform (Zᴴ * A * Z) c < t * qform (Zᴴ * Z) c) :
    finrank 𝕜 P ≤ eigCount hA t := by
  have hinj : Set.InjOn Z.mulVecLin P := by
    intro c hc d hd hcd
    by_contra hne
    have hcd' : Z *ᵥ (c - d) = 0 := by
      simp only [Matrix.mulVecLin_apply] at hcd
      rw [mulVec_sub, hcd, sub_self]
    have h := hP (c - d) (P.sub_mem hc hd) (sub_ne_zero.2 hne)
    rw [← qform_mulVec, ← nsq_mulVec, hcd', qform_zero_right, nsq_zero, mul_zero] at h
    exact lt_irrefl _ h
  have hrank : finrank 𝕜 (P.map Z.mulVecLin) = finrank 𝕜 P := by
    have hinj' : Function.Injective (Z.mulVecLin.domRestrict P) := by
      intro c d hcd
      exact Subtype.ext (hinj c.2 d.2 hcd)
    rw [← LinearMap.range_domRestrict, LinearMap.finrank_range_of_inj hinj']
  rw [← hrank]
  refine le_eigCount_of_forall_qform_lt hA _ fun x hx hx0 ↦ ?_
  obtain ⟨c, hc, rfl⟩ := Submodule.mem_map.1 hx
  simp only [Matrix.mulVecLin_apply] at hx0 ⊢
  have hc0 : c ≠ 0 := by rintro rfl; exact hx0 (mulVec_zero _)
  rw [qform_mulVec, nsq_mulVec]
  exact hP c hc hc0

/-- **Cauchy interlacing** in counting form: for `Z` with orthonormal columns (`Zᴴ Z = 1`) the
compression `Zᴴ A Z` has at most as many eigenvalues below `t` as `A`. -/
theorem eigCount_compression_le {k : Type*} [Fintype k] [DecidableEq k] (Z : Matrix n k 𝕜)
    (hZ : Zᴴ * Z = 1) (hM : (Zᴴ * A * Z).IsHermitian) (t : ℝ) :
    eigCount hM t ≤ eigCount hA t := by
  have h := le_eigCount_of_compression hA Z
    (eigSpan hM (Finset.univ.filter fun i ↦ hM.eigenvalues i < t)) fun c hc hc0 ↦ by
      rw [hZ, qform_one]
      exact qform_lt_of_mem_eigSpan hM (fun i hi ↦ (Finset.mem_filter.1 hi).2) hc hc0
  rwa [finrank_eigSpan] at h

/-- **Count through a compression.** Let `S` be a subspace and `Z : n × k` a matrix with
`range Z + S` the whole space. If the form of `A - t` is nonnegative on every vector `Z c + u`
with `c ∈ P` and `u ∈ S`, then `N(t) + dim P ≤ k`, that is `N(t) ≤ codim P`. -/
theorem eigCount_le_of_compression {k : Type*} [Fintype k] {t : ℝ} (Z : Matrix n k 𝕜)
    (S : Submodule 𝕜 (n → 𝕜)) (hspan : LinearMap.range Z.mulVecLin ⊔ S = ⊤)
    (P : Submodule 𝕜 (k → 𝕜))
    (hP : ∀ c ∈ P, ∀ u ∈ S, t * nsq (Z *ᵥ c + u) ≤ qform A (Z *ᵥ c + u)) :
    eigCount hA t + finrank 𝕜 P ≤ Fintype.card k := by
  obtain ⟨P', hP'⟩ := Submodule.exists_isCompl P
  have hPP' := Submodule.finrank_add_eq_of_isCompl hP'
  rw [Module.finrank_fintype_fun_eq_card] at hPP'
  set V := P.map Z.mulVecLin ⊔ S
  have htop : V ⊔ P'.map Z.mulVecLin = ⊤ := by
    rw [eq_top_iff, ← hspan, LinearMap.range_eq_map, ← hP'.sup_eq_top, Submodule.map_sup]
    exact sup_le (sup_le (le_sup_of_le_left le_sup_left) le_sup_right)
      (le_sup_of_le_left le_sup_right)
  have hdim : Fintype.card n ≤ finrank 𝕜 V + finrank 𝕜 P' := by
    have h1 := Submodule.finrank_add_le_finrank_add_finrank V (P'.map Z.mulVecLin)
    rw [htop, finrank_top, Module.finrank_fintype_fun_eq_card] at h1
    have h2 := Submodule.finrank_map_le Z.mulVecLin P'
    omega
  have hcount := eigCount_le_of_forall_le_qform hA V hdim fun x hx ↦ by
    obtain ⟨y, hy, u, hu, rfl⟩ := Submodule.mem_sup.1 hx
    obtain ⟨c, hc, rfl⟩ := Submodule.mem_map.1 hy
    exact hP c hc u hu
  omega

/-! ### Weyl monotonicity -/

/-- **Weyl monotonicity** in counting form: if `Re ⟪x, B x⟫ ≤ Re ⟪x, A x⟫` for all `x`, then
`A` has at most as many eigenvalues below any `t` as `B`. -/
theorem eigCount_le_eigCount_of_qform_le {B : Matrix n n 𝕜} (hB : B.IsHermitian)
    (h : ∀ x, qform B x ≤ qform A x) (t : ℝ) : eigCount hA t ≤ eigCount hB t := by
  set s := Finset.univ.filter fun i ↦ ¬hB.eigenvalues i < t
  have hs : s.card + eigCount hB t = Fintype.card n := by
    rw [add_comm]
    exact Finset.card_filter_add_card_filter_not _
  refine eigCount_le_of_forall_le_qform hA (eigSpan hB s) (by rw [finrank_eigSpan]; omega)
    fun x hx ↦ ?_
  exact (le_qform_of_mem_eigSpan hB (fun i hi ↦ not_lt.1 (Finset.mem_filter.1 hi).2) hx).trans
    (h x)

omit [DecidableEq n] in
/-- The form of a positive semidefinite matrix is nonnegative. -/
theorem qform_nonneg_of_posSemidef {M : Matrix n n 𝕜} (hM : M.PosSemidef) (x : n → 𝕜) :
    0 ≤ qform M x :=
  (RCLike.nonneg_iff.1 (hM.dotProduct_mulVec_nonneg x)).1

/-- **Weyl monotonicity** for the Loewner order `B ≤ A` (that is, `A - B` positive
semidefinite): `#{i : λᵢ(A) < t} ≤ #{i : λᵢ(B) < t}` for every `t`. -/
theorem eigCount_le_eigCount_of_posSemidef {B : Matrix n n 𝕜} (hB : B.IsHermitian)
    (h : (A - B).PosSemidef) (t : ℝ) : eigCount hA t ≤ eigCount hB t :=
  eigCount_le_eigCount_of_qform_le hA hB (fun x ↦ by
    have := qform_nonneg_of_posSemidef h x
    rw [qform_sub_left] at this
    linarith) t

/-! ### Sorted eigenvalues -/

omit [Fintype n] [DecidableEq n] in
/-- For two antitone (decreasingly sorted) sequences, domination of the counting functions,
`#{j : aⱼ < t} ≤ #{j : bⱼ < t}` for all `t`, is pointwise domination `bₖ ≤ aₖ`. -/
theorem le_of_countBelow_le_of_antitone {N : ℕ} {a b : Fin N → ℝ} (ha : Antitone a)
    (hb : Antitone b) (h : ∀ t, countBelow Finset.univ a t ≤ countBelow Finset.univ b t)
    (k : Fin N) : b k ≤ a k := by
  by_contra hlt
  push Not at hlt
  -- At `t = bₖ`: all `j ≥ k` have `aⱼ ≤ aₖ < t`, while `j ≤ k` have `bⱼ ≥ bₖ = t`.
  have h1 : (Finset.Ici k).card ≤ countBelow Finset.univ a (b k) :=
    Finset.card_le_card fun j hj ↦ Finset.mem_filter.2
      ⟨Finset.mem_univ _, (ha (Finset.mem_Ici.1 hj)).trans_lt hlt⟩
  have h2 : countBelow Finset.univ b (b k) ≤ (Finset.Ioi k).card :=
    Finset.card_le_card fun j hj ↦ Finset.mem_Ioi.2 (lt_of_not_ge fun hjk ↦
      (not_le.2 (Finset.mem_filter.1 hj).2) (hb hjk))
  have h3 := h (b k)
  rw [Fin.card_Ici] at h1
  rw [Fin.card_Ioi] at h2
  have := k.2
  omega

/-- The counting function is the same whether computed from `eigenvalues` or from the sorted
`eigenvalues₀`. -/
theorem eigCount_eq_countBelow_eigenvalues₀ (t : ℝ) :
    eigCount hA t = countBelow Finset.univ hA.eigenvalues₀ t := by
  refine Finset.card_equiv (Fintype.equivOfCardEq (Fintype.card_fin _)).symm fun i ↦ ?_
  rw [Finset.mem_filter, Finset.mem_filter]
  simp [IsHermitian.eigenvalues]

/-- **Weyl monotonicity for sorted eigenvalues.** If `B ≤ A` in the Loewner order then
`λₖ(B) ≤ λₖ(A)` for the decreasingly sorted eigenvalues `eigenvalues₀`. -/
theorem eigenvalues₀_le_of_posSemidef {B : Matrix n n 𝕜} (hB : B.IsHermitian)
    (h : (A - B).PosSemidef) (k : Fin (Fintype.card n)) :
    hB.eigenvalues₀ k ≤ hA.eigenvalues₀ k :=
  le_of_countBelow_le_of_antitone hA.eigenvalues₀_antitone hB.eigenvalues₀_antitone
    (fun t ↦ by
      rw [← eigCount_eq_countBelow_eigenvalues₀, ← eigCount_eq_countBelow_eigenvalues₀]
      exact eigCount_le_eigCount_of_posSemidef hA hB h t) k

/-! ### Residual bound -/

/-- **Residual bound.** For every vector `x` and real `θ`, some eigenvalue satisfies
`(λᵢ - θ)² ‖x‖² ≤ ‖A x - θ x‖²`. -/
theorem exists_eigenvalue_near [Nonempty n] (x : n → 𝕜) (θ : ℝ) :
    ∃ i, (hA.eigenvalues i - θ) ^ 2 * nsq x ≤ nsq (A *ᵥ x - (θ : 𝕜) • x) := by
  set y := (eigU hA)ᴴ *ᵥ x
  have hx : x = eigU hA *ᵥ y := (eigU_mulVec_star hA x).symm
  have hres : nsq (A *ᵥ x - (θ : 𝕜) • x) = qform ((A - (θ : 𝕜) • 1) ^ 2) x := by
    have hherm : (A - (θ : 𝕜) • 1)ᴴ = A - (θ : 𝕜) • 1 := by
      rw [conjTranspose_sub, hA.eq, conjTranspose_smul, conjTranspose_one, RCLike.star_def,
        RCLike.conj_ofReal]
    rw [show A *ᵥ x - (θ : 𝕜) • x = (A - (θ : 𝕜) • 1) *ᵥ x by
      rw [sub_mulVec, smul_mulVec, one_mulVec], nsq_mulVec, hherm, sq]
  have hsq : ((A - (θ : 𝕜) • 1) ^ 2).IsHermitian := by
    have : (A - (θ : 𝕜) • 1).IsHermitian := by
      refine IsHermitian.sub hA ?_
      rw [IsHermitian, conjTranspose_smul, conjTranspose_one, RCLike.star_def, RCLike.conj_ofReal]
    exact this.pow 2
  -- In the eigenbasis of `A`, `(A - θ)²` is diagonal with entries `(λᵢ - θ)²`.
  have hdiag : qform ((A - (θ : 𝕜) • 1) ^ 2) (eigU hA *ᵥ y) =
      ∑ i, (hA.eigenvalues i - θ) ^ 2 * ‖y i‖ ^ 2 := by
    have hconj : (eigU hA)ᴴ * (A - (θ : 𝕜) • 1) ^ 2 * eigU hA =
        diagonal fun i ↦ (((hA.eigenvalues i - θ) ^ 2 : ℝ) : 𝕜) := by
      have h1 : (eigU hA)ᴴ * (A - (θ : 𝕜) • 1) * eigU hA =
          diagonal fun i ↦ ((hA.eigenvalues i - θ : ℝ) : 𝕜) := by
        rw [Matrix.mul_sub, Matrix.sub_mul, eigU_conj, Matrix.mul_smul, Matrix.mul_one,
          Matrix.smul_mul, eigU_star_mul, ← diagonal_one, ← diagonal_smul, diagonal_sub]
        congr 1
        ext i
        simp
      have h2 : (eigU hA)ᴴ * (A - (θ : 𝕜) • 1) ^ 2 * eigU hA =
          ((eigU hA)ᴴ * (A - (θ : 𝕜) • 1) * eigU hA) *
            ((eigU hA)ᴴ * (A - (θ : 𝕜) • 1) * eigU hA) := by
        conv_rhs => rw [Matrix.mul_assoc, ← Matrix.mul_assoc (eigU hA),
          ← Matrix.mul_assoc (eigU hA), eigU_mul_star, Matrix.one_mul]
        simp only [sq, Matrix.mul_assoc]
      rw [h2, h1, diagonal_mul_diagonal]
      congr 1
      ext i
      push_cast
      ring
    rw [qform_mulVec, hconj]
    simp only [qform, mulVec_diagonal, dotProduct, Pi.star_apply, map_sum]
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    rw [mul_left_comm, RCLike.star_def, RCLike.conj_mul, ← RCLike.ofReal_pow,
      ← RCLike.ofReal_mul, RCLike.ofReal_re]
  obtain ⟨i, -, hi⟩ := Finset.exists_min_image Finset.univ
    (fun i ↦ (hA.eigenvalues i - θ) ^ 2) Finset.univ_nonempty
  refine ⟨i, ?_⟩
  rw [hres, hx, hdiag, nsq_eigU, nsq_eq_sum, Finset.mul_sum]
  exact Finset.sum_le_sum fun j _ ↦ mul_le_mul_of_nonneg_right (hi j (Finset.mem_univ j))
    (by positivity)

end CertifiedGaps.Spectral
