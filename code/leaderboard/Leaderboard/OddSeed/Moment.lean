/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import CertifiedGaps.ChainCheck
import CertifiedGaps.Transfer.MomentForm

/-!
# Odd seed lengths: the moment form at all ring lengths

`CertifiedGaps.HasMomentForm` expresses the even partition functions `Z L = ∑ᵢ μᵢ^L` with
`μᵢ ≥ 0` (the moduli of the eigenvalues of a self-adjoint transfer operator). When the transfer
operator is itself positive semidefinite (`X_β = K K* ⪰ 0`, as for the dimerized spin-1 chain,
whose intra-cell bond is also reflection positive with real factors), the eigenvalues are
nonnegative and the same representation `Z L = ∑ᵢ λᵢ^L` holds at **all** lengths `L ≥ n₀`, not
only even ones. This file defines the all-lengths variants

* `HasMomentFormAll`, `SpatialMomentBoundsAll`, `HasSpatialMomentBoundsAll`

together with the moment-inequality lemmas (`HasMomentFormAll.spatialMomentBoundsAll`, the limit
lemma `SpatialMomentBoundsAll.of_tendsto_hasMomentFormAll`) and
`hasMomentFormAll_of_trace_pow`, which builds `HasMomentFormAll` from a positive semidefinite
transfer matrix whose powers give every `Z L` (`Z L = Re tr (T ^ L)` for all `L ≥ 1`).

## Main results

* `HasMomentFormAll.spatialMomentBoundsAll`: the moment inequalities at every length.
* `SpatialMomentBoundsAll.toSpatialMomentBounds`: restriction to even lengths.
* `SpatialMomentBoundsAll.of_tendsto_hasMomentFormAll`: moment inequalities of a limit.
* `hasMomentFormAll_of_trace_pow`: the moment form from a positive semidefinite transfer matrix
  at all lengths.
-/

open Filter Finset Real
open scoped ComplexOrder

namespace CertifiedGaps

/-- The moment form at every length `L ≥ n₀` (not only even ones): there are `μᵢ ≥ 0` (the
nonnegative eigenvalues of a positive semidefinite transfer operator) with
`Z L = ∑ᵢ μᵢ^L` for every `L ≥ n₀`. -/
def HasMomentFormAll (Z : ℕ → ℝ) (n₀ : ℕ) : Prop :=
  ∃ μ : ℕ → ℝ, (∀ i, 0 ≤ μ i) ∧ ∀ L : ℕ, n₀ ≤ L → HasSum (fun i ↦ μ i ^ L) (Z L)

/-- The three moment inequalities at every length `n ≥ n₀` (no evenness hypothesis), the
all-lengths analogue of `SpatialMomentBounds`. -/
structure SpatialMomentBoundsAll (Z : ℕ → ℝ) (n₀ : ℕ) : Prop where
  /-- `Z_{2n} ≤ Z_n²`, that is `P_s(n) ≤ 1`, at every `n ≥ n₀`. -/
  sq_le : ∀ n, n₀ ≤ n → Z (2 * n) ≤ Z n ^ 2
  /-- Lemma 1 for the spatial spectrum, in power sums, at every `n ≥ n₀`. -/
  squaring : ∀ n, n₀ ≤ n →
    2 * (Z (2 * n) ^ 2 - Z (2 * (2 * n))) ≤ (Z n ^ 2 - Z (2 * n)) ^ 2
  /-- `ℓ^p` monotonicity: `Z_q ≤ Z_p^{q/p}` for all `n₀ ≤ p ≤ q`. -/
  le_rpow : ∀ p q, n₀ ≤ p → p ≤ q → Z q ≤ Z p ^ ((q : ℝ) / p)

/-- The all-lengths moment form implies the all-lengths moment inequalities. -/
theorem HasMomentFormAll.spatialMomentBoundsAll {Z : ℕ → ℝ} {n₀ : ℕ} (hn₀ : 0 < n₀)
    (h : HasMomentFormAll Z n₀) : SpatialMomentBoundsAll Z n₀ := by
  obtain ⟨μ, hμ, hsum⟩ := h
  have key : ∀ n, n₀ ≤ n → Z (2 * n) ≤ Z n ^ 2 ∧
      2 * (Z (2 * n) ^ 2 - Z (2 * (2 * n))) ≤ (Z n ^ 2 - Z (2 * n)) ^ 2 := fun n hn ↦ by
    have h₁ := hsum n hn
    have h₂ := hsum (2 * n) (by omega)
    have h₄ := hsum (2 * (2 * n)) (by omega)
    simp only [pow_mul'] at h₂ h₄
    have h₄' : HasSum (fun i ↦ (μ i ^ n) ^ 4) (Z (2 * (2 * n))) := by
      convert h₄ using 1; ext i; ring
    exact sq_sub_sum_sq_of_hasSum (fun i ↦ pow_nonneg (hμ i) n) h₁ h₂ h₄'
  exact ⟨fun n hn ↦ (key n hn).1, fun n hn ↦ (key n hn).2,
    fun p q hp hpq ↦ hasSum_pow_le_rpow hμ (by omega) hpq (hsum p hp) (hsum q (hp.trans hpq))⟩

/-- The all-lengths moment inequalities restrict to the even ones. -/
theorem SpatialMomentBoundsAll.toSpatialMomentBounds {Z : ℕ → ℝ} {n₀ : ℕ}
    (h : SpatialMomentBoundsAll Z n₀) : SpatialMomentBounds Z n₀ :=
  ⟨fun n hn _ ↦ h.sq_le n hn, fun n hn _ ↦ h.squaring n hn,
    fun p q hp hpq _ _ ↦ h.le_rpow p q hp hpq⟩

/-- The all-lengths moment inequalities pass to limits taken at each length `L ≥ n₀`. -/
theorem SpatialMomentBoundsAll.of_tendsto {α : Type*} {l : Filter α} [l.NeBot] {Zs : α → ℕ → ℝ}
    {Z : ℕ → ℝ} {n₀ : ℕ} (hn₀ : 0 < n₀) (hZs : ∀ᶠ a in l, SpatialMomentBoundsAll (Zs a) n₀)
    (hlim : ∀ L, n₀ ≤ L → Tendsto (fun a ↦ Zs a L) l (nhds (Z L))) :
    SpatialMomentBoundsAll Z n₀ where
  sq_le n hn := le_of_tendsto_of_tendsto (hlim _ (by omega)) ((hlim n hn).pow 2)
    (hZs.mono fun _ h ↦ h.sq_le n hn)
  squaring n hn := by
    have h₁ := hlim n hn
    have h₂ := hlim (2 * n) (by omega)
    have h₄ := hlim (2 * (2 * n)) (by omega)
    exact le_of_tendsto_of_tendsto (((h₂.pow 2).sub h₄).const_mul 2) (((h₁.pow 2).sub h₂).pow 2)
      (hZs.mono fun _ h ↦ h.squaring n hn)
  le_rpow p q hp hpq := by
    have hp0 : (0 : ℝ) < p := by exact_mod_cast hn₀.trans_le hp
    exact le_of_tendsto_of_tendsto (hlim q (hp.trans hpq))
      (Filter.Tendsto.rpow_const (hlim p hp) (Or.inr (by positivity : (0 : ℝ) ≤ (q : ℝ) / p)))
      (hZs.mono fun _ h ↦ h.le_rpow p q hp hpq)

/-- The all-lengths moment inequalities above `n₀` imply them above any `n₁ ≥ n₀`. -/
theorem SpatialMomentBoundsAll.mono {Z : ℕ → ℝ} {n₀ n₁ : ℕ} (h : SpatialMomentBoundsAll Z n₀)
    (hn : n₀ ≤ n₁) : SpatialMomentBoundsAll Z n₁ :=
  ⟨fun n hn' ↦ h.sq_le n (hn.trans hn'), fun n hn' ↦ h.squaring n (hn.trans hn'),
    fun p q hp hpq ↦ h.le_rpow p q (hn.trans hp) hpq⟩

/-- The route for approximants: if each `Zs a` has the all-lengths moment form and
`Zs a L → Z L` at every `L ≥ n₀`, then `Z` satisfies the all-lengths moment inequalities. -/
theorem SpatialMomentBoundsAll.of_tendsto_hasMomentFormAll {α : Type*} {l : Filter α} [l.NeBot]
    {Zs : α → ℕ → ℝ} {Z : ℕ → ℝ} {n₀ : ℕ} (hn₀ : 0 < n₀)
    (hZs : ∀ᶠ a in l, HasMomentFormAll (Zs a) n₀)
    (hlim : ∀ L, n₀ ≤ L → Tendsto (fun a ↦ Zs a L) l (nhds (Z L))) :
    SpatialMomentBoundsAll Z n₀ :=
  of_tendsto hn₀ (hZs.mono fun _ h ↦ h.spatialMomentBoundsAll hn₀) hlim

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

/-- The all-lengths moment inequalities (`SpatialMomentBoundsAll`) for the partition functions at
every `β > 0`. This is the spatial hypothesis of the doubling argument when the transfer operator
is positive semidefinite, which allows odd seed lengths. -/
def HasSpatialMomentBoundsAll (E : (L : ℕ) → ι L → ℝ) (n₀ : ℕ) : Prop :=
  ∀ β : ℝ, 0 < β → SpatialMomentBoundsAll (fun L ↦ partitionFunction E L β) n₀

/-- The all-lengths moment form of the spectrum at every `β > 0`. -/
def HasSpatialMomentFormAll (E : (L : ℕ) → ι L → ℝ) (n₀ : ℕ) : Prop :=
  ∀ β : ℝ, 0 < β → HasMomentFormAll (fun L ↦ partitionFunction E L β) n₀

omit [∀ L, Nonempty (ι L)] in
/-- The all-lengths moment form implies the all-lengths moment inequalities at every `β > 0`. -/
theorem HasSpatialMomentFormAll.hasSpatialMomentBoundsAll {E : (L : ℕ) → ι L → ℝ} {n₀ : ℕ}
    (hn₀ : 0 < n₀) (h : HasSpatialMomentFormAll E n₀) : HasSpatialMomentBoundsAll E n₀ :=
  fun β hβ ↦ (h β hβ).spatialMomentBoundsAll hn₀

/-- If `Z L = Re tr (T ^ L)` for all `L ≥ 1` with `T` positive semidefinite, then `Z` has the
all-lengths moment form from `L = 1` on, with weights the (nonnegative) eigenvalues of `T`. -/
theorem hasMomentFormAll_of_trace_pow {D : Type*} [Fintype D] [DecidableEq D]
    {T : Matrix D D ℂ} (hT : T.PosSemidef) {Z : ℕ → ℝ}
    (hZ : ∀ L : ℕ, 1 ≤ L → Z L = (Matrix.trace (T ^ L)).re) :
    HasMomentFormAll Z 1 := by
  set e := Fintype.equivFin D
  let μ : ℕ → ℝ := fun i ↦
    if h : i < Fintype.card D then hT.1.eigenvalues (e.symm ⟨i, h⟩) else 0
  refine ⟨μ, fun i ↦ ?_, fun L hL ↦ ?_⟩
  · simp only [μ]; split_ifs with h
    · exact hT.eigenvalues_nonneg _
    · exact le_refl 0
  · have hsum := hasSum_sum_of_ne_finset_zero (f := fun i ↦ μ i ^ L)
      (L := SummationFilter.unconditional ℕ) (s := Finset.range (Fintype.card D)) fun i hi ↦ by
        simp only [Finset.mem_range, not_lt] at hi
        simp [μ, show ¬ i < Fintype.card D by omega, show L ≠ 0 by omega]
    convert hsum using 1
    rw [hZ L hL, trace_pow_eq_sum_eigenvalues hT.1, Complex.ofReal_re, Finset.sum_range,
      ← e.symm.sum_comp]
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    simp only [μ, i.2, dite_true]

end CertifiedGaps
