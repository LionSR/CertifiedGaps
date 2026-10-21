/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Spectral.Counting
import Mathlib.Analysis.SpecialFunctions.Exp

/-!
# Bounds on `∑ᵢ f(λᵢ)` from eigenvalue counts

For an antitone weight `f` (in applications `f x = exp (-b x)` with `b ≥ 0`, so `∑ᵢ f(λᵢ)` is the
partition function `Tr e^{-bH}`), sums of `f` over a finite family are controlled by the counting
function `t ↦ #{i : aᵢ < t}` alone:

* `sum_le_sum_of_countBelow_le` (**comparison**): if every count of `a` is at most the
  corresponding count of `b` (that is, sorted `aₖ ≥ bₖ`, with possibly fewer `a`'s), then
  `∑ f(aᵢ) ≤ ∑ f(bⱼ)` for antitone `f ≥ 0`. Specializations: `sum_le_of_eigCount_le`
  (pointwise lower bounds `ℓ` on the sorted eigenvalues), `sum_le_sum_of_posSemidef`
  (`H₀ ≤ H` implies `Tr f(H) ≤ Tr f(H₀)`), and `sum_compression_le` (Rayleigh–Ritz:
  the Ritz values of an orthonormal compression bound the sum from below).
* `sum_le_binned` (**binned upper bound**, Abel summation): with thresholds
  `t₀ ≤ t₁ ≤ ⋯ ≤ tₘ`, all values `≥ t₀`, and counts `#{aᵢ < t_{k+1}} ≤ N_{k+1}`,
  `∑ f(aᵢ) ≤ ∑_{k<m} N_{k+1} (f(tₖ) - f(t_{k+1})) + |s| f(tₘ)`.
* `binned_le_sum` (**binned lower bound**): with counts `M_k ≤ #{aᵢ < tₖ}` and `f ≥ 0`,
  `∑_{k<m} M_k (f(tₖ) - f(t_{k+1})) ≤ ∑ f(aᵢ)`.
-/

open Matrix
open scoped ComplexOrder

namespace CertifiedGaps.Spectral

variable {ι κ : Type*}

theorem countBelow_mono_set {s s' : Finset ι} (h : s ⊆ s') (a : ι → ℝ) (t : ℝ) :
    countBelow s a t ≤ countBelow s' a t :=
  Finset.card_le_card (Finset.filter_subset_filter _ h)

/-- Removing a value below the cut lowers the count by one. -/
theorem countBelow_erase_add_one [DecidableEq ι] {s : Finset ι} {a : ι → ℝ} {t : ℝ} {i : ι}
    (hi : i ∈ s) (hit : a i < t) : countBelow (s.erase i) a t + 1 = countBelow s a t := by
  rw [countBelow, countBelow, Finset.filter_erase,
    Finset.card_erase_add_one (Finset.mem_filter.2 ⟨hi, hit⟩)]

/-- No value lies below a cut `t ≤ min a`. -/
theorem countBelow_eq_zero {s : Finset ι} {a : ι → ℝ} {t : ℝ} (h : ∀ i ∈ s, t ≤ a i) :
    countBelow s a t = 0 := by
  rw [countBelow, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  exact fun i hi ↦ not_lt.2 (h i hi)

/-- **Sum comparison by counts.** For an antitone nonnegative `f`, if
`#{i ∈ s : aᵢ < t} ≤ #{j ∈ u : bⱼ < t}` for every `t`, then `∑_{i ∈ s} f(aᵢ) ≤ ∑_{j ∈ u} f(bⱼ)`. -/
theorem sum_le_sum_of_countBelow_le [DecidableEq ι] [DecidableEq κ] {f : ℝ → ℝ}
    (hf : Antitone f) (hf0 : ∀ x, 0 ≤ f x) (s : Finset ι) (a : ι → ℝ) (u : Finset κ)
    (b : κ → ℝ) (h : ∀ t, countBelow s a t ≤ countBelow u b t) :
    ∑ i ∈ s, f (a i) ≤ ∑ j ∈ u, f (b j) := by
  induction hN : s.card generalizing s u with
  | zero =>
    rw [Finset.card_eq_zero.1 hN, Finset.sum_empty]
    exact Finset.sum_nonneg fun j _ ↦ hf0 _
  | succ N ih =>
    have hs : s.Nonempty := Finset.card_pos.1 (by omega)
    obtain ⟨i₀, hi₀, hmin⟩ := Finset.exists_min_image s a hs
    -- The smallest `a` forces a `b` at or below it.
    have hu : u.Nonempty := by
      have h1 : 0 < countBelow s a (a i₀ + 1) :=
        Finset.card_pos.2 ⟨i₀, Finset.mem_filter.2 ⟨hi₀, by linarith⟩⟩
      rw [← Finset.card_pos]
      exact (h1.trans_le (h _)).trans_le (Finset.card_filter_le _ _)
    obtain ⟨j₀, hj₀, hminb⟩ := Finset.exists_min_image u b hu
    have hle : b j₀ ≤ a i₀ := by
      by_contra hlt
      push Not at hlt
      have h1 := h (b j₀)
      rw [countBelow_eq_zero hminb] at h1
      have h2 : 0 < countBelow s a (b j₀) :=
        Finset.card_pos.2 ⟨i₀, Finset.mem_filter.2 ⟨hi₀, hlt⟩⟩
      omega
    rw [← Finset.add_sum_erase s _ hi₀, ← Finset.add_sum_erase u _ hj₀]
    refine add_le_add (hf hle) (ih (s.erase i₀) (u.erase j₀) (fun t ↦ ?_) ?_)
    · by_cases ht : a i₀ < t
      · have h1 := countBelow_erase_add_one hi₀ ht
        have h2 := countBelow_erase_add_one hj₀ (hle.trans_lt ht)
        have := h t
        omega
      · rw [countBelow_eq_zero fun i hi ↦ (not_lt.1 ht).trans
          (hmin i (Finset.mem_of_mem_erase hi))]
        exact Nat.zero_le _
    · rw [Finset.card_erase_of_mem hi₀, hN]; rfl

/-! ### Binned bounds -/

/-- Pointwise Abel bound behind `sum_le_binned`. -/
theorem le_binned_pointwise {f : ℝ → ℝ} (hf : Antitone f) {t : ℕ → ℝ} (ht : Monotone t)
    {x : ℝ} (hx : t 0 ≤ x) (m : ℕ) :
    f x ≤ ∑ k ∈ Finset.range m, (if x < t (k + 1) then f (t k) - f (t (k + 1)) else 0) +
      f (t m) := by
  have hΔ : ∀ k, 0 ≤ f (t k) - f (t (k + 1)) := fun k ↦ sub_nonneg.2 (hf (ht k.le_succ))
  induction m with
  | zero => simpa using hf hx
  | succ m ih =>
    rw [Finset.sum_range_succ]
    by_cases hxm : x < t (m + 1)
    · rw [ite_eq_left hxm]; linarith
    · have hsum : 0 ≤ ∑ k ∈ Finset.range m,
          (if x < t (k + 1) then f (t k) - f (t (k + 1)) else 0) :=
        Finset.sum_nonneg fun k _ ↦ by split_ifs <;> simp [hΔ k]
      rw [ite_eq_right hxm]
      have := hf (not_lt.1 hxm)
      linarith

/-- `∑_{i ∈ s} [aᵢ < t] c = #{aᵢ < t} c`. -/
theorem sum_ite_lt_eq (s : Finset ι) (a : ι → ℝ) (t c : ℝ) :
    ∑ i ∈ s, (if a i < t then c else 0) = countBelow s a t * c := by
  rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul, countBelow]

/-- **Binned upper bound.** Let `t₀ ≤ t₁ ≤ ⋯` be thresholds with every `aᵢ ≥ t₀`, and let
`#{i ∈ s : aᵢ < t_{k+1}} ≤ N_{k+1}` for `k < m`. Then for antitone `f`,
`∑ f(aᵢ) ≤ ∑_{k<m} N_{k+1} (f(tₖ) - f(t_{k+1})) + |s| f(tₘ)`. -/
theorem sum_le_binned {f : ℝ → ℝ} (hf : Antitone f) (s : Finset ι) (a : ι → ℝ) {t : ℕ → ℝ}
    (ht : Monotone t) (m : ℕ) (h0 : ∀ i ∈ s, t 0 ≤ a i) (N : ℕ → ℕ)
    (hN : ∀ k < m, countBelow s a (t (k + 1)) ≤ N (k + 1)) :
    ∑ i ∈ s, f (a i) ≤
      ∑ k ∈ Finset.range m, (N (k + 1) : ℝ) * (f (t k) - f (t (k + 1))) + s.card * f (t m) := by
  have hΔ : ∀ k, 0 ≤ f (t k) - f (t (k + 1)) := fun k ↦ sub_nonneg.2 (hf (ht k.le_succ))
  calc ∑ i ∈ s, f (a i)
      ≤ ∑ i ∈ s, (∑ k ∈ Finset.range m,
          (if a i < t (k + 1) then f (t k) - f (t (k + 1)) else 0) + f (t m)) :=
        Finset.sum_le_sum fun i hi ↦ le_binned_pointwise hf ht (h0 i hi) m
    _ = ∑ k ∈ Finset.range m, (countBelow s a (t (k + 1)) : ℝ) * (f (t k) - f (t (k + 1))) +
          s.card * f (t m) := by
        rw [Finset.sum_add_distrib, Finset.sum_comm, Finset.sum_const, nsmul_eq_mul]
        simp only [sum_ite_lt_eq]
    _ ≤ _ := by
        gcongr with k hk
        · exact hΔ k
        · exact_mod_cast hN k (Finset.mem_range.1 hk)

/-- Pointwise bound behind `binned_le_sum`: `∑_{k<m} [x < tₖ] (f(tₖ) - f(t_{k+1})) ≤ f x`. -/
theorem binned_le_pointwise {f : ℝ → ℝ} (hf : Antitone f) (hf0 : ∀ x, 0 ≤ f x) {t : ℕ → ℝ}
    (ht : Monotone t) (x : ℝ) (m : ℕ) :
    ∑ k ∈ Finset.range m, (if x < t k then f (t k) - f (t (k + 1)) else 0) +
      (if x < t m then f (t m) else 0) ≤ f x := by
  induction m with
  | zero =>
    simp only [Finset.range_zero, Finset.sum_empty, zero_add]
    split_ifs with h
    · exact hf h.le
    · exact hf0 x
  | succ m ih =>
    rw [Finset.sum_range_succ]
    by_cases hxm : x < t m
    · have hxm' : x < t (m + 1) := hxm.trans_le (ht m.le_succ)
      rw [ite_eq_left hxm, ite_eq_left hxm'] at *
      linarith
    · -- Past `tₘ` every earlier indicator vanishes, since `tₖ ≤ tₘ ≤ x`.
      have hzero : ∑ k ∈ Finset.range m,
          (if x < t k then f (t k) - f (t (k + 1)) else 0) = 0 :=
        Finset.sum_eq_zero fun k hk ↦ ite_eq_right fun hk' ↦ hxm (hk'.trans_le
          (ht (Finset.mem_range.1 hk).le))
      rw [hzero, ite_eq_right hxm, zero_add, zero_add]
      split_ifs with h
      · exact hf h.le
      · exact hf0 x

/-- **Binned lower bound.** For antitone `f ≥ 0`, thresholds `t₀ ≤ t₁ ≤ ⋯` and counts
`M_k ≤ #{i ∈ s : aᵢ < tₖ}` for `k < m`, `∑_{k<m} M_k (f(tₖ) - f(t_{k+1})) ≤ ∑ f(aᵢ)`. -/
theorem binned_le_sum {f : ℝ → ℝ} (hf : Antitone f) (hf0 : ∀ x, 0 ≤ f x) (s : Finset ι)
    (a : ι → ℝ) {t : ℕ → ℝ} (ht : Monotone t) (m : ℕ) (M : ℕ → ℕ)
    (hM : ∀ k < m, M k ≤ countBelow s a (t k)) :
    ∑ k ∈ Finset.range m, (M k : ℝ) * (f (t k) - f (t (k + 1))) ≤ ∑ i ∈ s, f (a i) := by
  have hΔ : ∀ k, 0 ≤ f (t k) - f (t (k + 1)) := fun k ↦ sub_nonneg.2 (hf (ht k.le_succ))
  calc ∑ k ∈ Finset.range m, (M k : ℝ) * (f (t k) - f (t (k + 1)))
      ≤ ∑ k ∈ Finset.range m, (countBelow s a (t k) : ℝ) * (f (t k) - f (t (k + 1))) := by
        gcongr with k hk
        · exact hΔ k
        · exact_mod_cast hM k (Finset.mem_range.1 hk)
    _ = ∑ i ∈ s, ∑ k ∈ Finset.range m,
          (if a i < t k then f (t k) - f (t (k + 1)) else 0) := by
        rw [Finset.sum_comm]
        simp only [sum_ite_lt_eq]
    _ ≤ ∑ i ∈ s, f (a i) := Finset.sum_le_sum fun i _ ↦ by
        have h := binned_le_pointwise hf hf0 ht (a i) m
        have : 0 ≤ (if a i < t m then f (t m) else 0) := by
          split_ifs
          · exact hf0 _
          · exact le_rfl
        linarith

/-! ### The Boltzmann weight -/

/-- The Boltzmann weight `x ↦ exp (-b x)` is antitone for `b ≥ 0`. -/
theorem antitone_exp_neg_mul {b : ℝ} (hb : 0 ≤ b) : Antitone fun x : ℝ ↦ Real.exp (-(b * x)) :=
  fun _ _ hxy ↦ Real.exp_le_exp.2 (neg_le_neg (mul_le_mul_of_nonneg_left hxy hb))

theorem exp_neg_mul_nonneg (b x : ℝ) : 0 ≤ Real.exp (-(b * x)) := (Real.exp_pos _).le

/-! ### Eigenvalue sums -/

variable {𝕜 : Type*} [RCLike 𝕜] {n : Type*} [Fintype n] [DecidableEq n]
  {A : Matrix n n 𝕜} (hA : A.IsHermitian)

/-- **Trace bound from lower bounds on sorted eigenvalues.** If the eigenvalue counts of `A` are
dominated by those of a family `ℓ` (equivalently, the sorted eigenvalues satisfy `λₖ ≥ ℓₖ`),
then `∑ᵢ f(λᵢ) ≤ ∑ⱼ f(ℓⱼ)` for antitone `f ≥ 0`. -/
theorem sum_le_of_eigCount_le [DecidableEq κ] {f : ℝ → ℝ} (hf : Antitone f)
    (hf0 : ∀ x, 0 ≤ f x) (u : Finset κ) (ℓ : κ → ℝ) (h : ∀ t, eigCount hA t ≤ countBelow u ℓ t) :
    ∑ i, f (hA.eigenvalues i) ≤ ∑ j ∈ u, f (ℓ j) :=
  sum_le_sum_of_countBelow_le hf hf0 _ _ u ℓ h

/-- **Trace comparison.** If `B ≤ A` in the Loewner order then `∑ f(λᵢ(A)) ≤ ∑ f(λᵢ(B))` for
antitone `f ≥ 0`; with `f x = exp (-b x)`, `Tr e^{-bA} ≤ Tr e^{-bB}`. -/
theorem sum_le_sum_of_posSemidef {f : ℝ → ℝ} (hf : Antitone f) (hf0 : ∀ x, 0 ≤ f x)
    {B : Matrix n n 𝕜} (hB : B.IsHermitian) (h : (A - B).PosSemidef) :
    ∑ i, f (hA.eigenvalues i) ≤ ∑ i, f (hB.eigenvalues i) :=
  sum_le_sum_of_countBelow_le hf hf0 _ _ _ _ fun t ↦
    eigCount_le_eigCount_of_posSemidef hA hB h t

/-- **Rayleigh–Ritz lower bound.** For `Z` with orthonormal columns, the Ritz values (the
eigenvalues of `Zᴴ A Z`) satisfy `∑ⱼ f(θⱼ) ≤ ∑ᵢ f(λᵢ(A))` for antitone `f ≥ 0`; with
`f x = exp (-b x)` this is the variational lower bound on the partition function. -/
theorem sum_compression_le {k : Type*} [Fintype k] [DecidableEq k] {f : ℝ → ℝ}
    (hf : Antitone f) (hf0 : ∀ x, 0 ≤ f x) (Z : Matrix n k 𝕜) (hZ : Zᴴ * Z = 1)
    (hM : (Zᴴ * A * Z).IsHermitian) :
    ∑ j, f (hM.eigenvalues j) ≤ ∑ i, f (hA.eigenvalues i) :=
  sum_le_sum_of_countBelow_le hf hf0 _ _ _ _ fun t ↦ eigCount_compression_le hA Z hZ hM t

end CertifiedGaps.Spectral
