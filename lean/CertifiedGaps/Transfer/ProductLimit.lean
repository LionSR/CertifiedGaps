/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib

/-!
# A first-order Lie–Trotter (Chernoff) product formula

We prove that for a list of complex square matrices `h₁, …, hₖ` and a real number `c`,
$$\Bigl(\prod_i \bigl(1 + \tfrac{c}{M} h_i\bigr)\Bigr)^M
  \longrightarrow \exp\Bigl(c \sum_i h_i\Bigr) \qquad (M \to \infty).$$

The proof goes through a general Chernoff-type lemma in a complete normed `ℂ`-algebra: if
`g : ℂ → 𝔸` satisfies `g 0 = 1` and `HasDerivAt g A 0`, then `g (1 / M) ^ M → exp A`.
It is obtained by comparing `g (1 / M)` with `exp ((1 / M) • A)`, whose `M`-th power is exactly
`exp A`, together with the telescoping estimate
`‖X ^ k - Y ^ k‖ ≤ k * K ^ k * ‖X - Y‖` for `‖X‖, ‖Y‖ ≤ K` and `1 ≤ K`.

## Main results

* `CertifiedGaps.norm_pow_sub_pow_le`: the telescoping estimate for powers.
* `CertifiedGaps.tendsto_pow_of_hasDerivAt`: the Chernoff lemma in a complete normed algebra.
* `CertifiedGaps.tendsto_list_prod_one_add_smul_pow`: the product formula for complex matrices
  (with the default product topology on `Matrix n n ℂ`).
-/

open Filter Topology

namespace CertifiedGaps

section General

variable {𝔸 : Type*} [NormedRing 𝔸]

/-- Telescoping estimate: if `‖X‖ ≤ K`, `‖Y‖ ≤ K` and `1 ≤ K`, then
`‖X ^ k - Y ^ k‖ ≤ k * K ^ k * ‖X - Y‖`. -/
theorem norm_pow_sub_pow_le [NormOneClass 𝔸] {X Y : 𝔸} {K : ℝ} (hK : 1 ≤ K) (hX : ‖X‖ ≤ K)
    (hY : ‖Y‖ ≤ K) (k : ℕ) : ‖X ^ k - Y ^ k‖ ≤ k * K ^ k * ‖X - Y‖ := by
  induction k with
  | zero => simp
  | succ k ih =>
    have hK0 : 0 ≤ K := zero_le_one.trans hK
    have hXk : ‖X ^ k‖ ≤ K ^ k :=
      (norm_pow_le X k).trans (pow_le_pow_left₀ (norm_nonneg _) hX k)
    have hpow : K ^ k ≤ K ^ (k + 1) := pow_le_pow_right₀ hK (Nat.le_succ k)
    have hd := norm_nonneg (X - Y)
    have heq : X ^ (k + 1) - Y ^ (k + 1) = X ^ k * (X - Y) + (X ^ k - Y ^ k) * Y := by
      rw [pow_succ, pow_succ]; noncomm_ring
    rw [heq]
    calc ‖X ^ k * (X - Y) + (X ^ k - Y ^ k) * Y‖
        ≤ ‖X ^ k‖ * ‖X - Y‖ + ‖X ^ k - Y ^ k‖ * ‖Y‖ :=
          (norm_add_le _ _).trans (add_le_add (norm_mul_le _ _) (norm_mul_le _ _))
      _ ≤ K ^ k * ‖X - Y‖ + (k * K ^ k * ‖X - Y‖) * K := by
          gcongr
      _ ≤ K ^ (k + 1) * ‖X - Y‖ + k * K ^ (k + 1) * ‖X - Y‖ := by
          rw [pow_succ]
          have : K ^ k * ‖X - Y‖ ≤ K ^ k * K * ‖X - Y‖ := by
            gcongr; rw [← pow_succ]; exact hpow
          nlinarith
      _ = ((k + 1 : ℕ) : ℝ) * K ^ (k + 1) * ‖X - Y‖ := by push_cast; ring

variable [NormedAlgebra ℂ 𝔸]

/-- If `g 0 = 1` and `g` has derivative `A` at `0`, then `M • (g (1 / M) - 1) → A`. -/
theorem tendsto_natCast_smul_sub_one {g : ℂ → 𝔸} {A : 𝔸} (hg0 : g 0 = 1)
    (hg : HasDerivAt g A 0) :
    Tendsto (fun M : ℕ ↦ (M : ℂ) • (g (M : ℂ)⁻¹ - 1)) atTop (𝓝 A) := by
  have hs := hg.tendsto_slope
  have hinv : Tendsto (fun M : ℕ ↦ (M : ℂ)⁻¹) atTop (𝓝[≠] 0) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _
      tendsto_inv_atTop_nhds_zero_nat ?_
    filter_upwards [eventually_ge_atTop 1] with M hM
    simp only [Set.mem_compl_iff, Set.mem_singleton_iff, inv_eq_zero, Nat.cast_eq_zero]
    omega
  refine (hs.comp hinv).congr' ?_
  filter_upwards with M
  simp [slope_def_module, hg0]

variable [CompleteSpace 𝔸] [NormedAlgebra ℚ 𝔸] [NormOneClass 𝔸]

/-- Chernoff lemma: if `X M` satisfies `M • (X M - 1) → A`, then `X M ^ M → exp A`. -/
theorem tendsto_pow_of_tendsto_natCast_smul {X : ℕ → 𝔸} {A : 𝔸}
    (hX : Tendsto (fun M : ℕ ↦ (M : ℂ) • (X M - 1)) atTop (𝓝 A)) :
    Tendsto (fun M : ℕ ↦ X M ^ M) atTop (𝓝 (NormedSpace.exp A)) := by
  set Y : ℕ → 𝔸 := fun M ↦ NormedSpace.exp (((M : ℂ)⁻¹) • A) with hYdef
  have hY : Tendsto (fun M : ℕ ↦ (M : ℂ) • (Y M - 1)) atTop (𝓝 A) := by
    have hd := hasDerivAt_exp_smul_const (𝕂 := ℂ) A 0
    simp only [zero_smul, NormedSpace.exp_zero, one_mul] at hd
    exact tendsto_natCast_smul_sub_one (g := fun u : ℂ ↦ NormedSpace.exp (u • A))
      (by simp) hd
  have hYpow : ∀ M : ℕ, 1 ≤ M → Y M ^ M = NormedSpace.exp A := by
    intro M hM
    rw [hYdef, ← NormedSpace.exp_nsmul, ← Nat.cast_smul_eq_nsmul ℂ, smul_smul,
      mul_inv_cancel₀ (by exact_mod_cast (show M ≠ 0 by omega)), one_smul]
  set C : ℝ := ‖A‖ + 1 with hC
  -- eventual norm bounds
  have hbound : ∀ Z : ℕ → 𝔸, Tendsto (fun M : ℕ ↦ (M : ℂ) • (Z M - 1)) atTop (𝓝 A) →
      ∀ᶠ M in atTop, ‖Z M‖ ≤ 1 + C / M := by
    intro Z hZ
    have h1 : ∀ᶠ M : ℕ in atTop, ‖(M : ℂ) • (Z M - 1)‖ < C :=
      (hZ.norm).eventually (gt_mem_nhds (by rw [hC]; linarith))
    filter_upwards [h1, eventually_ge_atTop 1] with M hM hM1
    have hMpos : (0 : ℝ) < M := by exact_mod_cast hM1
    rw [norm_smul, Complex.norm_natCast] at hM
    have : ‖Z M - 1‖ ≤ C / M := by
      rw [le_div_iff₀ hMpos]; linarith
    calc ‖Z M‖ = ‖1 + (Z M - 1)‖ := by congr 1; abel
      _ ≤ ‖(1 : 𝔸)‖ + ‖Z M - 1‖ := norm_add_le _ _
      _ ≤ 1 + C / M := by rw [norm_one]; linarith
  -- the difference
  have hdiff : Tendsto (fun M : ℕ ↦ (M : ℂ) • (X M - Y M)) atTop (𝓝 0) := by
    have := hX.sub hY
    rw [sub_self] at this
    refine this.congr' (Eventually.of_forall fun M ↦ ?_)
    simp only [← smul_sub]; congr 1; abel
  have hdiffn : Tendsto (fun M : ℕ ↦ Real.exp C * ‖(M : ℂ) • (X M - Y M)‖) atTop (𝓝 0) := by
    simpa using (hdiff.norm).const_mul (Real.exp C)
  have hmain : Tendsto (fun M : ℕ ↦ X M ^ M - Y M ^ M) atTop (𝓝 0) := by
    rw [tendsto_zero_iff_norm_tendsto_zero]
    refine squeeze_zero' (Eventually.of_forall fun _ ↦ norm_nonneg _) ?_ hdiffn
    filter_upwards [hbound X hX, hbound Y hY, eventually_ge_atTop 1] with M hXM hYM hM1
    have hMpos : (0 : ℝ) < M := by exact_mod_cast hM1
    have hC0 : 0 ≤ C := by rw [hC]; positivity
    have hK : 1 ≤ 1 + C / M := by have := div_nonneg hC0 hMpos.le; linarith
    have hKM : (1 + C / M) ^ M ≤ Real.exp C := by
      calc (1 + C / M) ^ M ≤ Real.exp (C / M) ^ M := by
            gcongr
            linarith [Real.add_one_le_exp (C / M)]
        _ = Real.exp C := by
            rw [← Real.exp_nat_mul]; congr 1; field_simp
    calc ‖X M ^ M - Y M ^ M‖ ≤ M * (1 + C / M) ^ M * ‖X M - Y M‖ :=
          norm_pow_sub_pow_le hK hXM hYM M
      _ = (1 + C / M) ^ M * ‖(M : ℂ) • (X M - Y M)‖ := by
          rw [norm_smul, Complex.norm_natCast]; ring
      _ ≤ Real.exp C * ‖(M : ℂ) • (X M - Y M)‖ := by gcongr
  have := hmain.add_const (NormedSpace.exp A)
  rw [zero_add] at this
  refine this.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with M hM
  rw [hYpow M hM, sub_add_cancel]

/-- Chernoff lemma, derivative form: if `g 0 = 1` and `HasDerivAt g A 0`, then
`g (1 / M) ^ M → exp A`. -/
theorem tendsto_pow_of_hasDerivAt {g : ℂ → 𝔸} {A : 𝔸} (hg0 : g 0 = 1)
    (hg : HasDerivAt g A 0) :
    Tendsto (fun M : ℕ ↦ g ((M : ℂ)⁻¹) ^ M) atTop (𝓝 (NormedSpace.exp A)) :=
  tendsto_pow_of_tendsto_natCast_smul (tendsto_natCast_smul_sub_one hg0 hg)

end General

section ListProd

variable {𝔸 : Type*} [NormedRing 𝔸] [NormedAlgebra ℂ 𝔸]

/-- The polynomial `τ ↦ ∏ᵢ (1 + τ • hᵢ)` has derivative `∑ᵢ hᵢ` at `τ = 0`. -/
theorem hasDerivAt_list_prod_one_add_smul (l : List 𝔸) :
    HasDerivAt (fun τ : ℂ ↦ (l.map fun h ↦ 1 + τ • h).prod) l.sum 0 := by
  induction l with
  | nil => simpa using hasDerivAt_const (0 : ℂ) (1 : 𝔸)
  | cons a l ih =>
    have ha : HasDerivAt (fun τ : ℂ ↦ 1 + τ • a) a 0 := by
      simpa using ((hasDerivAt_id (0 : ℂ)).smul_const a).const_add 1
    have := ha.fun_mul ih
    simpa [List.map_cons, List.prod_cons, List.sum_cons] using this

end ListProd

/-- First-order Lie–Trotter product formula for complex matrices:
`(∏ᵢ (1 + (c / M) • hᵢ)) ^ M → exp (c • ∑ᵢ hᵢ)` as `M → ∞`. -/
theorem tendsto_list_prod_one_add_smul_pow {n : Type*} [Fintype n] [DecidableEq n]
    (l : List (Matrix n n ℂ)) (c : ℝ) :
    Filter.Tendsto (fun M : ℕ ↦ ((l.map fun h ↦ 1 + (((c / M : ℝ)) : ℂ) • h).prod) ^ M)
      Filter.atTop (nhds (NormedSpace.exp ((c : ℂ) • l.sum))) := by
  rcases isEmpty_or_nonempty n with hn | hn
  · exact tendsto_const_nhds.congr fun _ ↦ Subsingleton.elim _ _
  let _ : NormedRing (Matrix n n ℂ) := Matrix.linftyOpNormedRing
  let _ : NormedAlgebra ℂ (Matrix n n ℂ) := Matrix.linftyOpNormedAlgebra
  let _ : NormedAlgebra ℚ (Matrix n n ℂ) := Matrix.linftyOpNormedAlgebra
  set f : ℂ → Matrix n n ℂ := fun τ ↦ (l.map fun h ↦ 1 + τ • h).prod with hf
  have := tendsto_pow_of_hasDerivAt (𝔸 := Matrix n n ℂ) (g := fun τ : ℂ ↦ f (τ * c))
    (A := (c : ℂ) • l.sum) (by simp [hf]) (by
      have := (hasDerivAt_list_prod_one_add_smul l).scomp_of_eq (0 : ℂ)
        ((hasDerivAt_id (0 : ℂ)).mul_const (c : ℂ)) (by simp)
      simpa [hf, Function.comp_def] using this)
  refine this.congr fun M ↦ ?_
  simp only [hf]
  congr 4
  ext h
  push_cast
  rw [div_eq_inv_mul]

end CertifiedGaps
