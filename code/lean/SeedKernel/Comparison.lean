/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Spectral

/-!
# Weighted comparison with the diagonal

A lower bound on the quadratic form of a Hermitian matrix `A` on the vectors that vanish on a set
`B` of basis vectors, for counting by comparison (`SeedKernel.Ising.card_eigenvalues_lt_le`):

* `le_re_quadratic_of_weights`: if positive weights `w` off `B` satisfy
  `∑_{τ ∉ B, τ ≠ σ} |A σ τ| w(τ) ≤ w(σ) (Re A σ σ - ρ)` for every `σ ∉ B`, then
  `Re ⟪x, A x⟫ ≥ ρ ‖x‖²` for every `x` vanishing on `B`.

This is the Collatz–Wielandt bound for the compression of `A` to the complement of `B`. With
`w = 1` it is the Gershgorin bound; a weight vector close to the Perron vector of the compressed
off-diagonal part makes it nearly sharp. It applies to any Hamiltonian written in a product
basis, and the weight condition is a finite check on integers once `A` is scaled to an integer
matrix.
-/

open Matrix Complex

namespace SeedKernel

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- `2 s t ‖a‖ ‖b‖ ≤ t² ‖a‖² + s² ‖b‖²`, from `(t ‖a‖ - s ‖b‖)² ≥ 0`. -/
theorem two_mul_norm_mul_le (a b s t : ℝ) :
    2 * s * t * (a * b) ≤ t ^ 2 * a ^ 2 + s ^ 2 * b ^ 2 := by
  nlinarith [sq_nonneg (t * a - s * b)]

/-- One off-diagonal term: `Re(x̄ a y) ≥ -(‖a‖ t/s ‖x‖² + ‖a‖ s/t ‖y‖²) / 2` for `s, t > 0`. -/
theorem neg_le_re_star_mul_mul (x a y : ℂ) {s t : ℝ} (hs : 0 < s) (ht : 0 < t) :
    -((‖a‖ * t / s) * normSq x + (‖a‖ * s / t) * normSq y) / 2 ≤ (star x * a * y).re := by
  have h1 : -‖star x * a * y‖ ≤ (star x * a * y).re := by
    have := Complex.abs_re_le_norm (star x * a * y)
    linarith [neg_abs_le (star x * a * y).re]
  have h2 : ‖star x * a * y‖ = ‖a‖ * (‖x‖ * ‖y‖) := by
    rw [norm_mul, norm_mul, norm_star]; ring
  have h3 := two_mul_norm_mul_le ‖x‖ ‖y‖ s t
  rw [← Complex.normSq_eq_norm_sq] at h3
  rw [← Complex.normSq_eq_norm_sq] at h3
  have ha := norm_nonneg a
  have h4 : ‖a‖ * (‖x‖ * ‖y‖) ≤ ((‖a‖ * t / s) * normSq x + (‖a‖ * s / t) * normSq y) / 2 := by
    rw [div_mul_eq_mul_div, div_mul_eq_mul_div, div_add_div _ _ hs.ne' ht.ne', div_div,
      le_div_iff₀ (by positivity)]
    nlinarith [mul_le_mul_of_nonneg_left h3 ha]
  linarith

/-- **Weighted comparison with the diagonal.** Let `A` be Hermitian, `x` vanish on `B`, and let
`w` be positive off `B` with `∑_{τ ∉ B, τ ≠ σ} ‖A σ τ‖ w(τ) ≤ w(σ) (Re A σ σ - ρ)` for `σ ∉ B`.
Then `Re ⟪x, A x⟫ ≥ ρ ‖x‖²`. -/
theorem le_re_quadratic_of_weights {A : Matrix n n ℂ} (hA : A.IsHermitian) (B : Finset n)
    (w : n → ℝ) (ρ : ℝ)
    (hw : ∀ σ, σ ∉ B → 0 < w σ ∧
      ∑ τ ∈ Finset.univ.filter (fun τ ↦ τ ∉ B ∧ τ ≠ σ), ‖A σ τ‖ * w τ ≤
        w σ * ((A σ σ).re - ρ))
    (x : n → ℂ) (hx : ∀ σ ∈ B, x σ = 0) :
    ρ * normSqSum x ≤ (star x ⬝ᵥ (A *ᵥ x)).re := by
  classical
  -- `c σ τ = ‖A σ τ‖ w τ / w σ` for distinct `σ, τ ∉ B`, and `0` otherwise
  set c : n → n → ℝ := fun σ τ ↦
    if σ ∈ B ∨ τ ∈ B ∨ σ = τ then 0 else ‖A σ τ‖ * w τ / w σ with hc
  have hsymm (σ τ : n) : ‖A τ σ‖ = ‖A σ τ‖ := by
    rw [← hA.apply σ τ, norm_star]
  -- every term is bounded below by its diagonal part minus the weighted off-diagonal parts
  have hterm (σ τ : n) :
      (if σ = τ then (A σ σ).re * normSq (x σ) else 0) -
          (c σ τ * normSq (x σ) + c τ σ * normSq (x τ)) / 2 ≤
        (star (x σ) * A σ τ * x τ).re := by
    by_cases hστ : σ = τ
    · subst hστ
      have : (star (x σ) * A σ σ * x σ).re = (A σ σ).re * normSq (x σ) := by
        rw [mul_comm (star (x σ)), mul_assoc, star_def, ← Complex.normSq_eq_conj_mul_self]
        simp
      rw [this]
      simp [hc]
    by_cases hσ : σ ∈ B
    · simp [hc, hσ, hx σ hσ, hστ]
    by_cases hτ : τ ∈ B
    · simp [hc, hτ, hx τ hτ, hστ]
    have hs := (hw σ hσ).1
    have ht := (hw τ hτ).1
    have := neg_le_re_star_mul_mul (x σ) (A σ τ) (x τ) hs ht
    simp only [hc, hσ, hτ, hστ, Ne.symm hστ, or_self, ite_false, hsymm σ τ, zero_sub]
    convert this using 2
    ring
  have hexp : (star x ⬝ᵥ (A *ᵥ x)).re = ∑ σ, ∑ τ, (star (x σ) * A σ τ * x τ).re := by
    simp only [dotProduct, mulVec, Pi.star_apply, Finset.mul_sum, Complex.re_sum, mul_assoc]
  have hsum := Finset.sum_le_sum fun σ (_ : σ ∈ Finset.univ) ↦
    Finset.sum_le_sum fun τ (_ : τ ∈ Finset.univ) ↦ hterm σ τ
  simp only [Finset.sum_sub_distrib, Finset.sum_ite_eq, Finset.mem_univ, ite_true,
    ← Finset.sum_div, Finset.sum_add_distrib] at hsum
  rw [Finset.sum_comm (f := fun σ τ ↦ c τ σ * normSq (x τ))] at hsum
  rw [hexp]
  refine le_trans ?_ hsum
  have hpt (σ : n) : ρ * normSq (x σ) ≤
      (A σ σ).re * normSq (x σ) - (∑ τ, c σ τ) * normSq (x σ) := by
    by_cases hσB : σ ∈ B
    · simp [hx σ hσB]
    obtain ⟨hwσ, hle⟩ := hw σ hσB
    have hcσ : ∑ τ, c σ τ =
        (∑ τ ∈ Finset.univ.filter (fun τ ↦ τ ∉ B ∧ τ ≠ σ), ‖A σ τ‖ * w τ) / w σ := by
      rw [Finset.sum_div, Finset.sum_filter]
      refine Finset.sum_congr rfl fun τ _ ↦ ?_
      by_cases hτ : τ ∈ B
      · simp [hc, hτ]
      by_cases hτσ : τ = σ
      · subst hτσ; simp [hc]
      · simp [hc, hσB, hτ, hτσ, Ne.symm hτσ]
    have hle' : ∑ τ, c σ τ ≤ (A σ σ).re - ρ := by
      rw [hcσ, div_le_iff₀ hwσ]; linarith
    nlinarith [normSq_nonneg (x σ)]
  calc ρ * normSqSum x = ∑ σ, ρ * normSq (x σ) := by
        rw [normSqSum, Finset.mul_sum]
    _ ≤ ∑ σ, ((A σ σ).re * normSq (x σ) - (∑ τ, c σ τ) * normSq (x σ)) :=
        Finset.sum_le_sum fun σ _ ↦ hpt σ
    _ = _ := by
        rw [Finset.sum_sub_distrib]
        simp only [Finset.sum_mul]
        ring

end SeedKernel
