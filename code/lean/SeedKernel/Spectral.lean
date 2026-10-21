/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Spectral.Counting

/-!
# Spectral bounds in the seed checker's squared-norm notation

The seed arithmetic uses `∑ σ, Complex.normSq (x σ)`. These adapters apply the shared
`CertifiedGaps.Spectral` residual and counting bounds in that notation.
-/

open Matrix Complex

namespace SeedKernel

variable {n : Type*} [Fintype n]

/-- The squared Euclidean norm of a complex vector. -/
noncomputable def normSqSum (x : n → ℂ) : ℝ := ∑ σ, normSq (x σ)

theorem normSqSum_eq_re_dotProduct (x : n → ℂ) : normSqSum x = (star x ⬝ᵥ x).re := by
  simp [normSqSum, dotProduct, Complex.re_sum, Complex.normSq_apply, mul_comm]

/-- The seed checker's sum agrees with the shared spectral squared norm. -/
theorem normSqSum_eq_nsq (x : n → ℂ) : normSqSum x = CertifiedGaps.Spectral.nsq x :=
  normSqSum_eq_re_dotProduct x

variable [DecidableEq n] {A : Matrix n n ℂ} (hA : A.IsHermitian)

/-- The shared residual bound in the seed checker's squared-norm notation. -/
theorem exists_eigenvalue_near [Nonempty n] (x : n → ℂ) (θ : ℝ) :
    ∃ i, (hA.eigenvalues i - θ) ^ 2 * normSqSum x ≤
      normSqSum (A *ᵥ x - (θ : ℂ) • x) := by
  simp_rw [normSqSum_eq_nsq]
  exact CertifiedGaps.Spectral.exists_eigenvalue_near hA x θ

/-- The shared counting bound in the seed checker's squared-norm notation. -/
theorem card_eigenvalues_lt_le (B : Finset n) (ρ : ℝ)
    (hB : ∀ x : n → ℂ, (∀ σ ∈ B, x σ = 0) →
      ρ * normSqSum x ≤ (star x ⬝ᵥ (A *ᵥ x)).re) :
    (Finset.univ.filter fun i ↦ hA.eigenvalues i < ρ).card ≤ B.card := by
  refine CertifiedGaps.Spectral.eigCount_le_card_of_vanishing hA B (t := ρ) fun x hx ↦ ?_
  rw [← normSqSum_eq_nsq]
  exact hB x hx

end SeedKernel
