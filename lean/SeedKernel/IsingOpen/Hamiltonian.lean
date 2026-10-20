/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.Model
import SeedKernel.Ising.Hamiltonian
import SeedKernel.Spectral

/-!
# The free Ising chain in the computational basis

The open chain of `IsingLongitudinal.model g h` with free ends (`RPRing.freeBoundary`) on
`k + 1` sites acts on a vector as

  `(H x)(σ) = E(σ) x(σ) - g ∑ⱼ x(flipⱼ σ)`,   `E(σ) = -∑_{j<k} zⱼ zⱼ₊₁ - h ∑ⱼ zⱼ`

(`openHamiltonian_mulVec_apply`). The crude comparison `Re ⟪x, H x⟫ ≥ ∑_σ (E(σ) - g n) |x(σ)|²`
of the ring is too weak on an open chain, whose end spins cost less to flip. This file proves
the weighted comparison instead (`le_re_quadratic_of_weights`): if `x` vanishes on a set `B` and
positive weights `w` off `B` satisfy `g ∑ⱼ w(flipⱼ σ) ≤ w(σ) (E(σ) - ρ)` (the sum over the
flips that leave `B`), then `Re ⟪x, H x⟫ ≥ ρ ‖x‖²`. It is the Collatz–Wielandt bound for the
compression of `H` to the complement of `B`.
-/

open Matrix Complex CertifiedGaps SeedKernel.Ising

namespace SeedKernel.IsingOpen

variable {k : ℕ}

/-- The classical energy `E(σ) = -∑_{j<k} zⱼ zⱼ₊₁ - h ∑ⱼ zⱼ` of a configuration of the open
chain of `k + 1` sites. -/
noncomputable def openE (h : ℝ) (σ : Fin (k + 1) → Fin 2) : ℝ :=
  -(∑ j : Fin k, zR (σ j.castSucc) * zR (σ j.succ)) - h * ∑ j, zR (σ j)

variable (g h : ℝ)

theorem model_bondOp_apply {L : ℕ} (j j' : Fin L) (σ τ : Fin L → Fin 2) :
    bondOp j j' (IsingLongitudinal.model g h).bondMatrix σ τ =
      if σ = τ then -((zR (σ j) * zR (σ j') : ℝ) : ℂ) else 0 := by
  have hD : (IsingLongitudinal.model g h).Dc 0 = (IsingLongitudinal.model g h).Cc 0 :=
    RPRing.Dc_eq_of_refl _ rfl 0
  refine bondOp_diag_apply _ _ _ (fun p ↦ -((zR p.1 * zR p.2 : ℝ) : ℂ)) (fun p q ↦ ?_) σ τ
  have hC : (IsingLongitudinal.model g h).Cc 0 = pauliZ.map (↑) := rfl
  simp only [RPRing.bondMatrix, Finset.univ_unique, Fin.default_eq_zero, Finset.sum_singleton,
    hD, hC, Matrix.neg_apply, kroneckerMap_apply, pauliZ_map_apply]
  obtain ⟨a, b⟩ := p
  obtain ⟨c, e⟩ := q
  by_cases hac : a = c <;> by_cases hbe : b = e <;> simp [hac, hbe]

/-- The free open chain on `k + 1` sites. -/
noncomputable abbrev freeH (k : ℕ) : Matrix (Fin (k + 1) → Fin 2) (Fin (k + 1) → Fin 2) ℂ :=
  (IsingLongitudinal.model g h).openHamiltonian (IsingLongitudinal.model g h).freeBoundary (k + 1)

/-- **The free chain on a vector.** `(H x)(σ) = E(σ) x(σ) - g ∑ⱼ x(flipⱼ σ)`. -/
theorem openHamiltonian_mulVec_apply (x : (Fin (k + 1) → Fin 2) → ℂ) (σ : Fin (k + 1) → Fin 2) :
    (freeH g h k *ᵥ x) σ = (openE h σ : ℂ) * x σ - (g : ℂ) * ∑ j, x (flipAt j σ) := by
  rw [freeH, RPRing.openHamiltonian_free, add_mulVec, sum_mulVec, sum_mulVec, Pi.add_apply,
    Finset.sum_apply, Finset.sum_apply]
  have hsite (i : Fin (k + 1)) :
      (siteOp i (IsingLongitudinal.model g h).oc *ᵥ x) σ =
        ((-(h * zR (σ i)) : ℝ) : ℂ) * x σ - (g : ℂ) * x (flipAt i σ) :=
    model_site_mulVec_apply g h i x σ
  have hbond (j : Fin k) :
      (bondOp j.castSucc j.succ (IsingLongitudinal.model g h).bondMatrix *ᵥ x) σ =
        -((zR (σ j.castSucc) * zR (σ j.succ) : ℝ) : ℂ) * x σ := by
    simp only [mulVec, dotProduct, model_bondOp_apply, ite_mul, zero_mul,
      Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  simp only [hsite, hbond, Finset.sum_sub_distrib, ← Finset.sum_mul, ← Finset.mul_sum, openE]
  push_cast
  simp only [Finset.sum_neg_distrib, ← Finset.mul_sum]
  ring

/-- A spin flip changes the configuration. -/
theorem flipAt_ne (j : Fin (k + 1)) (σ : Fin (k + 1) → Fin 2) : flipAt j σ ≠ σ := by
  intro h
  have := congrFun h j
  simp only [flipAt, Function.update_self] at this
  generalize σ j = a at this
  fin_cases a <;> simp at this

/-- **The matrix entries of the free chain.** `H σ τ = E(σ) [σ = τ] - g #{j : flipⱼ σ = τ}`. -/
theorem freeH_apply (σ τ : Fin (k + 1) → Fin 2) :
    freeH g h k σ τ = (if σ = τ then (openE h σ : ℂ) else 0) -
      (g : ℂ) * ∑ j, (if flipAt j σ = τ then (1 : ℂ) else 0) := by
  have h1 : (freeH g h k *ᵥ Pi.single τ 1) σ = freeH g h k σ τ := by
    simp [mulVec, dotProduct, Pi.single_apply]
  rw [← h1, openHamiltonian_mulVec_apply]
  simp [Pi.single_apply]

/-- **The weight condition of the free chain.** For `g ≥ 0`, positive weights `w` off `B` with
`g ∑ⱼ w(flipⱼ σ) ≤ w(σ) (E(σ) - ρ)` (the sum over the flips that leave `B`) satisfy the
hypothesis of the general weighted comparison `SeedKernel.le_re_quadratic_of_weights`. -/
theorem freeH_weights (hg : 0 ≤ g) (B : Finset (Fin (k + 1) → Fin 2))
    (w : (Fin (k + 1) → Fin 2) → ℝ) (ρ : ℝ)
    (hw : ∀ σ, σ ∉ B → 0 < w σ ∧
      g * ∑ j, (if flipAt j σ ∈ B then 0 else w (flipAt j σ)) ≤ w σ * (openE h σ - ρ))
    (σ : Fin (k + 1) → Fin 2) (hσ : σ ∉ B) :
    0 < w σ ∧ ∑ τ ∈ Finset.univ.filter (fun τ ↦ τ ∉ B ∧ τ ≠ σ), ‖freeH g h k σ τ‖ * w τ ≤
      w σ * ((freeH g h k σ σ).re - ρ) := by
  classical
  obtain ⟨hwσ, hle⟩ := hw σ hσ
  refine ⟨hwσ, ?_⟩
  have hdiag : (freeH g h k σ σ).re = openE h σ := by
    rw [freeH_apply]; simp [flipAt_ne]
  have hoff (τ : Fin (k + 1) → Fin 2) (hτ : τ ≠ σ) :
      ‖freeH g h k σ τ‖ = ∑ j, if flipAt j σ = τ then g else 0 := by
    have hc : freeH g h k σ τ = ((-(∑ j, if flipAt j σ = τ then g else 0) : ℝ) : ℂ) := by
      rw [freeH_apply, ite_eq_right_of_eq_false _ _ (eq_false (Ne.symm hτ)), zero_sub, Finset.mul_sum, Complex.ofReal_neg,
        Complex.ofReal_sum]
      congr 1
      refine Finset.sum_congr rfl fun j _ ↦ ?_
      split_ifs <;> simp
    rw [hc, Complex.norm_real, Real.norm_eq_abs, abs_neg,
      abs_of_nonneg (Finset.sum_nonneg fun j _ ↦ by split_ifs <;> linarith)]
  rw [hdiag]
  calc ∑ τ ∈ Finset.univ.filter (fun τ ↦ τ ∉ B ∧ τ ≠ σ), ‖freeH g h k σ τ‖ * w τ
      = ∑ τ ∈ Finset.univ.filter (fun τ ↦ τ ∉ B ∧ τ ≠ σ),
          ∑ j, if flipAt j σ = τ then g * w τ else 0 := by
        refine Finset.sum_congr rfl fun τ hτ ↦ ?_
        rw [hoff τ (Finset.mem_filter.1 hτ).2.2, Finset.sum_mul]
        exact Finset.sum_congr rfl fun j _ ↦ by split_ifs <;> simp
    _ = g * ∑ j, (if flipAt j σ ∈ B then 0 else w (flipAt j σ)) := by
        rw [Finset.sum_comm, Finset.mul_sum]
        refine Finset.sum_congr rfl fun j _ ↦ ?_
        rw [Finset.sum_ite_eq]
        by_cases hj : flipAt j σ ∈ B <;> simp [hj, flipAt_ne]
    _ ≤ w σ * (openE h σ - ρ) := hle

end SeedKernel.IsingOpen
