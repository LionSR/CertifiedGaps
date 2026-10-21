/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Models.IsingLongitudinal
import CertifiedGaps.Spectral.Schur

/-!
# The Ising ring Hamiltonian in the computational basis

`IsingLongitudinal.model g h` is built from `RPRing.bondOp` and `RPRing.siteOp`. This file
computes its action on a vector,

  `(H x)(σ) = E(σ) x(σ) - g ∑ⱼ x(flipₖ σ)`,   `E(σ) = -∑ⱼ zⱼ zⱼ₊₁ - h ∑ⱼ zⱼ`,

with `zⱼ = ±1` for `σ j = 0, 1` (`hamiltonian_mulVec_apply`), and the comparison of its quadratic
form with the classical energies, `Re ⟪x, H x⟫ ≥ ∑_σ (E(σ) - g L) |x(σ)|²` for `g ≥ 0`
(`re_quadratic_ge`). The entries of one-site and diagonal two-site operators come from
`Ring/Hamiltonian.lean` (`siteOp_mulVec_apply`, `bondOp_diag_apply`), the quadratic form and
Cauchy–Schwarz from `Spectral/` (`Spectral.qform`, `Spectral.re_dotProduct_sq_le`).
-/

open Matrix Complex CertifiedGaps CertifiedGaps.Spectral

namespace SeedKernel.Ising

variable {L : ℕ}

/-- The eigenvalue `±1` of `Z` on the basis state `a`. -/
def zR (a : Fin 2) : ℝ := if a = 0 then 1 else -1

/-- The classical energy `E(σ) = -∑ⱼ zⱼ zⱼ₊₁ - h ∑ⱼ zⱼ` of a configuration. -/
noncomputable def diagE (h : ℝ) (σ : Fin L → Fin 2) : ℝ :=
  -(∑ j, zR (σ j) * zR (σ (finRotate L j))) - h * ∑ j, zR (σ j)

/-- The configuration with spin `j` flipped. -/
def flipAt (j : Fin L) (σ : Fin L → Fin 2) : Fin L → Fin 2 := Function.update σ j (Fin.rev (σ j))

theorem flipAt_flipAt (j : Fin L) (σ : Fin L → Fin 2) : flipAt j (flipAt j σ) = σ := by
  ext i
  by_cases hi : i = j
  · subst hi; simp [flipAt]
  · simp [flipAt, hi]

/-- Flipping a spin as a permutation of the configurations. -/
def flipPerm (j : Fin L) : Equiv.Perm (Fin L → Fin 2) :=
  Function.Involutive.toPerm (flipAt j) (flipAt_flipAt j)

/-! ### The Ising Hamiltonian -/

variable (g h : ℝ)

theorem pauliZ_map_apply (a b : Fin 2) :
    (pauliZ.map (↑) : Matrix (Fin 2) (Fin 2) ℂ) a b = if a = b then (zR a : ℂ) else 0 := by
  fin_cases a <;> fin_cases b <;> simp [pauliZ, zR]

theorem model_bond_apply (j : Fin L) (σ τ : Fin L → Fin 2) :
    (IsingLongitudinal.model g h).bond L j σ τ =
      if σ = τ then -((zR (σ j) * zR (σ (finRotate L j)) : ℝ) : ℂ) else 0 := by
  have hD : (IsingLongitudinal.model g h).Dc 0 = (IsingLongitudinal.model g h).Cc 0 :=
    RPRing.Dc_eq_of_refl _ rfl 0
  refine bondOp_diag_apply _ _ _ (fun p ↦ -((zR p.1 * zR p.2 : ℝ) : ℂ)) (fun p q ↦ ?_) σ τ
  have hC : (IsingLongitudinal.model g h).Cc 0 = pauliZ.map (↑) := rfl
  simp only [RPRing.bondMatrix, Finset.univ_unique, Fin.default_eq_zero, Finset.sum_singleton,
    hD, hC, Matrix.neg_apply, kroneckerMap_apply, pauliZ_map_apply]
  obtain ⟨a, b⟩ := p
  obtain ⟨c, e⟩ := q
  by_cases hac : a = c <;> by_cases hbe : b = e <;> simp [hac, hbe]

theorem model_oc_apply (a b : Fin 2) :
    (IsingLongitudinal.model g h).oc a b =
      if b = a then ((-(h * zR a) : ℝ) : ℂ) else ((-g : ℝ) : ℂ) := by
  fin_cases a <;> fin_cases b <;>
    simp [RPRing.oc, IsingLongitudinal.model, pauliX, pauliZ, zR]

theorem model_site_mulVec_apply (j : Fin L) (x : (Fin L → Fin 2) → ℂ) (σ : Fin L → Fin 2) :
    ((IsingLongitudinal.model g h).site L j *ᵥ x) σ =
      ((-(h * zR (σ j)) : ℝ) : ℂ) * x σ - (g : ℂ) * x (flipAt j σ) := by
  rw [RPRing.site, siteOp_mulVec_apply, Fin.sum_univ_two, model_oc_apply, model_oc_apply,
    flipAt]
  have hσ : σ = Function.update σ j (σ j) := by simp
  rw [show x σ = x (Function.update σ j (σ j)) from congrArg x hσ]
  generalize σ j = s
  fin_cases s <;> simp [sub_eq_add_neg, add_comm]

/-- **The Hamiltonian on a vector.** `(H x)(σ) = E(σ) x(σ) - g ∑ⱼ x(flipⱼ σ)`. -/
theorem hamiltonian_mulVec_apply (x : (Fin L → Fin 2) → ℂ) (σ : Fin L → Fin 2) :
    ((IsingLongitudinal.model g h).hamiltonian L *ᵥ x) σ =
      (diagE h σ : ℂ) * x σ - (g : ℂ) * ∑ j, x (flipAt j σ) := by
  rw [RPRing.hamiltonian, add_mulVec, sum_mulVec, sum_mulVec, Pi.add_apply, Finset.sum_apply,
    Finset.sum_apply]
  simp only [model_site_mulVec_apply]
  have hbond (j : Fin L) : ((IsingLongitudinal.model g h).bond L j *ᵥ x) σ =
      -((zR (σ j) * zR (σ (finRotate L j)) : ℝ) : ℂ) * x σ := by
    simp only [mulVec, dotProduct, model_bond_apply, ite_mul, zero_mul,
      Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  simp only [hbond, Finset.sum_sub_distrib, ← Finset.sum_mul, ← Finset.mul_sum, diagE]
  push_cast
  simp only [Finset.sum_neg_distrib, ← Finset.mul_sum]
  ring

/-- **Comparison with the classical energies.** For `g ≥ 0`,
`Re ⟪x, H x⟫ ≥ ∑_σ (E(σ) - g L) |x(σ)|²`. -/
theorem re_quadratic_ge (hg : 0 ≤ g) (x : (Fin L → Fin 2) → ℂ) :
    ∑ σ, (diagE h σ - g * L) * ‖x σ‖ ^ 2 ≤
      qform ((IsingLongitudinal.model g h).hamiltonian L) x := by
  -- Each flip term is `Re ⟪x, x ∘ flipⱼ⟫`; flipping permutes the configurations, so `x ∘ flipⱼ`
  -- has the norm of `x` and Cauchy–Schwarz bounds the term by `‖x‖²`.
  have hflip (j : Fin L) : ∑ σ, (star (x σ) * x (flipAt j σ)).re ≤ nsq x := by
    have hperm : nsq (x ∘ flipAt j) = nsq x := by
      simp only [nsq_eq_sum, Function.comp_apply]
      exact Equiv.sum_comp (flipPerm j) fun σ ↦ ‖x σ‖ ^ 2
    have hcs := re_dotProduct_sq_le x (x ∘ flipAt j)
    rw [hperm, ← sq] at hcs
    have h := le_of_sq_le_sq hcs (nsq_nonneg x)
    simpa [dotProduct, Complex.re_sum] using h
  have hσ (σ : Fin L → Fin 2) :
      (star (x σ) * ((IsingLongitudinal.model g h).hamiltonian L *ᵥ x) σ).re =
        diagE h σ * ‖x σ‖ ^ 2 - g * ∑ j, (star (x σ) * x (flipAt j σ)).re := by
    rw [hamiltonian_mulVec_apply, mul_sub, sub_re, Finset.mul_sum, Finset.mul_sum, re_sum,
      Complex.sq_norm, normSq_apply]
    simp only [mul_re, mul_im, star_def, conj_re, conj_im, ofReal_re, ofReal_im, Finset.mul_sum]
    ring_nf
  have hexp : qform ((IsingLongitudinal.model g h).hamiltonian L) x =
      ∑ σ, diagE h σ * ‖x σ‖ ^ 2 - g * ∑ j, ∑ σ, (star (x σ) * x (flipAt j σ)).re := by
    rw [qform, RCLike.re_to_complex, dotProduct, re_sum]
    simp only [Pi.star_apply, hσ, Finset.sum_sub_distrib, ← Finset.mul_sum]
    rw [Finset.sum_comm]
  rw [hexp]
  have hsum := Finset.sum_le_sum (s := Finset.univ) fun j (_ : j ∈ Finset.univ) ↦ hflip j
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, nsq_eq_sum]
    at hsum
  simp only [sub_mul, Finset.sum_sub_distrib, mul_assoc, ← Finset.mul_sum]
  nlinarith [mul_le_mul_of_nonneg_left hsum hg]

end SeedKernel.Ising
