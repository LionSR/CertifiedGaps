/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Local.SpinOne

/-!
# Rotations of the spin-one matrices

In the Cartesian basis the spin-one site space is the vector representation of `SO(3)`: a
rotation `Q` acts on it as the real orthogonal matrix `Q` itself, and the spin matrices transform
as a vector,

  `Qᵀ S^α Q = ∑_β Q_{αβ} S^β`.

Entrywise this is `∑_{b'c'} ε_{α b' c'} Q_{b'b} Q_{c'c} = ∑_β Q_{αβ} ε_{βbc}`: `Q` preserves the
cross product (`PreservesCross`), which holds for every proper rotation. A model checks it for a
concrete `Q` by finite computation.

A two-site exchange twisted by a matrix `R`, `rotExchange R = ∑_{αβ} R_{αβ} S^α ⊗ S^β`
(`R = 1` is `S · S`), therefore transforms by `R ↦ Qᵀ R Q` under `Q ⊗ Q`
(`conj_rotExchange`). For a rotation `g`, the bond `∑_α S^α ⊗ g S^α g⁻¹` whose second site is
conjugated by `g` is `rotExchange gᵀ`; for an axis permutation `π`, the bond `∑_α S^α ⊗ S^{π α}`
is `rotExchange` of the permutation matrix `[β = π α]` (`rotExchange_permRot`).

## Main results

* `map_G_kronecker_G`: `G_α ⊗ G_β = -(S^α ⊗ S^β)`.
* `conj_S_of_preservesCross`: `Qᵀ S^α Q = ∑_β Q_{αβ} S^β`.
* `conj_rotExchange`: `(Q ⊗ Q)ᵀ (rotExchange R) (Q ⊗ Q) = rotExchange (Qᵀ R Q)`.
* `rotExchange_one`, `rotExchange_permRot`.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps.SpinOne

/-- `G_α ⊗ G_β = -(S^α ⊗ S^β)`, since `S^α = -i G_α`. -/
theorem map_G_kronecker_G (α β : Fin 3) :
    (G α).map (↑) ⊗ₖ (G β).map (↑) = -(S α ⊗ₖ S β) := by
  ext ⟨a, b⟩ ⟨c, e⟩
  simp only [G, S, Matrix.neg_apply, kroneckerMap_apply, map_apply, of_apply,
    Complex.ofReal_intCast]
  ring_nf
  rw [Complex.I_sq]
  ring

/-- The two-site exchange twisted by a real matrix `R`: `∑_{αβ} R_{αβ} S^α ⊗ S^β`. -/
noncomputable def rotExchange (R : Matrix (Fin 3) (Fin 3) ℝ) :
    Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ :=
  ∑ α, ∑ β, (R α β : ℂ) • (S α ⊗ₖ S β)

/-- The rotated exchange of a real matrix is Hermitian. -/
theorem rotExchange_isHermitian (R : Matrix (Fin 3) (Fin 3) ℝ) : (rotExchange R).IsHermitian :=
  isSelfAdjoint_sum _ fun α _ ↦ isSelfAdjoint_sum _ fun β _ ↦ by
    change ((R α β : ℂ) • (S α ⊗ₖ S β))ᴴ = (R α β : ℂ) • (S α ⊗ₖ S β)
    rw [conjTranspose_smul, conjTranspose_kronecker, (S_isHermitian α).eq, (S_isHermitian β).eq,
      Complex.star_def, Complex.conj_ofReal]

/-- The untwisted exchange is `S · S`. -/
theorem rotExchange_one : rotExchange 1 = exchange := by
  simp [rotExchange, one_apply, exchange, apply_ite Complex.ofReal]

/-- The permutation matrix `[β = π α]` of an axis permutation, as a rotation. -/
def permRot (π : Fin 3 → Fin 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  of fun α β ↦ if β = π α then 1 else 0

/-- The exchange twisted by an axis permutation is `∑_α S^α ⊗ S^{π α}`. -/
theorem rotExchange_permRot (π : Fin 3 → Fin 3) :
    rotExchange (permRot π) = ∑ α, S α ⊗ₖ S (π α) := by
  simp [rotExchange, permRot, apply_ite Complex.ofReal]

/-- A real `3 × 3` matrix **preserves the cross product**: `Q u × Q v = Q (u × v)`, entrywise
`∑_{b'c'} ε_{α b' c'} Q_{b'b} Q_{c'c} = ∑_β Q_{αβ} ε_{βbc}`. Every proper rotation does. -/
def PreservesCross (Q : Matrix (Fin 3) (Fin 3) ℝ) : Prop :=
  ∀ α b c, ∑ b', ∑ c', (leviCivita α b' c' : ℝ) * Q b' b * Q c' c =
    ∑ β, Q α β * (leviCivita β b c : ℝ)

/-- **The spin matrices are a vector.** For `Q` preserving the cross product,
`Qᵀ S^α Q = ∑_β Q_{αβ} S^β`. -/
theorem conj_S_of_preservesCross {Q : Matrix (Fin 3) (Fin 3) ℝ} (hQ : PreservesCross Q)
    (α : Fin 3) :
    (Q.map (↑))ᴴ * S α * Q.map (↑) = ∑ β, (Q α β : ℂ) • S β := by
  ext b c
  have h := congrArg (fun r : ℝ ↦ -Complex.I * (r : ℂ)) (hQ α b c)
  simp only [Fin.sum_univ_three, Complex.ofReal_add, Complex.ofReal_mul,
    Complex.ofReal_intCast] at h
  simp only [mul_apply, conjTranspose_apply, map_apply, S, of_apply, Matrix.sum_apply,
    Matrix.smul_apply, smul_eq_mul, Fin.sum_univ_three, Complex.star_def, Complex.conj_ofReal]
  linear_combination h

/-- The bilinear expansion of a tensor product of two linear combinations. -/
theorem sum_smul_kronecker_sum_smul (a b : Fin 3 → ℂ) (A B : Fin 3 → Matrix (Fin 3) (Fin 3) ℂ) :
    (∑ γ, a γ • A γ) ⊗ₖ (∑ δ, b δ • B δ) = ∑ γ, ∑ δ, (a γ * b δ) • (A γ ⊗ₖ B δ) := by
  simp only [Fin.sum_univ_three, add_kronecker, kronecker_add, smul_kronecker, kronecker_smul]
  module

/-- **Rotating a twisted exchange.** For an orthogonal `Q` preserving the cross product, the
two-site conjugation by `Q ⊗ Q` turns `rotExchange R` into `rotExchange (Qᵀ R Q)`. -/
theorem conj_rotExchange {Q : Matrix (Fin 3) (Fin 3) ℝ} (hQ : PreservesCross Q)
    (R : Matrix (Fin 3) (Fin 3) ℝ) :
    (Q.map (↑) ⊗ₖ Q.map (↑))ᴴ * rotExchange R * (Q.map (↑) ⊗ₖ Q.map (↑)) =
      rotExchange (Qᵀ * R * Q) := by
  set V : Matrix (Fin 3) (Fin 3) ℂ := Q.map (↑)
  have hS (α β : Fin 3) : (V ⊗ₖ V)ᴴ * (S α ⊗ₖ S β) * (V ⊗ₖ V) =
      ∑ γ, ∑ δ, ((Q α γ : ℂ) * Q β δ) • (S γ ⊗ₖ S δ) := by
    rw [conjTranspose_kronecker, ← mul_kronecker_mul, ← mul_kronecker_mul,
      conj_S_of_preservesCross hQ, conj_S_of_preservesCross hQ, sum_smul_kronecker_sum_smul]
  simp only [rotExchange, Finset.mul_sum, Finset.sum_mul, Matrix.mul_smul, Matrix.smul_mul, hS,
    Finset.smul_sum, smul_smul]
  simp only [Fin.sum_univ_three, mul_apply, transpose_apply, Complex.ofReal_add,
    Complex.ofReal_mul]
  module

/-- A real orthogonal matrix is unitary as a complex matrix. -/
theorem map_conjTranspose_mul_map {Q : Matrix (Fin 3) (Fin 3) ℝ} (hQ : Qᵀ * Q = 1) :
    (Q.map (↑) : Matrix (Fin 3) (Fin 3) ℂ)ᴴ * Q.map (↑) = 1 := by
  have h : (Q.map (↑) : Matrix (Fin 3) (Fin 3) ℂ)ᴴ = Complex.ofRealHom.mapMatrix Qᵀ := by
    ext i j
    simp [conjTranspose_apply]
  rw [h, show (Q.map (↑) : Matrix (Fin 3) (Fin 3) ℂ) = Complex.ofRealHom.mapMatrix Q from rfl,
    ← map_mul, hQ, map_one]

end CertifiedGaps.SpinOne
