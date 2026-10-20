/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas

/-!
# Two eigenvalue bounds for Hermitian matrices

Both come from the unitary diagonalisation `A = U D Uᴴ` of `Matrix.IsHermitian`:

* `exists_eigenvalue_near` (residual bound): for every vector `x` and real `θ` some eigenvalue
  satisfies `(λᵢ - θ)² ‖x‖² ≤ ‖A x - θ x‖²`.
* `card_eigenvalues_lt_le` (counting by comparison): if `Re ⟪x, A x⟫ ≥ ρ ‖x‖²` for every `x`
  vanishing on a set `B` of basis vectors, then at most `|B|` eigenvalues lie below `ρ`.

Squared norms are written `∑ σ, normSq (x σ)`, so that they can be evaluated exactly.
-/

open Matrix Complex

namespace SeedKernel.Ising

variable {n : Type*} [Fintype n]

/-- The squared Euclidean norm of a complex vector. -/
noncomputable def normSqSum (x : n → ℂ) : ℝ := ∑ σ, normSq (x σ)

theorem normSqSum_eq_re_dotProduct (x : n → ℂ) : normSqSum x = (star x ⬝ᵥ x).re := by
  simp [normSqSum, dotProduct, Complex.re_sum, Complex.normSq_apply, mul_comm]

theorem star_mulVec_dotProduct_mulVec {m : Type*} [Fintype m] (M : Matrix n m ℂ)
    (y z : m → ℂ) : star (M *ᵥ y) ⬝ᵥ (M *ᵥ z) = star y ⬝ᵥ ((Mᴴ * M) *ᵥ z) := by
  rw [star_mulVec, ← dotProduct_mulVec, mulVec_mulVec]

/-- A unitary matrix preserves the squared norm. -/
theorem normSqSum_unitary_mulVec [DecidableEq n] {U : Matrix n n ℂ} (hU : Uᴴ * U = 1) (y : n → ℂ) :
    normSqSum (U *ᵥ y) = normSqSum y := by
  rw [normSqSum_eq_re_dotProduct, normSqSum_eq_re_dotProduct, star_mulVec_dotProduct_mulVec,
    hU, one_mulVec]

variable [DecidableEq n] {A : Matrix n n ℂ} (hA : A.IsHermitian)

/-- The unitary of eigenvectors. -/
noncomputable abbrev eigU : Matrix n n ℂ := (hA.eigenvectorUnitary : Matrix n n ℂ)

theorem eigU_star_mul : (eigU hA)ᴴ * eigU hA = 1 := Unitary.coe_star_mul_self _

theorem eigU_mul_star : eigU hA * (eigU hA)ᴴ = 1 := Unitary.coe_mul_star_self _

/-- `Uᴴ A U` is the diagonal matrix of eigenvalues. -/
theorem eigU_conj : (eigU hA)ᴴ * A * eigU hA = diagonal (fun i ↦ (hA.eigenvalues i : ℂ)) := by
  have hspec : A = eigU hA * diagonal (fun i ↦ (hA.eigenvalues i : ℂ)) * (eigU hA)ᴴ := by
    conv_lhs => rw [hA.spectral_theorem, Unitary.conjStarAlgAut_apply]
    rfl
  generalize hU : eigU hA = U at hspec
  have h1 : Uᴴ * U = 1 := hU ▸ eigU_star_mul hA
  conv_lhs => rw [hspec]
  simp only [Matrix.mul_assoc]
  rw [h1, Matrix.mul_one, ← Matrix.mul_assoc, h1, Matrix.one_mul]

/-- `A (U y) = U (D y)`. -/
theorem mulVec_eigU (y : n → ℂ) :
    A *ᵥ (eigU hA *ᵥ y) = eigU hA *ᵥ fun i ↦ (hA.eigenvalues i : ℂ) * y i := by
  have h : A * eigU hA = eigU hA * diagonal (fun i ↦ (hA.eigenvalues i : ℂ)) := by
    rw [← eigU_conj, ← Matrix.mul_assoc, ← Matrix.mul_assoc, eigU_mul_star, Matrix.one_mul]
  rw [mulVec_mulVec, h, ← mulVec_mulVec]
  congr 1
  ext i
  simp [mulVec_diagonal]

/-- In the eigenbasis the quadratic form is `∑ λᵢ |yᵢ|²`. -/
theorem re_quadratic_eigU (y : n → ℂ) :
    (star (eigU hA *ᵥ y) ⬝ᵥ (A *ᵥ (eigU hA *ᵥ y))).re =
      ∑ i, hA.eigenvalues i * normSq (y i) := by
  rw [mulVec_eigU, star_mulVec_dotProduct_mulVec, eigU_star_mul, one_mulVec]
  simp only [dotProduct, Pi.star_apply, Complex.re_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp [Complex.normSq_apply]
  ring

/-- **Residual bound.** Some eigenvalue satisfies `(λᵢ - θ)² ‖x‖² ≤ ‖A x - θ x‖²`. -/
theorem exists_eigenvalue_near [Nonempty n] (x : n → ℂ) (θ : ℝ) :
    ∃ i, (hA.eigenvalues i - θ) ^ 2 * normSqSum x ≤ normSqSum (A *ᵥ x - (θ : ℂ) • x) := by
  set y := (eigU hA)ᴴ *ᵥ x
  have hx : x = eigU hA *ᵥ y := by
    rw [mulVec_mulVec, eigU_mul_star, one_mulVec]
  have hres : A *ᵥ x - (θ : ℂ) • x =
      eigU hA *ᵥ fun i ↦ ((hA.eigenvalues i - θ : ℝ) : ℂ) * y i := by
    rw [hx, mulVec_eigU, ← mulVec_smul, ← mulVec_sub]
    congr 1
    ext i
    simp [sub_mul]
  obtain ⟨i, -, hi⟩ := Finset.exists_min_image Finset.univ
    (fun i ↦ (hA.eigenvalues i - θ) ^ 2) Finset.univ_nonempty
  refine ⟨i, ?_⟩
  rw [hres, normSqSum_unitary_mulVec (eigU_star_mul hA), hx,
    normSqSum_unitary_mulVec (eigU_star_mul hA), normSqSum, normSqSum, Finset.mul_sum]
  refine Finset.sum_le_sum fun j _ ↦ ?_
  rw [map_mul, Complex.normSq_ofReal, ← sq]
  exact mul_le_mul_of_nonneg_right (hi j (Finset.mem_univ j)) (normSq_nonneg _)

/-- **Counting by comparison.** If `Re ⟪x, A x⟫ ≥ ρ ‖x‖²` for every `x` vanishing on `B`, then
fewer than `|B| + 1` eigenvalues lie below `ρ`. -/
theorem card_eigenvalues_lt_le (B : Finset n) (ρ : ℝ)
    (hB : ∀ x : n → ℂ, (∀ σ ∈ B, x σ = 0) →
      ρ * normSqSum x ≤ (star x ⬝ᵥ (A *ᵥ x)).re) :
    (Finset.univ.filter fun i ↦ hA.eigenvalues i < ρ).card ≤ B.card := by
  classical
  set S := Finset.univ.filter fun i ↦ hA.eigenvalues i < ρ
  by_contra hlt
  push Not at hlt
  let ext : (S → ℂ) →ₗ[ℂ] (n → ℂ) :=
    { toFun := fun c i ↦ if h : i ∈ S then c ⟨i, h⟩ else 0
      map_add' := fun c d ↦ by ext i; by_cases h : i ∈ S <;> simp [h]
      map_smul' := fun a c ↦ by ext i; by_cases h : i ∈ S <;> simp [h] }
  let f : (S → ℂ) →ₗ[ℂ] (B → ℂ) :=
    (LinearMap.pi fun b : B ↦ LinearMap.proj (R := ℂ) (φ := fun _ : n ↦ ℂ) b.1) ∘ₗ
      (Matrix.mulVecLin (eigU hA)) ∘ₗ ext
  have hker : LinearMap.ker f ≠ ⊥ := LinearMap.ker_ne_bot_of_finrank_lt (by
    simpa [Module.finrank_fintype_fun_eq_card] using hlt)
  obtain ⟨c, hc, hc0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hker
  set y := ext c
  have hyS : ∀ i, i ∉ S → y i = 0 := fun i hi ↦ by simp [y, ext, hi]
  obtain ⟨j, hj⟩ : ∃ j : S, c j ≠ 0 := by
    by_contra h
    push Not at h
    exact hc0 (funext h)
  have hyj : y j.1 ≠ 0 := by simpa [y, ext, j.2] using hj
  have hvan : ∀ σ ∈ B, (eigU hA *ᵥ y) σ = 0 := fun σ hσ ↦ by
    have := congrFun (LinearMap.mem_ker.1 hc) ⟨σ, hσ⟩
    simpa [f, y] using this
  have h := hB _ hvan
  rw [re_quadratic_eigU, normSqSum_unitary_mulVec (eigU_star_mul hA), normSqSum,
    Finset.mul_sum] at h
  refine absurd h (not_le.2 (Finset.sum_lt_sum (fun i _ ↦ ?_) ⟨j.1, Finset.mem_univ _, ?_⟩))
  · by_cases hi : i ∈ S
    · exact mul_le_mul_of_nonneg_right (Finset.mem_filter.1 hi).2.le (normSq_nonneg _)
    · simp [hyS i hi]
  · exact mul_lt_mul_of_pos_right (Finset.mem_filter.1 j.2).2 (normSq_pos.2 hyj)

end SeedKernel.Ising
