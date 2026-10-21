/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.MomentForm
import Mathlib.Data.Matrix.PEquiv

/-!
# `ℤ₃` sectors of a real symmetric matrix

Let `Y` be a real symmetric matrix on a finite set `D` and `σ` a permutation of `D` of order
dividing three that preserves `Y`, `Y (σ x) (σ y) = Y x y`. The permutation matrix `Π` of `σ`
commutes with `Y` and splits the space into its eigenspaces for the eigenvalues `1`, `ω`, `ω̄`
(`ω = e^{2πi/3}`). Complex conjugation exchanges the last two and preserves the real matrix `Y`,
so `Y` has the same spectrum on both. Hence there are real numbers `aᵢ` (the spectrum on the
fixed space) and `bᵢ` (the spectrum on the `ω` space) with

  `Tr Y^L = ∑ᵢ aᵢ^L + 2 ∑ᵢ bᵢ^L`,  `Tr (Π Y^L) = ∑ᵢ aᵢ^L - ∑ᵢ bᵢ^L`  (`L ≥ 1`).

The proof does not construct the eigenspaces. The projections onto them are the polynomials
`P_ζ = (1 + ζ̄ Π + ζ Π²) / 3` in `Π`, for `ζ = 1, ω`; each is Hermitian, idempotent and commutes
with `Y`, so `Tr (P_ζ Y^L) = Tr (P_ζ Y)^L` is the power sum of the real eigenvalues of the
Hermitian matrix `P_ζ Y` (zeros from the complement of the range do not contribute at `L ≥ 1`).
The two traces are `(s₀ + 2 s₁) / 3` and `(s₀ - s₁) / 3` with `s₀ = Tr Y^L` and
`s₁ = Tr (Π Y^L)`, because `Tr (Π² Y^L) = Tr (Π Y^L)` for real symmetric `Y` and
`ω + ω̄ = -1`.

## Main results

* `exists_trace_mul_pow_of_isIdempotent`: `Tr (P Y^L)` as a power sum of real numbers, for a
  Hermitian projection `P` commuting with a Hermitian `Y`.
* `exists_cyclic_sectors`: the `ℤ₃` sector decomposition above.
-/

open Matrix

namespace CertifiedGaps

variable {D : Type*} [Fintype D] [DecidableEq D]

/-- **Power sums on the range of a projection.** If a Hermitian projection `P` commutes with a
Hermitian matrix `Y`, then `Tr (P Y^L)` is, for every `L ≥ 1`, the `L`th power sum of real
numbers (the eigenvalues of the Hermitian matrix `P Y`). -/
theorem exists_trace_mul_pow_of_isIdempotent {P Y : Matrix D D ℂ} (hP : P.IsHermitian)
    (hPP : IsIdempotentElem P) (hY : Y.IsHermitian) (hc : Commute P Y) :
    ∃ μ : D → ℝ, ∀ L, 1 ≤ L → trace (P * Y ^ L) = ((∑ i, μ i ^ L : ℝ) : ℂ) := by
  have hM : (P * Y).IsHermitian := by
    rw [IsHermitian, conjTranspose_mul, hY.eq, hP.eq]
    exact hc.eq.symm
  refine ⟨hM.eigenvalues, fun L hL ↦ ?_⟩
  obtain ⟨k, rfl⟩ : ∃ k, L = k + 1 := ⟨L - 1, by omega⟩
  rw [← trace_pow_eq_sum_eigenvalues hM, hc.mul_pow, hPP.pow_succ_eq]

/-! ### The cyclic projections -/

section Cyclic

variable {A Y : Matrix D D ℂ}

/-- The projection `(1 + ζ̄ A + ζ A²) / 3` onto the `ζ` eigenspace of an `A` with `A³ = 1`. -/
noncomputable def cyclicProj (A : Matrix D D ℂ) (ζ : ℂ) : Matrix D D ℂ :=
  (1 / 3 : ℂ) • (1 + star ζ • A + ζ • A ^ 2)

/-- `A` acts on the range of the cyclic projection by `ζ`. -/
theorem mul_cyclicProj (hA : A ^ 3 = 1) {ζ : ℂ} (hζ : ζ * star ζ = 1) (hζ2 : ζ ^ 2 = star ζ) :
    A * cyclicProj A ζ = ζ • cyclicProj A ζ := by
  have h3 : A * A ^ 2 = 1 := by rw [← pow_succ', hA]
  simp only [cyclicProj, mul_smul_comm, mul_add, mul_one, h3, smul_add, smul_smul, ← sq]
  rw [show ζ * (1 / 3 * star ζ) = 1 / 3 by linear_combination (1 / 3 : ℂ) * hζ,
    show ζ * (1 / 3 * ζ) = 1 / 3 * star ζ by linear_combination (1 / 3 : ℂ) * hζ2]
  module

/-- The cyclic projection is idempotent. -/
theorem isIdempotentElem_cyclicProj (hA : A ^ 3 = 1) {ζ : ℂ} (hζ : ζ * star ζ = 1)
    (hζ2 : ζ ^ 2 = star ζ) : IsIdempotentElem (cyclicProj A ζ) := by
  have h1 := mul_cyclicProj hA hζ hζ2
  have h2 : A ^ 2 * cyclicProj A ζ = ζ ^ 2 • cyclicProj A ζ := by
    rw [sq A, Matrix.mul_assoc, h1, mul_smul_comm, h1, smul_smul, ← sq]
  have h3 : ζ * ζ ^ 2 = 1 := by rw [hζ2, hζ]
  have e : cyclicProj A ζ * cyclicProj A ζ = (1 / 3 : ℂ) • (cyclicProj A ζ +
      star ζ • (A * cyclicProj A ζ) + ζ • (A ^ 2 * cyclicProj A ζ)) := by
    conv_lhs => arg 1; rw [cyclicProj]
    simp only [Matrix.smul_mul, Matrix.add_mul, Matrix.one_mul]
  rw [IsIdempotentElem, e, h1, h2, smul_smul, smul_smul, mul_comm (star ζ), hζ, h3, one_smul]
  module

/-- The cyclic projection is Hermitian when `Aᴴ = A²`. -/
theorem cyclicProj_isHermitian (hA : A ^ 3 = 1) (hAH : Aᴴ = A ^ 2) (ζ : ℂ) :
    (cyclicProj A ζ).IsHermitian := by
  have hA4 : (A ^ 2) ^ 2 = A := by
    rw [← pow_mul, show 2 * 2 = 3 + 1 from rfl, pow_succ, hA, Matrix.one_mul]
  rw [IsHermitian, cyclicProj, conjTranspose_smul, conjTranspose_add, conjTranspose_add,
    conjTranspose_one, conjTranspose_smul, conjTranspose_smul, conjTranspose_pow, hAH, hA4,
    star_star, show star (1 / 3 : ℂ) = 1 / 3 by simp]
  abel_nf

/-- The cyclic projection commutes with every matrix that commutes with `A`. -/
theorem commute_cyclicProj (hc : Commute A Y) (ζ : ℂ) : Commute (cyclicProj A ζ) Y :=
  (((Commute.one_left Y).add_left (hc.smul_left _)).add_left
    ((hc.pow_left 2).smul_left _)).smul_left _

/-- The trace against the cyclic projection, in the traces against `1`, `A` and `A²`. -/
theorem trace_cyclicProj_mul (ζ : ℂ) (M : Matrix D D ℂ) :
    trace (cyclicProj A ζ * M) =
      (1 / 3) * (trace M + star ζ * trace (A * M) + ζ * trace (A ^ 2 * M)) := by
  simp only [cyclicProj, Matrix.smul_mul, Matrix.add_mul, Matrix.one_mul, trace_smul,
    trace_add, smul_eq_mul]

end Cyclic

/-- The primitive cube root of unity `ω = e^{2πi/3} = -1/2 + (√3/2) i`. -/
noncomputable def cubeRoot : ℂ := ⟨-1 / 2, √3 / 2⟩

/-- `ω ω̄ = 1`. -/
theorem cubeRoot_mul_star : cubeRoot * star cubeRoot = 1 := by
  have h := Real.sq_sqrt (show (0 : ℝ) ≤ 3 by norm_num)
  apply Complex.ext <;> simp [cubeRoot] <;> nlinarith

/-- `ω² = ω̄`. -/
theorem cubeRoot_sq : cubeRoot ^ 2 = star cubeRoot := by
  have h := Real.sq_sqrt (show (0 : ℝ) ≤ 3 by norm_num)
  apply Complex.ext <;> simp [cubeRoot, sq] <;> nlinarith

/-- `ω + ω̄ = -1`. -/
theorem cubeRoot_add_star : cubeRoot + star cubeRoot = -1 := by
  apply Complex.ext <;> simp [cubeRoot]

/-- **`ℤ₃` sectors.** Let `Y` be real symmetric and `σ` a permutation with `σ³ = 1` that
preserves `Y`. Then there are real numbers `aᵢ` and `bᵢ` with
`Tr Y^L = ∑ᵢ aᵢ^L + 2 ∑ᵢ bᵢ^L` and `∑ₓ (Y^L) (σ x) x = Tr (Π Y^L) = ∑ᵢ aᵢ^L - ∑ᵢ bᵢ^L` for
every `L ≥ 1`, with `Π` the permutation matrix of `σ`. -/
theorem exists_cyclic_sectors {Y : Matrix D D ℝ} (hY : Y.IsHermitian) (σ : D ≃ D)
    (hσ : ∀ x, σ (σ (σ x)) = x) (hYσ : ∀ x y, Y (σ x) (σ y) = Y x y) :
    ∃ a b : D → ℝ, ∀ L, 1 ≤ L →
      trace (Y ^ L) = ∑ i, a i ^ L + 2 * ∑ i, b i ^ L ∧
        ∑ x, (Y ^ L) (σ x) x = ∑ i, a i ^ L - ∑ i, b i ^ L := by
  set A : Matrix D D ℂ := σ.toPEquiv.toMatrix with hAdef
  set Yc : Matrix D D ℂ := Y.map (↑) with hYc
  have hAM : ∀ M : Matrix D D ℂ, A * M = M.submatrix σ id := PEquiv.toMatrix_toPEquiv_mul σ
  have hMA : ∀ M : Matrix D D ℂ, M * A = M.submatrix id σ.symm :=
    fun M ↦ PEquiv.mul_toMatrix_toPEquiv M σ
  have hA3 : A ^ 3 = 1 := by
    rw [show A ^ 3 = A * (A * (A * 1)) by simp only [Matrix.mul_one, pow_succ', pow_zero]]
    ext x y
    simp only [hAM]
    change (1 : Matrix D D ℂ) (σ (σ (σ x))) y = (1 : Matrix D D ℂ) x y
    rw [hσ]
  have hAH : Aᴴ = A ^ 2 := by
    rw [sq, ← Matrix.mul_one (A * A), Matrix.mul_assoc, hAM, hAM]
    ext x y
    simp only [conjTranspose_apply, hAdef, PEquiv.toMatrix_toPEquiv_apply, submatrix_apply, id,
      one_apply, Pi.single_apply]
    have : x = σ y ↔ σ (σ x) = y :=
      ⟨fun h ↦ by rw [h, hσ], fun h ↦ by rw [← h, hσ]⟩
    by_cases h : x = σ y
    · simp [h, hσ]
    · simp [h, this.not.mp h]
  have hcomm : Commute A Yc := by
    change A * Yc = Yc * A
    rw [hAM, hMA]
    ext x y
    simp only [submatrix_apply, id, hYc, map_apply]
    rw [← hYσ x (σ.symm y), Equiv.apply_symm_apply]
  have hYcH : Yc.IsHermitian := by
    ext i j
    simp [hYc, conjTranspose_apply, ← hY.apply i j]
  obtain ⟨a, ha⟩ := exists_trace_mul_pow_of_isIdempotent (cyclicProj_isHermitian hA3 hAH 1)
    (isIdempotentElem_cyclicProj hA3 (by simp) (by simp)) hYcH (commute_cyclicProj hcomm 1)
  obtain ⟨b, hb⟩ := exists_trace_mul_pow_of_isIdempotent
    (cyclicProj_isHermitian hA3 hAH cubeRoot)
    (isIdempotentElem_cyclicProj hA3 cubeRoot_mul_star cubeRoot_sq) hYcH
    (commute_cyclicProj hcomm cubeRoot)
  refine ⟨a, b, fun L hL ↦ ?_⟩
  set s₀ := trace (Y ^ L)
  set s₁ := ∑ x, (Y ^ L) (σ x) x
  have hpow : Yc ^ L = (Y ^ L).map (↑) := by
    change (Complex.ofRealHom.mapMatrix Y) ^ L = Complex.ofRealHom.mapMatrix (Y ^ L)
    rw [map_pow]
  have t0 : trace (Yc ^ L) = (s₀ : ℂ) := by
    rw [hpow]; simp [s₀, trace, Complex.ofReal_sum]
  have t1 : trace (A * Yc ^ L) = (s₁ : ℂ) := by
    rw [hAM, hpow]; simp [s₁, trace, Complex.ofReal_sum]
  -- `Tr (Π² Y^L) = Tr (Π Y^L)`: substitute `x = σ y` and use the symmetry of `Y^L`.
  have hs₂ : ∑ x, (Y ^ L) (σ (σ x)) x = s₁ := by
    rw [← Equiv.sum_comp σ]
    refine Finset.sum_congr rfl fun y _ ↦ ?_
    rw [hσ, ← (hY.pow L).apply (σ y) y, star_trivial]
  have t2 : trace (A ^ 2 * Yc ^ L) = (s₁ : ℂ) := by
    rw [sq, Matrix.mul_assoc, hAM, hAM, hpow, ← hs₂]
    simp [trace, Complex.ofReal_sum]
  have ea := ha L hL
  have eb := hb L hL
  rw [trace_cyclicProj_mul, t0, t1, t2] at ea eb
  have ea' : ∑ i, a i ^ L = (s₀ + 2 * s₁) / 3 := by
    apply Complex.ofReal_injective
    rw [← ea]; push_cast; simp; ring
  have eb' : ∑ i, b i ^ L = (s₀ - s₁) / 3 := by
    apply Complex.ofReal_injective
    rw [← eb]; push_cast
    linear_combination ((1 / 3) * (s₁ : ℂ)) * cubeRoot_add_star
  constructor <;> linarith

end CertifiedGaps
