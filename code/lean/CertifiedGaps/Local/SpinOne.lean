/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.LinearAlgebra.Matrix.Kronecker

/-!
# Spin-one building blocks

The one-site and two-site matrices shared by every spin-one model, in two bases.

**The Cartesian basis** `|x⟩, |y⟩, |z⟩`, in which `(S^α)_{s s'} = -i ε_{α s s'}`. There the
generators `G_α = i S^α = ε_α` are real and antisymmetric and `-S · S = ∑_α G_α ⊗ G_α`, which is
the reflection-positivity factorization of every spin-one model (`SpinOneLargeD`,
`SpinOneDimerized`, the Heisenberg ring). The `π` rotations about the three axes are diagonal
sign matrices (`rotX`, `rotY`, `rotZ`), with `R G_α R = r_α G_α`.

**The `S^z` basis**, ordered `m = -1, 0, 1`, in which the seed computations run because the
exchange has integer entries. The change of basis `szBasis` has as columns the Cartesian
components of the `S^z` eigenstates, with the Condon–Shortley phases

  `|-1⟩ = (|x⟩ - i|y⟩)/√2`, `|0⟩ = |z⟩`, `|1⟩ = -(|x⟩ + i|y⟩)/√2`,

for which `S⁺` has the nonnegative matrix elements `√2`.

## Main definitions

* `leviCivita`, `S α`, `G α`, `exchange = ∑_α S^α ⊗ S^α`.
* `rotX`, `rotY`, `rotZ`: the sign vectors of the `π` rotations about the axes.
* `szBasis`, `spinZ α`: the unitary change of basis and the spin matrices in the `S^z` basis.

## Main results

* `S_isHermitian`, `G_transpose`, `exchange_isHermitian`.
* `map_G_kronecker`: `G_α ⊗ G_α = -(S^α ⊗ S^α)`.
* `szBasis_conjTranspose_mul_self`, `conj_S`: `Uᴴ S^α U = spinZ α`; `conj_rotZ`.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps.SpinOne

/-! ### The Cartesian basis -/

/-- The Levi-Civita symbol on `Fin 3`: `ε₀₁₂ = 1`, totally antisymmetric. -/
def leviCivita : Fin 3 → Fin 3 → Fin 3 → ℤ :=
  ![![![0, 0, 0], ![0, 0, 1], ![0, -1, 0]],
    ![![0, 0, -1], ![0, 0, 0], ![1, 0, 0]],
    ![![0, 1, 0], ![-1, 0, 0], ![0, 0, 0]]]

/-- The Levi-Civita symbol is antisymmetric in its last two indices. -/
theorem leviCivita_swap (a b c : Fin 3) : leviCivita a c b = -leviCivita a b c := by
  fin_cases a <;> fin_cases b <;> fin_cases c <;> rfl

/-- The Cartesian spin-one matrices, `(S^α)_{s s'} = -i ε_{α s s'}`. -/
def S (α : Fin 3) : Matrix (Fin 3) (Fin 3) ℂ :=
  of fun s s' ↦ -Complex.I * (leviCivita α s s' : ℂ)

/-- Each spin-one matrix is Hermitian. -/
theorem S_isHermitian (α : Fin 3) : (S α).IsHermitian := by
  ext s s'
  simp only [conjTranspose_apply, S, of_apply, leviCivita_swap α s s', star_mul',
    star_neg, Complex.star_def, Complex.conj_I, Int.cast_neg, ← Complex.ofReal_intCast,
    Complex.conj_ofReal]
  push_cast
  ring

/-- The Cartesian spin matrices, entry by entry. -/
theorem S_eq (α : Fin 3) : S α = ![
    !![0, 0, 0; 0, 0, -Complex.I; 0, Complex.I, 0],
    !![0, 0, Complex.I; 0, 0, 0; -Complex.I, 0, 0],
    !![0, -Complex.I, 0; Complex.I, 0, 0; 0, 0, 0]] α := by
  ext i j
  fin_cases α <;> fin_cases i <;> fin_cases j <;> simp [S, leviCivita]

/-- The real generators `G_α = ε_α = i S^α`, with entries `ε_{α s s'}`. -/
def G (α : Fin 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  of fun s s' ↦ (leviCivita α s s' : ℝ)

/-- The generators `G_α` are antisymmetric. -/
theorem G_transpose (α : Fin 3) : (G α)ᵀ = -G α := by
  ext s s'; simp [G, leviCivita_swap α s s']

/-- `G_α ⊗ G_α = -(S^α ⊗ S^α)`, since `S^α = -i G_α`: the reflection-positivity factorization
`-S · S = ∑_α G_α ⊗ G_α` term by term. -/
theorem map_G_kronecker (α : Fin 3) :
    (G α).map (↑) ⊗ₖ (G α).map (↑) = -(S α ⊗ₖ S α) := by
  ext ⟨a, b⟩ ⟨c, e⟩
  simp only [G, S, Matrix.neg_apply, kroneckerMap_apply, map_apply, of_apply,
    Complex.ofReal_intCast]
  ring_nf
  rw [Complex.I_sq]
  ring

/-- The two-site exchange interaction `S · S = ∑_α S^α ⊗ S^α`. -/
def exchange : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ :=
  ∑ α, S α ⊗ₖ S α

/-- The exchange interaction is Hermitian. -/
theorem exchange_isHermitian : exchange.IsHermitian :=
  isSelfAdjoint_sum _ fun α _ ↦ by
    change (S α ⊗ₖ S α)ᴴ = S α ⊗ₖ S α
    rw [conjTranspose_kronecker, (S_isHermitian α).eq]

/-! ### The `π` rotations about the axes -/

/-- The signs of the `π` rotation about `x` on the Cartesian components. -/
def rotX : Fin 3 → ℝ := ![1, -1, -1]

/-- The signs of the `π` rotation about `y` on the Cartesian components. -/
def rotY : Fin 3 → ℝ := ![-1, 1, -1]

/-- The signs of the `π` rotation about `z` on the Cartesian components; the twist of the
twisted spin-one ring. -/
def rotZ : Fin 3 → ℝ := ![-1, -1, 1]

/-! ### The `S^z` basis -/

/-- The number `1/√2`, as a complex number. -/
noncomputable def invSqrtTwo : ℂ := ((Real.sqrt 2)⁻¹ : ℝ)

/-- `1/√2` is real. -/
theorem star_invSqrtTwo : star invSqrtTwo = invSqrtTwo := by
  simp [invSqrtTwo, Complex.conj_ofReal]

/-- `(1/√2)² = 1/2`: the only property of `√2` the basis change needs. -/
theorem invSqrtTwo_sq : invSqrtTwo ^ 2 = 1 / 2 := by
  rw [invSqrtTwo, ← Complex.ofReal_pow, inv_pow, Real.sq_sqrt (by norm_num)]
  push_cast; ring

/-- The change of basis from the `S^z` basis to the Cartesian basis: column `i` holds the
Cartesian components `(x, y, z)` of the `S^z` eigenstate with `m = i - 1`, namely
`(|x⟩ - i|y⟩)/√2`, `|z⟩` and `-(|x⟩ + i|y⟩)/√2`. -/
noncomputable def szBasis : Matrix (Fin 3) (Fin 3) ℂ :=
  !![invSqrtTwo, 0, -invSqrtTwo;
     -Complex.I * invSqrtTwo, 0, -Complex.I * invSqrtTwo;
     0, 1, 0]

/-- The conjugate transpose of `szBasis`, entry by entry. -/
theorem szBasis_conjTranspose : szBasisᴴ =
    !![invSqrtTwo, Complex.I * invSqrtTwo, 0;
       0, 0, 1;
       -invSqrtTwo, Complex.I * invSqrtTwo, 0] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [szBasis, star_invSqrtTwo]

/-- `szBasis` is unitary: `Uᴴ U = 1`. -/
theorem szBasis_conjTranspose_mul_self : szBasisᴴ * szBasis = 1 := by
  rw [szBasis_conjTranspose, szBasis, mul_fin_three]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [one_apply] <;> ring_nf <;> simp [invSqrtTwo_sq] <;>
    norm_num

/-- `szBasis` is unitary: `U Uᴴ = 1`. -/
theorem szBasis_mul_conjTranspose_self : szBasis * szBasisᴴ = 1 :=
  mul_eq_one_comm.mp szBasis_conjTranspose_mul_self

/-- The spin matrices `Sˣ = (S⁺ + S⁻)/2`, `Sʸ = (S⁺ - S⁻)/(2i)`, `Sᶻ = diag(-1, 0, 1)` in the
`S^z` basis, with `⟨m + 1| S⁺ |m⟩ = √2`. -/
noncomputable def spinZ : Fin 3 → Matrix (Fin 3) (Fin 3) ℂ := ![
  !![0, invSqrtTwo, 0; invSqrtTwo, 0, invSqrtTwo; 0, invSqrtTwo, 0],
  !![0, Complex.I * invSqrtTwo, 0; -Complex.I * invSqrtTwo, 0, Complex.I * invSqrtTwo;
     0, -Complex.I * invSqrtTwo, 0],
  !![-1, 0, 0; 0, 0, 0; 0, 0, 1]]

/-- The Cartesian spin matrices become `spinZ` in the `S^z` basis. -/
theorem conj_S (α : Fin 3) : szBasisᴴ * S α * szBasis = spinZ α := by
  rw [szBasis_conjTranspose, szBasis, S_eq]
  fin_cases α <;> simp only [spinZ, Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk,
    Matrix.cons_val, mul_fin_three] <;>
  · ext i j
    fin_cases i <;> fin_cases j <;> simp <;> ring_nf <;> simp [invSqrtTwo_sq]

/-- The `π` rotation about `z`, `diag(-1, -1, 1)` in the Cartesian basis, is
`e^{-iπSᶻ} = diag(-1, 1, -1)` in the `S^z` basis. -/
theorem conj_rotZ :
    szBasisᴴ * diagonal (fun α ↦ (rotZ α : ℂ)) * szBasis = diagonal ![-1, 1, -1] := by
  have hP : diagonal (fun α ↦ (rotZ α : ℂ)) = !![-1, 0, 0; 0, -1, 0; 0, 0, 1] := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp [rotZ]
  rw [hP, szBasis_conjTranspose, szBasis, mul_fin_three, mul_fin_three]
  ext i j
  fin_cases i <;> fin_cases j <;> simp <;> ring_nf <;> simp [invSqrtTwo_sq] <;> norm_num

end CertifiedGaps.SpinOne
