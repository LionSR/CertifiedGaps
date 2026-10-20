/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# The Eisenstein integers `ℤ[ω]`

The thermal table at `b = 21/2` of the spin-one paper (`thermal120.tex`) twists the ring by the
rotation `C` by `2π/3`. Its seam entries are `-ω` and `-ω⁻¹ = 1 + ω` with `ω = e^{2πi/3}`, so
the integer Horner recursion runs over `ℤ[ω]`: every coordinate is a pair `u + v ω` of integers,
multiplied with `ω² = -1 - ω` (the paper's `-ω(u + vω) = v + (v - u)ω`).

Mathlib has no Eisenstein integers (`GaussianInt` is `ℤ[i]`, `Zsqrtd (-3)` is `ℤ[√-3]`, an
index-two subring), so this file defines the minimal structure the recursion needs.

## Main definitions

* `EisensteinInt`: pairs `u + v ω`, with a computable commutative ring structure and decidable
  equality, so that integer identities between matrices over `ℤ[ω]` are checked by `decide`.
* `EisensteinInt.omegaC = -1/2 + (√3/2) i` and the ring hom `toComplex : ℤ[ω] →+* ℂ`.
* `EisensteinInt.conjHom`: the conjugation `ω ↦ ω² = -1 - ω`, which `toComplex` carries to complex
  conjugation (`toComplex_conj`).
* `EisensteinInt.normSq z = u² - uv + v²`, the squared modulus of `toComplex z`
  (`norm_toComplex_sq`).
* `EisensteinInt.fdiv z s = ⌊u/s⌋ + ⌊v/s⌋ ω`, the paper's coordinatewise floor; its error is at
  most `2` in modulus (`norm_toComplex_fdiv_sub_le`), the paper's bound `|δu| + |δv| ≤ 2`.
-/

namespace CertifiedGaps.Trace

/-- The **Eisenstein integers** `ℤ[ω]`, `ω = e^{2πi/3}`: the element `u + v ω` is the pair
`(u, v)`. -/
@[ext] structure EisensteinInt where
  /-- The coordinate along `1`. -/
  u : ℤ
  /-- The coordinate along `ω`. -/
  v : ℤ
deriving DecidableEq, Repr

namespace EisensteinInt

instance : Zero EisensteinInt := ⟨⟨0, 0⟩⟩
instance : One EisensteinInt := ⟨⟨1, 0⟩⟩
instance : Add EisensteinInt := ⟨fun z w ↦ ⟨z.u + w.u, z.v + w.v⟩⟩
instance : Neg EisensteinInt := ⟨fun z ↦ ⟨-z.u, -z.v⟩⟩
instance : Sub EisensteinInt := ⟨fun z w ↦ ⟨z.u - w.u, z.v - w.v⟩⟩

/-- `(a + bω)(c + dω) = (ac - bd) + (ad + bc - bd) ω`, from `ω² = -1 - ω`. -/
instance : Mul EisensteinInt :=
  ⟨fun z w ↦ ⟨z.u * w.u - z.v * w.v, z.u * w.v + z.v * w.u - z.v * w.v⟩⟩

instance : NatCast EisensteinInt := ⟨fun n ↦ ⟨n, 0⟩⟩
instance : IntCast EisensteinInt := ⟨fun n ↦ ⟨n, 0⟩⟩

/-- The generator `ω`. -/
def omega : EisensteinInt := ⟨0, 1⟩

@[simp] theorem zero_u : (0 : EisensteinInt).u = 0 := rfl
@[simp] theorem zero_v : (0 : EisensteinInt).v = 0 := rfl
@[simp] theorem one_u : (1 : EisensteinInt).u = 1 := rfl
@[simp] theorem one_v : (1 : EisensteinInt).v = 0 := rfl
@[simp] theorem add_u (z w : EisensteinInt) : (z + w).u = z.u + w.u := rfl
@[simp] theorem add_v (z w : EisensteinInt) : (z + w).v = z.v + w.v := rfl
@[simp] theorem neg_u (z : EisensteinInt) : (-z).u = -z.u := rfl
@[simp] theorem neg_v (z : EisensteinInt) : (-z).v = -z.v := rfl
@[simp] theorem sub_u (z w : EisensteinInt) : (z - w).u = z.u - w.u := rfl
@[simp] theorem sub_v (z w : EisensteinInt) : (z - w).v = z.v - w.v := rfl
@[simp] theorem mul_u (z w : EisensteinInt) : (z * w).u = z.u * w.u - z.v * w.v := rfl
@[simp] theorem mul_v (z w : EisensteinInt) : (z * w).v = z.u * w.v + z.v * w.u - z.v * w.v :=
  rfl
@[simp] theorem natCast_u (n : ℕ) : (n : EisensteinInt).u = n := rfl
@[simp] theorem natCast_v (n : ℕ) : (n : EisensteinInt).v = 0 := rfl
@[simp] theorem intCast_u (n : ℤ) : (n : EisensteinInt).u = n := rfl
@[simp] theorem intCast_v (n : ℤ) : (n : EisensteinInt).v = 0 := rfl
@[simp] theorem omega_u : omega.u = 0 := rfl
@[simp] theorem omega_v : omega.v = 1 := rfl

/-- The ring axioms, each checked coordinatewise by `ring`. The scalar actions and powers are
the recursive defaults, so that the kernel can still evaluate them. -/
instance : CommRing EisensteinInt where
  add_assoc _ _ _ := by ext <;> simp <;> ring
  zero_add _ := by ext <;> simp
  add_zero _ := by ext <;> simp
  nsmul := nsmulRec
  zsmul := zsmulRec
  neg_add_cancel _ := by ext <;> simp
  add_comm _ _ := by ext <;> simp <;> ring
  sub_eq_add_neg _ _ := by ext <;> simp <;> ring
  mul_assoc _ _ _ := by ext <;> simp <;> ring
  one_mul _ := by ext <;> simp
  mul_one _ := by ext <;> simp
  zero_mul _ := by ext <;> simp
  mul_zero _ := by ext <;> simp
  left_distrib _ _ _ := by ext <;> simp <;> ring
  right_distrib _ _ _ := by ext <;> simp <;> ring
  mul_comm _ _ := by ext <;> simp <;> ring
  natCast_zero := by ext <;> simp
  natCast_succ _ := by ext <;> simp
  intCast_ofNat _ := by ext <;> simp
  intCast_negSucc _ := by ext <;> simp [Int.negSucc_eq]

/-! ### The complex embedding -/

/-- `ω = e^{2πi/3} = -1/2 + (√3/2) i` in `ℂ`. -/
noncomputable def omegaC : ℂ := ⟨-1 / 2, √3 / 2⟩

/-- `ω² = -1 - ω` in `ℂ`. -/
theorem omegaC_sq : omegaC * omegaC = -1 - omegaC := by
  have h : √3 * √3 = 3 := Real.mul_self_sqrt (by norm_num)
  apply Complex.ext <;> simp [omegaC] <;> nlinarith

/-- `|ω| = 1`. -/
theorem norm_omegaC : ‖omegaC‖ = 1 := by
  have h : √3 * √3 = 3 := Real.mul_self_sqrt (by norm_num)
  rw [Complex.norm_def, Complex.normSq_apply, Real.sqrt_eq_one]
  simp only [omegaC]
  nlinarith

/-- The complex value `u + v ω` of an Eisenstein integer. -/
noncomputable def toC (z : EisensteinInt) : ℂ := z.u + z.v * omegaC

/-- **The complex embedding** `ℤ[ω] →+* ℂ`, `u + v ω ↦ u + v e^{2πi/3}`. -/
noncomputable def toComplex : EisensteinInt →+* ℂ where
  toFun := toC
  map_zero' := by simp [toC]
  map_one' := by simp [toC]
  map_add' z w := by simp only [toC, add_u, add_v]; push_cast; ring
  map_mul' z w := by
    simp only [toC, mul_u, mul_v]
    push_cast
    linear_combination -(z.v * w.v : ℂ) * omegaC_sq

theorem toComplex_apply (z : EisensteinInt) : toComplex z = z.u + z.v * omegaC := rfl

theorem toComplex_omega : toComplex omega = omegaC := by simp [toComplex_apply]

/-- An integer embeds as itself. -/
theorem toComplex_intCast (n : ℤ) : toComplex (n : EisensteinInt) = n := map_intCast _ n

/-! ### Conjugation -/

/-- The conjugate `u + v ω² = (u - v) - v ω` of `u + v ω`. -/
def conj (z : EisensteinInt) : EisensteinInt := ⟨z.u - z.v, -z.v⟩

/-- **Conjugation** `ω ↦ ω² = -1 - ω`, a ring automorphism of `ℤ[ω]`. -/
def conjHom : EisensteinInt →+* EisensteinInt where
  toFun := conj
  map_zero' := by ext <;> simp [conj]
  map_one' := by ext <;> simp [conj]
  map_add' _ _ := by ext <;> simp [conj] <;> ring
  map_mul' _ _ := by ext <;> simp [conj] <;> ring

theorem conjHom_apply (z : EisensteinInt) : conjHom z = ⟨z.u - z.v, -z.v⟩ := rfl

/-- The embedding carries the conjugation of `ℤ[ω]` to complex conjugation. -/
theorem toComplex_conj (z : EisensteinInt) : toComplex (conjHom z) = star (toComplex z) := by
  apply Complex.ext <;> simp [toComplex_apply, conjHom_apply, omegaC]
  ring

/-! ### Squared modulus -/

/-- The squared modulus `u² - uv + v²` of `u + v ω`, an integer. -/
def normSq (z : EisensteinInt) : ℤ := z.u ^ 2 - z.u * z.v + z.v ^ 2

/-- `|u + v ω|² = u² - uv + v²`. -/
theorem norm_toComplex_sq (z : EisensteinInt) : ‖toComplex z‖ ^ 2 = normSq z := by
  have h : √3 * √3 = 3 := Real.mul_self_sqrt (by norm_num)
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply]
  simp only [toComplex_apply, omegaC, normSq, Complex.add_re, Complex.intCast_re,
    Complex.mul_re, Complex.intCast_im, Complex.add_im, Complex.mul_im]
  push_cast
  nlinarith

/-! ### The coordinatewise floor -/

/-- The paper's coordinatewise floor `⌊u/s⌋ + ⌊v/s⌋ ω` of `(u + v ω)/s`. On `ℤ`, `/` is
`Int.ediv`, the floor for a positive divisor. -/
def fdiv (z : EisensteinInt) (s : ℕ) : EisensteinInt := ⟨z.u / (s : ℤ), z.v / (s : ℤ)⟩

/-- Integer division by a natural number is the floor of the real quotient. -/
theorem intDiv_eq_floor (x : ℤ) (s : ℕ) : x / (s : ℤ) = ⌊(x : ℝ) / (s : ℝ)⌋ := by
  rw [← Rat.floor_intCast_div_natCast, ← Rat.floor_cast (α := ℝ)]
  push_cast; rfl

/-- One coordinate of the floor is within `1` of the exact quotient. -/
theorem abs_intDiv_sub_le (x : ℤ) (s : ℕ) : |((x / (s : ℤ) : ℤ) : ℝ) - (x : ℝ) / s| ≤ 1 := by
  rw [intDiv_eq_floor]
  have h1 := Int.floor_le ((x : ℝ) / s)
  have h2 := Int.lt_floor_add_one ((x : ℝ) / s)
  rw [abs_le]; constructor <;> linarith

/-- **The floor error**: `|fl(z/s) - z/s| ≤ |δu| + |δv| |ω| ≤ 2`. -/
theorem norm_toComplex_fdiv_sub_le (z : EisensteinInt) (s : ℕ) :
    ‖toComplex (fdiv z s) - toComplex z / s‖ ≤ 2 := by
  have hdecomp : toComplex (fdiv z s) - toComplex z / s =
      ((((z.u / (s : ℤ) : ℤ) : ℝ) - (z.u : ℝ) / s : ℝ) : ℂ) +
        ((((z.v / (s : ℤ) : ℤ) : ℝ) - (z.v : ℝ) / s : ℝ) : ℂ) * omegaC := by
    simp only [toComplex_apply, fdiv]
    push_cast
    ring
  rw [hdecomp]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, norm_omegaC, mul_one, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs]
  linarith [abs_intDiv_sub_le z.u s, abs_intDiv_sub_le z.v s]

end EisensteinInt

end CertifiedGaps.Trace
