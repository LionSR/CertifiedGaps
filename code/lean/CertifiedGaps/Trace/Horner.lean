/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.InnerProductSpace.PiL2
import CertifiedGaps.Spectral.Basic

/-!
# Integer Horner columns with floor rounding

The spin-1 trace enclosures (OpenAI, `thermal72.tex`, §"Integer columns and their error")
evaluate a polynomial `∑_{j ≤ J} pⱼ Rʲ` of a scaled integer matrix `R = D / s` on a unit vector
`e_r` with integer arithmetic only: starting from `X = 0`, for `j = J, J - 1, …, 0`,

  `X ← ⌊D X / s⌋ + cⱼ e_r`   (coordinatewise floor),

with integer coefficients `cⱼ` approximating real targets `aⱼ` (`cⱼ = ⌊Q pⱼ⌋`, `aⱼ = Q pⱼ`).
This file proves the generic rounding bound used there, for any finite index type, any integer
matrix `D`, any positive scale `s` and any coefficients:

* `hornerFloor D s c r J`: the integer recursion, a computable definition;
* `hornerFloor_error`: if `R` is a contraction (`IsScaledContraction D s`) and `|cⱼ - aⱼ| ≤ 1`,
  then `‖X - (∑_{j ≤ J} aⱼ Rʲ) e_r‖₂ ≤ (J + 1)(δ + 1)` whenever `√(card ι) ≤ δ`;
  `hornerFloor_error_floor(_rat)` specialize to `cⱼ = ⌊Q pⱼ⌋`;
* `hornerFloor_error_nsq`: the same bound over `ℂ` with `Spectral.nsq`, for the complexified
  matrix `s⁻¹ • D`, so that it composes with matrix-level bounds stated for `Matrix ι ι ℂ`.

Each step adds at most `√(card ι)` (vector floor) plus `1` (coefficient error), and the
contraction does not enlarge the error already present; the bound is the sum over the `J + 1`
steps. Unlike the paper, self-adjointness of `R` is not needed, only `‖R v‖ ≤ ‖v‖`.

**Representation.** Vectors are functions `ι → ℤ` over a `Fintype` with `DecidableEq`, and the
recursion is a `List.foldr` over `[c₀, …, c_J]`. This keeps the statements index-generic and the
proofs a plain list induction, and the kernel evaluates it well (`decide +kernel`): its `whnf`
cache shares the repeated evaluations of the previous column. Compiled evaluation (`#eval`) has no
such sharing and is exponential in `J`; large blocks are meant to be computed by an external
program and only small ones checked here.
-/

open Matrix

namespace CertifiedGaps.Trace

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ### The integer recursion -/

/-- One integer Horner step `X ↦ ⌊D X / s⌋ + c e_r`. On `ℤ`, `/` is `Int.ediv`, which for a
positive divisor is floor division (`intDiv_eq_floor`), also for negative numerators. -/
def hornerStep (D : Matrix ι ι ℤ) (s : ℕ) (r : ι) (c : ℤ) (X : ι → ℤ) : ι → ℤ :=
  fun i ↦ (D *ᵥ X) i / (s : ℤ) + if i = r then c else 0

/-- The recursion over a coefficient list `[c₀, …, c_J]`: the last coefficient is processed
first, starting from `X = 0`. -/
def hornerList (D : Matrix ι ι ℤ) (s : ℕ) (r : ι) (l : List ℤ) : ι → ℤ :=
  l.foldr (hornerStep D s r) 0

/-- The integer Horner column: `X ← ⌊D X / s⌋ + cⱼ e_r` for `j = J, …, 0`, from `X = 0`. -/
def hornerFloor (D : Matrix ι ι ℤ) (s : ℕ) (c : ℕ → ℤ) (r : ι) (J : ℕ) : ι → ℤ :=
  hornerList D s r ((List.range (J + 1)).map c)

/-- `/` on `ℤ` is `Int.ediv`, the operation the kernel evaluates. -/
theorem intDiv_eq_ediv (x y : ℤ) : x / y = Int.ediv x y := rfl

/-- Integer division by a natural number is the floor of the real quotient. -/
theorem intDiv_eq_floor (x : ℤ) (s : ℕ) : x / (s : ℤ) = ⌊(x : ℝ) / (s : ℝ)⌋ := by
  rw [← Rat.floor_intCast_div_natCast, ← Rat.floor_cast (α := ℝ)]
  push_cast; rfl

/-- The recursion on the empty list is the zero vector. -/
@[simp] theorem hornerList_nil (D : Matrix ι ι ℤ) (s : ℕ) (r : ι) :
    hornerList D s r [] = 0 := rfl

/-- Unfolding one step of the recursion. -/
@[simp] theorem hornerList_cons (D : Matrix ι ι ℤ) (s : ℕ) (r : ι) (c : ℤ) (l : List ℤ) :
    hornerList D s r (c :: l) = hornerStep D s r c (hornerList D s r l) := rfl

/-! ### The real recursion -/

/-- The real matrix `R = D / s`. -/
noncomputable def scaledMatrix (D : Matrix ι ι ℤ) (s : ℕ) : Matrix ι ι ℝ :=
  (s : ℝ)⁻¹ • D.map ((↑) : ℤ → ℝ)

/-- `R = D / s` is a Euclidean contraction, stated without division:
`‖D v‖² ≤ s² ‖v‖²` for every real `v`. -/
def IsScaledContraction (D : Matrix ι ι ℤ) (s : ℕ) : Prop :=
  ∀ v : ι → ℝ, ∑ i, ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i ^ 2 ≤ (s : ℝ) ^ 2 * ∑ i, v i ^ 2

/-- The ideal recursion `Y ↦ R Y + aⱼ e_r` over a list of real coefficients. -/
noncomputable def idealList (R : Matrix ι ι ℝ) (r : ι) (a : List ℝ) : ι → ℝ :=
  a.foldr (fun aj Y ↦ R *ᵥ Y + Pi.single r aj) 0

/-- The ideal recursion on the empty list is the zero vector. -/
@[simp] theorem idealList_nil (R : Matrix ι ι ℝ) (r : ι) : idealList R r [] = 0 := rfl

/-- Unfolding one step of the ideal recursion. -/
@[simp] theorem idealList_cons (R : Matrix ι ι ℝ) (r : ι) (a : ℝ) (l : List ℝ) :
    idealList R r (a :: l) = R *ᵥ idealList R r l + Pi.single r a := rfl

/-- The ideal recursion is the polynomial `∑_{j ≤ J} aⱼ Rʲ` applied to `e_r` (Horner's rule). -/
theorem idealList_range_map (R : Matrix ι ι ℝ) (r : ι) (a : ℕ → ℝ) (J : ℕ) :
    idealList R r ((List.range (J + 1)).map a) =
      (∑ j ∈ Finset.range (J + 1), a j • R ^ j) *ᵥ Pi.single r 1 := by
  induction J generalizing a with
  | zero =>
    ext k
    simp [idealList, Pi.single_apply, Matrix.one_apply]
  | succ J ih =>
    rw [List.range_succ_eq_map, List.map_cons, List.map_map, idealList_cons, ih,
      Finset.sum_range_succ' _ (J + 1), add_mulVec, mulVec_mulVec, Finset.mul_sum]
    congr 1
    · simp only [Function.comp_apply, pow_succ', Matrix.mul_smul]
    · ext k
      simp [Pi.single_apply, Matrix.one_apply]

/-! ### Euclidean norm helpers -/

/-- Euclidean norm of a real vector, routed through `EuclideanSpace` to reuse its triangle
inequality; the public statements use the explicit `√(∑ i, v i ^ 2)` form. -/
private noncomputable abbrev l2 (v : ι → ℝ) : ℝ := ‖WithLp.toLp 2 v‖

omit [DecidableEq ι] in
private theorem l2_eq (v : ι → ℝ) : l2 v = Real.sqrt (∑ i, v i ^ 2) := by
  simp [l2, EuclideanSpace.norm_eq]

omit [DecidableEq ι] in
private theorem l2_add (u v : ι → ℝ) : l2 (u + v) ≤ l2 u + l2 v := by
  simpa [l2, WithLp.toLp_add] using norm_add_le (WithLp.toLp 2 u) (WithLp.toLp 2 v)

omit [DecidableEq ι] in
private theorem l2_le_sqrt_card (u : ι → ℝ) (hu : ∀ i, |u i| ≤ 1) :
    l2 u ≤ Real.sqrt (Fintype.card ι) := by
  rw [l2_eq]
  refine Real.sqrt_le_sqrt ?_
  calc ∑ i, u i ^ 2 ≤ ∑ _i : ι, (1 : ℝ) :=
        Finset.sum_le_sum fun i _ ↦ (sq_le_one_iff_abs_le_one _).2 (hu i)
    _ = Fintype.card ι := by simp

private theorem l2_single (r : ι) (x : ℝ) : l2 (Pi.single r x) = |x| := by
  rw [l2_eq, ← Real.sqrt_sq_eq_abs]
  congr 1
  simp [Pi.single_apply, ite_pow]

omit [DecidableEq ι] in
private theorem l2_scaled_le {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    (v : ι → ℝ) : l2 (scaledMatrix D s *ᵥ v) ≤ l2 v := by
  rw [l2_eq, l2_eq]
  refine Real.sqrt_le_sqrt ?_
  have hs' : (0 : ℝ) < s := Nat.cast_pos.2 hs
  have h := hD v
  simp only [scaledMatrix, smul_mulVec, Pi.smul_apply, smul_eq_mul, mul_pow, ← Finset.mul_sum]
  rw [inv_pow]
  calc ((s : ℝ) ^ 2)⁻¹ * ∑ i, ((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i ^ 2
      ≤ ((s : ℝ) ^ 2)⁻¹ * ((s : ℝ) ^ 2 * ∑ i, v i ^ 2) := by gcongr
    _ = ∑ i, v i ^ 2 := by field_simp

/-! ### One step and the accumulated error -/

omit [DecidableEq ι] in
/-- `R` applied to an integer vector is the exact integer product divided by `s`. -/
private theorem scaledMatrix_mulVec_intCast (D : Matrix ι ι ℤ) (s : ℕ) (X : ι → ℤ) (i : ι) :
    (scaledMatrix D s *ᵥ fun j ↦ (X j : ℝ)) i = ((D *ᵥ X) i : ℝ) / s := by
  simp [scaledMatrix, mulVec, dotProduct, div_eq_inv_mul, Finset.mul_sum, mul_assoc]

/-- One floored step moves the error by at most `√(card ι)` (vector floor) plus `1`
(coefficient error), on top of the propagated error, which `R` does not enlarge. -/
private theorem step_error {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    (r : ι) {c : ℤ} {a : ℝ} (hca : |(c : ℝ) - a| ≤ 1) (X : ι → ℤ) (Y : ι → ℝ) :
    l2 (fun i ↦ (hornerStep D s r c X i : ℝ) - (scaledMatrix D s *ᵥ Y + Pi.single r a : ι → ℝ) i) ≤
      l2 (fun i ↦ (X i : ℝ) - Y i) + Real.sqrt (Fintype.card ι) + 1 := by
  -- `u` is the floor error of the vector division
  obtain ⟨u, hu⟩ : ∃ u : ι → ℝ,
      u = fun i ↦ (((D *ᵥ X) i / (s : ℤ) : ℤ) : ℝ) - ((D *ᵥ X) i : ℝ) / s := ⟨_, rfl⟩
  have hdecomp : (fun i ↦ (hornerStep D s r c X i : ℝ) - (scaledMatrix D s *ᵥ Y +
      Pi.single r a : ι → ℝ) i) = u + scaledMatrix D s *ᵥ (fun i ↦ (X i : ℝ) - Y i) +
      Pi.single r ((c : ℝ) - a) := by
    funext i
    have hsub : (fun i ↦ (X i : ℝ) - Y i) = (fun i ↦ (X i : ℝ)) - Y := rfl
    rw [hsub, mulVec_sub, hu]
    simp only [hornerStep, Pi.add_apply, Pi.sub_apply, scaledMatrix_mulVec_intCast,
      Pi.single_apply]
    push_cast
    split_ifs <;> ring
  have hu1 : ∀ i, |u i| ≤ 1 := by
    intro i
    rw [hu]; dsimp only; rw [intDiv_eq_floor]
    have h1 := Int.floor_le (((D *ᵥ X) i : ℝ) / s)
    have h2 := Int.lt_floor_add_one (((D *ᵥ X) i : ℝ) / s)
    rw [abs_le]; constructor <;> linarith
  have h1 := l2_le_sqrt_card u hu1
  have h2 := l2_scaled_le hs hD (fun i ↦ (X i : ℝ) - Y i)
  have h3 := l2_single r ((c : ℝ) - a)
  have h4 := l2_add (u + scaledMatrix D s *ᵥ (fun i ↦ (X i : ℝ) - Y i))
    (Pi.single r ((c : ℝ) - a))
  have h5 := l2_add u (scaledMatrix D s *ᵥ (fun i ↦ (X i : ℝ) - Y i))
  rw [hdecomp]
  linarith

/-! ### Main error bounds -/

/-- **Rounding error of the integer Horner recursion (coefficient lists).** If `R = D / s` is a
contraction and every integer coefficient `cⱼ` is within `1` of its real target `aⱼ`, the
floored recursion stays within `(number of steps) · (δ + 1)` of the ideal recursion, for any
`δ ≥ √(card ι)`. -/
theorem hornerList_error {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι) {c : List ℤ} {a : List ℝ}
    (hca : List.Forall₂ (fun (cj : ℤ) (aj : ℝ) ↦ |(cj : ℝ) - aj| ≤ 1) c a) :
    Real.sqrt (∑ i, ((hornerList D s r c i : ℝ) - idealList (scaledMatrix D s) r a i) ^ 2) ≤
      c.length * (δ + 1) := by
  rw [← l2_eq]
  induction hca with
  | nil =>
    have h0 : (fun i ↦ (hornerList D s r [] i : ℝ) - idealList (scaledMatrix D s) r [] i) = 0 :=
      funext fun i ↦ by simp
    rw [h0]; simp [l2]
  | @cons cj aj c a hj _ ih =>
    rw [hornerList_cons, idealList_cons, List.length_cons]
    refine (step_error hs hD r hj _ _).trans ?_
    push_cast
    linarith

/-- **Rounding error of `hornerFloor`.** With `R = D / s` a contraction, `|cⱼ − aⱼ| ≤ 1` for
`j ≤ J` and `δ ≥ √(card ι)`, the integer column is within `(J + 1)(δ + 1)` (Euclidean norm) of
`(∑_{j ≤ J} aⱼ Rʲ) e_r`. -/
theorem hornerFloor_error {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι) {c : ℕ → ℤ} {a : ℕ → ℝ} (J : ℕ)
    (hca : ∀ j ≤ J, |(c j : ℝ) - a j| ≤ 1) :
    Real.sqrt (∑ i, ((hornerFloor D s c r J i : ℝ) -
        ((∑ j ∈ Finset.range (J + 1), a j • scaledMatrix D s ^ j) *ᵥ Pi.single r 1) i) ^ 2) ≤
      (J + 1) * (δ + 1) := by
  have h2 : List.Forall₂ (fun (cj : ℤ) (aj : ℝ) ↦ |(cj : ℝ) - aj| ≤ 1)
      ((List.range (J + 1)).map c) ((List.range (J + 1)).map a) := by
    rw [List.forall₂_map_left_iff, List.forall₂_map_right_iff, List.forall₂_same]
    intro j hj
    exact hca j (Nat.lt_succ_iff.1 (List.mem_range.1 hj))
  have := hornerList_error hs hD hδ r h2
  rw [idealList_range_map, List.length_map, List.length_range] at this
  simpa [hornerFloor] using this

/-- `hornerFloor_error` with the paper's constant `δ = √(card ι)`. -/
theorem hornerFloor_error_card {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    (r : ι) {c : ℕ → ℤ} {a : ℕ → ℝ} (J : ℕ) (hca : ∀ j ≤ J, |(c j : ℝ) - a j| ≤ 1) :
    Real.sqrt (∑ i, ((hornerFloor D s c r J i : ℝ) -
        ((∑ j ∈ Finset.range (J + 1), a j • scaledMatrix D s ^ j) *ᵥ Pi.single r 1) i) ^ 2) ≤
      (J + 1) * (Real.sqrt (Fintype.card ι) + 1) :=
  hornerFloor_error hs hD le_rfl r J hca

/-- Squared form of `hornerFloor_error`, convenient for integer or `nsq` comparisons. -/
theorem hornerFloor_error_sq {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι) {c : ℕ → ℤ} {a : ℕ → ℝ} (J : ℕ)
    (hca : ∀ j ≤ J, |(c j : ℝ) - a j| ≤ 1) :
    ∑ i, ((hornerFloor D s c r J i : ℝ) -
        ((∑ j ∈ Finset.range (J + 1), a j • scaledMatrix D s ^ j) *ᵥ Pi.single r 1) i) ^ 2 ≤
      ((J + 1) * (δ + 1)) ^ 2 := by
  have h := hornerFloor_error hs hD hδ r J hca
  have hδ0 : 0 ≤ δ := (Real.sqrt_nonneg _).trans hδ
  exact (Real.sqrt_le_left (by positivity)).1 h

/-- A floored coefficient is within `1` of its real value. -/
theorem abs_floor_sub_le_one (x : ℝ) : |((⌊x⌋ : ℤ) : ℝ) - x| ≤ 1 := by
  have h1 := Int.floor_le x
  have h2 := Int.lt_floor_add_one x
  rw [abs_le]; constructor <;> linarith

/-- **The paper's column bound before the Poisson truncation.** With coefficients
`cⱼ = ⌊Q pⱼ⌋` the integer column is within `(J + 1)(δ + 1)` of `Q (∑_{j ≤ J} pⱼ Rʲ) e_r`. -/
theorem hornerFloor_error_floor {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι) (Q : ℝ) (p : ℕ → ℝ) (J : ℕ) :
    Real.sqrt (∑ i, ((hornerFloor D s (fun j ↦ ⌊Q * p j⌋) r J i : ℝ) -
        ((∑ j ∈ Finset.range (J + 1), (Q * p j) • scaledMatrix D s ^ j) *ᵥ
          Pi.single r 1) i) ^ 2) ≤
      (J + 1) * (δ + 1) :=
  hornerFloor_error hs hD hδ r J fun _ _ ↦ abs_floor_sub_le_one _

/-- `hornerFloor_error_floor` for rational `Q` and `pⱼ`, with the floor taken in `ℚ` as the
certificate computes it. -/
theorem hornerFloor_error_floor_rat {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι) (Q : ℚ) (p : ℕ → ℚ) (J : ℕ) :
    Real.sqrt (∑ i, ((hornerFloor D s (fun j ↦ ⌊Q * p j⌋) r J i : ℝ) -
        ((∑ j ∈ Finset.range (J + 1), ((Q * p j : ℚ) : ℝ) • scaledMatrix D s ^ j) *ᵥ
          Pi.single r 1) i) ^ 2) ≤
      (J + 1) * (δ + 1) :=
  hornerFloor_error hs hD hδ r J fun j _ ↦ by
    rw [← Rat.floor_cast (α := ℝ)]; exact abs_floor_sub_le_one _

/-! ### Bridges to complex matrices and `Spectral.nsq` -/

omit [Fintype ι] [DecidableEq ι] in
/-- The complexification of `R = D / s` is `s⁻¹ • D` over `ℂ`. -/
theorem scaledMatrix_map_ofReal (D : Matrix ι ι ℤ) (s : ℕ) :
    (scaledMatrix D s).map ((↑) : ℝ → ℂ) = (s : ℂ)⁻¹ • D.map ((↑) : ℤ → ℂ) := by
  ext i j
  simp [scaledMatrix]

/-- The ideal column, complexified, is the complex polynomial in the complexified `R` applied to
`e_r`; this is the form a matrix-level (Poisson) bound over `ℂ` acts on. -/
theorem ofReal_poly_mulVec_single (R : Matrix ι ι ℝ) (a : ℕ → ℝ) (J : ℕ) (r : ι) :
    (fun i ↦ (((∑ j ∈ Finset.range (J + 1), a j • R ^ j) *ᵥ Pi.single r 1) i : ℂ)) =
      (∑ j ∈ Finset.range (J + 1), (a j : ℂ) • (R.map ((↑) : ℝ → ℂ)) ^ j) *ᵥ
        Pi.single r 1 := by
  funext i
  have h := RingHom.map_mulVec Complex.ofRealHom
    (∑ j ∈ Finset.range (J + 1), a j • R ^ j) (Pi.single r 1) i
  have hmap : (∑ j ∈ Finset.range (J + 1), a j • R ^ j).map Complex.ofRealHom =
      ∑ j ∈ Finset.range (J + 1), (a j : ℂ) • (R.map ((↑) : ℝ → ℂ)) ^ j := by
    rw [← RingHom.mapMatrix_apply, map_sum]
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    have hp : (R ^ j).map Complex.ofRealHom = (R.map Complex.ofRealHom) ^ j :=
      map_pow Complex.ofRealHom.mapMatrix R j
    rw [RingHom.mapMatrix_apply, Matrix.map_smul' _ _ _ fun x y ↦ map_mul _ x y, hp]
    rfl
  have hsingle : (Complex.ofRealHom ∘ Pi.single r (1 : ℝ)) = Pi.single r (1 : ℂ) := by
    funext k; simp [Pi.single_apply, apply_ite]
  rw [hmap, hsingle] at h
  simpa using h

omit [DecidableEq ι] in
/-- The contraction hypothesis follows from its complex form `‖D v‖² ≤ s² ‖v‖²` on `ι → ℂ`,
which is how operator bounds on `Matrix ι ι ℂ` are usually stated. -/
theorem IsScaledContraction.of_complex {D : Matrix ι ι ℤ} {s : ℕ}
    (h : ∀ v : ι → ℂ, Spectral.nsq (D.map ((↑) : ℤ → ℂ) *ᵥ v) ≤ (s : ℝ) ^ 2 * Spectral.nsq v) :
    IsScaledContraction D s := by
  intro v
  have hv := h (fun i ↦ (v i : ℂ))
  have hD : D.map ((↑) : ℤ → ℂ) = (D.map ((↑) : ℤ → ℝ)).map Complex.ofRealHom := by
    ext i j; simp
  have hmul : (D.map ((↑) : ℤ → ℂ) *ᵥ fun i ↦ (v i : ℂ)) =
      fun i ↦ (((D.map ((↑) : ℤ → ℝ)) *ᵥ v) i : ℂ) := by
    funext i
    rw [hD]
    exact (RingHom.map_mulVec Complex.ofRealHom _ v i).symm
  -- `Spectral.nsq_ofReal` at `𝕜 = ℂ`, with the coercion spelled `Complex.ofReal`.
  have hR : ∀ w : ι → ℝ, Spectral.nsq (fun i ↦ (w i : ℂ)) = ∑ i, w i ^ 2 :=
    Spectral.nsq_ofReal
  rwa [hmul, hR, hR] at hv

/-- **Complex form of the column bound**, ready to combine with a matrix bound over `ℂ` by the
triangle inequality for `√nsq`. -/
theorem hornerFloor_error_nsq {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : IsScaledContraction D s)
    {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι) {c : ℕ → ℤ} {a : ℕ → ℝ} (J : ℕ)
    (hca : ∀ j ≤ J, |(c j : ℝ) - a j| ≤ 1) :
    Spectral.nsq (fun i ↦ (hornerFloor D s c r J i : ℂ) -
        ((∑ j ∈ Finset.range (J + 1), (a j : ℂ) • ((s : ℂ)⁻¹ • D.map ((↑) : ℤ → ℂ)) ^ j) *ᵥ
          Pi.single r 1) i) ≤
      ((J + 1) * (δ + 1)) ^ 2 := by
  rw [← scaledMatrix_map_ofReal, ← ofReal_poly_mulVec_single]
  have h := hornerFloor_error_sq hs hD hδ r J hca
  have hR : ∀ w : ι → ℝ, Spectral.nsq (fun i ↦ (w i : ℂ)) = ∑ i, w i ^ 2 :=
    Spectral.nsq_ofReal
  rw [← hR] at h
  convert h using 3
  push_cast; rfl

/-! ### Kernel evaluation -/

/-- The recursion reduces in the kernel: three steps on a `2 × 2` block. The second coordinate
passes through `-21000 / 32 = -657` (floor), not `-656` (truncation). -/
example :
    hornerFloor !![20, -7; -7, 25] 32 (fun j ↦ 1000 * (j + 1)) (0 : Fin 2) 2 = ![3565, -1361] := by
  decide +kernel

end CertifiedGaps.Trace
