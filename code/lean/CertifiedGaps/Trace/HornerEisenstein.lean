/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.InnerProductSpace.PiL2
import CertifiedGaps.Spectral.Basic
import CertifiedGaps.Trace.Eisenstein
import CertifiedGaps.Trace.Horner

/-!
# Horner columns over `ℤ[ω]` with coordinatewise floors

The `b = 21/2` table of the spin-one paper (`thermal120.tex`, §"A rational polynomial and its
error") evaluates `∑_{j ≤ J} pⱼ Wʲ e_r`, `W = B / s`, for a matrix `B` with entries in the
Eisenstein integers `ℤ[ω]` (the ring twisted by the `2π/3` rotation `C`), with integer arithmetic
only: from `z = 0`, for `j = J, J - 1, …, 0`,

  `z ← fl(B z / s) + qⱼ e_r`   (eq. (thermal120-horner)),

where `fl(u + v ω) = ⌊u⌋ + ⌊v⌋ ω` floors both coordinates separately and `qⱼ ∈ ℤ`. This file is
the analogue of `Trace/Horner.lean` (real coordinates) for `ℤ[ω]`, with the error measured in the
complex Euclidean norm after the embedding `toComplex : ℤ[ω] → ℂ`:

* `hornerFloorE B s q r J`: the recursion, a computable definition;
* `IsScaledContractionE B s`: `‖B v‖ ≤ s ‖v‖` for every complex `v`, the hypothesis that the
  complex image `W = s⁻¹ B` is a contraction;
* `hornerFloorE_error_nsq`: if `|qⱼ - aⱼ| ≤ 1` then
  `‖z - (∑_{j ≤ J} aⱼ Wʲ) e_r‖ ≤ (J + 1)(2δ + 1)` whenever `δ ≥ √(card ι)`.

Each step adds the floor error, at most `2` per coordinate (`|δu| + |δv| |ω|`, the paper's
bound) and so at most `2 √(card ι)` in norm, plus `1` from the coefficient; `W` does not enlarge
the error already present. Integer vectors (`v = 0` throughout) are a special case, so the same
recursion serves the untwisted and `P`-twisted rings with the paper's constant `2√d + 1`.
-/

open Matrix
open CertifiedGaps.Spectral (nsq nsq_eq_sum nsq_nonneg)

namespace CertifiedGaps.Trace

open EisensteinInt

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ### The recursion over `ℤ[ω]` -/

/-- One Horner step `z ↦ fl(B z / s) + q e_r` over `ℤ[ω]`. -/
def hornerStepE (B : Matrix ι ι EisensteinInt) (s : ℕ) (r : ι) (q : ℤ)
    (z : ι → EisensteinInt) : ι → EisensteinInt :=
  fun i ↦ fdiv ((B *ᵥ z) i) s + if i = r then (q : EisensteinInt) else 0

/-- The recursion over a coefficient list `[q₀, …, q_J]`: the last coefficient is processed
first, from `z = 0`. -/
def hornerListE (B : Matrix ι ι EisensteinInt) (s : ℕ) (r : ι) (l : List ℤ) :
    ι → EisensteinInt :=
  l.foldr (hornerStepE B s r) 0

/-- **The Horner column over `ℤ[ω]`**: `z ← fl(B z / s) + qⱼ e_r` for `j = J, …, 0`, from
`z = 0`. -/
def hornerFloorE (B : Matrix ι ι EisensteinInt) (s : ℕ) (q : ℕ → ℤ) (r : ι) (J : ℕ) :
    ι → EisensteinInt :=
  hornerListE B s r ((List.range (J + 1)).map q)

@[simp] theorem hornerListE_nil (B : Matrix ι ι EisensteinInt) (s : ℕ) (r : ι) :
    hornerListE B s r [] = 0 := rfl

@[simp] theorem hornerListE_cons (B : Matrix ι ι EisensteinInt) (s : ℕ) (r : ι) (q : ℤ)
    (l : List ℤ) : hornerListE B s r (q :: l) = hornerStepE B s r q (hornerListE B s r l) := rfl

/-- **On an integer matrix the recursion over `ℤ[ω]` is the integer recursion**: the `ω`
coordinate stays zero, and `fl` floors the `1` coordinate as `hornerStep` does. -/
theorem hornerListE_map_intCast (D : Matrix ι ι ℤ) (s : ℕ) (r : ι) (l : List ℤ) :
    hornerListE (D.map (↑)) s r l = fun i ↦ ((hornerList D s r l i : ℤ) : EisensteinInt) := by
  induction l with
  | nil => funext i; ext <;> rfl
  | cons c l ih =>
    funext i
    have hmv : (D.map (↑) *ᵥ fun j ↦ ((hornerList D s r l j : ℤ) : EisensteinInt)) i =
        (((D *ᵥ hornerList D s r l) i : ℤ) : EisensteinInt) :=
      (RingHom.map_mulVec (Int.castRingHom EisensteinInt) D (hornerList D s r l) i).symm
    rw [hornerListE_cons, hornerList_cons, ih, hornerStepE, hmv, hornerStep]
    split_ifs <;> ext <;> simp [fdiv]

/-- `hornerFloorE` of an integer matrix. -/
theorem hornerFloorE_map_intCast (D : Matrix ι ι ℤ) (s : ℕ) (q : ℕ → ℤ) (r : ι) (J : ℕ) :
    hornerFloorE (D.map (↑)) s q r J = fun i ↦ ((hornerFloor D s q r J i : ℤ) : EisensteinInt) :=
  hornerListE_map_intCast D s r _

/-! ### The complex recursion -/

/-- `B` is `s` times a contraction over `ℂ`: `‖B v‖² ≤ s² ‖v‖²` for every complex `v`. -/
def IsScaledContractionE (B : Matrix ι ι EisensteinInt) (s : ℕ) : Prop :=
  ∀ v : ι → ℂ, nsq (B.map toComplex *ᵥ v) ≤ (s : ℝ) ^ 2 * nsq v

/-- The ideal recursion `Y ↦ W Y + aⱼ e_r` over a list of complex coefficients. -/
noncomputable def idealListC (W : Matrix ι ι ℂ) (r : ι) (a : List ℂ) : ι → ℂ :=
  a.foldr (fun aj Y ↦ W *ᵥ Y + Pi.single r aj) 0

@[simp] theorem idealListC_nil (W : Matrix ι ι ℂ) (r : ι) : idealListC W r [] = 0 := rfl

@[simp] theorem idealListC_cons (W : Matrix ι ι ℂ) (r : ι) (a : ℂ) (l : List ℂ) :
    idealListC W r (a :: l) = W *ᵥ idealListC W r l + Pi.single r a := rfl

/-- The ideal recursion is the polynomial `∑_{j ≤ J} aⱼ Wʲ` applied to `e_r` (Horner's rule). -/
theorem idealListC_range_map (W : Matrix ι ι ℂ) (r : ι) (a : ℕ → ℂ) (J : ℕ) :
    idealListC W r ((List.range (J + 1)).map a) =
      (∑ j ∈ Finset.range (J + 1), a j • W ^ j) *ᵥ Pi.single r 1 := by
  induction J generalizing a with
  | zero =>
    ext k
    simp [idealListC, Pi.single_apply, Matrix.one_apply]
  | succ J ih =>
    rw [List.range_succ_eq_map, List.map_cons, List.map_map, idealListC_cons, ih,
      Finset.sum_range_succ' _ (J + 1), add_mulVec, mulVec_mulVec, Finset.mul_sum]
    congr 1
    · simp only [Function.comp_apply, pow_succ', Matrix.mul_smul]
    · ext k
      simp [Pi.single_apply, Matrix.one_apply]

/-! ### Euclidean norm helpers -/

/-- The Euclidean norm of a complex vector, routed through `EuclideanSpace` for its triangle
inequality. -/
private noncomputable abbrev l2 (v : ι → ℂ) : ℝ := ‖WithLp.toLp 2 v‖

omit [DecidableEq ι] in
private theorem l2_sq (v : ι → ℂ) : l2 v ^ 2 = nsq v := by
  rw [l2, EuclideanSpace.norm_eq, Real.sq_sqrt (Finset.sum_nonneg fun _ _ ↦ by positivity),
    nsq_eq_sum]

omit [DecidableEq ι] in
private theorem l2_add (u v : ι → ℂ) : l2 (u + v) ≤ l2 u + l2 v := by
  simpa [l2, WithLp.toLp_add] using norm_add_le (WithLp.toLp 2 u) (WithLp.toLp 2 v)

omit [DecidableEq ι] in
private theorem l2_le_of_nsq_le {v : ι → ℂ} {a : ℝ} (ha : 0 ≤ a) (h : nsq v ≤ a ^ 2) :
    l2 v ≤ a := by
  rw [← l2_sq] at h
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) ha two_ne_zero).1 h

omit [DecidableEq ι] in
/-- A vector with coordinates of modulus at most `2` has norm at most `2 √(card ι)`. -/
private theorem l2_le_two_sqrt_card (u : ι → ℂ) (hu : ∀ i, ‖u i‖ ≤ 2) :
    l2 u ≤ 2 * Real.sqrt (Fintype.card ι) := by
  refine l2_le_of_nsq_le (by positivity) ?_
  rw [nsq_eq_sum, mul_pow, Real.sq_sqrt (Nat.cast_nonneg _)]
  calc ∑ i, ‖u i‖ ^ 2 ≤ ∑ _i : ι, (2 : ℝ) ^ 2 :=
        Finset.sum_le_sum fun i _ ↦ pow_le_pow_left₀ (norm_nonneg _) (hu i) 2
    _ = 2 ^ 2 * Fintype.card ι := by simp [mul_comm]

private theorem l2_single (r : ι) (x : ℂ) : l2 (Pi.single r x) = ‖x‖ := by
  refine le_antisymm (l2_le_of_nsq_le (norm_nonneg _) ?_) ?_
  · rw [nsq_eq_sum, Finset.sum_eq_single r (fun b _ hb ↦ by simp [hb]) (by simp)]
    simp
  · rw [l2, EuclideanSpace.norm_eq]
    refine Real.le_sqrt_of_sq_le ?_
    exact Finset.single_le_sum (f := fun i ↦ ‖(WithLp.toLp 2 (Pi.single r x) :
      EuclideanSpace ℂ ι) i‖ ^ 2) (fun _ _ ↦ by positivity) (Finset.mem_univ r) |>.trans_eq'
      (by simp)

omit [DecidableEq ι] in
private theorem l2_scaled_le {B : Matrix ι ι EisensteinInt} {s : ℕ} (hs : 0 < s)
    (hB : IsScaledContractionE B s) (v : ι → ℂ) :
    l2 (((s : ℂ)⁻¹ • B.map toComplex) *ᵥ v) ≤ l2 v := by
  have hs' : (0 : ℝ) < s := Nat.cast_pos.2 hs
  refine l2_le_of_nsq_le (norm_nonneg _) ?_
  rw [l2_sq, Matrix.smul_mulVec, Spectral.nsq_smul, norm_inv, Complex.norm_natCast, inv_pow]
  calc ((s : ℝ) ^ 2)⁻¹ * nsq (B.map toComplex *ᵥ v) ≤ ((s : ℝ) ^ 2)⁻¹ * ((s : ℝ) ^ 2 * nsq v) :=
        mul_le_mul_of_nonneg_left (hB v) (by positivity)
    _ = nsq v := by field_simp

/-! ### One step and the accumulated error -/

omit [DecidableEq ι] in
/-- The complex image of `B z` is the complex matrix applied to the complex image of `z`. -/
private theorem toComplex_mulVec (B : Matrix ι ι EisensteinInt) (z : ι → EisensteinInt)
    (i : ι) : toComplex ((B *ᵥ z) i) = (B.map toComplex *ᵥ fun j ↦ toComplex (z j)) i :=
  RingHom.map_mulVec toComplex B z i

/-- One floored step moves the error by at most `2 √(card ι)` (coordinatewise floor) plus `1`
(coefficient error), on top of the propagated error, which `W` does not enlarge. -/
private theorem step_error {B : Matrix ι ι EisensteinInt} {s : ℕ} (hs : 0 < s)
    (hB : IsScaledContractionE B s) (r : ι) {q : ℤ} {a : ℂ} (hqa : ‖(q : ℂ) - a‖ ≤ 1)
    (z : ι → EisensteinInt) (Y : ι → ℂ) :
    l2 (fun i ↦ toComplex (hornerStepE B s r q z i) -
        ((((s : ℂ)⁻¹ • B.map toComplex) *ᵥ Y + Pi.single r a : ι → ℂ) i)) ≤
      l2 (fun i ↦ toComplex (z i) - Y i) + 2 * Real.sqrt (Fintype.card ι) + 1 := by
  set W := (s : ℂ)⁻¹ • B.map toComplex with hW
  -- `u` is the floor error of the vector division.
  set u : ι → ℂ := fun i ↦ toComplex (fdiv ((B *ᵥ z) i) s) - toComplex ((B *ᵥ z) i) / s with hu
  have hdecomp : (fun i ↦ toComplex (hornerStepE B s r q z i) -
      ((W *ᵥ Y + Pi.single r a : ι → ℂ) i)) =
      u + W *ᵥ (fun i ↦ toComplex (z i) - Y i) + Pi.single r ((q : ℂ) - a) := by
    funext i
    have hsub : (fun i ↦ toComplex (z i) - Y i) = (fun i ↦ toComplex (z i)) - Y := rfl
    rw [hsub, mulVec_sub]
    simp only [hornerStepE, hu, Pi.add_apply, Pi.sub_apply, map_add, Pi.single_apply, hW,
      Matrix.smul_mulVec, Pi.smul_apply, smul_eq_mul, toComplex_mulVec]
    split_ifs
    · rw [toComplex_intCast]; ring
    · rw [map_zero]; ring
  have h1 := l2_le_two_sqrt_card u fun i ↦ norm_toComplex_fdiv_sub_le _ s
  have h2 := l2_scaled_le hs hB (fun i ↦ toComplex (z i) - Y i)
  have h3 := l2_single r ((q : ℂ) - a)
  have h4 := l2_add (u + W *ᵥ (fun i ↦ toComplex (z i) - Y i)) (Pi.single r ((q : ℂ) - a))
  have h5 := l2_add u (W *ᵥ (fun i ↦ toComplex (z i) - Y i))
  rw [hdecomp]
  linarith

/-! ### Main error bounds -/

/-- **Rounding error of the recursion over `ℤ[ω]` (coefficient lists).** If `W = B / s` is a
contraction and every integer coefficient is within `1` of its complex target, the floored
recursion stays within `(number of steps) · (2δ + 1)` of the ideal recursion, for any
`δ ≥ √(card ι)`. -/
theorem hornerListE_error {B : Matrix ι ι EisensteinInt} {s : ℕ} (hs : 0 < s)
    (hB : IsScaledContractionE B s) {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι)
    {q : List ℤ} {a : List ℂ}
    (hqa : List.Forall₂ (fun (qj : ℤ) (aj : ℂ) ↦ ‖(qj : ℂ) - aj‖ ≤ 1) q a) :
    nsq (fun i ↦ toComplex (hornerListE B s r q i) -
        idealListC ((s : ℂ)⁻¹ • B.map toComplex) r a i) ≤ (q.length * (2 * δ + 1)) ^ 2 := by
  rw [← l2_sq]
  have hδ0 : 0 ≤ δ := (Real.sqrt_nonneg _).trans hδ
  refine pow_le_pow_left₀ (norm_nonneg _) ?_ 2
  change l2 _ ≤ _
  induction hqa with
  | nil =>
    have h0 : (fun i ↦ toComplex (hornerListE B s r [] i) -
        idealListC ((s : ℂ)⁻¹ • B.map toComplex) r [] i) = 0 :=
      funext fun i ↦ by simp
    rw [h0]; simp [l2]
  | @cons qj aj q a hj _ ih =>
    rw [hornerListE_cons, idealListC_cons, List.length_cons]
    refine (step_error hs hB r hj _ _).trans ?_
    push_cast
    nlinarith

/-- **Rounding error of `hornerFloorE`** (eq. (thermal120-column-error) before the Poisson
truncation). With `W = B / s` a contraction over `ℂ`, `|qⱼ - aⱼ| ≤ 1` for `j ≤ J` and
`δ ≥ √(card ι)`, the column over `ℤ[ω]`, embedded in `ℂ`, is within `(J + 1)(2δ + 1)` of
`(∑_{j ≤ J} aⱼ Wʲ) e_r`. -/
theorem hornerFloorE_error_nsq {B : Matrix ι ι EisensteinInt} {s : ℕ} (hs : 0 < s)
    (hB : IsScaledContractionE B s) {δ : ℝ} (hδ : Real.sqrt (Fintype.card ι) ≤ δ) (r : ι)
    {q : ℕ → ℤ} {a : ℕ → ℝ} (J : ℕ) (hqa : ∀ j ≤ J, |(q j : ℝ) - a j| ≤ 1) :
    nsq (fun i ↦ toComplex (hornerFloorE B s q r J i) -
        ((∑ j ∈ Finset.range (J + 1), (a j : ℂ) • ((s : ℂ)⁻¹ • B.map toComplex) ^ j) *ᵥ
          Pi.single r 1) i) ≤
      ((J + 1) * (2 * δ + 1)) ^ 2 := by
  have h2 : List.Forall₂ (fun (qj : ℤ) (aj : ℂ) ↦ ‖(qj : ℂ) - aj‖ ≤ 1)
      ((List.range (J + 1)).map q) ((List.range (J + 1)).map fun j ↦ (a j : ℂ)) := by
    rw [List.forall₂_map_left_iff, List.forall₂_map_right_iff, List.forall₂_same]
    intro j hj
    have h := hqa j (Nat.lt_succ_iff.1 (List.mem_range.1 hj))
    rwa [show (q j : ℂ) - (a j : ℂ) = (((q j : ℝ) - a j : ℝ) : ℂ) by push_cast; rfl,
      Complex.norm_real, Real.norm_eq_abs]
  have := hornerListE_error hs hB hδ r h2
  rw [idealListC_range_map, List.length_map, List.length_range] at this
  simpa [hornerFloorE] using this

/-! ### Kernel evaluation -/

/-- The recursion reduces in the kernel: two steps on a `1 × 1` block with the seam entry
`-ω` and coefficients `q₀ = 1000`, `q₁ = 7`. The first step gives `7`; the second floors
`-7ω / 16` to `-ω` (floor, not truncation) and adds `1000`. -/
example :
    hornerFloorE !![-EisensteinInt.omega] 16 (fun j ↦ 1000 * (j + 1) - 1993 * j) (0 : Fin 1) 1 =
      ![⟨1000, -1⟩] := by
  decide +kernel

end CertifiedGaps.Trace
