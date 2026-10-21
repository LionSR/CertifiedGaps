/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Hamiltonian
import CertifiedGaps.Spectral.Basic
import Mathlib.Analysis.InnerProductSpace.PiL2

/-!
# Gibbs traces from representative columns

The thermal seeds of the spin-1 certificate (OpenAI, `thermal72.tex`, Sections "Representative
columns" and "Exact totals and rational enclosures") never form the Gibbs operator. They write

  `Z = Tr e^{-bH} = e^{bc} ‖F‖_F²`,  `F = exp (-(b/2) (H + c I))`,

compute one integer column `X_ρ ≈ Q F e_ρ` per symmetry orbit of basis vectors, and turn the
integer `W = ∑_ρ ν_ρ ‖X_ρ‖²` into a rational enclosure of `Z`. This file proves that chain for
an arbitrary Hermitian matrix on an arbitrary finite index type, so that each concrete
enclosure needs only the per-column error `E` and the integer `W` as data.

* `gibbsHalf H b c = exp (-(b/2) • (H + c • 1))` is the half-step factor `F`.
* `gibbsTrace_eq_exp_mul_sum_nsq`: `Z = e^{bc} ∑_r ‖F e_r‖²`, from `Fᴴ F = e^{-bc} e^{-bH}` and
  `Spectral.sum_nsq_mulVec_single` (`∑_r ‖M e_r‖² = Re Tr (Mᴴ M)`);
  `sum_nsq_gibbsHalf_eq_sum`: `∑_r ‖F e_r‖² = ∑_i e^{-b(λᵢ + c)}`; a scalar shift `H + a`
  multiplies the trace by `e^{-ba}` (`gibbsTrace_add_smul_one`) and adds to `c`
  (`gibbsHalf_add_smul_one`).
* `nsq_gibbsHalf_single_eq`: a unitary `U` that commutes with `H` and maps `e_r` to a phase
  times `e_{r'}` gives the columns `r` and `r'` of `F` equal norms.
* `abs_sqrt_wsum_sub_le`: Minkowski's inequality in the weighted `ℓ²` direct sum,
  `|√(∑ ν ‖X‖²) - √(∑ ν ‖Y‖²)| ≤ √(∑ ν) E` when every `‖X_ρ - Y_ρ‖ ≤ E`; and the triangle
  inequality `nsq_add_le_sq` for squared norms.
* `gibbsTrace_enclosure`: the rational enclosure
  `s max(0, r - ℰ)² / Q² ≤ Z ≤ (s + t)(r + 1 + ℰ)² / Q²`, with the variants
  `gibbsTrace_enclosure_of_rep` (representatives from a map `rep`, fiber sizes `ν`, `∑ ν = dim`;
  the regrouping is Mathlib's `Finset.sum_comp`) and
  `gibbsTrace_enclosure_of_intColumns` (integer columns, `r = Nat.sqrt W`).

The paper states `|√W - Q ‖F‖_F| ≤ √(3^n) E < ℰ`; only `≤ ℰ` is used, so the hypotheses below
ask for `(∑ ν) E² ≤ ℰ²`, a rational inequality.
-/

open Matrix NormedSpace
open scoped Finset
open CertifiedGaps.Spectral (nsq nsq_eq_sum nsq_nonneg nsq_smul nsq_intCast nsq_mulVec_of_unitary
  sum_nsq_mulVec_single)

namespace CertifiedGaps.Trace

/-! ### The half-step Gibbs factor -/

section Gibbs

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The half-step Gibbs factor `F = exp (-(b/2) (H + c I))` of the shifted Hamiltonian, whose
squared Frobenius norm is `e^{-bc} Tr e^{-bH}`. -/
noncomputable def gibbsHalf (H : Matrix n n ℂ) (b c : ℝ) : Matrix n n ℂ :=
  exp (((-(b / 2) : ℝ) : ℂ) • (H + (c : ℂ) • 1))

/-- `Fᴴ F = e^{-bc} e^{-bH}`: the two half steps recombine into one full Gibbs step. -/
theorem conjTranspose_gibbsHalf_mul_self {H : Matrix n n ℂ} (hH : H.IsHermitian) (b c : ℝ) :
    (gibbsHalf H b c)ᴴ * gibbsHalf H b c =
      (Real.exp (-(b * c)) : ℂ) • exp ((-(b : ℂ)) • H) := by
  have hA : (H + (c : ℂ) • 1).IsHermitian :=
    hH.add (isHermitian_one.smul (Complex.conj_ofReal c))
  have hF : (gibbsHalf H b c).IsHermitian := (hA.smul (Complex.conj_ofReal (-(b / 2)))).exp
  rw [hF.eq, gibbsHalf, ← exp_add_of_commute _ _ (Commute.refl _), ← add_smul]
  -- `-(b/2) - (b/2) = -b` splits into a scalar part and the full step `-bH`.
  have hsplit : (((-(b / 2) : ℝ) : ℂ) + ((-(b / 2) : ℝ) : ℂ)) • (H + (c : ℂ) • 1) =
      (-(b : ℂ)) • H + ((-(b * c) : ℝ) : ℂ) • (1 : Matrix n n ℂ) := by
    rw [smul_add, smul_smul]; push_cast; congr 1 <;> [congr 1; congr 1] <;> ring
  rw [hsplit, exp_add_of_commute _ _ ((Commute.one_right _).smul_right _), exp_smul_one,
    Matrix.mul_smul, Matrix.mul_one, ← Complex.exp_eq_exp_ℂ, Complex.ofReal_exp]

/-- **Frobenius form of the Gibbs trace.** `Z = Tr e^{-bH} = e^{bc} ∑_r ‖F e_r‖²`. -/
theorem gibbsTrace_eq_exp_mul_sum_nsq {H : Matrix n n ℂ} (hH : H.IsHermitian) (b c : ℝ) :
    gibbsTrace H b = Real.exp (b * c) * ∑ r, nsq (gibbsHalf H b c *ᵥ Pi.single r 1) := by
  rw [sum_nsq_mulVec_single, conjTranspose_gibbsHalf_mul_self hH, trace_smul, smul_eq_mul,
    RCLike.re_to_complex, Complex.re_ofReal_mul, gibbsTrace, ← mul_assoc, ← Real.exp_add]
  simp

/-- The squared Frobenius norm of `F` in terms of the spectrum:
`∑_r ‖F e_r‖² = ∑ᵢ e^{-b(λᵢ + c)}`. -/
theorem sum_nsq_gibbsHalf_eq_sum {H : Matrix n n ℂ} (hH : H.IsHermitian) (b c : ℝ) :
    ∑ r, nsq (gibbsHalf H b c *ᵥ Pi.single r 1) =
      ∑ i, Real.exp (-(b * (hH.eigenvalues i + c))) := by
  have h := gibbsTrace_eq_exp_mul_sum_nsq hH b c
  rw [gibbsTrace_eq_sum hH] at h
  rw [eq_comm, ← mul_right_inj' (Real.exp_pos (b * c)).ne', ← h, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [← Real.exp_add]; ring_nf

/-- A scalar shift `H + a` multiplies the Gibbs trace by `e^{-ba}`. -/
theorem gibbsTrace_add_smul_one (H : Matrix n n ℂ) (a b : ℝ) :
    gibbsTrace (H + (a : ℂ) • 1) b = Real.exp (-(b * a)) * gibbsTrace H b := by
  have hsplit : (-(b : ℂ)) • (H + (a : ℂ) • 1) =
      (-(b : ℂ)) • H + ((-(b * a) : ℝ) : ℂ) • (1 : Matrix n n ℂ) := by
    rw [smul_add, smul_smul]; push_cast; ring_nf
  rw [gibbsTrace, hsplit, exp_add_of_commute _ _ ((Commute.one_right _).smul_right _),
    exp_smul_one, Matrix.mul_smul, Matrix.mul_one, trace_smul, smul_eq_mul,
    ← Complex.exp_eq_exp_ℂ, ← Complex.ofReal_exp, Complex.re_ofReal_mul, gibbsTrace]

/-- The half-step factor of a shifted Hamiltonian: the shifts add up. -/
theorem gibbsHalf_add_smul_one (H : Matrix n n ℂ) (a b c : ℝ) :
    gibbsHalf (H + (a : ℂ) • 1) b c = gibbsHalf H b (a + c) := by
  rw [gibbsHalf, gibbsHalf, add_assoc, ← add_smul, Complex.ofReal_add]

end Gibbs

/-! ### Symmetries -/

section Symmetry

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- A matrix commuting with `H` commutes with the half-step factor `F`. -/
theorem commute_gibbsHalf {U H : Matrix n n ℂ} (hUH : Commute U H) (b c : ℝ) :
    Commute U (gibbsHalf H b c) :=
  ((hUH.add_right ((Commute.one_right U).smul_right _)).smul_right _).exp_right

/-- A unitary symmetry `U` of `H` can be moved through `F` without changing norms:
`‖F U v‖ = ‖U F v‖ = ‖F v‖`. -/
theorem nsq_gibbsHalf_mulVec_unitary {U H : Matrix n n ℂ} (hU : Uᴴ * U = 1)
    (hUH : Commute U H) (b c : ℝ) (v : n → ℂ) :
    nsq (gibbsHalf H b c *ᵥ (U *ᵥ v)) = nsq (gibbsHalf H b c *ᵥ v) := by
  rw [mulVec_mulVec, ← (commute_gibbsHalf hUH b c).eq, ← mulVec_mulVec,
    nsq_mulVec_of_unitary hU]

/-- **Signed and phase permutations.** If a unitary symmetry `U` of `H` maps the basis vector
`e_r` to `u e_{r'}` with `|u| = 1`, the columns `r` and `r'` of `F` have equal norms. -/
theorem nsq_gibbsHalf_single_eq {U H : Matrix n n ℂ} (hU : Uᴴ * U = 1) (hUH : Commute U H)
    (b c : ℝ) {r r' : n} {u : ℂ} (hu : ‖u‖ = 1)
    (hUr : U *ᵥ Pi.single r 1 = u • Pi.single r' 1) :
    nsq (gibbsHalf H b c *ᵥ Pi.single r' 1) = nsq (gibbsHalf H b c *ᵥ Pi.single r 1) := by
  rw [← nsq_gibbsHalf_mulVec_unitary hU hUH b c (Pi.single r 1), hUr, mulVec_smul, nsq_smul,
    hu, one_pow, one_mul]

end Symmetry

/-! ### Minkowski's inequality in the weighted direct sum -/

section Minkowski

variable {ι m : Type*} [Fintype m]

/-- The family `(Z_ρ)_{ρ ∈ R}` as one vector of the direct sum that repeats `Z_ρ` with weight
`w_ρ`: the coordinate `(ρ, i)` is `√w_ρ · Z_ρ(i)`. -/
private noncomputable def wEmb (R : Finset ι) (w : ι → ℝ) (Z : ι → m → ℂ) :
    EuclideanSpace ℂ (R × m) :=
  WithLp.toLp 2 fun p ↦ (√(w p.1) : ℂ) * Z p.1 p.2

private theorem norm_wEmb (R : Finset ι) {w : ι → ℝ} (hw : ∀ ρ, 0 ≤ w ρ) (Z : ι → m → ℂ) :
    ‖wEmb R w Z‖ = √(∑ ρ ∈ R, w ρ * nsq (Z ρ)) := by
  rw [EuclideanSpace.norm_eq]
  congr 1
  simp only [wEmb, Fintype.sum_prod_type, norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (Real.sqrt_nonneg _), mul_pow, Real.sq_sqrt (hw _), ← Finset.mul_sum,
    nsq_eq_sum]
  exact Finset.sum_coe_sort R fun ρ ↦ w ρ * ∑ i, ‖Z ρ i‖ ^ 2

omit [Fintype m] in
private theorem wEmb_sub (R : Finset ι) (w : ι → ℝ) (X Y : ι → m → ℂ) :
    wEmb R w X - wEmb R w Y = wEmb R w (X - Y) := by
  ext p
  simp [wEmb, mul_sub]

/-- **Minkowski in the weighted direct sum.** If every column error is at most `E`, the
weighted `ℓ²` norms of the two families differ by at most `√(∑ w) E`. With integer weights
`w_ρ = ν_ρ` this is the triangle inequality in the direct sum that repeats each column `ν_ρ`
times. -/
theorem abs_sqrt_wsum_sub_le (R : Finset ι) {w : ι → ℝ} (hw : ∀ ρ, 0 ≤ w ρ)
    (X Y : ι → m → ℂ) {E : ℝ} (hE : 0 ≤ E) (hXY : ∀ ρ ∈ R, nsq (X ρ - Y ρ) ≤ E ^ 2) :
    |√(∑ ρ ∈ R, w ρ * nsq (X ρ)) - √(∑ ρ ∈ R, w ρ * nsq (Y ρ))| ≤ √(∑ ρ ∈ R, w ρ) * E := by
  rw [← norm_wEmb R hw, ← norm_wEmb R hw]
  refine (abs_norm_sub_norm_le _ _).trans ?_
  rw [wEmb_sub, norm_wEmb R hw, ← Real.sqrt_sq hE, ← Real.sqrt_mul (Finset.sum_nonneg
    fun ρ _ ↦ hw ρ), Finset.sum_mul]
  exact Real.sqrt_le_sqrt (Finset.sum_le_sum fun ρ hρ ↦
    mul_le_mul_of_nonneg_left (hXY ρ hρ) (hw ρ))

/-- **The triangle inequality for squared norms**: `‖x‖ ≤ a` and `‖y‖ ≤ b` give
`‖x + y‖ ≤ a + b`. Used to add the rounding error of an integer column to the truncation error
of the polynomial it approximates. -/
theorem nsq_add_le_sq [Fintype ι] {x y : ι → ℂ} {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hx : nsq x ≤ a ^ 2) (hy : nsq y ≤ b ^ 2) : nsq (x + y) ≤ (a + b) ^ 2 := by
  have hn : ∀ z : ι → ℂ, nsq z = ‖WithLp.toLp 2 z‖ ^ 2 := fun z ↦ by
    rw [EuclideanSpace.norm_eq, Real.sq_sqrt (Finset.sum_nonneg fun _ _ ↦ by positivity),
      nsq_eq_sum]
  rw [hn] at hx hy ⊢
  have hx' := (pow_le_pow_iff_left₀ (norm_nonneg _) ha two_ne_zero).1 hx
  have hy' := (pow_le_pow_iff_left₀ (norm_nonneg _) hb two_ne_zero).1 hy
  rw [WithLp.toLp_add]
  exact pow_le_pow_left₀ (norm_nonneg _) ((norm_add_le _ _).trans (add_le_add hx' hy')) 2

end Minkowski

/-! ### The rational enclosure -/

section Enclosure

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Rational enclosure of a Gibbs trace from representative columns.**

Data: representatives `R` with weights `ν` that reproduce the squared Frobenius norm of
`F = gibbsHalf H b c` (`horbit`), computed columns `X_ρ` within `E` of `Q F e_ρ` (`hcol`),
their weighted total `W` (`hW`), an integer `r` with `r² ≤ W < (r + 1)²`, a bound
`(∑ ν) E² ≤ ℰ²`, and `s ≤ e^{bc} ≤ s + t`. Conclusion:
`s max(0, r - ℰ)² / Q² ≤ Tr e^{-bH} ≤ (s + t)(r + 1 + ℰ)² / Q²`. -/
theorem gibbsTrace_enclosure {H : Matrix n n ℂ} (hH : H.IsHermitian) (b c : ℝ)
    (R : Finset n) (ν : n → ℕ)
    (horbit : ∑ r, nsq (gibbsHalf H b c *ᵥ Pi.single r 1) =
      ∑ ρ ∈ R, (ν ρ : ℝ) * nsq (gibbsHalf H b c *ᵥ Pi.single ρ 1))
    (X : n → n → ℂ) {Q E ℰ : ℝ} (hQ : 0 < Q) (hE : 0 ≤ E)
    (hcol : ∀ ρ ∈ R, nsq (X ρ - (Q : ℂ) • (gibbsHalf H b c *ᵥ Pi.single ρ 1)) ≤ E ^ 2)
    (hℰ0 : 0 ≤ ℰ) (hℰ : (∑ ρ ∈ R, (ν ρ : ℝ)) * E ^ 2 ≤ ℰ ^ 2)
    {W r : ℕ} (hW : ∑ ρ ∈ R, (ν ρ : ℝ) * nsq (X ρ) = W) (hrW : r ^ 2 ≤ W)
    (hWr : W < (r + 1) ^ 2) {s t : ℝ} (hs : s ≤ Real.exp (b * c))
    (ht : Real.exp (b * c) ≤ s + t) :
    s * max 0 ((r : ℝ) - ℰ) ^ 2 / Q ^ 2 ≤ gibbsTrace H b ∧
      gibbsTrace H b ≤ (s + t) * ((r : ℝ) + 1 + ℰ) ^ 2 / Q ^ 2 := by
  set F := gibbsHalf H b c
  set N2 := ∑ r, nsq (F *ᵥ Pi.single r 1) with hN2
  have hN2nn : 0 ≤ N2 := Finset.sum_nonneg fun _ _ ↦ nsq_nonneg _
  have hZ : gibbsTrace H b = Real.exp (b * c) * N2 := gibbsTrace_eq_exp_mul_sum_nsq hH b c
  -- The exact direct sum has squared norm `Q² ‖F‖_F²`, the computed one `W`.
  have hY : ∑ ρ ∈ R, (ν ρ : ℝ) * nsq ((Q : ℂ) • (F *ᵥ Pi.single ρ 1)) = (Q * √N2) ^ 2 := by
    simp_rw [nsq_smul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hQ, mul_pow,
      Real.sq_sqrt hN2nn, horbit, Finset.mul_sum]
    exact Finset.sum_congr rfl fun ρ _ ↦ by ring
  have hmink := abs_sqrt_wsum_sub_le R (fun ρ ↦ Nat.cast_nonneg (ν ρ)) X
    (fun ρ ↦ (Q : ℂ) • (F *ᵥ Pi.single ρ 1)) hE hcol
  have hQN : 0 ≤ Q * √N2 := mul_nonneg hQ.le (Real.sqrt_nonneg _)
  rw [hW, hY, Real.sqrt_sq hQN] at hmink
  have hℰ' : √(∑ ρ ∈ R, (ν ρ : ℝ)) * E ≤ ℰ := by
    rw [← Real.sqrt_sq hE, ← Real.sqrt_mul (Finset.sum_nonneg fun ρ _ ↦ Nat.cast_nonneg _),
      ← Real.sqrt_sq hℰ0]
    exact Real.sqrt_le_sqrt hℰ
  have hdiff := (abs_le.1 (hmink.trans hℰ'))
  -- `r ≤ √W < r + 1`.
  have hr1 : (r : ℝ) ≤ √(W : ℝ) := Real.le_sqrt_of_sq_le (by exact_mod_cast hrW)
  have hr2 : √(W : ℝ) < r + 1 :=
    (Real.sqrt_lt' (by positivity)).2 (by exact_mod_cast hWr)
  have hQ2 : 0 < Q ^ 2 := by positivity
  have hsq : (Q * √N2) ^ 2 = Q ^ 2 * N2 := by rw [mul_pow, Real.sq_sqrt hN2nn]
  constructor
  · -- Lower bound: `Q ‖F‖_F ≥ max(0, r - ℰ)`.
    have hm : max 0 ((r : ℝ) - ℰ) ≤ Q * √N2 := max_le hQN (by linarith)
    have hm2 : max 0 ((r : ℝ) - ℰ) ^ 2 / Q ^ 2 ≤ N2 := by
      rw [div_le_iff₀ hQ2, mul_comm, ← hsq]
      exact pow_le_pow_left₀ (le_max_left _ _) hm 2
    rw [hZ, mul_div_assoc]
    exact (mul_le_mul_of_nonneg_right hs (by positivity)).trans
      (mul_le_mul_of_nonneg_left hm2 (Real.exp_pos _).le)
  · -- Upper bound: `Q ‖F‖_F ≤ √W + ℰ < r + 1 + ℰ`.
    have hM : Q * √N2 ≤ (r : ℝ) + 1 + ℰ := by linarith
    have hM2 : N2 ≤ ((r : ℝ) + 1 + ℰ) ^ 2 / Q ^ 2 := by
      rw [le_div_iff₀ hQ2, mul_comm, ← hsq]
      exact pow_le_pow_left₀ hQN hM 2
    rw [hZ, mul_div_assoc]
    exact (mul_le_mul_of_nonneg_right ht hN2nn).trans
      (mul_le_mul_of_nonneg_left hM2 ((Real.exp_pos _).le.trans ht))

/-- `gibbsTrace_enclosure` with the representatives produced by a map `rep` that is constant
on column norms: `ν_ρ` is the size of the fiber of `ρ`, so `∑ ν = |n|` and the error budget
reads `|n| E² ≤ ℰ²`. -/
theorem gibbsTrace_enclosure_of_rep {H : Matrix n n ℂ} (hH : H.IsHermitian) (b c : ℝ)
    (rep : n → n)
    (hrep : ∀ r, nsq (gibbsHalf H b c *ᵥ Pi.single r 1) =
      nsq (gibbsHalf H b c *ᵥ Pi.single (rep r) 1))
    (X : n → n → ℂ) {Q E ℰ : ℝ} (hQ : 0 < Q) (hE : 0 ≤ E)
    (hcol : ∀ ρ ∈ Finset.univ.image rep,
      nsq (X ρ - (Q : ℂ) • (gibbsHalf H b c *ᵥ Pi.single ρ 1)) ≤ E ^ 2)
    (hℰ0 : 0 ≤ ℰ) (hℰ : (Fintype.card n : ℝ) * E ^ 2 ≤ ℰ ^ 2) {W r : ℕ}
    (hW : ∑ ρ ∈ Finset.univ.image rep, (#{x | rep x = ρ} : ℝ) * nsq (X ρ) = W)
    (hrW : r ^ 2 ≤ W) (hWr : W < (r + 1) ^ 2) {s t : ℝ} (hs : s ≤ Real.exp (b * c))
    (ht : Real.exp (b * c) ≤ s + t) :
    s * max 0 ((r : ℝ) - ℰ) ^ 2 / Q ^ 2 ≤ gibbsTrace H b ∧
      gibbsTrace H b ≤ (s + t) * ((r : ℝ) + 1 + ℰ) ^ 2 / Q ^ 2 := by
  -- Group the columns by their representative (`Finset.sum_comp`); the fibers partition `n`.
  have horbit : ∑ r, nsq (gibbsHalf H b c *ᵥ Pi.single r 1) = ∑ ρ ∈ Finset.univ.image rep,
      (#{x | rep x = ρ} : ℝ) * nsq (gibbsHalf H b c *ᵥ Pi.single ρ 1) := by
    rw [Finset.sum_congr rfl fun r _ ↦ hrep r,
      Finset.sum_comp (fun ρ ↦ nsq (gibbsHalf H b c *ᵥ Pi.single ρ 1)) rep]
    simp_rw [nsmul_eq_mul]
  refine gibbsTrace_enclosure hH b c (Finset.univ.image rep) (fun ρ ↦ #{x | rep x = ρ})
    horbit X hQ hE hcol hℰ0 ?_ hW hrW hWr hs ht
  rw [← Nat.cast_sum, ← Finset.card_eq_sum_card_image, Finset.card_univ]
  exact hℰ

/-- `gibbsTrace_enclosure` for integer columns `X_ρ ∈ ℤⁿ` and `r = ⌊√W⌋`, where the integer
total is `W = ∑_ρ ν_ρ ∑_w X_ρ(w)²`. -/
theorem gibbsTrace_enclosure_of_intColumns {H : Matrix n n ℂ} (hH : H.IsHermitian) (b c : ℝ)
    (R : Finset n) (ν : n → ℕ)
    (horbit : ∑ r, nsq (gibbsHalf H b c *ᵥ Pi.single r 1) =
      ∑ ρ ∈ R, (ν ρ : ℝ) * nsq (gibbsHalf H b c *ᵥ Pi.single ρ 1))
    (X : n → n → ℤ) {Q E ℰ : ℝ} (hQ : 0 < Q) (hE : 0 ≤ E)
    (hcol : ∀ ρ ∈ R,
      nsq ((fun w ↦ (X ρ w : ℂ)) - (Q : ℂ) • (gibbsHalf H b c *ᵥ Pi.single ρ 1)) ≤ E ^ 2)
    (hℰ0 : 0 ≤ ℰ) (hℰ : (∑ ρ ∈ R, (ν ρ : ℝ)) * E ^ 2 ≤ ℰ ^ 2) {W : ℕ}
    (hW : ∑ ρ ∈ R, (ν ρ : ℤ) * ∑ w, X ρ w ^ 2 = W) {s t : ℝ} (hs : s ≤ Real.exp (b * c))
    (ht : Real.exp (b * c) ≤ s + t) :
    s * max 0 ((Nat.sqrt W : ℝ) - ℰ) ^ 2 / Q ^ 2 ≤ gibbsTrace H b ∧
      gibbsTrace H b ≤ (s + t) * ((Nat.sqrt W : ℝ) + 1 + ℰ) ^ 2 / Q ^ 2 := by
  refine gibbsTrace_enclosure hH b c R ν horbit (fun ρ w ↦ (X ρ w : ℂ)) hQ hE hcol hℰ0 hℰ
    ?_ (Nat.sqrt_le' W) (Nat.lt_succ_sqrt' W) hs ht
  simp_rw [nsq_intCast]
  exact_mod_cast hW

end Enclosure

end CertifiedGaps.Trace
