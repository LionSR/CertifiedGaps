/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Trace.Frobenius
import CertifiedGaps.Trace.Poisson
import CertifiedGaps.Trace.WeightedHorner

/-!
# Gibbs traces from integer Horner columns of a similar matrix

The thermal appendix of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026)
encloses a Gibbs trace `Z = Tr e^{-bH}` through the half-step factor

  `F = gibbsHalf H b c = exp (m (A - 1))`,  `A` a Hermitian contraction,

approximated by integer Horner columns of the integer matrix `D = s J A J⁻¹`, where `J` is a
positive diagonal matrix. Column `r` of the integer recursion approximates `Q J P(A) J⁻¹ e_r`,
with `P` the conditioned Poisson polynomial; undoing the similarity gives the matrix
`Ê = J⁻¹ X J / Q`, whose column `r` is `(j_r / j_w) X_r(w) / Q`. This file proves, for any
finite index type,

* `nsq_column_sub_le`: `‖(j_r/j_w · X_r(w))_w - Q F e_r‖ ≤ j_r (K + 1)(δ + γ) + Q η(m, K)`,
  from the weighted Horner bound (`hornerFloor_error_weighted`) and the Poisson estimate
  (`nsq_exp_sub_poisson_mulVec_le_of_qform`);
* `gibbsTrace_enclosure_weighted`: with `W = ∑_{r,w} (j_r/j_w)² X_r(w)²` and
  `ε ≥ √|ι| (j_max (K + 1)(δ + γ) + Q η)`,

    `e^{bc} max(0, √W - ε)² / Q² ≤ Tr e^{-bH} ≤ e^{bc} (√W + ε)² / Q²`.

The paper bounds the rounding error after undoing the similarity, through
`‖J⁻¹‖ ‖J‖ ≤ √2`; here the rounding error is bounded in the weighted norm `‖J⁻¹ ·‖`, in which
`D / s` is a contraction, and the similarity is undone once per column (the factor `j_r`).
-/

open Matrix NormedSpace
open CertifiedGaps.Spectral (nsq nsq_eq_sum nsq_nonneg nsq_smul qform)

namespace CertifiedGaps.Trace

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Powers of a similarity: `(J A J⁻¹)ᵏ = J Aᵏ J⁻¹`. -/
theorem pow_similar {J Ji A : Matrix ι ι ℂ} (hJ : Ji * J = 1) (hJ' : J * Ji = 1) (k : ℕ) :
    (J * A * Ji) ^ k = J * A ^ k * Ji := by
  induction k with
  | zero => rw [pow_zero, pow_zero, Matrix.mul_one, hJ']
  | succ k ih =>
    rw [pow_succ, ih]
    calc J * A ^ k * Ji * (J * A * Ji) = J * A ^ k * (Ji * J) * A * Ji := by
          simp only [Matrix.mul_assoc]
      _ = J * A ^ (k + 1) * Ji := by
          rw [hJ, Matrix.mul_one, pow_succ, Matrix.mul_assoc J (A ^ k) A]

/-- The complexified ideal Horner column of `R = D / s`, for `D = s J A J⁻¹`, is
`J (∑ aₖ Aᵏ) J⁻¹ e_r`. -/
theorem ideal_column_eq {A : Matrix ι ι ℂ} {j : ι → ℝ} (hj : ∀ i, 0 < j i) {D : Matrix ι ι ℤ}
    {s : ℕ} (hs : 0 < s)
    (hD : D.map ((↑) : ℤ → ℂ) = (s : ℂ) • (diagonal (fun i ↦ (j i : ℂ)) * A *
      diagonal (fun i ↦ ((j i)⁻¹ : ℂ)))) (a : ℕ → ℝ) (K : ℕ) (r w : ι) :
    ((((∑ k ∈ Finset.range (K + 1), a k • scaledMatrix D s ^ k) *ᵥ Pi.single r 1) w : ℝ) : ℂ)
      = (j w : ℂ) * ((j r)⁻¹ : ℂ) *
        ((∑ k ∈ Finset.range (K + 1), (a k : ℂ) • A ^ k) *ᵥ Pi.single r 1) w := by
  have hcol := congrFun (ofReal_poly_mulVec_single (scaledMatrix D s) a K r) w
  rw [hcol, scaledMatrix_map_ofReal, hD, smul_smul, inv_mul_cancel₀ (by exact_mod_cast hs.ne'),
    one_smul]
  set J : Matrix ι ι ℂ := diagonal fun i ↦ (j i : ℂ)
  set Ji : Matrix ι ι ℂ := diagonal fun i ↦ ((j i)⁻¹ : ℂ)
  have hJ : Ji * J = 1 := by
    rw [diagonal_mul_diagonal, ← diagonal_one]
    congr 1; funext i
    have : (j i : ℂ) ≠ 0 := by exact_mod_cast (hj i).ne'
    field_simp
  have hJ' : J * Ji = 1 := mul_eq_one_comm.mp hJ
  have hpow : ∀ k, (J * A * Ji) ^ k = J * A ^ k * Ji := pow_similar hJ hJ'
  simp_rw [hpow]
  have hsum : ∑ k ∈ Finset.range (K + 1), (a k : ℂ) • (J * A ^ k * Ji) =
      J * (∑ k ∈ Finset.range (K + 1), (a k : ℂ) • A ^ k) * Ji := by
    rw [Matrix.mul_sum, Matrix.sum_mul]
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [Matrix.mul_smul, Matrix.smul_mul]
  rw [hsum, ← mulVec_mulVec, ← mulVec_mulVec]
  have hJi : Ji *ᵥ Pi.single r 1 = ((j r)⁻¹ : ℂ) • Pi.single r 1 := by
    funext i
    rw [mulVec_diagonal]
    by_cases h : i = r
    · subst h; simp
    · simp [h]
  rw [hJi, mulVec_smul, mulVec_diagonal, Pi.smul_apply, smul_eq_mul]
  ring

/-- **One column.** With `F = exp (m (A - 1))` for a Hermitian contraction `A`,
`D = s J A J⁻¹`, integer coefficients within `1` of `Q pₖ` (the conditioned Poisson weights of
mean `m`), `δ ≥ √(∑ j⁻²)` and `γ ≥ j_r⁻¹`, the rescaled integer column
`(j_r / j_w) X_r(w)` is within `j_r (K + 1)(δ + γ) + Q η(m, K)` of `Q F e_r`. -/
theorem nsq_column_sub_le {A : Matrix ι ι ℂ} (hA : A.IsHermitian)
    (hcon : ∀ x, |qform A x| ≤ nsq x) {j : ι → ℝ} (hj : ∀ i, 0 < j i) {D : Matrix ι ι ℤ}
    {s : ℕ} (hs : 0 < s)
    (hD : D.map ((↑) : ℤ → ℂ) = (s : ℂ) • (diagonal (fun i ↦ (j i : ℂ)) * A *
      diagonal (fun i ↦ ((j i)⁻¹ : ℂ)))) {m : ℚ} (hm : 0 ≤ m) (K : ℕ) {Q : ℝ} (hQ : 0 ≤ Q)
    (cf : ℕ → ℤ) (hcf : ∀ k ≤ K, |(cf k : ℝ) - Q * poissonWeight m K k| ≤ 1) {δ γ : ℝ}
    (hδ : Real.sqrt (∑ i, (j i)⁻¹ ^ 2) ≤ δ) (r : ι) (hγ : (j r)⁻¹ ≤ γ) :
    nsq ((fun w ↦ (((j r) * (j w)⁻¹ * hornerFloor D s cf r K w : ℝ) : ℂ)) -
        (Q : ℂ) • (exp ((m : ℂ) • (A - 1)) *ᵥ Pi.single r 1)) ≤
      (j r * ((K + 1) * (δ + γ)) + Q * poissonError m K) ^ 2 := by
  set a : ℕ → ℝ := fun k ↦ Q * poissonWeight m K k
  set P : Matrix ι ι ℂ := ∑ k ∈ Finset.range (K + 1), (poissonWeight m K k : ℂ) • A ^ k
  set F : Matrix ι ι ℂ := exp ((m : ℂ) • (A - 1))
  set Y : ι → ℝ := (∑ k ∈ Finset.range (K + 1), a k • scaledMatrix D s ^ k) *ᵥ Pi.single r 1
  have hω : ∀ i, 0 ≤ (j i)⁻¹ ^ 2 := fun i ↦ sq_nonneg _
  have hW := isWeightedContraction_of_similar hA hcon hj hD
  have hγ' : Real.sqrt ((j r)⁻¹ ^ 2) ≤ γ := by
    rw [Real.sqrt_sq (inv_nonneg.2 (hj r).le)]; exact hγ
  have hH := hornerFloor_error_weighted hω hs hW hδ r hγ' K hcf
  -- The ideal column: `(j_r / j_w) Y(w) = Q (P e_r)(w)`.
  have hY : ∀ w, (((j r) * (j w)⁻¹ * Y w : ℝ) : ℂ) = (Q : ℂ) * (P *ᵥ Pi.single r 1) w := by
    intro w
    have h := ideal_column_eq hj hs hD a K r w
    have hP : (∑ k ∈ Finset.range (K + 1), (a k : ℂ) • A ^ k) = (Q : ℂ) • P := by
      simp only [P, Finset.smul_sum, smul_smul, a]
      push_cast; rfl
    rw [hP, smul_mulVec, Pi.smul_apply, smul_eq_mul] at h
    push_cast
    rw [h]
    have hjr : (j r : ℂ) ≠ 0 := by exact_mod_cast (hj r).ne'
    have hjw : (j w : ℂ) ≠ 0 := by exact_mod_cast (hj w).ne'
    field_simp
  -- Split the error into the rounding part and the Poisson part.
  have hsplit : (fun w ↦ (((j r) * (j w)⁻¹ * hornerFloor D s cf r K w : ℝ) : ℂ)) -
      (Q : ℂ) • (F *ᵥ Pi.single r 1) =
      (fun w ↦ (((j r) * ((j w)⁻¹ * ((hornerFloor D s cf r K w : ℝ) - Y w)) : ℝ) : ℂ)) +
        (Q : ℂ) • ((P - F) *ᵥ Pi.single r 1) := by
    funext w
    have := hY w
    simp only [Pi.sub_apply, Pi.add_apply, Pi.smul_apply, smul_eq_mul, sub_mulVec] at this ⊢
    push_cast at this ⊢
    linear_combination this
  have hδ0 : 0 ≤ δ := (Real.sqrt_nonneg _).trans hδ
  have hγ0 : 0 ≤ γ := (inv_nonneg.2 (hj r).le).trans hγ
  have hR : ∀ v : ι → ℝ, nsq (fun i ↦ (v i : ℂ)) = ∑ i, v i ^ 2 := Spectral.nsq_ofReal
  rw [hsplit]
  refine nsq_add_le_sq (by have := hj r; positivity) (by
    have := poissonError_nonneg (m := (m : ℝ)) (by exact_mod_cast hm) K; positivity) ?_ ?_
  · -- The rounding part, in the weighted norm.
    rw [hR]
    have hsq : ∑ w, ((j r) * ((j w)⁻¹ * ((hornerFloor D s cf r K w : ℝ) - Y w))) ^ 2 =
        (j r) ^ 2 * ∑ w, (j w)⁻¹ ^ 2 * ((hornerFloor D s cf r K w : ℝ) - Y w) ^ 2 := by
      rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun w _ ↦ by ring
    rw [hsq, mul_pow]
    have h0 : 0 ≤ ∑ w, (j w)⁻¹ ^ 2 * ((hornerFloor D s cf r K w : ℝ) - Y w) ^ 2 :=
      Finset.sum_nonneg fun w _ ↦ by positivity
    have h1 := (Real.sqrt_le_left (by positivity)).1 hH
    exact mul_le_mul_of_nonneg_left h1 (sq_nonneg _)
  · -- The Poisson part.
    rw [nsq_smul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hQ, mul_pow]
    have h := nsq_exp_sub_poisson_mulVec_le_of_qform hA hcon hm K (Pi.single r 1)
    have h1 : nsq (Pi.single r (1 : ℂ)) = 1 := by
      rw [nsq_eq_sum, Finset.sum_eq_single r (fun i _ hi ↦ by simp [hi]) (by simp)]
      simp
    rw [h1, mul_one] at h
    have hPF : (P - F) *ᵥ Pi.single r 1 = -((F - P) *ᵥ Pi.single r 1) := by
      rw [← neg_mulVec, neg_sub]
    have hneg : ∀ v : ι → ℂ, nsq (-v) = nsq v := fun v ↦ by simp [nsq_eq_sum]
    rw [hPF, hneg]
    exact mul_le_mul_of_nonneg_left h (sq_nonneg _)

/-- **Gibbs-trace enclosure from integer Horner columns of a similar matrix.** Under the
hypotheses of `nsq_column_sub_le`, if moreover `gibbsHalf H b c = exp (m (A - 1))`,
`j ≤ j_max` and `ε ≥ √|ι| (j_max (K + 1)(δ + γ) + Q η)`, then with
`W = ∑_{r,w} (j_r/j_w)² X_r(w)²`,
`e^{bc} max(0, √W - ε)² / Q² ≤ Tr e^{-bH} ≤ e^{bc} (√W + ε)² / Q²`. -/
theorem gibbsTrace_enclosure_weighted {H A : Matrix ι ι ℂ} (hH : H.IsHermitian) (b c : ℝ)
    (hA : A.IsHermitian) (hcon : ∀ x, |qform A x| ≤ nsq x) {m : ℚ} (hm : 0 ≤ m)
    (hF : gibbsHalf H b c = exp ((m : ℂ) • (A - 1))) {j : ι → ℝ} (hj : ∀ i, 0 < j i)
    {jmax : ℝ} (hjmax0 : 0 ≤ jmax) (hjmax : ∀ i, j i ≤ jmax) {D : Matrix ι ι ℤ} {s : ℕ} (hs : 0 < s)
    (hD : D.map ((↑) : ℤ → ℂ) = (s : ℂ) • (diagonal (fun i ↦ (j i : ℂ)) * A *
      diagonal (fun i ↦ ((j i)⁻¹ : ℂ)))) (K : ℕ) {Q : ℝ} (hQ : 0 < Q) (cf : ℕ → ℤ)
    (hcf : ∀ k ≤ K, |(cf k : ℝ) - Q * poissonWeight m K k| ≤ 1) {δ γ : ℝ}
    (hδ : Real.sqrt (∑ i, (j i)⁻¹ ^ 2) ≤ δ) (hγ0 : 0 ≤ γ) (hγ : ∀ i, (j i)⁻¹ ≤ γ) {ε : ℝ}
    (hε : Real.sqrt (Fintype.card ι) *
      (jmax * ((K + 1) * (δ + γ)) + Q * poissonError m K) ≤ ε) :
    Real.exp (b * c) * max 0 (Real.sqrt (∑ r, ∑ w, (j r / j w) ^ 2 *
        (hornerFloor D s cf r K w : ℝ) ^ 2) - ε) ^ 2 / Q ^ 2 ≤ gibbsTrace H b ∧
      gibbsTrace H b ≤ Real.exp (b * c) * (Real.sqrt (∑ r, ∑ w, (j r / j w) ^ 2 *
        (hornerFloor D s cf r K w : ℝ) ^ 2) + ε) ^ 2 / Q ^ 2 := by
  set W := ∑ r, ∑ w, (j r / j w) ^ 2 * (hornerFloor D s cf r K w : ℝ) ^ 2
  set F := gibbsHalf H b c
  set N2 := ∑ r, nsq (F *ᵥ Pi.single r 1) with hN2
  have hN2nn : 0 ≤ N2 := Finset.sum_nonneg fun _ _ ↦ nsq_nonneg _
  have hZ : gibbsTrace H b = Real.exp (b * c) * N2 := gibbsTrace_eq_exp_mul_sum_nsq hH b c
  have hη : 0 ≤ poissonError (m : ℝ) K := poissonError_nonneg (by exact_mod_cast hm) K
  have hδ0 : 0 ≤ δ := (Real.sqrt_nonneg _).trans hδ
  set E := jmax * ((K + 1) * (δ + γ)) + Q * poissonError m K
  -- Each column error is at most the uniform `E`.
  set X : ι → ι → ℂ := fun r w ↦ (((j r) * (j w)⁻¹ * hornerFloor D s cf r K w : ℝ) : ℂ)
  set Y : ι → ι → ℂ := fun r ↦ (Q : ℂ) • (F *ᵥ Pi.single r 1)
  have hcol : ∀ r ∈ (Finset.univ : Finset ι), nsq (X r - Y r) ≤ E ^ 2 := by
    intro r _
    have hγr := hγ r
    have hγ0 : 0 ≤ γ := (inv_nonneg.2 (hj r).le).trans hγr
    have h := nsq_column_sub_le hA hcon hj hs hD hm K hQ.le cf hcf hδ r hγr
    simp only [X, Y, hF]
    refine h.trans (pow_le_pow_left₀ (by have := hj r; positivity) ?_ 2)
    have : j r * ((K + 1) * (δ + γ)) ≤ jmax * ((K + 1) * (δ + γ)) :=
      mul_le_mul_of_nonneg_right (hjmax r) (by positivity)
    linarith
  have hE : 0 ≤ E := by positivity
  have hmink := abs_sqrt_wsum_sub_le Finset.univ (w := fun _ ↦ 1) (fun _ ↦ zero_le_one) X Y hE
    hcol
  have hR : ∀ v : ι → ℝ, nsq (fun i ↦ (v i : ℂ)) = ∑ i, v i ^ 2 := Spectral.nsq_ofReal
  have hXW : ∑ r, (1 : ℝ) * nsq (X r) = W := by
    refine Finset.sum_congr rfl fun r _ ↦ ?_
    rw [one_mul, hR]
    refine Finset.sum_congr rfl fun w _ ↦ ?_
    rw [div_eq_mul_inv]; ring
  have hYN : ∑ r, (1 : ℝ) * nsq (Y r) = (Q * Real.sqrt N2) ^ 2 := by
    rw [mul_pow, Real.sq_sqrt hN2nn, hN2, Finset.mul_sum]
    refine Finset.sum_congr rfl fun r _ ↦ ?_
    simp only [Y, one_mul, nsq_smul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hQ]
  have hQN : 0 ≤ Q * Real.sqrt N2 := mul_nonneg hQ.le (Real.sqrt_nonneg _)
  rw [hXW, hYN, Real.sqrt_sq hQN, Finset.sum_const, Finset.card_univ, nsmul_eq_mul,
    mul_one] at hmink
  have hdiff := abs_le.1 (hmink.trans hε)
  have hQ2 : 0 < Q ^ 2 := by positivity
  have hsq : (Q * Real.sqrt N2) ^ 2 = Q ^ 2 * N2 := by rw [mul_pow, Real.sq_sqrt hN2nn]
  have hexp := Real.exp_pos (b * c)
  constructor
  · -- Lower bound: `Q √N2 ≥ max(0, √W - ε)`.
    have hm' : max 0 (Real.sqrt W - ε) ≤ Q * Real.sqrt N2 := max_le hQN (by linarith)
    have hm2 : max 0 (Real.sqrt W - ε) ^ 2 / Q ^ 2 ≤ N2 := by
      rw [div_le_iff₀ hQ2, mul_comm, ← hsq]
      exact pow_le_pow_left₀ (le_max_left _ _) hm' 2
    rw [hZ, mul_div_assoc]
    exact mul_le_mul_of_nonneg_left hm2 hexp.le
  · -- Upper bound: `Q √N2 ≤ √W + ε`.
    have hM : Q * Real.sqrt N2 ≤ Real.sqrt W + ε := by linarith
    have hM2 : N2 ≤ (Real.sqrt W + ε) ^ 2 / Q ^ 2 := by
      rw [le_div_iff₀ hQ2, mul_comm, ← hsq]
      exact pow_le_pow_left₀ hQN hM 2
    rw [hZ, mul_div_assoc]
    exact mul_le_mul_of_nonneg_left hM2 hexp.le

end CertifiedGaps.Trace
