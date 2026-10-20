/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.MomentBounds
import CertifiedGaps.Spectral.FunctionalCalculus
import CertifiedGaps.Spectral.Schur
import CertifiedGaps.Transfer.MomentForm

/-!
# Two-boundary pairings with an isolated leading line

The boundary moment inequalities (`OpenMomentWitness`) are proved in `Open/MomentBounds.lean`
from a single pairing `⟨v, Xⁿ v⟩` and the full eigendecomposition of `X`. The thermal appendix
of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026) uses three further
facts, all about pairings `⟨v, Xⁿ W v⟩` whose right boundary vector is moved by an operator `W`
commuting with `X`:

* the same-field chain at odd lengths is `⟨v, Xⁿ U v⟩` for a unitary `U` (the `π` rotation
  about `x`);
* the averaged trace `I(n) = ⟨v, Xⁿ Π v⟩` uses the orthogonal projection `Π` onto the vectors
  fixed by a rotation group;
* both `U` and `Π` fix the leading eigenvector `ψ` of `X`.

This file isolates the linear algebra. Suppose `X` is Hermitian with a unit eigenvector `ψ`,
`X ψ = λ₀ ψ`, and that the rest `R = X - λ₀ |ψ⟩⟨ψ|` has quadratic form at most `r ≤ |λ₀|` in
modulus (`IsLeadingLine X ψ λ₀ r`; the eigenvalue is named `e₀` in the code). Then
`Xⁿ = λ₀ⁿ |ψ⟩⟨ψ| + Rⁿ` for `n ≥ 1`, and in the eigenbasis of `R`:

* the ring traces, the even pairings `⟨v, X²ᴺ v⟩` and the odd pairings `⟨v, X²ᴺ⁺¹ U v⟩` satisfy
  the boundary moment inequalities with `μ = |λ₀|`, leading weight `c = |⟨ψ, v⟩|²` and a ratio
  `ρ` with `ρ |λ₀| ≤ r` (`IsLeadingLine.openMomentWitness`). The odd pairing differs from
  `s c μ²ᴺ⁺¹` by `⟨Rᴷ v, R^{2N+1-K} U v⟩`, which Cauchy–Schwarz bounds by the even remainder;
* for every polynomial `f(x) = ∑ₙ aₙ xⁿ` with `f ≥ 0` and no constant term,
  `c f(λ₀) ≤ ∑ₙ aₙ ⟨v, Xⁿ Π v⟩` (`IsLeadingLine.mul_poly_le`): the rest contributes
  `⟨Π v, f(R) Π v⟩ ≥ 0`, and `Π ψ = ψ` makes the leading weight of `Π v` equal to `c`.

The second item is the paper's eq. (appth:polynomial) with `f = P_*²`; no positivity of odd
powers of `X` is used.

## Main results

* `OpenMomentWitness.of_split`: the inequalities from a moment form with an explicit leading term
  and an explicit bound on the odd remainders.
* `IsLeadingLine.pow_eq`: `Xⁿ = λ₀ⁿ |ψ⟩⟨ψ| + Rⁿ`.
* `IsLeadingLine.eq_smul_of_eigen`: for `r < |λ₀|` the eigenspace of `λ₀` is the line of `ψ`.
* `IsLeadingLine.openMomentWitness`, `IsLeadingLine.mul_poly_le`.
-/

open Finset Matrix

namespace CertifiedGaps

/-! ### A moment form with an explicit leading term -/

/-- For `0 ≤ x ≤ y` and `K ≤ N`, `y^K x^N ≤ y^N x^K`. -/
theorem pow_mul_pow_le_pow_mul_pow {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) {K N : ℕ}
    (hKN : K ≤ N) : y ^ K * x ^ N ≤ y ^ N * x ^ K := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hKN
  have h0 : 0 ≤ y ^ K * x ^ K := mul_nonneg (pow_nonneg (hx.trans hxy) K) (pow_nonneg hx K)
  calc y ^ K * x ^ (K + d) = y ^ K * x ^ K * x ^ d := by ring
    _ ≤ y ^ K * x ^ K * y ^ d := mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hx hxy d) h0
    _ = y ^ (K + d) * x ^ K := by ring

/-- **The boundary moment inequalities from a split moment form.** If the ring traces are
`Z(2N) = μ²ᴺ + ∑ᵢ rᵢ²ᴺ`, the even open-chain traces are `G(2N) = c μ²ᴺ + ∑ᵢ wᵢ rᵢ²ᴺ` with
`wᵢ ≥ 0`, `|rᵢ| ≤ μ`, and every odd trace lies within `μ^{2N+1-2K} ∑ᵢ wᵢ rᵢ²ᴷ` of `s c μ²ᴺ⁺¹`,
then `(μ, c, ρ, s)` is a witness for a ratio `ρ` with `ρ μ ≤ b` for every bound `b ≥ 0` on the
`|rᵢ|`. -/
theorem OpenMomentWitness.of_split {κ : Type*} [Fintype κ] {μ c s : ℝ} (r w : κ → ℝ)
    (hμ : 0 ≤ μ) (hc : 0 ≤ c) (hs : s ^ 2 = 1) (hw : ∀ i, 0 ≤ w i) (hr : ∀ i, |r i| ≤ μ)
    {Z G : ℕ → ℝ} (hZ : ∀ N, 1 ≤ N → Z (2 * N) = μ ^ (2 * N) + ∑ i, r i ^ (2 * N))
    (hGe : ∀ N, 1 ≤ N → G (2 * N) = c * μ ^ (2 * N) + ∑ i, w i * r i ^ (2 * N))
    (hGo : ∀ K N, 1 ≤ K → K ≤ N → |G (2 * N + 1) - s * c * μ ^ (2 * N + 1)| ≤
      μ ^ (2 * N + 1 - 2 * K) * ∑ i, w i * r i ^ (2 * K)) :
    ∃ ρ, OpenMomentWitness Z G μ c ρ s ∧ ∀ b, 0 ≤ b → (∀ i, |r i| ≤ b) → ρ * μ ≤ b := by
  classical
  -- The largest modulus `m` of the rest bounds every `|rᵢ|` and is attained (or is `0`).
  obtain ⟨m, hm0, hmμ, hmT, hm⟩ : ∃ m, 0 ≤ m ∧ m ≤ μ ∧ (∀ i, |r i| ≤ m) ∧
      (m = 0 ∨ ∃ i, |r i| = m) := by
    rcases isEmpty_or_nonempty κ with hκ | hκ
    · exact ⟨0, le_rfl, hμ, fun i ↦ (IsEmpty.false i).elim, Or.inl rfl⟩
    · obtain ⟨i₁, -, hmax⟩ := exists_max_image univ (fun i ↦ |r i|) univ_nonempty
      exact ⟨|r i₁|, abs_nonneg _, hr i₁, fun i ↦ hmax i (mem_univ _), Or.inr ⟨i₁, rfl⟩⟩
  have hρμ : m / μ * μ = m := by
    rcases hμ.eq_or_lt with h | h
    · rw [← h, mul_zero]; linarith [h ▸ hmμ]
    · exact div_mul_cancel₀ m h.ne'
  refine ⟨m / μ, ?_, fun b hb hrb ↦ ?_⟩
  swap
  · rw [hρμ]
    rcases hm with h | ⟨i, hi⟩
    · rw [h]; exact hb
    · rw [← hi]; exact hrb i
  have hp0 : ∀ i, 0 ≤ r i ^ 2 := fun i ↦ sq_nonneg _
  have hpm : ∀ i, r i ^ 2 ≤ μ ^ 2 := fun i ↦ by
    rw [← sq_abs (r i)]; exact pow_le_pow_left₀ (abs_nonneg _) (hr i) 2
  have hpr : ∀ i, r i ^ 2 ≤ (m / μ) ^ 2 * μ ^ 2 := fun i ↦ by
    rw [← mul_pow, hρμ, ← sq_abs (r i)]; exact pow_le_pow_left₀ (abs_nonneg _) (hmT i) 2
  have hZ' : ∀ N, 1 ≤ N → Z (2 * N) = μ ^ (2 * N) + ∑ i, (r i ^ 2) ^ N := fun N hN ↦ by
    rw [hZ N hN]; simp only [pow_mul]
  have hGe' : ∀ N, 1 ≤ N → G (2 * N) = c * μ ^ (2 * N) + ∑ i, w i * (r i ^ 2) ^ N :=
    fun N hN ↦ by rw [hGe N hN]; simp only [pow_mul]
  have hTZ : ∀ N, 0 ≤ ∑ i, (r i ^ 2) ^ N := fun N ↦ sum_nonneg fun i _ ↦ pow_nonneg (hp0 i) N
  have hTG : ∀ N, 0 ≤ ∑ i, w i * (r i ^ 2) ^ N := fun N ↦
    sum_nonneg fun i _ ↦ mul_nonneg (hw i) (pow_nonneg (hp0 i) N)
  exact
    { mu_nonneg := hμ
      c_nonneg := hc
      rho_nonneg := div_nonneg hm0 hμ
      rho_le_one := div_le_one_of_le₀ hmμ hμ
      s_sq := hs
      pow_le_Z := fun N hN ↦ by rw [hZ' N hN]; linarith [hTZ N]
      Z_anti := fun K N hK hKN ↦ by
        rw [hZ' N (hK.trans hKN), hZ' K hK]
        have h : μ ^ (2 * K) * ∑ i, (r i ^ 2) ^ N ≤ μ ^ (2 * N) * ∑ i, (r i ^ 2) ^ K := by
          rw [mul_sum, mul_sum, pow_mul, pow_mul]
          exact sum_le_sum fun i _ ↦ pow_mul_pow_le_pow_mul_pow (hp0 i) (hpm i) hKN
        linear_combination h
      Z_two_mul_le := fun N hN ↦ by
        rw [hZ' _ (by omega), hZ' N hN]
        have h : ∑ i, (r i ^ 2) ^ (2 * N) ≤ μ ^ (2 * N) * ∑ i, (r i ^ 2) ^ N := by
          rw [mul_sum, pow_mul]
          refine sum_le_sum fun i _ ↦ ?_
          rw [two_mul, pow_add]
          exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (hp0 i) (hpm i) N)
            (pow_nonneg (hp0 i) N)
        linear_combination h
      Z_two_mul_le_sq := fun N hN ↦ by
        rw [hZ' _ (by omega), hZ' N hN, add_sub_cancel_left]
        have h : ∑ i, (r i ^ 2) ^ (2 * N) ≤ (∑ i, (r i ^ 2) ^ N) ^ 2 :=
          calc ∑ i, (r i ^ 2) ^ (2 * N) = ∑ i, ((r i ^ 2) ^ N) ^ 2 :=
                sum_congr rfl fun i _ ↦ by ring
            _ ≤ _ := (sq_sub_sum_sq_finset univ _ fun i _ ↦ pow_nonneg (hp0 i) N).1
        linarith
      tail := fun N hN ↦ by
        rw [hZ' N hN, add_sub_cancel_left, hρμ]
        rcases hm with rfl | ⟨i₁, rfl⟩
        · rw [zero_pow (by omega)]; exact hTZ N
        · rw [pow_mul, sq_abs]
          exact single_le_sum (f := fun i ↦ (r i ^ 2) ^ N)
            (fun i _ ↦ pow_nonneg (hp0 i) N) (mem_univ i₁)
      le_G := fun N hN ↦ by rw [hGe' N hN]; linarith [hTG N]
      G_decay := fun K N hK hKN ↦ by
        obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hKN
        rw [hGe' _ (hK.trans hKN), hGe' K hK, show 2 * (K + d) - 2 * K = 2 * d by omega,
          add_sub_cancel_left, add_sub_cancel_left, mul_sum, mul_sum]
        refine sum_le_sum fun i _ ↦ ?_
        have h0 : 0 ≤ (μ ^ 2) ^ K * w i * (r i ^ 2) ^ K :=
          mul_nonneg (mul_nonneg (pow_nonneg (sq_nonneg μ) K) (hw i)) (pow_nonneg (hp0 i) K)
        calc μ ^ (2 * K) * (w i * (r i ^ 2) ^ (K + d))
            = (μ ^ 2) ^ K * w i * (r i ^ 2) ^ K * (r i ^ 2) ^ d := by ring
          _ ≤ (μ ^ 2) ^ K * w i * (r i ^ 2) ^ K * ((m / μ) ^ 2 * μ ^ 2) ^ d :=
            mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (hp0 i) (hpr i) d) h0
          _ = (m / μ) ^ (2 * d) * μ ^ (2 * (K + d)) * (w i * (r i ^ 2) ^ K) := by ring
      G_odd := fun K N hK hKN ↦ by
        rw [hGe K hK, add_sub_cancel_left]
        exact hGo K N hK hKN }

/-! ### An isolated leading line -/

section LeadingLine

open Spectral

variable {D : Type*} [Fintype D]

/-- The rank-one projector `|ψ⟩⟨ψ|`. -/
noncomputable def lineProj (ψ : D → ℂ) : Matrix D D ℂ := vecMulVec ψ (star ψ)

/-- **An isolated leading line.** `X` is Hermitian, `ψ` is a unit eigenvector with eigenvalue
`λ₀`, and the rest `R = X - λ₀ |ψ⟩⟨ψ|` has quadratic form at most `r` in modulus. -/
structure IsLeadingLine (X : Matrix D D ℂ) (ψ : D → ℂ) (e₀ r : ℝ) : Prop where
  hermitian : X.IsHermitian
  unit : star ψ ⬝ᵥ ψ = 1
  eigen : X *ᵥ ψ = (e₀ : ℂ) • ψ
  rest : ∀ y, |qform (X - (e₀ : ℂ) • lineProj ψ) y| ≤ r * nsq y

namespace IsLeadingLine

variable {X : Matrix D D ℂ} {ψ : D → ℂ} {e₀ r : ℝ} (h : IsLeadingLine X ψ e₀ r)
include h

/-- The projector is idempotent. -/
theorem proj_mul_proj : lineProj ψ * lineProj ψ = lineProj ψ := by
  rw [lineProj, vecMulVec_mul_vecMulVec, h.unit, one_smul]

/-- `X |ψ⟩⟨ψ| = λ₀ |ψ⟩⟨ψ|`. -/
theorem mul_proj : X * lineProj ψ = (e₀ : ℂ) • lineProj ψ := by
  rw [lineProj, mul_vecMulVec, h.eigen]
  ext i j; simp [vecMulVec_apply, mul_assoc]

/-- `|ψ⟩⟨ψ| X = λ₀ |ψ⟩⟨ψ|`, since `X` is Hermitian and `λ₀` real. -/
theorem proj_mul : lineProj ψ * X = (e₀ : ℂ) • lineProj ψ := by
  have hv : star ψ ᵥ* X = (e₀ : ℂ) • star ψ := by
    have := congrArg star h.eigen
    rw [star_mulVec, h.hermitian.eq, star_smul, Complex.star_def, Complex.conj_ofReal] at this
    exact this
  rw [lineProj, vecMulVec_mul, hv]
  ext i j; simp [vecMulVec_apply, mul_left_comm]

omit h [Fintype D] in
/-- The projector is Hermitian. -/
theorem proj_isHermitian (ψ : D → ℂ) : (lineProj ψ).IsHermitian := by
  rw [IsHermitian, lineProj, conjTranspose_vecMulVec, star_star]

/-- The rest is Hermitian. -/
theorem rest_isHermitian : (X - (e₀ : ℂ) • lineProj ψ).IsHermitian := by
  refine h.hermitian.sub ?_
  rw [IsHermitian, conjTranspose_smul, (proj_isHermitian ψ).eq, Complex.star_def,
    Complex.conj_ofReal]

/-- `R |ψ⟩⟨ψ| = 0`. -/
theorem rest_mul_proj : (X - (e₀ : ℂ) • lineProj ψ) * lineProj ψ = 0 := by
  rw [Matrix.sub_mul, h.mul_proj, smul_mul_assoc, h.proj_mul_proj, sub_self]

/-- `|ψ⟩⟨ψ| R = 0`. -/
theorem proj_mul_rest : lineProj ψ * (X - (e₀ : ℂ) • lineProj ψ) = 0 := by
  rw [Matrix.mul_sub, h.proj_mul, mul_smul_comm, h.proj_mul_proj, sub_self]

/-- **The leading line splits off every power**: `Xⁿ = λ₀ⁿ |ψ⟩⟨ψ| + Rⁿ` for `n ≥ 1`. -/
theorem pow_eq [DecidableEq D] {n : ℕ} (hn : 1 ≤ n) :
    X ^ n = ((e₀ ^ n : ℝ) : ℂ) • lineProj ψ + (X - (e₀ : ℂ) • lineProj ψ) ^ n := by
  induction n, hn using Nat.le_induction with
  | base => simp
  | succ n hn ih =>
    set R := X - (e₀ : ℂ) • lineProj ψ
    have hRn : R ^ n * lineProj ψ = 0 := by
      obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
      rw [pow_succ, Matrix.mul_assoc, h.rest_mul_proj, Matrix.mul_zero]
    have hX : X = (e₀ : ℂ) • lineProj ψ + R := by simp [R]
    rw [pow_succ, ih]
    conv_lhs => rw [hX]
    rw [Matrix.add_mul, Matrix.mul_add, Matrix.mul_add, smul_mul_smul_comm, h.proj_mul_proj,
      smul_mul_assoc, h.proj_mul_rest, smul_zero, mul_smul_comm, hRn, smul_zero, ← pow_succ]
    push_cast
    simp only [zero_add, add_zero, ← pow_succ]

/-- The bound `r` is nonnegative: test the form of the rest on `ψ`, where it vanishes. -/
theorem rest_nonneg : 0 ≤ r := by
  have hx := h.rest ψ
  have h1 : nsq ψ = 1 := by rw [nsq, h.unit]; simp
  rw [h1, mul_one] at hx
  exact (abs_nonneg _).trans hx

/-- Every eigenvalue of the rest is at most `r` in modulus. -/
theorem abs_eigenvalues_rest_le [DecidableEq D] (i : D) :
    |h.rest_isHermitian.eigenvalues i| ≤ r := by
  have hx := h.rest (eigU h.rest_isHermitian *ᵥ Pi.single i 1)
  rw [qform_eigU, nsq_eigU, nsq_eq_sum] at hx
  simpa [Pi.single_apply, apply_ite] using hx

omit h [Fintype D] in
/-- `|ψ⟩⟨ψ| w = ⟨ψ, w⟩ ψ`. -/
theorem lineProj_mulVec (ψ w : D → ℂ) [Fintype D] :
    lineProj ψ *ᵥ w = (star ψ ⬝ᵥ w) • ψ := by
  ext i
  simp only [lineProj, mulVec, dotProduct, vecMulVec_apply, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul]
  exact Finset.sum_congr rfl fun j _ ↦ by ring

omit h in
/-- `⟨u, |ψ⟩⟨ψ| w⟩ = ⟨u, ψ⟩ ⟨ψ, w⟩`. -/
theorem dotProduct_lineProj_mulVec (ψ u w : D → ℂ) :
    star u ⬝ᵥ (lineProj ψ *ᵥ w) = (star u ⬝ᵥ ψ) * (star ψ ⬝ᵥ w) := by
  rw [lineProj_mulVec, dotProduct_smul, smul_eq_mul, mul_comm]

omit h in
/-- `⟨v, ψ⟩ ⟨ψ, v⟩ = |⟨ψ, v⟩|²`. -/
theorem dotProduct_mul_dotProduct_self (ψ v : D → ℂ) :
    (star v ⬝ᵥ ψ) * (star ψ ⬝ᵥ v) = (Complex.normSq (star ψ ⬝ᵥ v) : ℂ) := by
  have h1 : star v ⬝ᵥ ψ = star (star ψ ⬝ᵥ v) := by
    rw [star_dotProduct]
  rw [h1, Complex.star_def, Complex.normSq_eq_conj_mul_self]

omit h in
/-- An operator `W` that fixes `ψ` and has `Wᴴ ψ = ψ` commutes with `|ψ⟩⟨ψ|` and absorbs it. -/
theorem lineProj_mul_of_fix {W : Matrix D D ℂ} (hW : W *ᵥ ψ = ψ) (hW' : Wᴴ *ᵥ ψ = ψ) :
    W * lineProj ψ = lineProj ψ ∧ lineProj ψ * W = lineProj ψ := by
  refine ⟨by rw [lineProj, mul_vecMulVec, hW], ?_⟩
  have : star ψ ᵥ* W = star ψ := by
    have := congrArg star hW'
    rwa [star_mulVec, conjTranspose_conjTranspose] at this
  rw [lineProj, vecMulVec_mul, this]

omit h in
/-- An operator commuting with `X` and fixing `ψ` (also under `Wᴴ`) commutes with the rest. -/
theorem commute_rest {W : Matrix D D ℂ} (hWX : W * X = X * W) (hW : W *ᵥ ψ = ψ)
    (hW' : Wᴴ *ᵥ ψ = ψ) :
    W * (X - (e₀ : ℂ) • lineProj ψ) = (X - (e₀ : ℂ) • lineProj ψ) * W := by
  obtain ⟨h1, h2⟩ := lineProj_mul_of_fix hW hW'
  rw [Matrix.mul_sub, Matrix.sub_mul, hWX, mul_smul_comm, smul_mul_assoc, h1, h2]

/-- **The leading line is the whole eigenspace of `λ₀`.** If `r < |λ₀|`, every eigenvector of `X`
with eigenvalue `λ₀` is a multiple of `ψ`: its part orthogonal to `ψ` is an eigenvector of the
rest with eigenvalue `λ₀`, which the bound `r` excludes. -/
theorem eq_smul_of_eigen (hr : r < |e₀|) {φ : D → ℂ} (hφ : X *ᵥ φ = (e₀ : ℂ) • φ) :
    φ = (star ψ ⬝ᵥ φ) • ψ := by
  set a := star ψ ⬝ᵥ φ
  set χ := φ - a • ψ
  have hψχ : star ψ ⬝ᵥ χ = 0 := by
    simp only [χ, dotProduct_sub, dotProduct_smul, h.unit, smul_eq_mul, mul_one, a, sub_self]
  have hRχ : (X - (e₀ : ℂ) • lineProj ψ) *ᵥ χ = (e₀ : ℂ) • χ := by
    rw [sub_mulVec, smul_mulVec, lineProj_mulVec, hψχ, zero_smul, smul_zero, sub_zero]
    simp only [χ, mulVec_sub, mulVec_smul, hφ, h.eigen, smul_sub, smul_comm a]
  have hq : qform (X - (e₀ : ℂ) • lineProj ψ) χ = e₀ * nsq χ := by
    rw [qform, hRχ, dotProduct_smul, smul_eq_mul, nsq, RCLike.re_eq_complex_re,
      Complex.re_ofReal_mul]
  have hb := h.rest χ
  rw [hq, abs_mul] at hb
  have h0 : nsq χ = 0 := by
    by_contra hne
    have hpos : 0 < nsq χ := lt_of_le_of_ne (nsq_nonneg χ) (Ne.symm hne)
    rw [abs_of_pos hpos] at hb
    exact absurd (le_of_mul_le_mul_right hb hpos) (not_le.2 hr)
  have := nsq_eq_zero_iff.1 h0
  exact (sub_eq_zero.1 this)


variable [DecidableEq D]

omit h in
/-- In the eigenbasis of a Hermitian `A`, `Re ⟨x, Aⁿ x⟩ = ∑ᵢ λᵢⁿ |(Uᴴ x)ᵢ|²`. -/
theorem qform_pow_eq_sum {A : Matrix D D ℂ} (hA : A.IsHermitian) (n : ℕ) (x : D → ℂ) :
    qform (A ^ n) x = ∑ i, hA.eigenvalues i ^ n * ‖((eigU hA)ᴴ *ᵥ x) i‖ ^ 2 := by
  rw [hermFun_pow hA, hermFun_eq, qform_LDL]

omit h in
/-- `Re ⟨x, A y⟩` for Hermitian `A` is the real dot product of `A x` and `y`. -/
theorem re_dotProduct_pow_mulVec {A : Matrix D D ℂ} (hA : A.IsHermitian) (k : ℕ)
    (x y : D → ℂ) : star x ⬝ᵥ (A ^ k *ᵥ y) = star (A ^ k *ᵥ x) ⬝ᵥ y := by
  rw [star_mulVec, ← dotProduct_mulVec, (hA.pow k).eq]

/-- The traces: `Re Tr Xⁿ = λ₀ⁿ + ∑ᵢ ρᵢⁿ` for `n ≥ 1`, with `ρᵢ` the eigenvalues of the rest. -/
theorem re_trace_pow {n : ℕ} (hn : 1 ≤ n) :
    (trace (X ^ n)).re = e₀ ^ n + ∑ i, h.rest_isHermitian.eigenvalues i ^ n := by
  have htr : trace (lineProj ψ) = 1 := by
    rw [lineProj, trace_vecMulVec, dotProduct_comm, h.unit]
  rw [h.pow_eq hn, trace_add, trace_smul, htr, trace_pow_eq_sum_eigenvalues h.rest_isHermitian]
  simp only [smul_eq_mul, mul_one, Complex.add_re, Complex.ofReal_re]

/-- The pairings: `Re ⟨u, Xⁿ w⟩ = λ₀ⁿ Re (⟨u, ψ⟩⟨ψ, w⟩) + Re ⟨u, Rⁿ w⟩` for `n ≥ 1`. -/
theorem re_dotProduct_pow {n : ℕ} (hn : 1 ≤ n) (u w : D → ℂ) :
    (star u ⬝ᵥ (X ^ n *ᵥ w)).re = e₀ ^ n * ((star u ⬝ᵥ ψ) * (star ψ ⬝ᵥ w)).re +
      (star u ⬝ᵥ ((X - (e₀ : ℂ) • lineProj ψ) ^ n *ᵥ w)).re := by
  rw [h.pow_eq hn, add_mulVec, dotProduct_add, smul_mulVec, dotProduct_smul,
    dotProduct_lineProj_mulVec, Complex.add_re, smul_eq_mul, Complex.re_ofReal_mul]

/-- **The boundary moment inequalities of two-boundary pairings with a leading line.** Let
`U` be unitary, commute with `X` and fix `ψ`. The ring traces `Z(2N) = Re Tr X²ᴺ`, the even
pairings `G(2N) = Re ⟨v, X²ᴺ v⟩` and the odd pairings `G(2N+1) = Re ⟨v, X²ᴺ⁺¹ U v⟩` admit the
witness `(|λ₀|, |⟨ψ, v⟩|², ρ, sign λ₀)` with `ρ |λ₀| ≤ r`. -/
theorem openMomentWitness (hr : r ≤ |e₀|) (v : D → ℂ) {U : Matrix D D ℂ} (hU : Uᴴ * U = 1)
    (hUX : U * X = X * U) (hUψ : U *ᵥ ψ = ψ) {Z G : ℕ → ℝ}
    (hZ : ∀ N, 1 ≤ N → Z (2 * N) = (trace (X ^ (2 * N))).re)
    (hGe : ∀ N, 1 ≤ N → G (2 * N) = (star v ⬝ᵥ (X ^ (2 * N) *ᵥ v)).re)
    (hGo : ∀ N, 1 ≤ N → G (2 * N + 1) = (star v ⬝ᵥ (X ^ (2 * N + 1) *ᵥ (U *ᵥ v))).re) :
    ∃ ρ, OpenMomentWitness Z G |e₀| (Complex.normSq (star ψ ⬝ᵥ v)) ρ
      (if 0 ≤ e₀ then 1 else -1) ∧ ρ * |e₀| ≤ r := by
  set hR := h.rest_isHermitian
  set R := X - (e₀ : ℂ) • lineProj ψ
  set ρs := hR.eigenvalues
  set w : D → ℝ := fun i ↦ ‖((eigU hR)ᴴ *ᵥ v) i‖ ^ 2
  set c := Complex.normSq (star ψ ⬝ᵥ v)
  set s : ℝ := if 0 ≤ e₀ then 1 else -1
  have hs : s ^ 2 = 1 := by simp only [s]; split_ifs <;> norm_num
  have he₀ : e₀ = s * |e₀| := by
    simp only [s]
    split_ifs with he
    · rw [abs_of_nonneg he, one_mul]
    · rw [abs_of_neg (not_le.1 he)]; ring
  have hρ : ∀ i, |ρs i| ≤ |e₀| := fun i ↦ (h.abs_eigenvalues_rest_le i).trans hr
  have hr0 : 0 ≤ r := h.rest_nonneg
  -- The even powers do not see the sign.
  have heven : ∀ N, e₀ ^ (2 * N) = |e₀| ^ (2 * N) := fun N ↦ by
    rw [pow_mul, pow_mul, sq_abs]
  have hpair : ∀ n, qform (R ^ n) v = ∑ i, ρs i ^ n * w i := fun n ↦
    qform_pow_eq_sum hR n v
  have hUψ' : Uᴴ *ᵥ ψ = ψ := by
    conv_lhs => rw [← hUψ]
    rw [mulVec_mulVec, hU, one_mulVec]
  have hUR : U * R = R * U := commute_rest hUX hUψ hUψ'
  -- The squared norms `‖Rᵏ v‖² = ∑ᵢ ρᵢ²ᵏ wᵢ`.
  have hnsq : ∀ k, nsq (R ^ k *ᵥ v) = ∑ i, w i * ρs i ^ (2 * k) := fun k ↦ by
    rw [nsq_mulVec, (hR.pow k).eq, ← pow_add, ← two_mul, hpair]
    exact sum_congr rfl fun i _ ↦ mul_comm _ _
  refine OpenMomentWitness.of_split ρs w (abs_nonneg _) (Complex.normSq_nonneg _) hs
    (fun i ↦ by positivity) hρ (fun N hN ↦ ?_) (fun N hN ↦ ?_) (fun K N hK hKN ↦ ?_)
    |>.imp fun ρ hρ' ↦ ⟨hρ'.1, ?_⟩
  · rw [hZ N hN, h.re_trace_pow (by omega), heven]
  · rw [hGe N hN, h.re_dotProduct_pow (by omega), dotProduct_mul_dotProduct_self,
      Complex.ofReal_re, heven, mul_comm (|e₀| ^ _)]
    congr 1
    rw [← RCLike.re_eq_complex_re]
    exact (hpair _).trans (sum_congr rfl fun i _ ↦ mul_comm _ _)
  · -- The odd pairing minus its leading term is `Re ⟨Rᴷ v, R^{2N+1-K} U v⟩`.
    have hlead : star ψ ⬝ᵥ (U *ᵥ v) = star ψ ⬝ᵥ v := by
      rw [dotProduct_mulVec, ← star_star (star ψ ᵥ* U), star_vecMul, star_star, hUψ']
    have hodd : e₀ ^ (2 * N + 1) = s * |e₀| ^ (2 * N + 1) := by
      conv_lhs => rw [he₀]
      rw [mul_pow, pow_succ s, pow_mul, hs, one_pow, one_mul]
    have hpow : R ^ K *ᵥ (R ^ (2 * N + 1 - K) *ᵥ (U *ᵥ v)) = R ^ (2 * N + 1) *ᵥ (U *ᵥ v) := by
      rw [mulVec_mulVec, ← pow_add, show K + (2 * N + 1 - K) = 2 * N + 1 by omega]
    have hsplit : G (2 * N + 1) - s * c * |e₀| ^ (2 * N + 1) =
        RCLike.re (star (R ^ K *ᵥ v) ⬝ᵥ (R ^ (2 * N + 1 - K) *ᵥ (U *ᵥ v))) := by
      rw [hGo N (hK.trans hKN), h.re_dotProduct_pow (by omega), hlead,
        dotProduct_mul_dotProduct_self, Complex.ofReal_re, ← re_dotProduct_pow_mulVec hR, hpow,
        RCLike.re_eq_complex_re, hodd]
      ring
    set S := ∑ i, w i * ρs i ^ (2 * K)
    have hS : 0 ≤ S := sum_nonneg fun i _ ↦ mul_nonneg (by positivity) (by
      rw [pow_mul]; positivity)
    -- `‖R^b U v‖² = ‖R^b v‖² ≤ |λ₀|^{2(b-K)} ‖Rᴷ v‖²`.
    have hb : nsq (R ^ (2 * N + 1 - K) *ᵥ (U *ᵥ v)) ≤ (|e₀| ^ (2 * N + 1 - 2 * K)) ^ 2 * S := by
      have hcomm : R ^ (2 * N + 1 - K) * U = U * R ^ (2 * N + 1 - K) :=
        ((Commute.pow_left hUR.symm _).eq)
      rw [mulVec_mulVec, hcomm, ← mulVec_mulVec, nsq_mulVec_of_unitary hU, hnsq, mul_sum]
      refine sum_le_sum fun i _ ↦ ?_
      have hpi : ρs i ^ 2 ≤ |e₀| ^ 2 := by
        rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) (hρ i) 2
      have e1 : ρs i ^ (2 * (2 * N + 1 - K)) =
          (ρs i ^ 2) ^ (2 * N + 1 - 2 * K) * ρs i ^ (2 * K) := by
        rw [← pow_mul, ← pow_add]; congr 1; omega
      have e2 : (|e₀| ^ (2 * N + 1 - 2 * K)) ^ 2 = (|e₀| ^ 2) ^ (2 * N + 1 - 2 * K) := by
        rw [← pow_mul, ← pow_mul, mul_comm]
      rw [e1, e2]
      have : 0 ≤ ρs i ^ (2 * K) := by rw [pow_mul]; positivity
      have hw0 : 0 ≤ w i := by positivity
      calc w i * ((ρs i ^ 2) ^ (2 * N + 1 - 2 * K) * ρs i ^ (2 * K))
          = (ρs i ^ 2) ^ (2 * N + 1 - 2 * K) * (w i * ρs i ^ (2 * K)) := by ring
        _ ≤ (|e₀| ^ 2) ^ (2 * N + 1 - 2 * K) * (w i * ρs i ^ (2 * K)) :=
          mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (sq_nonneg _) hpi _)
            (mul_nonneg hw0 this)
    have hcs := re_dotProduct_sq_le (R ^ K *ᵥ v) (R ^ (2 * N + 1 - K) *ᵥ (U *ᵥ v))
    rw [hnsq] at hcs
    rw [hsplit]
    refine abs_le_of_sq_le_sq ?_ (by positivity)
    calc _ ≤ S * nsq (R ^ (2 * N + 1 - K) *ᵥ (U *ᵥ v)) := hcs
      _ ≤ S * ((|e₀| ^ (2 * N + 1 - 2 * K)) ^ 2 * S) := mul_le_mul_of_nonneg_left hb hS
      _ = (|e₀| ^ (2 * N + 1 - 2 * K) * S) ^ 2 := by ring
  · exact hρ'.2 r hr0 fun i ↦ h.abs_eigenvalues_rest_le i

/-- **The leading weight under an averaging projection** (eq. (appth:polynomial) of the paper).
Let `Π` be an orthogonal projection commuting with `X` and fixing `ψ`, and let
`f(x) = ∑_{n ∈ S} aₙ xⁿ` be nonnegative with `n ≥ 1` on `S`. Then
`|⟨ψ, v⟩|² f(λ₀) ≤ ∑_{n ∈ S} aₙ Re ⟨v, Xⁿ Π v⟩`. -/
theorem mul_poly_le (v : D → ℂ) {P : Matrix D D ℂ} (hP : Pᴴ = P) (hPP : P * P = P)
    (hPX : P * X = X * P) (hPψ : P *ᵥ ψ = ψ) (S : Finset ℕ) (a : ℕ → ℝ)
    (hS : ∀ n ∈ S, 1 ≤ n) (hf : ∀ x : ℝ, 0 ≤ ∑ n ∈ S, a n * x ^ n) :
    Complex.normSq (star ψ ⬝ᵥ v) * ∑ n ∈ S, a n * e₀ ^ n ≤
      ∑ n ∈ S, a n * (star v ⬝ᵥ (X ^ n *ᵥ (P *ᵥ v))).re := by
  set hR := h.rest_isHermitian
  set R := X - (e₀ : ℂ) • lineProj ψ
  have hPψ' : Pᴴ *ᵥ ψ = ψ := by rw [hP, hPψ]
  have hPR : P * R = R * P := commute_rest hPX hPψ hPψ'
  have hlead : star ψ ⬝ᵥ (P *ᵥ v) = star ψ ⬝ᵥ v := by
    rw [dotProduct_mulVec, ← star_star (star ψ ᵥ* P), star_vecMul, star_star, hPψ']
  -- The rest contributes `Re ⟨P v, Rⁿ P v⟩`.
  have hrest : ∀ n, (star v ⬝ᵥ (R ^ n *ᵥ (P *ᵥ v))).re = qform (R ^ n) (P *ᵥ v) := fun n ↦ by
    have hc : R ^ n * P = P * R ^ n := (Commute.pow_left hPR.symm n).eq
    have hPRP : P * (R ^ n * P) = R ^ n * P := by rw [hc, ← Matrix.mul_assoc, hPP]
    rw [qform, star_mulVec, ← dotProduct_mulVec, hP, mulVec_mulVec, mulVec_mulVec, hPRP,
      RCLike.re_eq_complex_re]
  have hterm : ∀ n ∈ S, (star v ⬝ᵥ (X ^ n *ᵥ (P *ᵥ v))).re =
      e₀ ^ n * Complex.normSq (star ψ ⬝ᵥ v) +
        ∑ i, hR.eigenvalues i ^ n * ‖((eigU hR)ᴴ *ᵥ (P *ᵥ v)) i‖ ^ 2 := fun n hn ↦ by
    rw [h.re_dotProduct_pow (hS n hn), hlead, dotProduct_mul_dotProduct_self, Complex.ofReal_re,
      hrest, qform_pow_eq_sum]
  have hsum : ∑ n ∈ S, a n * (star v ⬝ᵥ (X ^ n *ᵥ (P *ᵥ v))).re =
      ∑ n ∈ S, a n * (e₀ ^ n * Complex.normSq (star ψ ⬝ᵥ v) +
        ∑ i, hR.eigenvalues i ^ n * ‖((eigU hR)ᴴ *ᵥ (P *ᵥ v)) i‖ ^ 2) :=
    sum_congr rfl fun n hn ↦ by rw [hterm n hn]
  have hlead' : Complex.normSq (star ψ ⬝ᵥ v) * ∑ n ∈ S, a n * e₀ ^ n =
      ∑ n ∈ S, a n * (e₀ ^ n * Complex.normSq (star ψ ⬝ᵥ v)) := by
    rw [mul_sum]; exact sum_congr rfl fun n _ ↦ by ring
  rw [hsum, hlead']
  simp only [mul_add, sum_add_distrib]
  refine le_add_of_nonneg_right ?_
  simp only [mul_sum]
  rw [sum_comm]
  refine sum_nonneg fun i _ ↦ ?_
  simp only [← mul_assoc, ← sum_mul]
  exact mul_nonneg (hf _) (by positivity)

end IsLeadingLine

end LeadingLine

end CertifiedGaps
