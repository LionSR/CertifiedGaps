/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity
import Mathlib.Data.Nat.Log

/-!
# A two-purity doubling criterion for a uniform spectral gap

This file formalizes the abstract part of the note *A thermal-MPO derivation of the spin-1 gap
bound* (dated October 7, 2026; cited below as *the thermal-MPO note*): a uniform lower bound on
the spectral gap of a family of finite-dimensional Hamiltonians, derived from two finite-size
purities, for any family whose partition functions are moment sums of a spatial transfer
operator.

## Setting

For each length `L`, the Hamiltonian is represented by its spectrum with multiplicity
`E L : ι L → ℝ`, over a finite nonempty index type. Then
* `partitionFunction E L β = ∑ₘ exp (-β Eₘ(L))` is `Z_L(β) = Tr e^{-β H_L}` (note, eq. (2));
* `thermalPurity E L β = Z_L(2β) / Z_L(β)²` is `P_τ(L, β)` (note, eq. (3));
* `spatialPurity E n β = Z_{2n}(β) / Z_n(β)²` is `P_s(n, β)` (note, eq. (6));
* `HasSpatialMomentForm E n₀` is the moment form `Z_L(β) = ∑ᵢ |λᵢ(β)|^L` (note, eq. (5)) for
  every even `L ≥ n₀`;
* `HasSpatialMomentBounds E n₀` is the weaker hypothesis actually used: the three moment
  inequalities `SpatialMomentBounds` (`P_s ≤ 1`, squaring, `ℓ^p` monotonicity) at every `β > 0`.
  They follow from the moment form and pass to limits, so approximants with the moment form,
  such as finite Trotter transfer matrices, also give them;
* `GapGT e c` says that the spectrum `e` has a nondegenerate lowest level with every other level
  more than `c` above it, which for ordered eigenvalues `E₀ ≤ E₁ ≤ ⋯` is `E₁ - E₀ > c`
  (`GapGT.lt_sub_of_monotone`).

## Main results

* `one_sub_squaredPurity_le_squaringDefect`: Lemma 1 (Squaring) of the note.
* `spatialPurity_two_mul`, `thermalPurity_four_mul`: Lemma 2 (Coupling), eqs. (9) and (10).
* `exp_neg_mul_sub_le_of_isMin`: the purity-to-gap bound, eq. (4).
* `doubling`: Lemma 3 (Doubling).
* `gapGT_of_seeds`: Proposition 1 with general seed defect `ε₀ ≤ 1/100`, giving
  `Δ_L > log (1/(5ε₀)) / β₀`, under the moment inequalities.
* `gapGT_of_purities`: Proposition 1 (Two-purity criterion) as stated in the note.

## Scope

The criterion is abstract: it needs only the spectra `E L` and the moment inequalities. They are
proved for every reflection-positive ring in `CertifiedGaps.Ring.Gap`, and Theorem 1 of the note
(the spin-one Heisenberg chain) is stated for the Lean Hamiltonian in `CertifiedGaps.Heisenberg`.
The numerical seeds `P_s(n₀, β₀)`, `P_τ(2n₀, β₀)` are the only inputs left as hypotheses.

The moment form is stated with the eigenvalue moduli `|λᵢ(β)|` indexed by `ℕ` (a finite spectrum
is padded by zeros) and only for even lengths `L ≥ n₀`, which is all that the proofs use; this is
implied by the note's eq. (5).
-/

open Filter Finset Real

namespace CertifiedGaps

/-! ### Lemma 1: squaring -/

/-- Finite form of the squaring inequality, in homogeneous power sums: for nonnegative `pᵢ`, with
`S_k = ∑ pᵢ^k`, one has `S₂ ≤ S₁²` and `2 (S₂² - S₄) ≤ (S₁² - S₂)²`. -/
theorem sq_sub_sum_sq_finset {κ : Type*} (s : Finset κ) (p : κ → ℝ) (hp : ∀ i ∈ s, 0 ≤ p i) :
    ∑ i ∈ s, p i ^ 2 ≤ (∑ i ∈ s, p i) ^ 2 ∧
      2 * ((∑ i ∈ s, p i ^ 2) ^ 2 - ∑ i ∈ s, p i ^ 4) ≤
        ((∑ i ∈ s, p i) ^ 2 - ∑ i ∈ s, p i ^ 2) ^ 2 := by
  classical
  suffices h : 0 ≤ ∑ i ∈ s, p i ∧ ∑ i ∈ s, p i ^ 2 ≤ (∑ i ∈ s, p i) ^ 2 ∧
      2 * ((∑ i ∈ s, p i ^ 2) ^ 2 - ∑ i ∈ s, p i ^ 4) ≤
        ((∑ i ∈ s, p i) ^ 2 - ∑ i ∈ s, p i ^ 2) ^ 2 from h.2
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    obtain ⟨h1, h2, h3⟩ := ih fun i hi ↦ hp i (mem_insert_of_mem hi)
    have hq := hp a (mem_insert_self a s)
    rw [sum_insert ha, sum_insert ha, sum_insert ha]
    have hT : 0 ≤ (∑ i ∈ s, p i) ^ 2 - ∑ i ∈ s, p i ^ 2 := by linarith
    refine ⟨by positivity, by nlinarith [mul_nonneg hq h1], ?_⟩
    nlinarith [mul_nonneg (mul_nonneg hq h1) hT, mul_nonneg (sq_nonneg (p a)) hT]

/-- The squaring inequality for summable nonnegative families, in homogeneous power sums. -/
theorem sq_sub_sum_sq_of_hasSum {κ : Type*} {p : κ → ℝ} (hp : ∀ i, 0 ≤ p i) {s₁ s₂ s₄ : ℝ}
    (h₁ : HasSum p s₁) (h₂ : HasSum (fun i ↦ p i ^ 2) s₂) (h₄ : HasSum (fun i ↦ p i ^ 4) s₄) :
    s₂ ≤ s₁ ^ 2 ∧ 2 * (s₂ ^ 2 - s₄) ≤ (s₁ ^ 2 - s₂) ^ 2 := by
  have hfin := fun s : Finset κ ↦ sq_sub_sum_sq_finset s p fun i _ ↦ hp i
  constructor
  · refine le_of_tendsto_of_tendsto' h₂ (h₁.pow 2) fun s ↦ (hfin s).1
  · refine le_of_tendsto_of_tendsto' ((((h₂.pow 2).sub h₄).const_mul 2))
      (((h₁.pow 2).sub h₂).pow 2) fun s ↦ (hfin s).2

/-- The function `f(z) = z² / (2 (1 - z)²)` of the note's Lemma 1, eq. (8), which bounds the
purity defect after squaring a spectrum. -/
noncomputable def squaringDefect (z : ℝ) : ℝ := z ^ 2 / (2 * (1 - z) ^ 2)

/-- `f` is increasing on `[0, 1)` (the thermal-MPO note, after Lemma 1). -/
theorem squaringDefect_mono {z w : ℝ} (hz : 0 ≤ z) (hzw : z ≤ w) (hw : w < 1) :
    squaringDefect z ≤ squaringDefect w := by
  unfold squaringDefect
  rw [div_le_div_iff₀ (by nlinarith) (by nlinarith)]
  have h1 : z * (1 - w) ≤ w * (1 - z) := by nlinarith
  have h2 : 0 ≤ z * (1 - w) := mul_nonneg hz (by linarith)
  nlinarith [mul_le_mul h1 h1 h2 (by nlinarith)]

/-- Homogeneous form of Lemma 1: `1 - S₄/S₂² ≤ f(1 - S₂/S₁²)`. -/
theorem one_sub_div_le_squaringDefect {s₁ s₂ s₄ : ℝ} (hs₁ : 0 < s₁) (hs₂ : 0 < s₂)
    (h : 2 * (s₂ ^ 2 - s₄) ≤ (s₁ ^ 2 - s₂) ^ 2) :
    1 - s₄ / s₂ ^ 2 ≤ squaringDefect (1 - s₂ / s₁ ^ 2) := by
  unfold squaringDefect
  rw [sub_sub_cancel]
  have e1 : (1 - s₂ / s₁ ^ 2) ^ 2 / (2 * (s₂ / s₁ ^ 2) ^ 2) = (s₁ ^ 2 - s₂) ^ 2 / (2 * s₂ ^ 2) := by
    field_simp
  have e2 : 1 - s₄ / s₂ ^ 2 = (s₂ ^ 2 - s₄) / s₂ ^ 2 := by field_simp
  rw [e1, e2, div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left h (sq_nonneg s₂)]

/-- **Lemma 1 (Squaring)** of the thermal-MPO note, eq. (8). Let `pᵢ ≥ 0` with `∑ pᵢ = 1` and
purity `P = ∑ pᵢ²`. The squared spectrum `pᵢ² / P` has purity `P'` with `1 - P' ≤ f(1 - P)`. -/
theorem one_sub_squaredPurity_le_squaringDefect {κ : Type*} {p : κ → ℝ} (hp : ∀ i, 0 ≤ p i)
    (h₁ : HasSum p 1) :
    1 - ∑' i, (p i ^ 2 / ∑' j, p j ^ 2) ^ 2 ≤ squaringDefect (1 - ∑' i, p i ^ 2) := by
  have hle : ∀ i, p i ≤ 1 := fun i ↦ le_hasSum h₁ i fun j _ ↦ hp j
  have hs₂ : Summable fun i ↦ p i ^ 2 :=
    h₁.summable.of_nonneg_of_le (fun i ↦ by positivity) fun i ↦ by
      nlinarith [hp i, hle i]
  have hs₄ : Summable fun i ↦ p i ^ 4 :=
    h₁.summable.of_nonneg_of_le (fun i ↦ by positivity) fun i ↦ by
      have := pow_le_one₀ (n := 3) (hp i) (hle i)
      nlinarith [hp i, pow_nonneg (hp i) 3]
  set P := ∑' i, p i ^ 2
  have hP : 0 < P := by
    refine lt_of_le_of_ne (tsum_nonneg fun i ↦ by positivity) fun h ↦ ?_
    have hz : (fun i ↦ p i ^ 2) = 0 :=
      (hasSum_zero_iff_of_nonneg fun i ↦ by positivity).1 (by rw [h]; exact hs₂.hasSum)
    have hp0 : p = fun _ ↦ 0 := funext fun i ↦
      pow_eq_zero_iff (n := 2) (by norm_num) |>.1 (congrFun hz i)
    have : (1 : ℝ) = 0 := h₁.unique (by rw [hp0]; exact hasSum_zero)
    norm_num at this
  have e : ∑' i, (p i ^ 2 / P) ^ 2 = (∑' i, p i ^ 4) / P ^ 2 := by
    rw [← tsum_div_const]
    exact tsum_congr fun i ↦ by ring
  rw [e]
  have h := (sq_sub_sum_sq_of_hasSum hp h₁ hs₂.hasSum hs₄.hasSum).2
  simpa using one_sub_div_le_squaringDefect one_pos hP h

/-- The numerical step in eq. (12) of the note: `f(1 - (1 - ε)³) ≤ 5 ε²` for `0 ≤ ε ≤ 1/100`. -/
theorem squaringDefect_le {ε : ℝ} (hε₀ : 0 ≤ ε) (hε₁ : ε ≤ 1 / 100) :
    squaringDefect (1 - (1 - ε) ^ 3) ≤ 5 * ε ^ 2 := by
  unfold squaringDefect
  have hu1 : (1 - ε) ^ 3 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
  have hu0 : 1 - 3 * ε ≤ (1 - ε) ^ 3 := by
    nlinarith [mul_nonneg (sq_nonneg ε) (by linarith : (0 : ℝ) ≤ 3 - ε)]
  rw [sub_sub_cancel, div_le_iff₀ (by nlinarith)]
  generalize (1 - ε) ^ 3 = u at hu0 hu1 ⊢
  have h1 : (1 - u) ^ 2 ≤ (3 * ε) ^ 2 := pow_le_pow_left₀ (by linarith) (by linarith) 2
  have h2 : 9 ≤ 10 * u ^ 2 := by nlinarith
  nlinarith [mul_le_mul_of_nonneg_left h2 (sq_nonneg ε)]

/-- The squaring step used twice in Lemma 3, in power sums `S₁, S₂, S₄` of a nonnegative family:
if `S₂ ≤ S₁²`, `2 (S₂² - S₄) ≤ (S₁² - S₂)²` and the purity is `S₂/S₁² ≥ (1 - ε)³`, then the
squared family has purity `S₄/S₂² ≥ 1 - 5ε²`. -/
theorem squaring_step_of_le {s₁ s₂ s₄ : ℝ} (hs₁ : 0 < s₁) (hs₂ : 0 < s₂) (h21 : s₂ ≤ s₁ ^ 2)
    (hmain : 2 * (s₂ ^ 2 - s₄) ≤ (s₁ ^ 2 - s₂) ^ 2) {ε : ℝ} (hε₀ : 0 ≤ ε) (hε₁ : ε ≤ 1 / 100)
    (h : (1 - ε) ^ 3 ≤ s₂ / s₁ ^ 2) : 1 - 5 * ε ^ 2 ≤ s₄ / s₂ ^ 2 := by
  have hP1 : s₂ / s₁ ^ 2 ≤ 1 := (div_le_one (by positivity)).2 h21
  have hstep := one_sub_div_le_squaringDefect hs₁ hs₂ hmain
  have hu : 0 < (1 - ε) ^ 3 := pow_pos (by linarith) 3
  have hmono : squaringDefect (1 - s₂ / s₁ ^ 2) ≤ squaringDefect (1 - (1 - ε) ^ 3) :=
    squaringDefect_mono (by linarith) (by linarith) (by linarith)
  linarith [squaringDefect_le hε₀ hε₁]

/-- The squaring step for a nonnegative summable family, `squaring_step_of_le` with the two
power-sum inequalities supplied by `sq_sub_sum_sq_of_hasSum`. -/
theorem squaring_step {κ : Type*} {p : κ → ℝ} (hp : ∀ i, 0 ≤ p i) {s₁ s₂ s₄ : ℝ}
    (h₁ : HasSum p s₁) (h₂ : HasSum (fun i ↦ p i ^ 2) s₂) (h₄ : HasSum (fun i ↦ p i ^ 4) s₄)
    (hs₂ : 0 < s₂) {ε : ℝ} (hε₀ : 0 ≤ ε) (hε₁ : ε ≤ 1 / 100) (h : (1 - ε) ^ 3 ≤ s₂ / s₁ ^ 2) :
    1 - 5 * ε ^ 2 ≤ s₄ / s₂ ^ 2 := by
  obtain ⟨h21, hmain⟩ := sq_sub_sum_sq_of_hasSum hp h₁ h₂ h₄
  have hs₁ : 0 < s₁ := by
    have := hasSum_le (fun i ↦ hp i) hasSum_zero h₁
    rcases this.lt_or_eq with h' | h'
    · exact h'
    · subst h'; nlinarith
  exact squaring_step_of_le hs₁ hs₂ h21 hmain hε₀ hε₁ h

/-- Monotonicity of `ℓ^p` norms, as used in eq. (14): for nonnegative `μ` and `0 < p ≤ q`,
`∑ μᵢ^q ≤ (∑ μᵢ^p)^{q/p}`. -/
theorem hasSum_pow_le_rpow {μ : ℕ → ℝ} (hμ : ∀ i, 0 ≤ μ i) {p q : ℕ} (hp : 0 < p) (hpq : p ≤ q)
    {Zp Zq : ℝ} (hZp : HasSum (fun i ↦ μ i ^ p) Zp) (hZq : HasSum (fun i ↦ μ i ^ q) Zq) :
    Zq ≤ Zp ^ ((q : ℝ) / p) := by
  have hZp0 : 0 ≤ Zp := hasSum_le (fun i ↦ pow_nonneg (hμ i) p) hasSum_zero hZp
  set r : ℝ := (q : ℝ) / p
  have hp' : (0 : ℝ) < p := by exact_mod_cast hp
  have hr1 : 1 ≤ r := (one_le_div hp').2 (by exact_mod_cast hpq)
  have hterm : ∀ i, μ i ^ q = (μ i ^ p) ^ r := fun i ↦ by
    rw [← Real.rpow_natCast, ← Real.rpow_natCast, ← Real.rpow_mul (hμ i)]
    congr 1; simp only [r]; field_simp
  have hle : ∀ i, μ i ^ q ≤ Zp ^ (r - 1) * μ i ^ p := fun i ↦ by
    have h0 : 0 ≤ μ i ^ p := pow_nonneg (hμ i) p
    rw [hterm, show r = (r - 1) + 1 by ring, Real.rpow_add' h0 (by linarith), Real.rpow_one,
      add_sub_cancel_right]
    exact mul_le_mul_of_nonneg_right (Real.rpow_le_rpow h0
      (le_hasSum hZp i fun j _ ↦ pow_nonneg (hμ j) p) (by linarith)) h0
  calc Zq ≤ Zp ^ (r - 1) * Zp := hasSum_le hle hZq (hZp.mul_left _)
    _ = Zp ^ r := by
      rcases hZp0.lt_or_eq with h | h
      · rw [← Real.rpow_add_one h.ne']; ring_nf
      · rw [← h, mul_zero, Real.zero_rpow (by linarith)]

/-! ### The spatial moment interface

The doubling argument uses the spatial transfer operator only through three inequalities between
the moments `Z_L(β)` at one inverse temperature. `SpatialMomentBounds` records exactly these.
They follow from the moment form (`HasMomentForm.spatialMomentBounds`) and, being closed
conditions, pass to pointwise limits (`SpatialMomentBounds.of_tendsto`). The second fact lets a
moment form for approximants (for example finite Trotter transfer matrices) prove the bounds for
the exact partition functions without a moment form for the limit. -/

/-- The moment form at one inverse temperature: there are `μᵢ ≥ 0` (the moduli `|λᵢ|` of the
eigenvalues of a self-adjoint transfer operator) with `Z L = ∑ᵢ μᵢ^L` for every even `L ≥ n₀`. -/
def HasMomentForm (Z : ℕ → ℝ) (n₀ : ℕ) : Prop :=
  ∃ μ : ℕ → ℝ, (∀ i, 0 ≤ μ i) ∧ ∀ L : ℕ, n₀ ≤ L → Even L → HasSum (fun i ↦ μ i ^ L) (Z L)

/-- The three moment inequalities used by the doubling argument, at one inverse temperature:
`P_s ≤ 1`, the squaring inequality of Lemma 1 in power-sum form, and `ℓ^p` monotonicity. -/
structure SpatialMomentBounds (Z : ℕ → ℝ) (n₀ : ℕ) : Prop where
  /-- `Z_{2n} ≤ Z_n²`, that is `P_s(n) ≤ 1`. -/
  sq_le : ∀ n, n₀ ≤ n → Even n → Z (2 * n) ≤ Z n ^ 2
  /-- Lemma 1 for the spatial spectrum `σ_n`, in power sums: `2 (Z_{2n}² - Z_{4n}) ≤
  (Z_n² - Z_{2n})²`. -/
  squaring : ∀ n, n₀ ≤ n → Even n →
    2 * (Z (2 * n) ^ 2 - Z (2 * (2 * n))) ≤ (Z n ^ 2 - Z (2 * n)) ^ 2
  /-- `ℓ^p` monotonicity: `Z_q ≤ Z_p^{q/p}` for even `n₀ ≤ p ≤ q`. -/
  le_rpow : ∀ p q, n₀ ≤ p → p ≤ q → Even p → Even q → Z q ≤ Z p ^ ((q : ℝ) / p)

/-- The moment form implies the moment inequalities. -/
theorem HasMomentForm.spatialMomentBounds {Z : ℕ → ℝ} {n₀ : ℕ} (hn₀ : 0 < n₀)
    (h : HasMomentForm Z n₀) : SpatialMomentBounds Z n₀ := by
  obtain ⟨μ, hμ, hsum⟩ := h
  have key : ∀ n, n₀ ≤ n → Even n → Z (2 * n) ≤ Z n ^ 2 ∧
      2 * (Z (2 * n) ^ 2 - Z (2 * (2 * n))) ≤ (Z n ^ 2 - Z (2 * n)) ^ 2 := fun n hn hne ↦ by
    have h₁ := hsum n hn hne
    have h₂ := hsum (2 * n) (by omega) (by simp)
    have h₄ := hsum (2 * (2 * n)) (by omega) (by simp)
    simp only [pow_mul'] at h₂ h₄
    have h₄' : HasSum (fun i ↦ (μ i ^ n) ^ 4) (Z (2 * (2 * n))) := by
      convert h₄ using 1; ext i; ring
    exact sq_sub_sum_sq_of_hasSum (fun i ↦ pow_nonneg (hμ i) n) h₁ h₂ h₄'
  exact ⟨fun n hn hne ↦ (key n hn hne).1, fun n hn hne ↦ (key n hn hne).2,
    fun p q hp hpq hpe hqe ↦ hasSum_pow_le_rpow hμ (by omega) hpq (hsum p hp hpe)
      (hsum q (hp.trans hpq) hqe)⟩

/-- The moment inequalities pass to limits taken at each even length `L ≥ n₀`. -/
theorem SpatialMomentBounds.of_tendsto {α : Type*} {l : Filter α} [l.NeBot] {Zs : α → ℕ → ℝ}
    {Z : ℕ → ℝ} {n₀ : ℕ} (hn₀ : 0 < n₀) (hZs : ∀ᶠ a in l, SpatialMomentBounds (Zs a) n₀)
    (hlim : ∀ L, n₀ ≤ L → Even L → Tendsto (fun a ↦ Zs a L) l (nhds (Z L))) :
    SpatialMomentBounds Z n₀ where
  sq_le n hn hne := le_of_tendsto_of_tendsto (hlim _ (by omega) (by simp))
    ((hlim n hn hne).pow 2) (hZs.mono fun _ h ↦ h.sq_le n hn hne)
  squaring n hn hne := by
    have h₁ := hlim n hn hne
    have h₂ := hlim (2 * n) (by omega) (by simp)
    have h₄ := hlim (2 * (2 * n)) (by omega) (by simp)
    exact le_of_tendsto_of_tendsto (((h₂.pow 2).sub h₄).const_mul 2) (((h₁.pow 2).sub h₂).pow 2)
      (hZs.mono fun _ h ↦ h.squaring n hn hne)
  le_rpow p q hp hpq hpe hqe := by
    have hp0 : (0 : ℝ) < p := by exact_mod_cast hn₀.trans_le hp
    exact le_of_tendsto_of_tendsto (hlim q (hp.trans hpq) hqe)
      (Filter.Tendsto.rpow_const (hlim p hp hpe) (Or.inr (by positivity : (0 : ℝ) ≤ (q : ℝ) / p)))
      (hZs.mono fun _ h ↦ h.le_rpow p q hp hpq hpe hqe)

/-- The moment inequalities above `n₀` imply them above any `n₁ ≥ n₀`. -/
theorem SpatialMomentBounds.mono {Z : ℕ → ℝ} {n₀ n₁ : ℕ} (h : SpatialMomentBounds Z n₀)
    (hn : n₀ ≤ n₁) : SpatialMomentBounds Z n₁ :=
  ⟨fun n hn' ↦ h.sq_le n (hn.trans hn'), fun n hn' ↦ h.squaring n (hn.trans hn'),
    fun p q hp ↦ h.le_rpow p q (hn.trans hp)⟩

/-- The route for approximants: if each `Zs a` has the moment form (for instance a trace of
powers of a finite positive semidefinite transfer matrix) and `Zs a L → Z L` at every even
`L ≥ n₀`, then `Z` satisfies the moment inequalities. -/
theorem SpatialMomentBounds.of_tendsto_hasMomentForm {α : Type*} {l : Filter α} [l.NeBot]
    {Zs : α → ℕ → ℝ} {Z : ℕ → ℝ} {n₀ : ℕ} (hn₀ : 0 < n₀)
    (hZs : ∀ᶠ a in l, HasMomentForm (Zs a) n₀)
    (hlim : ∀ L, n₀ ≤ L → Even L → Tendsto (fun a ↦ Zs a L) l (nhds (Z L))) :
    SpatialMomentBounds Z n₀ :=
  of_tendsto hn₀ (hZs.mono fun _ h ↦ h.spatialMomentBounds hn₀) hlim

/-! ### Thermal and spatial purities -/

/-- The spatial purity `P_s(n, β) = Z_{2n}(β) / Z_n(β)²` of any family of partition functions
`Z L β`. Every model's spatial purity is this one for its partition function. -/
noncomputable def spatialPurityOf (Z : ℕ → ℝ → ℝ) (n : ℕ) (β : ℝ) : ℝ :=
  Z (2 * n) β / Z n β ^ 2

/-- The thermal purity `P_τ(L, β) = Z_L(2β) / Z_L(β)²` of any family of partition functions. -/
noncomputable def thermalPurityOf (Z : ℕ → ℝ → ℝ) (L : ℕ) (β : ℝ) : ℝ :=
  Z L (2 * β) / Z L β ^ 2

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

/-- The partition function `Z_L(β) = Tr e^{-β H_L} = ∑ₘ e^{-β Eₘ(L)}` of a Hamiltonian with
spectrum (with multiplicity) `E L` (the thermal-MPO note, eq. (2)). -/
noncomputable def partitionFunction (E : (L : ℕ) → ι L → ℝ) (L : ℕ) (β : ℝ) : ℝ :=
  ∑ m, exp (-(β * E L m))

/-- The thermal purity `P_τ(L, β) = Tr ρ_L(β)² = Z_L(2β) / Z_L(β)²` (the thermal-MPO note,
eq. (3)). -/
noncomputable def thermalPurity (E : (L : ℕ) → ι L → ℝ) (L : ℕ) (β : ℝ) : ℝ :=
  partitionFunction E L (2 * β) / partitionFunction E L β ^ 2

/-- The spatial purity `P_s(n, β) = Z_{2n}(β) / Z_n(β)²` (the thermal-MPO note, eq. (6)). -/
noncomputable def spatialPurity (E : (L : ℕ) → ι L → ℝ) (n : ℕ) (β : ℝ) : ℝ :=
  partitionFunction E (2 * n) β / partitionFunction E n β ^ 2

/-- The purities of a spectrum are those of its partition function. -/
theorem spatialPurity_eq_spatialPurityOf (E : (L : ℕ) → ι L → ℝ) (n : ℕ) (β : ℝ) :
    spatialPurity E n β = spatialPurityOf (partitionFunction E) n β := rfl

theorem thermalPurity_eq_thermalPurityOf (E : (L : ℕ) → ι L → ℝ) (L : ℕ) (β : ℝ) :
    thermalPurity E L β = thermalPurityOf (partitionFunction E) L β := rfl

/-- The moment form of the thermal-MPO note, eq. (5): for every `β > 0` there are nonnegative
numbers `μᵢ = |λᵢ(β)|`, the moduli of the eigenvalues of the spatial transfer operator `X_β`,
with `Z_L(β) = ∑ᵢ μᵢ^L` for every even length `L ≥ n₀`. -/
def HasSpatialMomentForm (E : (L : ℕ) → ι L → ℝ) (n₀ : ℕ) : Prop :=
  ∀ β : ℝ, 0 < β → HasMomentForm (fun L ↦ partitionFunction E L β) n₀

/-- The moment inequalities (`SpatialMomentBounds`) for the partition functions at every
`β > 0`. This is the spatial hypothesis of the doubling argument. -/
def HasSpatialMomentBounds (E : (L : ℕ) → ι L → ℝ) (n₀ : ℕ) : Prop :=
  ∀ β : ℝ, 0 < β → SpatialMomentBounds (fun L ↦ partitionFunction E L β) n₀

omit [∀ L, Nonempty (ι L)] in
/-- The moment form implies the moment inequalities at every `β > 0`. -/
theorem HasSpatialMomentForm.hasSpatialMomentBounds {E : (L : ℕ) → ι L → ℝ} {n₀ : ℕ}
    (hn₀ : 0 < n₀) (h : HasSpatialMomentForm E n₀) : HasSpatialMomentBounds E n₀ :=
  fun β hβ ↦ (h β hβ).spatialMomentBounds hn₀

/-- The spectrum `e` has a nondegenerate lowest level `e m₀` and every other level lies more than
`c` above it. For ordered eigenvalues `E₀ ≤ E₁ ≤ ⋯` this is the gap bound `E₁ - E₀ > c`. -/
def GapGT {κ : Type*} (e : κ → ℝ) (c : ℝ) : Prop :=
  ∃ m₀, ∀ m, m ≠ m₀ → e m₀ + c < e m

/-- For an ordered spectrum `E₀ ≤ E₁ ≤ ⋯`, `GapGT` gives `Δ = E₁ - E₀ > c`. -/
theorem GapGT.lt_sub_of_monotone {N : ℕ} {e : Fin (N + 2) → ℝ} (he : Monotone e) {c : ℝ}
    (h : GapGT e c) : c < e 1 - e 0 := by
  obtain ⟨m₀, hm₀⟩ := h
  by_cases h0 : m₀ = 0
  · subst h0; linarith [hm₀ 1 (by simp)]
  · have := hm₀ 0 (Ne.symm h0)
    have := he (Fin.zero_le m₀)
    have := he (Fin.zero_le (1 : Fin (N + 2)))
    linarith

section Basic

variable (E : (L : ℕ) → ι L → ℝ)

/-- `Z_L(β) > 0`. -/
theorem partitionFunction_pos (L : ℕ) (β : ℝ) : 0 < partitionFunction E L β :=
  Finset.sum_pos (fun _ _ ↦ exp_pos _) univ_nonempty

/-- `e^{-2βe} = (e^{-βe})²`. -/
theorem exp_neg_mul_two_mul (β e : ℝ) : exp (-(2 * β * e)) = exp (-(β * e)) ^ 2 := by
  rw [← exp_nat_mul]; ring_nf

omit [∀ L, Nonempty (ι L)] in
/-- `Z_L(2β) = ∑ₘ (e^{-βEₘ})²`. -/
theorem partitionFunction_two_mul (L : ℕ) (β : ℝ) :
    partitionFunction E L (2 * β) = ∑ m, exp (-(β * E L m)) ^ 2 := by
  simp only [partitionFunction, exp_neg_mul_two_mul]

omit [∀ L, Nonempty (ι L)] in
/-- `Z_L(4β) = ∑ₘ (e^{-βEₘ})⁴`. -/
theorem partitionFunction_four_mul (L : ℕ) (β : ℝ) :
    partitionFunction E L (2 * (2 * β)) = ∑ m, exp (-(β * E L m)) ^ 4 := by
  simp only [partitionFunction, exp_neg_mul_two_mul, ← pow_mul]

/-- `P_τ(L, β) > 0`. -/
theorem thermalPurity_pos (L : ℕ) (β : ℝ) : 0 < thermalPurity E L β :=
  div_pos (partitionFunction_pos E _ _) (pow_pos (partitionFunction_pos E _ _) 2)

/-- `P_s(n, β) > 0`. -/
theorem spatialPurity_pos (n : ℕ) (β : ℝ) : 0 < spatialPurity E n β :=
  div_pos (partitionFunction_pos E _ _) (pow_pos (partitionFunction_pos E _ _) 2)

/-- `P_τ(L, β) ≤ 1`. -/
theorem thermalPurity_le_one (L : ℕ) (β : ℝ) : thermalPurity E L β ≤ 1 := by
  rw [thermalPurity, div_le_one (pow_pos (partitionFunction_pos E _ _) 2),
    partitionFunction_two_mul]
  exact (sq_sub_sum_sq_finset univ _ fun m _ ↦ (exp_pos _).le).1

/-- Lemma 1 applied to the Gibbs spectrum: squaring `ρ_L(β)` gives `ρ_L(2β)`. -/
theorem thermalPurity_two_mul_ge (L : ℕ) (β : ℝ) {ε : ℝ} (hε₀ : 0 ≤ ε) (hε₁ : ε ≤ 1 / 100)
    (h : (1 - ε) ^ 3 ≤ thermalPurity E L β) : 1 - 5 * ε ^ 2 ≤ thermalPurity E L (2 * β) := by
  have hp : ∀ m, 0 ≤ exp (-(β * E L m)) := fun m ↦ (exp_pos _).le
  have h₁ := hasSum_fintype fun m ↦ exp (-(β * E L m))
  have h₂ := hasSum_fintype fun m ↦ exp (-(β * E L m)) ^ 2
  have h₄ := hasSum_fintype fun m ↦ exp (-(β * E L m)) ^ 4
  rw [← partitionFunction_two_mul] at h₂
  rw [← partitionFunction_four_mul] at h₄
  exact squaring_step hp h₁ h₂ h₄ (partitionFunction_pos E _ _) hε₀ hε₁ h

/-- **Lemma 2 (Coupling)** of the thermal-MPO note, eq. (9):
`P_s(n, 2β) = P_s(n, β)² P_τ(2n, β) / P_τ(n, β)²`. -/
theorem spatialPurity_two_mul (n : ℕ) (β : ℝ) :
    spatialPurity E n (2 * β) =
      spatialPurity E n β ^ 2 * thermalPurity E (2 * n) β / thermalPurity E n β ^ 2 := by
  have := partitionFunction_pos E n β
  have := partitionFunction_pos E n (2 * β)
  have := partitionFunction_pos E (2 * n) β
  simp only [spatialPurity, thermalPurity]
  field_simp

/-- **Lemma 2 (Coupling)** of the thermal-MPO note, eq. (10):
`P_τ(4n, β) = P_τ(2n, β)² P_s(2n, 2β) / P_s(2n, β)²`. -/
theorem thermalPurity_four_mul (n : ℕ) (β : ℝ) :
    thermalPurity E (4 * n) β =
      thermalPurity E (2 * n) β ^ 2 * spatialPurity E (2 * n) (2 * β) /
        spatialPurity E (2 * n) β ^ 2 := by
  have e : 4 * n = 2 * (2 * n) := by ring
  have := partitionFunction_pos E (2 * n) β
  have := partitionFunction_pos E (2 * n) (2 * β)
  have := partitionFunction_pos E (2 * (2 * n)) β
  simp only [spatialPurity, thermalPurity, e]
  field_simp

/-- Under the moment inequalities, `P_s(n, β) ≤ 1` for even `n ≥ n₀`. -/
theorem spatialPurity_le_one {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ} (hn : n₀ ≤ n)
    (hne : Even n) {β : ℝ} (hβ : 0 < β) : spatialPurity E n β ≤ 1 := by
  rw [spatialPurity, div_le_one (pow_pos (partitionFunction_pos E _ _) 2)]
  exact (hZ β hβ).sq_le n hn hne

/-- Lemma 1 applied to the spatial spectrum: squaring `σ_n(β)` gives `σ_{2n}(β)`. -/
theorem spatialPurity_two_mul_ge {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ}
    (hn : n₀ ≤ n) (hne : Even n) {β : ℝ} (hβ : 0 < β) {ε : ℝ} (hε₀ : 0 ≤ ε) (hε₁ : ε ≤ 1 / 100)
    (h : (1 - ε) ^ 3 ≤ spatialPurity E n β) : 1 - 5 * ε ^ 2 ≤ spatialPurity E (2 * n) β :=
  squaring_step_of_le (partitionFunction_pos E _ _) (partitionFunction_pos E _ _)
    ((hZ β hβ).sq_le n hn hne) ((hZ β hβ).squaring n hn hne) hε₀ hε₁ h

/-- The purity-to-gap bound, eq. (4) of the thermal-MPO note: if `m₀` is a lowest level and
`m ≠ m₀`, then `e^{-β (E_m - E_{m₀})} ≤ (1 - P_τ) / P_τ`. -/
theorem exp_neg_mul_sub_le_of_isMin (L : ℕ) {β : ℝ} (hβ : 0 ≤ β) {m₀ m : ι L}
    (hmin : ∀ m', E L m₀ ≤ E L m') (hne : m ≠ m₀) :
    exp (-(β * (E L m - E L m₀))) ≤ (1 - thermalPurity E L β) / thermalPurity E L β := by
  classical
  set a : ι L → ℝ := fun m ↦ exp (-(β * E L m))
  have ha : ∀ m, 0 < a m := fun m ↦ exp_pos _
  set Z := ∑ m, a m
  set S := ∑ m, a m ^ 2
  have hZ : Z = partitionFunction E L β := rfl
  have hS : S = partitionFunction E L (2 * β) := (partitionFunction_two_mul E L β).symm
  have hZpos : 0 < Z := partitionFunction_pos E L β
  have hSpos : 0 < S := hS ▸ partitionFunction_pos E L (2 * β)
  have hmax : ∀ m', a m' ≤ a m₀ := fun m' ↦
    exp_le_exp.2 (by nlinarith [hmin m'])
  have hSle : S ≤ a m₀ * Z := by
    rw [mul_sum]
    exact sum_le_sum fun m' _ ↦ by nlinarith [hmax m', ha m']
  have hpair : a m + a m₀ ≤ Z := by
    rw [← sum_pair hne]
    exact sum_le_sum_of_subset_of_nonneg (subset_univ _) fun i _ _ ↦ (ha i).le
  have hexp : exp (-(β * (E L m - E L m₀))) = a m / a m₀ := by
    simp only [a, ← exp_sub]; ring_nf
  have hP : thermalPurity E L β = S / Z ^ 2 := by rw [thermalPurity, hS, hZ]
  rw [hexp, hP, show (1 - S / Z ^ 2) / (S / Z ^ 2) = (Z ^ 2 - S) / S by field_simp,
    div_le_div_iff₀ (ha m₀) hSpos]
  nlinarith [mul_le_mul_of_nonneg_right hpair hSpos.le, mul_le_mul_of_nonneg_left hSle hZpos.le]

end Basic

/-! ### One step of the exact recurrence, interpolation, and the gap

The doubling step for arbitrary defects (`doubling_step`), the interpolation of the thermal purity
across a dyadic window (`thermalPurity_ge_of_window`), and the bound from a thermal purity to a gap
(`gapGT_of_thermalPurity`). Lemma 3 and Proposition 1 below are their equal-defect cases, and
`CertifiedGaps.Chain` runs them for finitely many steps with unequal defects. -/

/-- The coupling defect `z(a, b) = 1 - (1 - b)(1 - a)²`: if `P_s ≥ 1 - a` and `P_τ ≥ 1 - b` then
the coupled purity is at least `1 - z(a, b)`. -/
def couplingDefect (a b : ℝ) : ℝ := 1 - (1 - b) * (1 - a) ^ 2

/-- `f(1 - P) ≤ f(z)` from `0 ≤ 1 - P ≤ z < 1`, as an upper bound on a defect. -/
theorem squaringDefect_le_of_le {P z : ℝ} (hP : P ≤ 1) (hPz : 1 - z ≤ P) (hz : z < 1) :
    squaringDefect (1 - P) ≤ squaringDefect z :=
  squaringDefect_mono (by linarith) (by linarith) hz

section Steps

variable (E : (L : ℕ) → ι L → ℝ)

/-- Lemma 1 applied to the Gibbs spectrum: `1 - P_τ(L, 2β) ≤ f(1 - P_τ(L, β))`. -/
theorem one_sub_thermalPurity_two_mul_le (L : ℕ) (β : ℝ) :
    1 - thermalPurity E L (2 * β) ≤ squaringDefect (1 - thermalPurity E L β) := by
  have h₁ := hasSum_fintype fun m ↦ exp (-(β * E L m))
  have h₂ := hasSum_fintype fun m ↦ exp (-(β * E L m)) ^ 2
  have h₄ := hasSum_fintype fun m ↦ exp (-(β * E L m)) ^ 4
  rw [← partitionFunction_two_mul] at h₂
  rw [← partitionFunction_four_mul] at h₄
  obtain ⟨-, hmain⟩ := sq_sub_sum_sq_of_hasSum (fun m ↦ (exp_pos _).le) h₁ h₂ h₄
  exact one_sub_div_le_squaringDefect (partitionFunction_pos E L β)
    (partitionFunction_pos E L (2 * β)) hmain

/-- Lemma 1 applied to the spatial spectrum: `1 - P_s(2n, β) ≤ f(1 - P_s(n, β))`. -/
theorem one_sub_spatialPurity_two_mul_le {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ}
    (hn : n₀ ≤ n) (hne : Even n) {β : ℝ} (hβ : 0 < β) :
    1 - spatialPurity E (2 * n) β ≤ squaringDefect (1 - spatialPurity E n β) :=
  one_sub_div_le_squaringDefect (partitionFunction_pos E n β) (partitionFunction_pos E _ β)
    ((hZ β hβ).squaring n hn hne)

/-- **One step of the exact doubling recurrence.** If `P_s(n, β) ≥ 1 - s` and
`P_τ(2n, β) ≥ 1 - t`, then `P_s(2n, 2β) ≥ 1 - s'` and `P_τ(4n, 2β) ≥ 1 - t'` for any
`s' ≥ f(z(s, t))` and `t' ≥ f(z(t, s'))`, provided both arguments of `f` are below `1`. -/
theorem doubling_step {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ} (hn : n₀ ≤ n)
    (hne : Even n) {β : ℝ} (hβ : 0 < β) {s t s' t' : ℝ} (hs1 : s ≤ 1) (ht1 : t ≤ 1)
    (hs'1 : s' ≤ 1) (hz₁ : couplingDefect s t < 1)
    (hs' : squaringDefect (couplingDefect s t) ≤ s') (hz₂ : couplingDefect t s' < 1)
    (ht' : squaringDefect (couplingDefect t s') ≤ t') (hs : 1 - s ≤ spatialPurity E n β)
    (ht : 1 - t ≤ thermalPurity E (2 * n) β) :
    1 - s' ≤ spatialPurity E (2 * n) (2 * β) ∧ 1 - t' ≤ thermalPurity E (4 * n) (2 * β) := by
  have h1s : 0 ≤ 1 - s := by linarith
  have h1t : 0 ≤ 1 - t := by linarith
  have h1s' : 0 ≤ 1 - s' := by linarith
  -- Eq. (9) and `P_τ(n, β) ≤ 1`: `P_s(n, 2β) ≥ (1 - s)² (1 - t)`.
  have hsn : 1 - couplingDefect s t ≤ spatialPurity E n (2 * β) := by
    have hτn := thermalPurity_le_one E n β
    have hτn0 := thermalPurity_pos E n β
    rw [couplingDefect, sub_sub_cancel, spatialPurity_two_mul, le_div_iff₀ (by positivity)]
    have hsq : (1 - s) ^ 2 ≤ spatialPurity E n β ^ 2 := pow_le_pow_left₀ h1s hs 2
    have hτsq : thermalPurity E n β ^ 2 ≤ 1 := pow_le_one₀ hτn0.le hτn
    have hprod := mul_le_mul ht hsq (by positivity) (thermalPurity_pos E (2 * n) β).le
    have : 0 ≤ (1 - t) * (1 - s) ^ 2 := by positivity
    nlinarith [mul_le_mul_of_nonneg_left hτsq this]
  have hsn1 := spatialPurity_le_one E hZ hn hne (by positivity : (0 : ℝ) < 2 * β)
  have hsp : 1 - s' ≤ spatialPurity E (2 * n) (2 * β) := by
    have := one_sub_spatialPurity_two_mul_le E hZ hn hne (by positivity : (0 : ℝ) < 2 * β)
    have := squaringDefect_le_of_le hsn1 hsn hz₁
    linarith
  refine ⟨hsp, ?_⟩
  -- Eq. (10) and `P_s(2n, β) ≤ 1`: `P_τ(4n, β) ≥ (1 - t)² (1 - s')`.
  have hτ4 : 1 - couplingDefect t s' ≤ thermalPurity E (4 * n) β := by
    rw [couplingDefect, sub_sub_cancel, thermalPurity_four_mul,
      le_div_iff₀ (pow_pos (spatialPurity_pos E _ _) 2)]
    have h2n := spatialPurity_le_one E hZ (by omega : n₀ ≤ 2 * n) (by simp) hβ
    have hs2 : spatialPurity E (2 * n) β ^ 2 ≤ 1 :=
      pow_le_one₀ (spatialPurity_pos E _ _).le h2n
    have htsq : (1 - t) ^ 2 ≤ thermalPurity E (2 * n) β ^ 2 := pow_le_pow_left₀ h1t ht 2
    have hprod := mul_le_mul htsq hsp h1s' (by positivity)
    have : 0 ≤ (1 - s') * (1 - t) ^ 2 := by positivity
    nlinarith [mul_le_mul_of_nonneg_left hs2 this]
  have := one_sub_thermalPurity_two_mul_le E (4 * n) β
  have := squaringDefect_le_of_le (thermalPurity_le_one E (4 * n) β) hτ4 hz₂
  linarith

/-- **Interpolation**, eqs. (14) and (15) of the thermal-MPO note, for unequal defects: if
`P_s(n, β) ≥ 1 - s` and `P_τ(2n, β) ≥ 1 - t`, then every even `L` with `n ≤ L ≤ 2n` has
`P_τ(L, β) ≥ 1 - z(s, t)`. -/
theorem thermalPurity_ge_of_window {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ}
    (hn₀n : n₀ ≤ n) (hne : Even n) (hn0 : 0 < n) {β : ℝ} (hβ : 0 < β) {s t : ℝ} (hs1 : s ≤ 1)
    (hs : 1 - s ≤ spatialPurity E n β) (ht : 1 - t ≤ thermalPurity E (2 * n) β) {L : ℕ}
    (hnL : n ≤ L) (hL2 : L ≤ 2 * n) (hLe : Even L) :
    1 - couplingDefect s t ≤ thermalPurity E L β := by
  set A := partitionFunction E L β
  set B := partitionFunction E n β
  set C := partitionFunction E L (2 * β)
  set D := partitionFunction E (2 * n) (2 * β)
  have hA := partitionFunction_pos E L β
  have hB := partitionFunction_pos E n β
  have hC := partitionFunction_pos E L (2 * β)
  have hD := partitionFunction_pos E (2 * n) (2 * β)
  set r : ℝ := (L : ℝ) / (2 * n)
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn0
  have hLpos : (0 : ℝ) < L := by exact_mod_cast lt_of_lt_of_le hn0 hnL
  have hr0 : 0 < r := by positivity
  have hr1 : r ≤ 1 := (div_le_one (by positivity)).2 (by exact_mod_cast hL2)
  have hAB : A ≤ B ^ (2 * r) := by
    have := (hZ β hβ).le_rpow n L hn₀n hnL hne hLe
    convert this using 2; simp only [r]; field_simp
  have hDC : D ≤ C ^ r⁻¹ := by
    have := (hZ (2 * β) (by positivity)).le_rpow L (2 * n) (hn₀n.trans hnL) hL2 hLe (by simp)
    convert this using 2; simp [r]
  have hDr : D ^ r ≤ C := by
    calc D ^ r ≤ (C ^ r⁻¹) ^ r := Real.rpow_le_rpow hD.le hDC hr0.le
      _ = C := Real.rpow_inv_rpow hC.le hr0.ne'
  have hA2 : A ^ 2 ≤ (B ^ 4) ^ r := by
    calc A ^ 2 ≤ (B ^ (2 * r)) ^ 2 := pow_le_pow_left₀ hA.le hAB 2
      _ = (B ^ 4) ^ r := by
        rw [← Real.rpow_mul_natCast hB.le, ← Real.rpow_natCast_mul hB.le]; ring_nf
  set x := D / B ^ 4
  have hx : x = thermalPurity E (2 * n) β * spatialPurity E n β ^ 2 := by
    have := partitionFunction_pos E (2 * n) β
    simp only [x, D, B, thermalPurity, spatialPurity]
    field_simp
  have hx0 : 0 < x := by positivity
  have hx1 : x ≤ 1 := by
    rw [hx]
    have h1 := thermalPurity_le_one E (2 * n) β
    have h2 := spatialPurity_le_one E hZ hn₀n hne hβ
    have h3 := thermalPurity_pos E (2 * n) β
    nlinarith [pow_le_one₀ (n := 2) (spatialPurity_pos E n β).le h2]
  have hxz : 1 - couplingDefect s t ≤ x := by
    rw [hx, couplingDefect, sub_sub_cancel]
    have h1s : 0 ≤ 1 - s := by linarith
    have : (1 - s) ^ 2 ≤ spatialPurity E n β ^ 2 := pow_le_pow_left₀ h1s hs 2
    rcases le_or_gt 0 (1 - t) with h1t | h1t
    · nlinarith [mul_le_mul ht this (by positivity) (thermalPurity_pos E (2 * n) β).le]
    · have : 0 ≤ thermalPurity E (2 * n) β * spatialPurity E n β ^ 2 := by
        have := thermalPurity_pos E (2 * n) β; positivity
      nlinarith [sq_nonneg (1 - s)]
  calc 1 - couplingDefect s t ≤ x := hxz
    _ = x ^ (1 : ℝ) := (Real.rpow_one x).symm
    _ ≤ x ^ r := Real.rpow_le_rpow_of_exponent_ge hx0 hx1 hr1
    _ = D ^ r / (B ^ 4) ^ r := Real.div_rpow hD.le (by positivity) r
    _ ≤ C / A ^ 2 := div_le_div₀ hC.le hDr (by positivity) hA2
    _ = thermalPurity E L β := rfl

/-- **From a thermal purity to a gap**, eq. (4) of the thermal-MPO note: if
`P_τ(L, β) ≥ 1 - z` with `z < 1` and `e^{βc} z < 1 - z`, then `Δ_L > c`. -/
theorem gapGT_of_thermalPurity (L : ℕ) {β : ℝ} (hβ : 0 < β) {z c : ℝ} (hz : z < 1)
    (hP : 1 - z ≤ thermalPurity E L β) (hc : exp (β * c) * z < 1 - z) : GapGT (E L) c := by
  obtain ⟨m₀, -, hm₀⟩ := exists_min_image univ (E L) univ_nonempty
  refine ⟨m₀, fun m hm ↦ ?_⟩
  have h4 := exp_neg_mul_sub_le_of_isMin E L hβ.le (fun m' ↦ hm₀ m' (mem_univ _)) hm
  have hPpos := thermalPurity_pos E L β
  have h1z : 0 < 1 - z := by linarith
  -- `(1 - P)/P ≤ z/(1 - z) < e^{-βc}`.
  have hle : (1 - thermalPurity E L β) / thermalPurity E L β ≤ z / (1 - z) := by
    rw [div_le_div_iff₀ hPpos h1z]; nlinarith
  have hlt : z / (1 - z) < exp (-(β * c)) := by
    rw [div_lt_iff₀ h1z, exp_neg, ← div_eq_inv_mul, lt_div_iff₀ (exp_pos _)]
    linarith [mul_comm z (exp (β * c))]
  have := lt_of_le_of_lt (h4.trans hle) hlt
  rw [exp_lt_exp] at this
  nlinarith

end Steps

/-! ### Lemma 3 and Proposition 1

The doubling step and the criterion are proved for an arbitrary coefficient `C` bounding
`f(1 - (1 - u)³) = G(u) u²` (Proposition 4.3 of OpenAI's *The periodic spin-one Haldane gap*),
and specialized to `C = 5`, `u ≤ 1/100` (Lemma 3 and Proposition 1 of the note). -/

/-- The doubling coefficient `G(u) = ½ ((1 - u)⁻¹ + (1 - u)⁻² + (1 - u)⁻³)²`, with
`f(1 - (1 - u)³) = G(u) u²` (`squaringDefect_couplingDefect_self`). -/
noncomputable def doublingCoefficient (u : ℝ) : ℝ :=
  ((1 - u)⁻¹ + (1 - u)⁻¹ ^ 2 + (1 - u)⁻¹ ^ 3) ^ 2 / 2

/-- `G(u) > 0` for `u < 1`. -/
theorem doublingCoefficient_pos {u : ℝ} (hu : u < 1) : 0 < doublingCoefficient u := by
  have : 0 < (1 - u)⁻¹ := inv_pos.2 (by linarith)
  unfold doublingCoefficient; positivity

/-- `f(z(u, u)) = G(u) u²`. -/
theorem squaringDefect_couplingDefect_self {u : ℝ} (hu : u < 1) :
    squaringDefect (couplingDefect u u) = doublingCoefficient u * u ^ 2 := by
  have h : 1 - u ≠ 0 := by linarith
  unfold squaringDefect couplingDefect doublingCoefficient
  have h3 : 1 - (1 - (1 - u) * (1 - u) ^ 2) = (1 - u) ^ 3 := by ring
  rw [h3]
  field_simp
  ring

/-- `G` is increasing below `1`. -/
theorem doublingCoefficient_mono {u v : ℝ} (huv : u ≤ v) (hv : v < 1) :
    doublingCoefficient u ≤ doublingCoefficient v := by
  have ha : (1 - u)⁻¹ ≤ (1 - v)⁻¹ := inv_anti₀ (by linarith) (by linarith)
  have ha0 : 0 ≤ (1 - u)⁻¹ := inv_nonneg.2 (by linarith)
  unfold doublingCoefficient
  gcongr

/-- **The doubling step with coefficient `C`.** If `P_s(n, β), P_τ(2n, β) ≥ 1 - u` with `u < 1`,
`G(u) ≤ C` and `C u ≤ 1`, then `P_s(2n, 2β), P_τ(4n, 2β) ≥ 1 - C u²`. -/
theorem doubling_of_coefficient (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ}
    (hZ : HasSpatialMomentBounds E n₀) {n : ℕ} (hn : n₀ ≤ n) (hne : Even n) {β : ℝ}
    (hβ : 0 < β) {u C : ℝ} (hu0 : 0 ≤ u) (hu1 : u < 1) (hG : doublingCoefficient u ≤ C)
    (hCu : C * u ≤ 1) (hs : 1 - u ≤ spatialPurity E n β)
    (ht : 1 - u ≤ thermalPurity E (2 * n) β) :
    1 - C * u ^ 2 ≤ spatialPurity E (2 * n) (2 * β) ∧
      1 - C * u ^ 2 ≤ thermalPurity E (4 * n) (2 * β) := by
  -- `doubling_step` with `s = t = u` and `s' = t' = C u²`, where `z(u, u) = 1 - (1 - u)³`.
  have hCu2 : C * u ^ 2 ≤ u := by nlinarith
  have hCu20 : 0 ≤ C * u ^ 2 := by
    have : 0 ≤ doublingCoefficient u := by unfold doublingCoefficient; positivity
    nlinarith [sq_nonneg u]
  have hz : couplingDefect u u = 1 - (1 - u) ^ 3 := by unfold couplingDefect; ring
  have hz1 : couplingDefect u u < 1 := by
    rw [hz]; linarith [pow_pos (by linarith : (0 : ℝ) < 1 - u) 3]
  have hf : squaringDefect (couplingDefect u u) ≤ C * u ^ 2 := by
    rw [squaringDefect_couplingDefect_self hu1]
    exact mul_le_mul_of_nonneg_right hG (sq_nonneg u)
  have hsq : 0 ≤ (1 - u) ^ 2 := sq_nonneg _
  have hz2 : couplingDefect u (C * u ^ 2) ≤ couplingDefect u u := by
    unfold couplingDefect; nlinarith [mul_le_mul_of_nonneg_right hCu2 hsq]
  have hz20 : 0 ≤ couplingDefect u (C * u ^ 2) := by
    have hsq1 : (1 - u) ^ 2 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    unfold couplingDefect; nlinarith [mul_nonneg hCu20 hsq]
  exact doubling_step E hZ hn hne hβ hu1.le hu1.le (by linarith) hz1 hf (hz2.trans_lt hz1)
    ((squaringDefect_mono hz20 hz2 hz1).trans hf) hs ht

/-- `G(1/100) ≤ 5`, the coefficient of the note's eq. (12). -/
theorem doublingCoefficient_le_five {ε : ℝ} (hε₁ : ε ≤ 1 / 100) :
    doublingCoefficient ε ≤ 5 :=
  (doublingCoefficient_mono hε₁ (by norm_num)).trans (by
    unfold doublingCoefficient; norm_num)

/-- **Lemma 3 (Doubling)** of the thermal-MPO note. Assume the moment form (5), let `n` be even
and `0 ≤ ε ≤ 1/100`. If `P_s(n, β), P_τ(2n, β) ≥ 1 - ε`, then
`P_s(2n, 2β), P_τ(4n, 2β) ≥ 1 - 5ε²`. -/
theorem doubling (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ}
    (hn : n₀ ≤ n) (hne : Even n) {β : ℝ} (hβ : 0 < β) {ε : ℝ} (hε₀ : 0 ≤ ε)
    (hε₁ : ε ≤ 1 / 100) (hs : 1 - ε ≤ spatialPurity E n β)
    (ht : 1 - ε ≤ thermalPurity E (2 * n) β) :
    1 - 5 * ε ^ 2 ≤ spatialPurity E (2 * n) (2 * β) ∧
      1 - 5 * ε ^ 2 ≤ thermalPurity E (4 * n) (2 * β) :=
  doubling_of_coefficient E hZ hn hne hβ hε₀ (by linarith) (doublingCoefficient_le_five hε₁)
    (by linarith) hs ht

/-- Every even length `L ≥ n₀ > 0` lies in a dyadic window `n₀ 2^j ≤ L ≤ 2 n₀ 2^j`. -/
theorem exists_dyadic_window {n₀ L : ℕ} (hn₀ : 0 < n₀) (hL : n₀ ≤ L) :
    ∃ j, n₀ * 2 ^ j ≤ L ∧ L ≤ 2 * (n₀ * 2 ^ j) := by
  have hq : L / n₀ ≠ 0 := (Nat.div_pos hL hn₀).ne'
  refine ⟨Nat.log 2 (L / n₀), ?_, ?_⟩
  · calc n₀ * 2 ^ Nat.log 2 (L / n₀) ≤ n₀ * (L / n₀) :=
          Nat.mul_le_mul_left _ (Nat.pow_log_le_self 2 hq)
      _ ≤ L := Nat.mul_div_le L n₀
  · have h1 := Nat.lt_pow_succ_log_self one_lt_two (L / n₀)
    have h2 : L < n₀ * (L / n₀ + 1) := by
      rw [mul_comm]; exact Nat.lt_mul_of_div_lt (Nat.lt_succ_self _) hn₀
    calc L ≤ n₀ * (L / n₀ + 1) := h2.le
      _ ≤ n₀ * 2 ^ (Nat.log 2 (L / n₀)).succ := Nat.mul_le_mul_left _ h1
      _ = 2 * (n₀ * 2 ^ Nat.log 2 (L / n₀)) := by rw [pow_succ]; ring

/-- The defects `u_j = (C u₀)^{2^j} / C` of the doubling induction with coefficient `C`, so that
`u₀` is the seed defect and `u_{j+1} = C u_j²`. -/
noncomputable def coefficientDefect (C u₀ : ℝ) (j : ℕ) : ℝ := (C * u₀) ^ 2 ^ j / C

/-- `u_{j+1} = C u_j²`. -/
theorem coefficientDefect_succ {C : ℝ} (hC : C ≠ 0) (u₀ : ℝ) (j : ℕ) :
    coefficientDefect C u₀ (j + 1) = C * coefficientDefect C u₀ j ^ 2 := by
  simp only [coefficientDefect, pow_succ 2 j, pow_mul]
  field_simp

/-- `0 < u_j ≤ u₀` when `0 < C u₀ ≤ 1`. -/
theorem coefficientDefect_mem {C u₀ : ℝ} (hC : 0 < C) (hu₀ : 0 < u₀) (hr : C * u₀ ≤ 1)
    (j : ℕ) : 0 < coefficientDefect C u₀ j ∧ coefficientDefect C u₀ j ≤ u₀ := by
  refine ⟨by unfold coefficientDefect; positivity, ?_⟩
  unfold coefficientDefect
  rw [div_le_iff₀ hC, mul_comm u₀ C]
  exact pow_le_of_le_one (by positivity) hr (pow_pos two_pos j).ne'

/-- The induction on scale with coefficient `C`: with `n_j = n₀ 2^j`, `β_j = β₀ 2^j`,
`P_s(n_j, β_j), P_τ(2n_j, β_j) ≥ 1 - u_j`. -/
theorem purities_ge_coefficientDefect (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : Even n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀) {u₀ C : ℝ} (hu₀ : 0 < u₀)
    (hu₁ : u₀ < 1) (hG : doublingCoefficient u₀ ≤ C) (hr : C * u₀ ≤ 1)
    (hs : 1 - u₀ ≤ spatialPurity E n₀ β₀) (ht : 1 - u₀ ≤ thermalPurity E (2 * n₀) β₀)
    (j : ℕ) :
    1 - coefficientDefect C u₀ j ≤ spatialPurity E (n₀ * 2 ^ j) (β₀ * 2 ^ j) ∧
      1 - coefficientDefect C u₀ j ≤ thermalPurity E (2 * (n₀ * 2 ^ j)) (β₀ * 2 ^ j) := by
  have hC : 0 < C := lt_of_lt_of_le (doublingCoefficient_pos (by linarith)) hG
  induction j with
  | zero =>
    simp only [coefficientDefect, pow_zero, pow_one, mul_one]
    rw [mul_div_cancel_left₀ _ hC.ne']
    exact ⟨hs, ht⟩
  | succ j ih =>
    obtain ⟨hu0, hu⟩ := coefficientDefect_mem hC hu₀ hr j
    have h := doubling_of_coefficient E hZ (n := n₀ * 2 ^ j)
      (Nat.le_mul_of_pos_right _ (by positivity)) (hn₀.mul_right _) (by positivity) hu0.le
      (hu.trans_lt hu₁) ((doublingCoefficient_mono hu hu₁).trans hG)
      ((mul_le_mul_of_nonneg_left hu hC.le).trans hr) ih.1 ih.2
    rw [coefficientDefect_succ hC.ne']
    have e1 : n₀ * 2 ^ (j + 1) = 2 * (n₀ * 2 ^ j) := by ring
    have e2 : 2 * (n₀ * 2 ^ (j + 1)) = 4 * (n₀ * 2 ^ j) := by ring
    have e3 : β₀ * 2 ^ (j + 1) = 2 * (β₀ * 2 ^ j) := by ring
    rw [e1, e3]; rw [e1] at e2; rw [e2]
    exact h

/-- **The two-purity criterion with coefficient `C`** (Proposition 4.3 of OpenAI's *The
periodic spin-one Haldane gap*). Assume the moment inequalities above an even `n₀ > 0`, and let
`0 < u₀ < 1/6`, `G(u₀) ≤ C`, `r = C u₀ < 1` and `C (1 - 3u₀) > 3`. If
`P_s(n₀, β₀), P_τ(2n₀, β₀) ≥ 1 - u₀`, then every even `L ≥ n₀` has a nondegenerate ground level
and gap `Δ_L > log (1/r) / β₀`. -/
theorem gapGT_of_coefficient (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀) (hn₀e : Even n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀) {u₀ C : ℝ} (hu₀ : 0 < u₀)
    (hu₁ : u₀ < 1 / 6) (hG : doublingCoefficient u₀ ≤ C) (hr : C * u₀ < 1)
    (hC3 : 3 < C * (1 - 3 * u₀)) (hs : 1 - u₀ ≤ spatialPurity E n₀ β₀)
    (ht : 1 - u₀ ≤ thermalPurity E (2 * n₀) β₀) {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) :
    GapGT (E L) (log (1 / (C * u₀)) / β₀) := by
  have hC : 0 < C := lt_of_lt_of_le (doublingCoefficient_pos (by linarith)) hG
  -- The dyadic window `n_j ≤ L ≤ 2 n_j` and the purities `≥ 1 - u_j` at scale `j`.
  obtain ⟨j, hjL, hLj⟩ := exists_dyadic_window hn₀ hL
  obtain ⟨hsj, htj⟩ := purities_ge_coefficientDefect E hn₀e hβ₀ hZ hu₀ (by linarith) hG hr.le
    hs ht j
  set u := coefficientDefect C u₀ j
  obtain ⟨hu0, hu⟩ := coefficientDefect_mem hC hu₀ hr.le j
  -- Interpolation: `P_τ(L, β_j) ≥ 1 - z(u_j, u_j)` with `0 ≤ z(u, u) ≤ 3u`.
  have hP := thermalPurity_ge_of_window E hZ (n := n₀ * 2 ^ j)
    (Nat.le_mul_of_pos_right _ (by positivity)) (hn₀e.mul_right _) (by positivity)
    (by positivity : (0 : ℝ) < β₀ * 2 ^ j) (by linarith) hsj htj hjL hLj hLe
  have hz3 : couplingDefect u u ≤ 3 * u := by
    unfold couplingDefect; nlinarith [sq_nonneg u, mul_nonneg hu0.le (sq_nonneg u)]
  have hz0 : 0 ≤ couplingDefect u u := by
    unfold couplingDefect
    nlinarith [pow_le_one₀ (n := 2) (by linarith : (0 : ℝ) ≤ 1 - u) (by linarith),
      sq_nonneg (1 - u)]
  refine gapGT_of_thermalPurity E L (by positivity) (by linarith) hP ?_
  -- From purity to a gap: `e^{β_j c} = (C u₀)^{-2^j} = 1/(C u_j)`, and `3/C < 1 - 3u_j`.
  have hexp : exp (β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀)) * (C * u) = 1 := by
    rw [show β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀) = (2 ^ j : ℕ) * log (1 / (C * u₀)) by
        field_simp; push_cast; ring, exp_nat_mul, exp_log (by positivity)]
    simp only [u, coefficientDefect]
    rw [mul_div_cancel₀ _ hC.ne', ← mul_pow]
    rw [one_div_mul_cancel (by positivity), one_pow]
  have hE := exp_pos (β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀))
  set e := exp (β₀ * 2 ^ j * (log (1 / (C * u₀)) / β₀))
  have h1 : e * couplingDefect u u ≤ 3 / C := by
    rw [le_div_iff₀ hC]
    nlinarith [mul_le_mul_of_nonneg_left hz3 hE.le]
  have h2 : 3 / C < 1 - 3 * u := by
    rw [div_lt_iff₀ hC]; nlinarith
  linarith

/-- The defect sequence `ε_j = (5 ε₀)^{2^j} / 5` of the doubling induction. For `ε₀ = 1/100`
it is `ε_j = 20^{-2^j} / 5`, eq. (13) of the thermal-MPO note. -/
noncomputable def defect (ε₀ : ℝ) (j : ℕ) : ℝ := (1 / 5) * (5 * ε₀) ^ 2 ^ j

/-- `ε_0 = ε₀`. -/
@[simp] theorem defect_zero (ε₀ : ℝ) : defect ε₀ 0 = ε₀ := by simp [defect]

/-- `ε_{j+1} = 5 ε_j²`. -/
theorem defect_succ (ε₀ : ℝ) (j : ℕ) : defect ε₀ (j + 1) = 5 * defect ε₀ j ^ 2 := by
  simp only [defect, pow_succ 2 j, pow_mul]; ring

/-- `ε_j > 0` for `ε₀ > 0`. -/
theorem defect_pos {ε₀ : ℝ} (hε₀ : 0 < ε₀) (j : ℕ) : 0 < defect ε₀ j := by
  unfold defect; positivity

/-- `(5 ε₀)^{2^j} ≤ 5 ε₀` for `0 ≤ ε₀ ≤ 1/5`. -/
theorem five_mul_pow_le {ε₀ : ℝ} (hε₀ : 0 ≤ ε₀) (hε₁ : ε₀ ≤ 1 / 5) (j : ℕ) :
    (5 * ε₀) ^ 2 ^ j ≤ 5 * ε₀ :=
  pow_le_of_le_one (by positivity) (by linarith) (pow_pos two_pos j).ne'

/-- `ε_j ≤ 1/100` for `0 ≤ ε₀ ≤ 1/100`. -/
theorem defect_le {ε₀ : ℝ} (hε₀ : 0 ≤ ε₀) (hε₁ : ε₀ ≤ 1 / 100) (j : ℕ) :
    defect ε₀ j ≤ 1 / 100 := by
  have := five_mul_pow_le hε₀ (by linarith) j
  unfold defect; linarith

/-- The defects with coefficient `5` are the `ε_j`. -/
theorem coefficientDefect_five (ε₀ : ℝ) (j : ℕ) : coefficientDefect 5 ε₀ j = defect ε₀ j := by
  simp only [coefficientDefect, defect]; ring

/-- The induction on scale in the proof of Proposition 1: with `n_j = n₀ 2^j`, `β_j = β₀ 2^j`,
`P_s(n_j, β_j), P_τ(2n_j, β_j) ≥ 1 - ε_j`. -/
theorem purities_ge_defect (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : Even n₀) {β₀ : ℝ}
    (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀) {ε₀ : ℝ} (hε₀ : 0 < ε₀)
    (hε₁ : ε₀ ≤ 1 / 100) (hs : 1 - ε₀ ≤ spatialPurity E n₀ β₀)
    (ht : 1 - ε₀ ≤ thermalPurity E (2 * n₀) β₀) (j : ℕ) :
    1 - defect ε₀ j ≤ spatialPurity E (n₀ * 2 ^ j) (β₀ * 2 ^ j) ∧
      1 - defect ε₀ j ≤ thermalPurity E (2 * (n₀ * 2 ^ j)) (β₀ * 2 ^ j) := by
  rw [← coefficientDefect_five]
  exact purities_ge_coefficientDefect E hn₀ hβ₀ hZ hε₀ (by linarith)
    (doublingCoefficient_le_five hε₁) (by linarith) hs ht j

/-- **Proposition 1 (Two-purity criterion)** of the thermal-MPO note, in the general form used
for all later refinements. Assume the moment inequalities `HasSpatialMomentBounds E n₀`. Let
`n₀ > 0` be even, `β₀ > 0` and `0 < ε₀ ≤ 1/100`, with seeds `P_s(n₀, β₀), P_τ(2n₀, β₀) ≥ 1 - ε₀`.
Then every even `L ≥ n₀` has a nondegenerate ground level and gap `Δ_L > log (1/(5ε₀)) / β₀`. -/
theorem gapGT_of_seeds (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀) (hn₀e : Even n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀) {ε₀ : ℝ} (hε₀ : 0 < ε₀)
    (hε₁ : ε₀ ≤ 1 / 100) (hs : 1 - ε₀ ≤ spatialPurity E n₀ β₀)
    (ht : 1 - ε₀ ≤ thermalPurity E (2 * n₀) β₀) {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) :
    GapGT (E L) (log (1 / (5 * ε₀)) / β₀) :=
  gapGT_of_coefficient E hn₀ hn₀e hβ₀ hZ hε₀ (by linarith)
    (doublingCoefficient_le_five hε₁) (by linarith) (by linarith) hs ht hL hLe


/-- **Proposition 1 (Two-purity criterion)** of the thermal-MPO note, as stated there. Assume the
moment form (5). Let `n₀` be even and `β₀ > 0`, with `P_s(n₀, β₀), P_τ(2n₀, β₀) ≥ 0.99`. Then the
spectral gap satisfies `Δ_L > (log 20) / β₀` for every even `L ≥ n₀`.

The note assumes `n₀ ≥ 4`; the proof uses only `n₀ > 0`. -/
theorem gapGT_of_purities (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀) (hn₀e : Even n₀)
    {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentForm E n₀)
    (hs : 99 / 100 ≤ spatialPurity E n₀ β₀) (ht : 99 / 100 ≤ thermalPurity E (2 * n₀) β₀)
    {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) : GapGT (E L) (log 20 / β₀) := by
  have h := gapGT_of_seeds E hn₀ hn₀e hβ₀ (hZ.hasSpatialMomentBounds hn₀) (ε₀ := 1 / 100)
    (by norm_num) le_rfl (by linarith) (by linarith) hL hLe
  norm_num at h
  exact h

end CertifiedGaps
