/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.PurityCriterion
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Topology.MetricSpace.Sequences

/-!
# Boundary moment inequalities for open chains

An open chain of `n` sites with boundary terms at both ends has the partition function
`G n = ⟨v, Xⁿ v⟩`, where `X` is the spatial transfer operator of the bulk (the same one whose
traces give the ring partition functions `Z L = Tr X^L`) and `v` is a boundary vector. Writing
`X = ∑ᵢ λᵢ |ψᵢ⟩⟨ψᵢ|`, this is the *boundary moment form*

  `Z (2N) = ∑ᵢ λᵢ^{2N}`,   `G n = ∑ᵢ wᵢ λᵢⁿ`,   `wᵢ = |⟨ψᵢ, v⟩|² ≥ 0`.

The open-chain argument of the note (Appendix A) uses this form only through finitely many
closed inequalities involving a leading eigenvalue modulus `μ = |λ₀|`, its boundary weight
`c = w₀`, the sign `s` of `λ₀`, and the ratio `ρ ≤ 1` of the next modulus to `μ`.
`OpenMomentWitness Z G μ c ρ s` records them, and `OpenMomentBounds Z G` says that some
witness exists. As in `CertifiedGaps.SpatialMomentBounds` for rings, these inequalities follow
from the moment form of a finite transfer matrix (`OpenMomentBounds.of_hermitian`) and pass to
limits (`OpenMomentBounds.of_tendsto`, by compactness of the witnesses). Finite Trotter
approximants thus prove them for the exact partition functions.

With `B m = G m / μ^m`, the inequalities say that `B m = c + D m` at even `m`, where the
remainder `D m ≥ 0` decays like `ρ^m` (`G_decay`), and that at odd `n` the value `B n` lies
within `D M` of `s c` for every even `M ≤ n` (`G_odd`).

## Main definitions

* `OpenMomentWitness Z G μ c ρ s`, `OpenMomentBounds Z G`.
* `leadingRoot Z = inf_N Z(2N)^{1/(2N)}`: the leading modulus `μ`, in terms of ring traces.
* `leadingOverlap Z G = inf_N G(2N) / μ^{2N}`: the leading boundary weight `c`.
* `boundaryDefect Z G k = 1 - c μ^k / G k`: the boundary defect `e_k` of the note, eq. (A13).

## Main results

* `OpenMomentBounds.of_eigen`, `OpenMomentBounds.of_hermitian`: the bounds from a finite
  boundary moment form.
* `OpenMomentBounds.of_tendsto`: the bounds pass to limits.
* `OpenMomentWitness.leadingRoot_eq`: `μ = leadingRoot Z` for every witness.
* `OpenMomentWitness.leadingOverlap_eq`: `c = leadingOverlap Z G` for every witness with `ρ < 1`.
-/

open Filter Finset Real Topology

namespace CertifiedGaps

/-- The closed inequalities of the boundary moment form at one inverse temperature, for the ring
traces `Z` and the open-chain traces `G`, with leading modulus `μ`, leading boundary weight `c`,
residual ratio `ρ` and leading sign `s`. Even lengths are written `2N`, odd ones `2N + 1`. -/
structure OpenMomentWitness (Z G : ℕ → ℝ) (μ c ρ s : ℝ) : Prop where
  mu_nonneg : 0 ≤ μ
  c_nonneg : 0 ≤ c
  rho_nonneg : 0 ≤ ρ
  rho_le_one : ρ ≤ 1
  /-- The sign is `±1`. -/
  s_sq : s ^ 2 = 1
  /-- The leading term of a ring trace: `μ^{2N} ≤ Z(2N)`. -/
  pow_le_Z : ∀ N, 1 ≤ N → μ ^ (2 * N) ≤ Z (2 * N)
  /-- `Z(2N)/μ^{2N}` decreases in `N`. -/
  Z_anti : ∀ K N, 1 ≤ K → K ≤ N → μ ^ (2 * K) * Z (2 * N) ≤ μ ^ (2 * N) * Z (2 * K)
  /-- The purity of `σ_{2N}` is at most its largest weight: `Z(4N) ≤ μ^{2N} Z(2N)`. -/
  Z_two_mul_le : ∀ N, 1 ≤ N → Z (2 * (2 * N)) ≤ μ ^ (2 * N) * Z (2 * N)
  /-- The purity of `σ_{2N}` is at most `p₀² + (1 - p₀)²`. -/
  Z_two_mul_le_sq : ∀ N, 1 ≤ N →
    Z (2 * (2 * N)) ≤ μ ^ (2 * (2 * N)) + (Z (2 * N) - μ ^ (2 * N)) ^ 2
  /-- The next modulus `ρ μ` appears among the remaining terms: `(ρ μ)^{2N} ≤ Z(2N) - μ^{2N}`. -/
  tail : ∀ N, 1 ≤ N → (ρ * μ) ^ (2 * N) ≤ Z (2 * N) - μ ^ (2 * N)
  /-- The leading term of an even open-chain trace: `c μ^{2N} ≤ G(2N)`. -/
  le_G : ∀ N, 1 ≤ N → c * μ ^ (2 * N) ≤ G (2 * N)
  /-- The remainder `D(2N) = G(2N)/μ^{2N} - c` decays: `D(2N) ≤ ρ^{2N - 2K} D(2K)`. -/
  G_decay : ∀ K N, 1 ≤ K → K ≤ N →
    μ ^ (2 * K) * (G (2 * N) - c * μ ^ (2 * N)) ≤
      ρ ^ (2 * N - 2 * K) * μ ^ (2 * N) * (G (2 * K) - c * μ ^ (2 * K))
  /-- An odd open-chain trace is `s c μ^{2N+1}` up to the even remainder `D(2K)`. -/
  G_odd : ∀ K N, 1 ≤ K → K ≤ N →
    |G (2 * N + 1) - s * c * μ ^ (2 * N + 1)| ≤
      μ ^ (2 * N + 1 - 2 * K) * (G (2 * K) - c * μ ^ (2 * K))

/-- The boundary moment inequalities: some witness `(μ, c, ρ, s)` exists. -/
def OpenMomentBounds (Z G : ℕ → ℝ) : Prop :=
  ∃ μ c ρ s : ℝ, OpenMomentWitness Z G μ c ρ s

/-! ### From a finite moment form -/

/-- For `0 ≤ x ≤ y` and `K ≤ N`, `y^K x^N ≤ y^N x^K`: raising to the extra `N - K` powers gains
at least as much for the larger base. -/
private theorem pow_mul_pow_le {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) {K N : ℕ} (hKN : K ≤ N) :
    y ^ K * x ^ N ≤ y ^ N * x ^ K := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hKN
  have h0 : 0 ≤ y ^ K * x ^ K := mul_nonneg (pow_nonneg (hx.trans hxy) K) (pow_nonneg hx K)
  calc y ^ K * x ^ (K + d) = y ^ K * x ^ K * x ^ d := by ring
    _ ≤ y ^ K * x ^ K * y ^ d := mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hx hxy d) h0
    _ = y ^ (K + d) * x ^ K := by ring

/-- The bounds from a finite boundary moment form `Z(2N) = ∑ᵢ λᵢ^{2N}`, `G n = ∑ᵢ wᵢ λᵢⁿ` with
`wᵢ ≥ 0`, for `N ≥ 1` and `n ≥ 2`. -/
theorem OpenMomentBounds.of_eigen {κ : Type*} [Fintype κ] [Nonempty κ] (ev w : κ → ℝ)
    (hw : ∀ i, 0 ≤ w i) {Z G : ℕ → ℝ} (hZ : ∀ N, 1 ≤ N → Z (2 * N) = ∑ i, ev i ^ (2 * N))
    (hG : ∀ n, 2 ≤ n → G n = ∑ i, w i * ev i ^ n) : OpenMomentBounds Z G := by
  classical
  -- The leading index `i₀` has the largest modulus `μ`; the others form `T`.
  obtain ⟨i₀, -, hi₀⟩ := exists_max_image univ (fun i ↦ |ev i|) univ_nonempty
  set T := univ.erase i₀
  set μ := |ev i₀|
  have hμ0 : 0 ≤ μ := abs_nonneg _
  -- The next modulus `r` bounds every `|λᵢ|` with `i ∈ T` and is attained there (or is `0`).
  obtain ⟨r, hr0, hrμ, hrT, hr⟩ : ∃ r, 0 ≤ r ∧ r ≤ μ ∧ (∀ i ∈ T, |ev i| ≤ r) ∧
      (r = 0 ∨ ∃ i ∈ T, |ev i| = r) := by
    rcases T.eq_empty_or_nonempty with hTe | hTne
    · exact ⟨0, le_rfl, hμ0, by simp [hTe], Or.inl rfl⟩
    · obtain ⟨i₁, hi₁, hmax⟩ := exists_max_image T (fun i ↦ |ev i|) hTne
      exact ⟨|ev i₁|, abs_nonneg _, hi₀ i₁ (mem_univ _), hmax, Or.inr ⟨i₁, hi₁, rfl⟩⟩
  set s : ℝ := if 0 ≤ ev i₀ then 1 else -1
  have hs : s ^ 2 = 1 := by simp only [s]; split_ifs <;> norm_num
  have hev₀ : ev i₀ = s * μ := by
    simp only [s, μ]
    split_ifs with h
    · rw [abs_of_nonneg h, one_mul]
    · rw [abs_of_neg (not_le.1 h)]; ring
  -- `ρ = r / μ` (which is `0` when `μ = 0`, since then `r = 0` too).
  have hρμ : r / μ * μ = r := by
    rcases hμ0.eq_or_lt with h | h
    · rw [← h, mul_zero]; linarith [h ▸ hrμ]
    · exact div_mul_cancel₀ r h.ne'
  refine ⟨μ, w i₀, r / μ, s, ?_⟩
  -- Everything is a sum of powers of `pᵢ = λᵢ² ≤ μ²`, with `pᵢ ≤ ρ² μ²` off `i₀`.
  have hp0 : ∀ i, 0 ≤ ev i ^ 2 := fun i ↦ sq_nonneg _
  have hpm : ∀ i, ev i ^ 2 ≤ μ ^ 2 := fun i ↦ by
    rw [← sq_abs (ev i)]; exact pow_le_pow_left₀ (abs_nonneg _) (hi₀ i (mem_univ _)) 2
  have hpr : ∀ i ∈ T, ev i ^ 2 ≤ (r / μ) ^ 2 * μ ^ 2 := fun i hi ↦ by
    rw [← mul_pow, hρμ, ← sq_abs (ev i)]; exact pow_le_pow_left₀ (abs_nonneg _) (hrT i hi) 2
  have hsplit : ∀ f : κ → ℝ, ∑ i, f i = f i₀ + ∑ i ∈ T, f i := fun f ↦
    (add_sum_erase _ _ (mem_univ i₀)).symm
  have hZ' : ∀ N, 1 ≤ N → Z (2 * N) = μ ^ (2 * N) + ∑ i ∈ T, (ev i ^ 2) ^ N := fun N hN ↦ by
    rw [hZ N hN, hsplit, hev₀, mul_pow, pow_mul s, hs, one_pow, one_mul]
    simp only [pow_mul]
  have hGe : ∀ N, 1 ≤ N →
      G (2 * N) = w i₀ * μ ^ (2 * N) + ∑ i ∈ T, w i * (ev i ^ 2) ^ N := fun N hN ↦ by
    rw [hG _ (by omega), hsplit, hev₀, mul_pow, pow_mul s, hs, one_pow, one_mul]
    simp only [pow_mul]
  have hGo : ∀ N, 1 ≤ N → G (2 * N + 1) =
      s * w i₀ * μ ^ (2 * N + 1) + ∑ i ∈ T, w i * ev i * (ev i ^ 2) ^ N := fun N hN ↦ by
    rw [hG _ (by omega), hsplit]
    congr 1
    · rw [hev₀, mul_pow, pow_succ, pow_mul, hs, one_pow]; ring
    · exact sum_congr rfl fun i _ ↦ by rw [pow_succ, pow_mul]; ring
  have hTZ : ∀ N, 0 ≤ ∑ i ∈ T, (ev i ^ 2) ^ N := fun N ↦
    sum_nonneg fun i _ ↦ pow_nonneg (hp0 i) N
  have hTG : ∀ N, 0 ≤ ∑ i ∈ T, w i * (ev i ^ 2) ^ N := fun N ↦
    sum_nonneg fun i _ ↦ mul_nonneg (hw i) (pow_nonneg (hp0 i) N)
  refine
    { mu_nonneg := hμ0
      c_nonneg := hw i₀
      rho_nonneg := div_nonneg hr0 hμ0
      rho_le_one := div_le_one_of_le₀ hrμ hμ0
      s_sq := hs
      pow_le_Z := fun N hN ↦ by rw [hZ' N hN]; linarith [hTZ N]
      Z_anti := fun K N hK hKN ↦ by
        rw [hZ' N (hK.trans hKN), hZ' K hK]
        have h : μ ^ (2 * K) * ∑ i ∈ T, (ev i ^ 2) ^ N ≤
            μ ^ (2 * N) * ∑ i ∈ T, (ev i ^ 2) ^ K := by
          rw [mul_sum, mul_sum, pow_mul, pow_mul]
          exact sum_le_sum fun i _ ↦ pow_mul_pow_le (hp0 i) (hpm i) hKN
        linear_combination h
      Z_two_mul_le := fun N hN ↦ by
        rw [hZ' _ (by omega), hZ' N hN]
        have h : ∑ i ∈ T, (ev i ^ 2) ^ (2 * N) ≤ μ ^ (2 * N) * ∑ i ∈ T, (ev i ^ 2) ^ N := by
          rw [mul_sum, pow_mul]
          refine sum_le_sum fun i _ ↦ ?_
          rw [two_mul, pow_add]
          exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (hp0 i) (hpm i) N)
            (pow_nonneg (hp0 i) N)
        linear_combination h
      Z_two_mul_le_sq := fun N hN ↦ by
        rw [hZ' _ (by omega), hZ' N hN, add_sub_cancel_left]
        have h : ∑ i ∈ T, (ev i ^ 2) ^ (2 * N) ≤ (∑ i ∈ T, (ev i ^ 2) ^ N) ^ 2 :=
          calc ∑ i ∈ T, (ev i ^ 2) ^ (2 * N) = ∑ i ∈ T, ((ev i ^ 2) ^ N) ^ 2 :=
                sum_congr rfl fun i _ ↦ by ring
            _ ≤ _ := (sq_sub_sum_sq_finset T _ fun i _ ↦ pow_nonneg (hp0 i) N).1
        linear_combination h
      tail := fun N hN ↦ by
        rw [hZ' N hN, add_sub_cancel_left, hρμ]
        rcases hr with rfl | ⟨i₁, hi₁, rfl⟩
        · rw [zero_pow (by omega)]; exact hTZ N
        · rw [pow_mul, sq_abs]
          exact single_le_sum (f := fun i ↦ (ev i ^ 2) ^ N)
            (fun i _ ↦ pow_nonneg (hp0 i) N) hi₁
      le_G := fun N hN ↦ by rw [hGe N hN]; linarith [hTG N]
      G_decay := fun K N hK hKN ↦ by
        obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hKN
        rw [hGe _ (hK.trans hKN), hGe K hK, show 2 * (K + d) - 2 * K = 2 * d by omega,
          add_sub_cancel_left, add_sub_cancel_left, mul_sum, mul_sum]
        refine sum_le_sum fun i hi ↦ ?_
        have h0 : 0 ≤ (μ ^ 2) ^ K * w i * (ev i ^ 2) ^ K :=
          mul_nonneg (mul_nonneg (pow_nonneg (sq_nonneg μ) K) (hw i)) (pow_nonneg (hp0 i) K)
        calc μ ^ (2 * K) * (w i * (ev i ^ 2) ^ (K + d))
            = (μ ^ 2) ^ K * w i * (ev i ^ 2) ^ K * (ev i ^ 2) ^ d := by ring
          _ ≤ (μ ^ 2) ^ K * w i * (ev i ^ 2) ^ K * ((r / μ) ^ 2 * μ ^ 2) ^ d :=
            mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (hp0 i) (hpr i hi) d) h0
          _ = (r / μ) ^ (2 * d) * μ ^ (2 * (K + d)) * (w i * (ev i ^ 2) ^ K) := by ring
      G_odd := fun K N hK hKN ↦ by
        obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hKN
        rw [hGo _ (hK.trans hKN), hGe K hK, add_sub_cancel_left, add_sub_cancel_left,
          show 2 * (K + d) + 1 - 2 * K = 2 * d + 1 by omega, mul_sum]
        refine (abs_sum_le_sum_abs _ _).trans (sum_le_sum fun i _ ↦ ?_)
        rw [abs_mul, abs_mul, abs_of_nonneg (hw i), abs_of_nonneg (pow_nonneg (hp0 i) _)]
        have h0 : 0 ≤ w i * (ev i ^ 2) ^ K := mul_nonneg (hw i) (pow_nonneg (hp0 i) K)
        calc w i * |ev i| * (ev i ^ 2) ^ (K + d)
            = (w i * (ev i ^ 2) ^ K) * (|ev i| * (ev i ^ 2) ^ d) := by ring
          _ ≤ (w i * (ev i ^ 2) ^ K) * (μ * (μ ^ 2) ^ d) :=
            mul_le_mul_of_nonneg_left (mul_le_mul (hi₀ i (mem_univ _))
              (pow_le_pow_left₀ (hp0 i) (hpm i) d) (pow_nonneg (hp0 i) d) hμ0) h0
          _ = μ ^ (2 * d + 1) * (w i * (ev i ^ 2) ^ K) := by ring }

/-- The bounds for the traces of a Hermitian transfer matrix `X` and a boundary vector `v`:
`Z(2N) = Re Tr X^{2N}` and `G n = Re ⟨v, Xⁿ v⟩`. -/
theorem OpenMomentBounds.of_hermitian {D : Type*} [Fintype D] [DecidableEq D] [Nonempty D]
    {X : Matrix D D ℂ} (hX : X.IsHermitian) (v : D → ℂ) {Z G : ℕ → ℝ}
    (hZ : ∀ N, 1 ≤ N → Z (2 * N) = (Matrix.trace (X ^ (2 * N))).re)
    (hG : ∀ n, 2 ≤ n → G n = (star v ⬝ᵥ (X ^ n).mulVec v).re) : OpenMomentBounds Z G := by
  -- In the eigenbasis `U` of `X`, `Xⁿ = U diag(λⁿ) U*`; both traces are then sums over `i`,
  -- with boundary weights `wᵢ = |(U* v)ᵢ|²`.
  have hpow : ∀ n : ℕ, X ^ n = (hX.eigenvectorUnitary : Matrix D D ℂ) *
      Matrix.diagonal (fun i ↦ ((hX.eigenvalues i ^ n : ℝ) : ℂ)) *
        star (hX.eigenvectorUnitary : Matrix D D ℂ) := fun n ↦ by
    have h : (RCLike.ofReal ∘ hX.eigenvalues) ^ n =
        fun i ↦ ((hX.eigenvalues i ^ n : ℝ) : ℂ) := by
      funext i; simp
    conv_lhs => rw [hX.spectral_theorem, ← map_pow, Matrix.diagonal_pow,
      Unitary.conjStarAlgAut_apply, h]
  refine OpenMomentBounds.of_eigen hX.eigenvalues
    (fun i ↦ Complex.normSq ((star (hX.eigenvectorUnitary : Matrix D D ℂ)).mulVec v i))
    (fun i ↦ Complex.normSq_nonneg _) (fun N hN ↦ ?_) (fun n hn ↦ ?_)
  · rw [hZ N hN, hpow, Matrix.trace_mul_cycle, Unitary.coe_star_mul_self, one_mul,
      Matrix.trace_diagonal, Complex.re_sum]
    simp only [Complex.ofReal_re]
  · rw [hG n hn, hpow, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
      show Matrix.vecMul (star v) (hX.eigenvectorUnitary : Matrix D D ℂ) =
        star ((star (hX.eigenvectorUnitary : Matrix D D ℂ)).mulVec v) by
      rw [Matrix.star_mulVec, Matrix.star_eq_conjTranspose, Matrix.conjTranspose_conjTranspose]]
    simp only [dotProduct, Matrix.mulVec_diagonal, Complex.re_sum]
    refine sum_congr rfl fun i _ ↦ ?_
    rw [mul_left_comm, Pi.star_apply, Complex.star_def, ← Complex.normSq_eq_conj_mul_self,
      ← Complex.ofReal_mul, Complex.ofReal_re, mul_comm]

/-! ### Limits -/

/-- A compact box containing every witness for traces close to `Z 2`, `Z 4`, `G 2`: the bounds
`μ² ≤ Z 2`, `Z 4 ≤ μ² Z 2` and `c μ² ≤ G 2` confine `μ` and `c`, and `ρ, s ∈ [-1, 1]`. -/
private noncomputable def witnessBox (z₂ z₄ g₂ : ℝ) : Set (ℝ × ℝ × ℝ × ℝ) :=
  Set.Icc 0 (1 + 2 * z₂) ×ˢ Set.Icc 0 ((g₂ + 1) / (z₄ / (4 * z₂))) ×ˢ Set.Icc 0 1 ×ˢ
    Set.Icc (-1) 1

/-- The witness box is compact, as a product of closed intervals. -/
private theorem isCompact_witnessBox (z₂ z₄ g₂ : ℝ) : IsCompact (witnessBox z₂ z₄ g₂) :=
  isCompact_Icc.prod (isCompact_Icc.prod (isCompact_Icc.prod isCompact_Icc))

/-- A witness for traces `Z' 2 < 2 z₂`, `Z' 4 > z₄ / 2`, `G' 2 < g₂ + 1` lies in the box. -/
private theorem OpenMomentWitness.mem_witnessBox {Z' G' : ℕ → ℝ} {μ c ρ s z₂ z₄ g₂ : ℝ}
    (h : OpenMomentWitness Z' G' μ c ρ s) (hz₂ : 0 < z₂) (hz₄ : 0 < z₄) (h2 : Z' 2 < 2 * z₂)
    (h4 : z₄ / 2 < Z' 4) (hg : G' 2 < g₂ + 1) :
    (μ, c, ρ, s) ∈ witnessBox z₂ z₄ g₂ := by
  have hμ2 : μ ^ 2 ≤ Z' 2 := by simpa using h.pow_le_Z 1 le_rfl
  have hμ4 : Z' 4 ≤ μ ^ 2 * Z' 2 := by simpa using h.Z_two_mul_le 1 le_rfl
  have hc : c * μ ^ 2 ≤ G' 2 := by simpa using h.le_G 1 le_rfl
  have hs := h.s_sq
  have hμ0 := h.mu_nonneg
  have hc0 := h.c_nonneg
  -- `z₄ < 2 Z' 4 ≤ 2 μ² Z' 2 < 4 μ² z₂`, a lower bound on `μ²` that then bounds `c`.
  have hμa : z₄ / (4 * z₂) < μ ^ 2 := by
    rw [div_lt_iff₀ (by positivity)]
    nlinarith [mul_le_mul_of_nonneg_left h2.le (sq_nonneg μ)]
  refine ⟨⟨hμ0, by nlinarith [sq_nonneg (μ - 1)]⟩, ⟨hc0, ?_⟩, ⟨h.rho_nonneg, h.rho_le_one⟩,
    ⟨by nlinarith, by nlinarith⟩⟩
  rw [le_div_iff₀ (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hμa.le hc0]

/-- The bounds pass to limits: if the approximants `Zs M`, `Gs M` satisfy them and converge at
every length `n ≥ 2` to `Z`, `G` with `Z 2, Z 4 > 0`, then `Z`, `G` satisfy them. -/
theorem OpenMomentBounds.of_tendsto {Zs Gs : ℕ → ℕ → ℝ} {Z G : ℕ → ℝ}
    (h : ∀ᶠ M in atTop, OpenMomentBounds (Zs M) (Gs M))
    (hZ : ∀ n, 2 ≤ n → Tendsto (fun M ↦ Zs M n) atTop (𝓝 (Z n)))
    (hG : ∀ n, 2 ≤ n → Tendsto (fun M ↦ Gs M n) atTop (𝓝 (G n)))
    (hZ2 : 0 < Z 2) (hZ4 : 0 < Z 4) : OpenMomentBounds Z G := by
  -- Choose a witness `x M = (μ, c, ρ, s)` for every large `M` (and junk below `M₀`).
  obtain ⟨M₀, hM₀⟩ := eventually_atTop.1 h
  have hx : ∀ M, ∃ x : ℝ × ℝ × ℝ × ℝ, M₀ ≤ M →
      OpenMomentWitness (Zs M) (Gs M) x.1 x.2.1 x.2.2.1 x.2.2.2 := fun M ↦ by
    by_cases hM : M₀ ≤ M
    · obtain ⟨μ, c, ρ, s, hw⟩ := hM₀ M hM
      exact ⟨(μ, c, ρ, s), fun _ ↦ hw⟩
    · exact ⟨0, fun h ↦ absurd h hM⟩
  choose x hx using hx
  have hw : ∀ᶠ M in atTop,
      OpenMomentWitness (Zs M) (Gs M) (x M).1 (x M).2.1 (x M).2.2.1 (x M).2.2.2 :=
    eventually_atTop.2 ⟨M₀, hx⟩
  -- The witnesses eventually lie in one compact box, since `Zs M 2`, `Zs M 4`, `Gs M 2` are
  -- eventually close to their limits, and `Z 2`, `Z 4` are positive.
  have hbox : ∀ᶠ M in atTop, x M ∈ witnessBox (Z 2) (Z 4) (G 2) := by
    filter_upwards [hw, (hZ 2 le_rfl).eventually (eventually_lt_nhds (by linarith : Z 2 < 2 * Z 2)),
      (hZ 4 (by norm_num)).eventually (eventually_gt_nhds (by linarith : Z 4 / 2 < Z 4)),
      (hG 2 le_rfl).eventually (eventually_lt_nhds (by linarith : G 2 < G 2 + 1))]
      with M hM h2 h4 hg
    exact hM.mem_witnessBox hZ2 hZ4 h2 h4 hg
  obtain ⟨y, -, φ, hφ, hlim⟩ := (isCompact_witnessBox _ _ _).tendsto_subseq' hbox.frequently
  -- Along the subsequence `φ`, every field is a closed inequality between convergent sequences.
  have hw' := hφ.tendsto_atTop.eventually hw
  have tμ : Tendsto (fun k ↦ (x (φ k)).1) atTop (𝓝 y.1) := (continuous_fst.tendsto y).comp hlim
  have tc : Tendsto (fun k ↦ (x (φ k)).2.1) atTop (𝓝 y.2.1) :=
    ((continuous_fst.comp continuous_snd).tendsto y).comp hlim
  have tρ : Tendsto (fun k ↦ (x (φ k)).2.2.1) atTop (𝓝 y.2.2.1) :=
    ((continuous_fst.comp (continuous_snd.comp continuous_snd)).tendsto y).comp hlim
  have ts : Tendsto (fun k ↦ (x (φ k)).2.2.2) atTop (𝓝 y.2.2.2) :=
    ((continuous_snd.comp (continuous_snd.comp continuous_snd)).tendsto y).comp hlim
  have tZ : ∀ n, 2 ≤ n → Tendsto (fun k ↦ Zs (φ k) n) atTop (𝓝 (Z n)) := fun n hn ↦
    (hZ n hn).comp hφ.tendsto_atTop
  have tG : ∀ n, 2 ≤ n → Tendsto (fun k ↦ Gs (φ k) n) atTop (𝓝 (G n)) := fun n hn ↦
    (hG n hn).comp hφ.tendsto_atTop
  exact ⟨y.1, y.2.1, y.2.2.1, y.2.2.2,
    { mu_nonneg := ge_of_tendsto tμ (hw'.mono fun _ h ↦ h.mu_nonneg)
      c_nonneg := ge_of_tendsto tc (hw'.mono fun _ h ↦ h.c_nonneg)
      rho_nonneg := ge_of_tendsto tρ (hw'.mono fun _ h ↦ h.rho_nonneg)
      rho_le_one := le_of_tendsto tρ (hw'.mono fun _ h ↦ h.rho_le_one)
      s_sq := tendsto_nhds_unique (ts.pow 2)
        (tendsto_const_nhds.congr' (hw'.mono fun _ h ↦ h.s_sq.symm))
      pow_le_Z := fun N hN ↦ le_of_tendsto_of_tendsto (tμ.pow _) (tZ _ (by omega))
        (hw'.mono fun _ h ↦ h.pow_le_Z N hN)
      Z_anti := fun K N hK hKN ↦ le_of_tendsto_of_tendsto ((tμ.pow _).mul (tZ _ (by omega)))
        ((tμ.pow _).mul (tZ _ (by omega))) (hw'.mono fun _ h ↦ h.Z_anti K N hK hKN)
      Z_two_mul_le := fun N hN ↦ le_of_tendsto_of_tendsto (tZ _ (by omega))
        ((tμ.pow _).mul (tZ _ (by omega))) (hw'.mono fun _ h ↦ h.Z_two_mul_le N hN)
      Z_two_mul_le_sq := fun N hN ↦ le_of_tendsto_of_tendsto (tZ _ (by omega))
        ((tμ.pow _).add (((tZ _ (by omega)).sub (tμ.pow _)).pow 2))
        (hw'.mono fun _ h ↦ h.Z_two_mul_le_sq N hN)
      tail := fun N hN ↦ le_of_tendsto_of_tendsto ((tρ.mul tμ).pow _)
        ((tZ _ (by omega)).sub (tμ.pow _)) (hw'.mono fun _ h ↦ h.tail N hN)
      le_G := fun N hN ↦ le_of_tendsto_of_tendsto (tc.mul (tμ.pow _)) (tG _ (by omega))
        (hw'.mono fun _ h ↦ h.le_G N hN)
      G_decay := fun K N hK hKN ↦ le_of_tendsto_of_tendsto
        ((tμ.pow _).mul ((tG _ (by omega)).sub (tc.mul (tμ.pow _))))
        (((tρ.pow _).mul (tμ.pow _)).mul ((tG _ (by omega)).sub (tc.mul (tμ.pow _))))
        (hw'.mono fun _ h ↦ h.G_decay K N hK hKN)
      G_odd := fun K N hK hKN ↦ le_of_tendsto_of_tendsto
        ((tG _ (by omega)).sub ((ts.mul tc).mul (tμ.pow _))).abs
        ((tμ.pow _).mul ((tG _ (by omega)).sub (tc.mul (tμ.pow _))))
        (hw'.mono fun _ h ↦ h.G_odd K N hK hKN) }⟩

/-! ### The leading modulus and weight in terms of traces -/

/-- The leading modulus `μ = inf_N Z(2N)^{1/(2N)}` of the transfer operator, from ring traces. -/
noncomputable def leadingRoot (Z : ℕ → ℝ) : ℝ :=
  ⨅ N : ℕ, Z (2 * (N + 1)) ^ ((1 : ℝ) / (2 * (N + 1) : ℕ))

/-- The leading boundary weight `c = inf_N G(2N) / μ^{2N}`. -/
noncomputable def leadingOverlap (Z G : ℕ → ℝ) : ℝ :=
  ⨅ N : ℕ, G (2 * (N + 1)) / leadingRoot Z ^ (2 * (N + 1))

/-- The boundary defect `e_k = 1 - c μ^k / G k`: the fraction of the open-chain trace at length
`k` not carried by the leading eigenvector (the note, eq. (A13)). -/
noncomputable def boundaryDefect (Z G : ℕ → ℝ) (k : ℕ) : ℝ :=
  1 - leadingOverlap Z G * leadingRoot Z ^ k / G k

namespace OpenMomentWitness

variable {Z G : ℕ → ℝ} {μ c ρ s : ℝ}

/-- `μ > 0` once `Z 4 > 0`. -/
theorem mu_pos (h : OpenMomentWitness Z G μ c ρ s) (hZ4 : 0 < Z 4) : 0 < μ := by
  -- `Z 4 ≤ μ² Z 2` would force `Z 4 ≤ 0` at `μ = 0`.
  refine h.mu_nonneg.lt_of_ne fun h0 ↦ ?_
  have h4 := h.Z_two_mul_le 1 le_rfl
  rw [← h0] at h4
  norm_num at h4
  linarith

/-- `1 / (2 (N + 1)) → 0`, the exponent of the roots in `leadingRoot`. -/
private theorem tendsto_one_div_two_mul_succ :
    Tendsto (fun N : ℕ ↦ (1 : ℝ) / (2 * (N + 1) : ℕ)) atTop (𝓝 0) :=
  tendsto_const_nhds.div_atTop <| tendsto_natCast_atTop_atTop.comp <|
    tendsto_atTop_mono (fun N ↦ show N ≤ 2 * (N + 1) by omega) tendsto_id

/-- The witness modulus is the leading root of the ring traces. -/
theorem leadingRoot_eq (h : OpenMomentWitness Z G μ c ρ s) (hZ4 : 0 < Z 4) :
    leadingRoot Z = μ := by
  have hμ := h.mu_pos hZ4
  have hZ2 : μ ^ 2 ≤ Z 2 := by simpa using h.pow_le_Z 1 le_rfl
  -- Each root `Z(2N)^{1/(2N)}` lies between `μ` and `μ K^{1/(2N)}`, with `K = Z 2 / μ²`; the
  -- upper bound tends to `μ`.
  set K := Z 2 / μ ^ 2
  have hK : 0 < K := div_pos ((by positivity : (0 : ℝ) < μ ^ 2).trans_le hZ2) (by positivity)
  have hlo : ∀ N : ℕ, μ ≤ Z (2 * (N + 1)) ^ ((1 : ℝ) / (2 * (N + 1) : ℕ)) := fun N ↦ by
    calc μ = (μ ^ (2 * (N + 1))) ^ ((1 : ℝ) / (2 * (N + 1) : ℕ)) := by
          rw [one_div, Real.pow_rpow_inv_natCast hμ.le (by omega)]
      _ ≤ _ := Real.rpow_le_rpow (by positivity) (h.pow_le_Z _ (by omega)) (by positivity)
  have hhi : ∀ N : ℕ, Z (2 * (N + 1)) ^ ((1 : ℝ) / (2 * (N + 1) : ℕ)) ≤
      μ * K ^ ((1 : ℝ) / (2 * (N + 1) : ℕ)) := fun N ↦ by
    have hle : Z (2 * (N + 1)) ≤ μ ^ (2 * (N + 1)) * K := by
      have h1 := h.Z_anti 1 (N + 1) le_rfl (by omega)
      rw [mul_one] at h1
      rw [mul_div_assoc', le_div_iff₀ (by positivity)]
      linear_combination h1
    calc _ ≤ (μ ^ (2 * (N + 1)) * K) ^ ((1 : ℝ) / (2 * (N + 1) : ℕ)) :=
          Real.rpow_le_rpow ((pow_nonneg hμ.le _).trans (h.pow_le_Z _ (by omega))) hle
            (by positivity)
      _ = _ := by
        rw [Real.mul_rpow (by positivity) hK.le, one_div,
          Real.pow_rpow_inv_natCast hμ.le (by omega)]
  have hlim : Tendsto (fun N : ℕ ↦ μ * K ^ ((1 : ℝ) / (2 * (N + 1) : ℕ))) atTop (𝓝 μ) := by
    have := (tendsto_const_nhds.rpow tendsto_one_div_two_mul_succ (Or.inl hK.ne')).const_mul μ
    rwa [Real.rpow_zero, mul_one] at this
  exact le_antisymm
    (ge_of_tendsto' hlim fun N ↦ (ciInf_le ⟨μ, Set.forall_mem_range.2 hlo⟩ N).trans (hhi N))
    (le_ciInf hlo)

/-- With `ρ < 1`, the witness weight is the leading overlap of the open-chain traces. -/
theorem leadingOverlap_eq (h : OpenMomentWitness Z G μ c ρ s) (hZ4 : 0 < Z 4) (hρ : ρ < 1) :
    leadingOverlap Z G = c := by
  have hμ := h.mu_pos hZ4
  rw [leadingOverlap, h.leadingRoot_eq hZ4]
  -- Each ratio `B(2N) = G(2N)/μ^{2N}` lies between `c` and `c + ρ^{2N-2} (B 2 - c)`, and the
  -- upper bound tends to `c` since `ρ < 1`.
  have hlo : ∀ N : ℕ, c ≤ G (2 * (N + 1)) / μ ^ (2 * (N + 1)) := fun N ↦
    (le_div_iff₀ (by positivity)).2 (h.le_G _ (by omega))
  have hhi : ∀ N : ℕ, G (2 * (N + 1)) / μ ^ (2 * (N + 1)) ≤
      c + (ρ ^ 2) ^ N * (G 2 / μ ^ 2 - c) := fun N ↦ by
    have h1 := h.G_decay 1 (N + 1) le_rfl (by omega)
    rw [show 2 * (N + 1) - 2 * 1 = 2 * N by omega, mul_one, pow_mul] at h1
    have hm : (0 : ℝ) < μ ^ (2 * (N + 1)) := by positivity
    have hm2 : (0 : ℝ) < μ ^ 2 := by positivity
    rw [div_le_iff₀ hm, show (c + (ρ ^ 2) ^ N * (G 2 / μ ^ 2 - c)) * μ ^ (2 * (N + 1)) =
      c * μ ^ (2 * (N + 1)) + (ρ ^ 2) ^ N * μ ^ (2 * (N + 1)) * (G 2 - c * μ ^ 2) / μ ^ 2 by
        field_simp, ← sub_le_iff_le_add', le_div_iff₀ hm2]
    linear_combination h1
  have hlim : Tendsto (fun N : ℕ ↦ c + (ρ ^ 2) ^ N * (G 2 / μ ^ 2 - c)) atTop (𝓝 c) := by
    have := ((tendsto_pow_atTop_nhds_zero_of_lt_one (sq_nonneg ρ)
      (pow_lt_one₀ h.rho_nonneg hρ two_ne_zero)).mul_const (G 2 / μ ^ 2 - c)).const_add c
    rwa [zero_mul, add_zero] at this
  exact le_antisymm
    (ge_of_tendsto' hlim fun N ↦ (ciInf_le ⟨c, Set.forall_mem_range.2 hlo⟩ N).trans (hhi N))
    (le_ciInf hlo)

end OpenMomentWitness

end CertifiedGaps
