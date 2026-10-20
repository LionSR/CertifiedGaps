/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.MomentBounds

/-!
# Boundary witnesses with a leading-weight bound, and their limits

The initialization of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026)
needs, besides the boundary moment inequalities at `b₀`, three closed facts about the same
witness `(μ, c, ρ, s)`: a lower bound `ℓ ≤ μ`, an upper bound `ρ ≤ ρmax`, and the leading-weight
bound `c f(s μ) ≤ ∑_{n ∈ S} aₙ I(n)` for a polynomial `f(x) = ∑_{n ∈ S} aₙ xⁿ` and the averaged
traces `I(n)` (eq. (appth:polynomial)). `LeadingWitness` bundles them. All the conditions are
closed, so they pass to limits of approximants along the same convergent subsequence of
witnesses (`LeadingWitness.of_tendsto`), exactly as in `OpenMomentBounds.of_tendsto`.

The leading eigenvalue `s μ` of the approximants carries a sign `s`. In the limit the sign is
`+1` as soon as `c > 0`, `ρ < 1` and the odd open-chain traces are positive
(`OpenMomentWitness.eq_zero_or_sign_eq_one`): a negative leading eigenvalue would make the odd
traces `≈ -c μ²ᴺ⁺¹` negative, as in the paper's proof that `λⱼ > 0` (Proposition
per:concentration). The leading-weight bound then holds at `μ` itself
(`LeadingWitness.mul_poly_le`).

## Main results

* `LeadingWitness`: the witness with the three closed facts.
* `LeadingWitness.of_tendsto`: they pass to limits.
* `OpenMomentWitness.eq_zero_or_sign_eq_one`, `LeadingWitness.mul_poly_le`.
-/

open Filter Finset Topology

namespace CertifiedGaps

/-- **A boundary witness with a window and a leading-weight bound**: a witness `(μ, c, ρ, s)` of
the boundary moment inequalities for `Z`, `G` with `ℓ ≤ μ`, `ρ ≤ ρmax`, and
`c ∑_{n ∈ S} aₙ (s μ)ⁿ ≤ ∑_{n ∈ S} aₙ I(n)` for the traces `I`. -/
structure LeadingWitness (Z G I : ℕ → ℝ) (S : Finset ℕ) (a : ℕ → ℝ) (ℓ ρmax : ℝ)
    (μ c ρ s : ℝ) : Prop extends OpenMomentWitness Z G μ c ρ s where
  ell_le : ℓ ≤ μ
  rho_le : ρ ≤ ρmax
  poly : c * ∑ n ∈ S, a n * (s * μ) ^ n ≤ ∑ n ∈ S, a n * I n

/-- A compact box containing every witness for traces close to `Z 2`, `Z 4`, `G 2` (as in
`Open/MomentBounds.lean`). -/
private noncomputable def leadingBox (z₂ z₄ g₂ : ℝ) : Set (ℝ × ℝ × ℝ × ℝ) :=
  Set.Icc 0 (1 + 2 * z₂) ×ˢ Set.Icc 0 ((g₂ + 1) / (z₄ / (4 * z₂))) ×ˢ Set.Icc 0 1 ×ˢ
    Set.Icc (-1) 1

private theorem isCompact_leadingBox (z₂ z₄ g₂ : ℝ) : IsCompact (leadingBox z₂ z₄ g₂) :=
  isCompact_Icc.prod (isCompact_Icc.prod (isCompact_Icc.prod isCompact_Icc))

/-- A witness for traces `Z' 2 < 2 z₂`, `Z' 4 > z₄ / 2`, `G' 2 < g₂ + 1` lies in the box. -/
private theorem mem_leadingBox {Z' G' : ℕ → ℝ} {μ c ρ s z₂ z₄ g₂ : ℝ}
    (h : OpenMomentWitness Z' G' μ c ρ s) (hz₂ : 0 < z₂) (hz₄ : 0 < z₄) (h2 : Z' 2 < 2 * z₂)
    (h4 : z₄ / 2 < Z' 4) (hg : G' 2 < g₂ + 1) :
    (μ, c, ρ, s) ∈ leadingBox z₂ z₄ g₂ := by
  have hμ2 : μ ^ 2 ≤ Z' 2 := by simpa using h.pow_le_Z 1 le_rfl
  have hμ4 : Z' 4 ≤ μ ^ 2 * Z' 2 := by simpa using h.Z_two_mul_le 1 le_rfl
  have hc : c * μ ^ 2 ≤ G' 2 := by simpa using h.le_G 1 le_rfl
  have hs := h.s_sq
  have hμ0 := h.mu_nonneg
  have hc0 := h.c_nonneg
  have hμa : z₄ / (4 * z₂) < μ ^ 2 := by
    rw [div_lt_iff₀ (by positivity)]
    nlinarith [mul_le_mul_of_nonneg_left h2.le (sq_nonneg μ)]
  refine ⟨⟨hμ0, by nlinarith [sq_nonneg (μ - 1)]⟩, ⟨hc0, ?_⟩, ⟨h.rho_nonneg, h.rho_le_one⟩,
    ⟨by nlinarith, by nlinarith⟩⟩
  rw [le_div_iff₀ (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hμa.le hc0]

/-- **The bundled facts pass to limits.** If the approximants `Zs M`, `Gs M`, `Is M` have
witnesses for every large `M` and converge to `Z`, `G` (at every length `n ≥ 2`) and `I` (on
`S`) with `Z 2, Z 4 > 0`, then `Z`, `G`, `I` have a witness. -/
theorem LeadingWitness.of_tendsto {Zs Gs Is : ℕ → ℕ → ℝ} {Z G I : ℕ → ℝ} {S : Finset ℕ}
    {a : ℕ → ℝ} {ℓ ρmax : ℝ}
    (h : ∀ᶠ M in atTop, ∃ μ c ρ s, LeadingWitness (Zs M) (Gs M) (Is M) S a ℓ ρmax μ c ρ s)
    (hZ : ∀ n, 2 ≤ n → Tendsto (fun M ↦ Zs M n) atTop (𝓝 (Z n)))
    (hG : ∀ n, 2 ≤ n → Tendsto (fun M ↦ Gs M n) atTop (𝓝 (G n)))
    (hI : ∀ n ∈ S, Tendsto (fun M ↦ Is M n) atTop (𝓝 (I n)))
    (hZ2 : 0 < Z 2) (hZ4 : 0 < Z 4) : ∃ μ c ρ s, LeadingWitness Z G I S a ℓ ρmax μ c ρ s := by
  obtain ⟨M₀, hM₀⟩ := eventually_atTop.1 h
  have hx : ∀ M, ∃ x : ℝ × ℝ × ℝ × ℝ, M₀ ≤ M →
      LeadingWitness (Zs M) (Gs M) (Is M) S a ℓ ρmax x.1 x.2.1 x.2.2.1 x.2.2.2 := fun M ↦ by
    by_cases hM : M₀ ≤ M
    · obtain ⟨μ, c, ρ, s, hw⟩ := hM₀ M hM
      exact ⟨(μ, c, ρ, s), fun _ ↦ hw⟩
    · exact ⟨0, fun h ↦ absurd h hM⟩
  choose x hx using hx
  have hw : ∀ᶠ M in atTop,
      LeadingWitness (Zs M) (Gs M) (Is M) S a ℓ ρmax (x M).1 (x M).2.1 (x M).2.2.1 (x M).2.2.2 :=
    eventually_atTop.2 ⟨M₀, hx⟩
  have hbox : ∀ᶠ M in atTop, x M ∈ leadingBox (Z 2) (Z 4) (G 2) := by
    filter_upwards [hw, (hZ 2 le_rfl).eventually (eventually_lt_nhds (by linarith : Z 2 < 2 * Z 2)),
      (hZ 4 (by norm_num)).eventually (eventually_gt_nhds (by linarith : Z 4 / 2 < Z 4)),
      (hG 2 le_rfl).eventually (eventually_lt_nhds (by linarith : G 2 < G 2 + 1))]
      with M hM h2 h4 hg
    exact mem_leadingBox hM.toOpenMomentWitness hZ2 hZ4 h2 h4 hg
  obtain ⟨y, -, φ, hφ, hlim⟩ := (isCompact_leadingBox _ _ _).tendsto_subseq' hbox.frequently
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
  have tI : Tendsto (fun k ↦ ∑ n ∈ S, a n * Is (φ k) n) atTop (𝓝 (∑ n ∈ S, a n * I n)) :=
    tendsto_finsetSum _ fun n hn ↦ ((hI n hn).comp hφ.tendsto_atTop).const_mul (a n)
  have tP : Tendsto (fun k ↦ (x (φ k)).2.1 * ∑ n ∈ S, a n * ((x (φ k)).2.2.2 * (x (φ k)).1) ^ n)
      atTop (𝓝 (y.2.1 * ∑ n ∈ S, a n * (y.2.2.2 * y.1) ^ n)) :=
    tc.mul (tendsto_finsetSum _ fun n _ ↦ ((ts.mul tμ).pow n).const_mul (a n))
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
        (hw'.mono fun _ h ↦ h.G_odd K N hK hKN)
      ell_le := ge_of_tendsto tμ (hw'.mono fun _ h ↦ h.ell_le)
      rho_le := le_of_tendsto tρ (hw'.mono fun _ h ↦ h.rho_le)
      poly := le_of_tendsto_of_tendsto tP tI (hw'.mono fun _ h ↦ h.poly) }⟩

/-- **The sign of the leading eigenvalue.** If `ρ < 1`, `μ > 0` and the odd open-chain traces
are positive, then `c = 0` or `s = 1`: with `s = -1` the odd traces would be
`G(2N+1) ≈ -c μ²ᴺ⁺¹ < 0`. -/
theorem OpenMomentWitness.eq_zero_or_sign_eq_one {Z G : ℕ → ℝ} {μ c ρ s : ℝ}
    (h : OpenMomentWitness Z G μ c ρ s) (hρ : ρ < 1) (hμ : 0 < μ)
    (hG : ∀ N, 1 ≤ N → 0 < G (2 * N + 1)) : c = 0 ∨ s = 1 := by
  have hs : s = 1 ∨ s = -1 := by
    have e : (s - 1) * (s + 1) = 0 := by linear_combination h.s_sq
    rcases mul_eq_zero.1 e with h1 | h1
    · left; linarith
    · right; linarith
  rcases hs with hs | hs
  · exact Or.inr hs
  left
  by_contra hc
  have hc0 : 0 < c := lt_of_le_of_ne h.c_nonneg (Ne.symm hc)
  -- The remainder `D = G 2 / μ² - c ≥ 0` controls `c < ρ^{2N-2} D` for every `N ≥ 1`.
  set D := G 2 / μ ^ 2 - c
  have hD : 0 ≤ D := by
    have := h.le_G 1 le_rfl
    simp only [mul_one] at this
    rw [sub_nonneg, le_div_iff₀ (by positivity)]
    simpa using this
  have hbound : ∀ N, 1 ≤ N → c ≤ (ρ ^ 2) ^ (N - 1) * D := fun N hN ↦ by
    have hodd := h.G_odd N N hN le_rfl
    rw [hs, show 2 * N + 1 - 2 * N = 1 by omega, pow_one] at hodd
    have hpos := hG N hN
    have h1 : c * μ ^ (2 * N + 1) ≤ μ * (G (2 * N) - c * μ ^ (2 * N)) := by
      have := (abs_le.1 hodd).2
      nlinarith
    have hdec := h.G_decay 1 N le_rfl hN
    rw [show 2 * N - 2 * 1 = 2 * (N - 1) by omega, mul_one] at hdec
    have hμN : (0 : ℝ) < μ ^ (2 * N) := by positivity
    have h2 : c * μ ^ (2 * N) ≤ G (2 * N) - c * μ ^ (2 * N) := by
      have e : c * μ ^ (2 * N + 1) = μ * (c * μ ^ (2 * N)) := by ring
      rw [e] at h1
      exact le_of_mul_le_mul_left h1 hμ
    have h3 : μ ^ 2 * (c * μ ^ (2 * N)) ≤ ρ ^ (2 * (N - 1)) * μ ^ (2 * N) * (G 2 - c * μ ^ 2) :=
      (mul_le_mul_of_nonneg_left h2 (by positivity)).trans hdec
    have hG2 : G 2 - c * μ ^ 2 = μ ^ 2 * D := by
      simp only [D]; field_simp
    rw [hG2, pow_mul ρ 2 (N - 1)] at h3
    have : μ ^ 2 * μ ^ (2 * N) * c ≤ μ ^ 2 * μ ^ (2 * N) * ((ρ ^ 2) ^ (N - 1) * D) :=
      calc μ ^ 2 * μ ^ (2 * N) * c = μ ^ 2 * (c * μ ^ (2 * N)) := by ring
        _ ≤ _ := h3
        _ = _ := by ring
    exact le_of_mul_le_mul_left this (by positivity)
  -- As `N → ∞` the bound tends to `0`.
  have hρ2 : ρ ^ 2 < 1 := by nlinarith [h.rho_nonneg]
  obtain ⟨k, hk⟩ := exists_pow_lt_of_lt_one (div_pos hc0 (by linarith : (0 : ℝ) < D + 1)) hρ2
  have := hbound (k + 1) (by omega)
  rw [Nat.add_sub_cancel] at this
  have hk' : (ρ ^ 2) ^ k * (D + 1) < c := by
    rwa [lt_div_iff₀ (by linarith)] at hk
  nlinarith [pow_nonneg (sq_nonneg ρ) k]

/-- **The leading-weight bound at `μ`.** For a witness with `ρmax < 1`, `ℓ > 0` and positive odd
traces, `∑ aₙ I(n) ≤ B` with `B ≥ 0` gives `c ∑ aₙ μⁿ ≤ B`. -/
theorem LeadingWitness.mul_poly_le {Z G I : ℕ → ℝ} {S : Finset ℕ} {a : ℕ → ℝ} {ℓ ρmax : ℝ}
    {μ c ρ s : ℝ} (h : LeadingWitness Z G I S a ℓ ρmax μ c ρ s) (hρmax : ρmax < 1) (hℓ : 0 < ℓ)
    (hG : ∀ N, 1 ≤ N → 0 < G (2 * N + 1)) {B : ℝ} (hB : 0 ≤ B)
    (hI : ∑ n ∈ S, a n * I n ≤ B) : c * ∑ n ∈ S, a n * μ ^ n ≤ B := by
  rcases h.eq_zero_or_sign_eq_one (h.rho_le.trans_lt hρmax) (hℓ.trans_le h.ell_le) hG with
    hc | hs
  · rw [hc, zero_mul]; exact hB
  · have := h.poly
    rw [hs, one_mul] at this
    exact this.trans hI

end CertifiedGaps
