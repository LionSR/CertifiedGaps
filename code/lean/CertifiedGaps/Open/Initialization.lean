/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.MomentBounds
import CertifiedGaps.Seed.Lift
import CertifiedGaps.Transfer.RatioMonotone

/-!
# Starting the boundary recurrence from short-chain data

The boundary recurrence (`move_one`) needs, at its first temperature `β₀`, a bound on the
boundary defect `e_m = 1 - c μᵐ / G(m)` and lower bounds on the open-chain thermal purities
`A(m) = G_m(2β) / G_m(β)²`. Lemma in:seedbounds of OpenAI's *A boundary-field gap for the
spin-one Heisenberg chain* (2026) derives both from four kinds of finite data at `β₀`:

* a window `ℓ ≤ μ ≤ U` for the leading modulus and a residual ratio `ρ ≤ ρ̄`;
* a short boundary moment `G(2K) ≤ Γ` (the paper's `Y₆`, `K = 3`);
* a bound `c P(μ)² ≤ C` on the leading weight through a polynomial `P`;
* a ground weight `g ≤ G(2N)`, from a variational energy.

The residual part of `G(2N)` is at most `(ρ̄ U)^{2N-2K} Γ` (`residual_le_of_bounds`); dividing by
the ground weight bounds the defect (`one_sub_mul_le_of_ground`). If `x^{2N} / P(x)²` increases on
`[ℓ, U]` (`PolynomialFilter.RatioCheck`), the leading part is at most `C U^{2N} / P(U)²`
(`mul_pow_le_of_ratioCheck`), so `G(2N) ≤ g (1 + R)` for an explicit `R`. One chosen ground
Boltzmann factor then carries all but a fraction `R` of the open-chain trace; at `2^J β₀` the
fraction becomes `R^{2^J}`, which bounds the thermal purity (`thermalPurity_ge_of_ground`).

The data are naturally stated for energy-shifted traces `tⁿ G(n)`; `OpenMomentWitness.scale`
moves a witness between the two normalizations.

## Main results

* `OpenMomentWitness.scale`: the witness inequalities are invariant under `G n ↦ tⁿ G n`.
* `OpenMomentWitness.residual_le_of_bounds`: `G(2N) - c μ^{2N} ≤ (ρ̄ U)^{2N-2K} Γ`.
* `OpenMomentWitness.one_sub_mul_le_of_ground`: the first inequality of Lemma in:seedbounds,
  `(1 - ε) G(2N) ≤ c μ^{2N}`.
* `mul_pow_le_of_ratioCheck`: the leading part `c μ^{2N} ≤ C U^{2N} / P(U)²`.
* `OpenMomentWitness.le_add_of_bounds`: `G(2N) ≤ C U^{2N} / P(U)² + (ρ̄ U)^{2N-2K} Γ`.
* `sum_pow_sq_le_of_ground`, `thermalPurity_ge_of_ground`: the second inequality of Lemma
  in:seedbounds, `A(m, J) ≥ (1 + R^{2^J})⁻²`.
* `OpenMomentWitness.boundaryDefect_le`: the boundary defect from `(1 - E) G(k) ≤ c μᵏ`.
-/

open Finset Real

namespace CertifiedGaps

open PolynomialFilter

namespace OpenMomentWitness

variable {Z G : ℕ → ℝ} {μ c ρ s : ℝ}

/-- **Energy shifts.** The witness inequalities are homogeneous: a witness for `(Z, G)` with
modulus `μ` is one for `(tᴸ Z(L), tⁿ G(n))` with modulus `t μ`, for every `t > 0`. A shift of the
energy by `κ` per site at inverse temperature `β` is `t = e^{-βκ}`. -/
theorem scale (h : OpenMomentWitness Z G μ c ρ s) {t : ℝ} (ht : 0 < t) :
    OpenMomentWitness (fun L ↦ t ^ L * Z L) (fun n ↦ t ^ n * G n) (t * μ) c ρ s where
  mu_nonneg := mul_nonneg ht.le h.mu_nonneg
  c_nonneg := h.c_nonneg
  rho_nonneg := h.rho_nonneg
  rho_le_one := h.rho_le_one
  s_sq := h.s_sq
  pow_le_Z N hN := by
    rw [mul_pow]; exact mul_le_mul_of_nonneg_left (h.pow_le_Z N hN) (by positivity)
  Z_anti K N hK hKN := by
    have h1 := mul_le_mul_of_nonneg_left (h.Z_anti K N hK hKN)
      (by positivity : (0 : ℝ) ≤ t ^ (2 * K) * t ^ (2 * N))
    calc (t * μ) ^ (2 * K) * (t ^ (2 * N) * Z (2 * N))
        = t ^ (2 * K) * t ^ (2 * N) * (μ ^ (2 * K) * Z (2 * N)) := by ring
      _ ≤ t ^ (2 * K) * t ^ (2 * N) * (μ ^ (2 * N) * Z (2 * K)) := h1
      _ = (t * μ) ^ (2 * N) * (t ^ (2 * K) * Z (2 * K)) := by ring
  Z_two_mul_le N hN := by
    have h1 := mul_le_mul_of_nonneg_left (h.Z_two_mul_le N hN)
      (by positivity : (0 : ℝ) ≤ t ^ (2 * N) * t ^ (2 * N))
    calc t ^ (2 * (2 * N)) * Z (2 * (2 * N))
        = t ^ (2 * N) * t ^ (2 * N) * Z (2 * (2 * N)) := by ring
      _ ≤ t ^ (2 * N) * t ^ (2 * N) * (μ ^ (2 * N) * Z (2 * N)) := h1
      _ = (t * μ) ^ (2 * N) * (t ^ (2 * N) * Z (2 * N)) := by ring
  Z_two_mul_le_sq N hN := by
    have h1 := mul_le_mul_of_nonneg_left (h.Z_two_mul_le_sq N hN)
      (by positivity : (0 : ℝ) ≤ t ^ (2 * N) * t ^ (2 * N))
    calc t ^ (2 * (2 * N)) * Z (2 * (2 * N))
        = t ^ (2 * N) * t ^ (2 * N) * Z (2 * (2 * N)) := by ring
      _ ≤ t ^ (2 * N) * t ^ (2 * N) *
          (μ ^ (2 * (2 * N)) + (Z (2 * N) - μ ^ (2 * N)) ^ 2) := h1
      _ = (t * μ) ^ (2 * (2 * N)) + (t ^ (2 * N) * Z (2 * N) - (t * μ) ^ (2 * N)) ^ 2 := by
        ring
  tail N hN := by
    have h1 := mul_le_mul_of_nonneg_left (h.tail N hN) (by positivity : (0 : ℝ) ≤ t ^ (2 * N))
    calc (ρ * (t * μ)) ^ (2 * N) = t ^ (2 * N) * (ρ * μ) ^ (2 * N) := by ring
      _ ≤ t ^ (2 * N) * (Z (2 * N) - μ ^ (2 * N)) := h1
      _ = t ^ (2 * N) * Z (2 * N) - (t * μ) ^ (2 * N) := by ring
  le_G N hN := by
    have h1 := mul_le_mul_of_nonneg_left (h.le_G N hN) (by positivity : (0 : ℝ) ≤ t ^ (2 * N))
    calc c * (t * μ) ^ (2 * N) = t ^ (2 * N) * (c * μ ^ (2 * N)) := by ring
      _ ≤ t ^ (2 * N) * G (2 * N) := h1
  G_decay K N hK hKN := by
    have h1 := mul_le_mul_of_nonneg_left (h.G_decay K N hK hKN)
      (by positivity : (0 : ℝ) ≤ t ^ (2 * K) * t ^ (2 * N))
    calc (t * μ) ^ (2 * K) * (t ^ (2 * N) * G (2 * N) - c * (t * μ) ^ (2 * N))
        = t ^ (2 * K) * t ^ (2 * N) * (μ ^ (2 * K) * (G (2 * N) - c * μ ^ (2 * N))) := by ring
      _ ≤ t ^ (2 * K) * t ^ (2 * N) *
          (ρ ^ (2 * N - 2 * K) * μ ^ (2 * N) * (G (2 * K) - c * μ ^ (2 * K))) := h1
      _ = ρ ^ (2 * N - 2 * K) * (t * μ) ^ (2 * N) *
          (t ^ (2 * K) * G (2 * K) - c * (t * μ) ^ (2 * K)) := by ring
  G_odd K N hK hKN := by
    have e : t ^ (2 * N + 1 - 2 * K) * t ^ (2 * K) = t ^ (2 * N + 1) := by
      rw [← pow_add, Nat.sub_add_cancel (by omega)]
    have h1 := mul_le_mul_of_nonneg_left (h.G_odd K N hK hKN)
      (by positivity : (0 : ℝ) ≤ t ^ (2 * N + 1))
    calc |t ^ (2 * N + 1) * G (2 * N + 1) - s * c * (t * μ) ^ (2 * N + 1)|
        = t ^ (2 * N + 1) * |G (2 * N + 1) - s * c * μ ^ (2 * N + 1)| := by
          rw [← abs_of_pos (pow_pos ht (2 * N + 1)), ← abs_mul, abs_of_pos (pow_pos ht _)]
          ring_nf
      _ ≤ t ^ (2 * N + 1) * (μ ^ (2 * N + 1 - 2 * K) * (G (2 * K) - c * μ ^ (2 * K))) := h1
      _ = t ^ (2 * N + 1 - 2 * K) * t ^ (2 * K) *
          (μ ^ (2 * N + 1 - 2 * K) * (G (2 * K) - c * μ ^ (2 * K))) := by rw [e]
      _ = (t * μ) ^ (2 * N + 1 - 2 * K) *
          (t ^ (2 * K) * G (2 * K) - c * (t * μ) ^ (2 * K)) := by ring

/-- The residual part of an even boundary moment: `G(2N) - c μ^{2N} ≤ (ρ μ)^{2N-2K} G(2K)` for
`1 ≤ K ≤ N`, the paper's eq. (be:residual) with the bound `D(2K) ≤ B(2K)`. -/
theorem residual_le (h : OpenMomentWitness Z G μ c ρ s) (hμ : 0 < μ) {K N : ℕ} (hK : 1 ≤ K)
    (hKN : K ≤ N) : G (2 * N) - c * μ ^ (2 * N) ≤ (ρ * μ) ^ (2 * N - 2 * K) * G (2 * K) := by
  have hsplit : μ ^ (2 * N) = μ ^ (2 * K) * μ ^ (2 * N - 2 * K) := by
    rw [← pow_add]; congr 1; omega
  have hd := h.G_decay K N hK hKN
  have hc : 0 ≤ c * μ ^ (2 * K) := mul_nonneg h.c_nonneg (pow_nonneg hμ.le _)
  have h1 : μ ^ (2 * K) * (G (2 * N) - c * μ ^ (2 * N)) ≤
      μ ^ (2 * K) * ((ρ * μ) ^ (2 * N - 2 * K) * G (2 * K)) := by
    refine hd.trans ?_
    have hρa : 0 ≤ ρ ^ (2 * N - 2 * K) * μ ^ (2 * N) :=
      mul_nonneg (pow_nonneg h.rho_nonneg _) (pow_nonneg hμ.le _)
    calc ρ ^ (2 * N - 2 * K) * μ ^ (2 * N) * (G (2 * K) - c * μ ^ (2 * K))
        ≤ ρ ^ (2 * N - 2 * K) * μ ^ (2 * N) * G (2 * K) :=
          mul_le_mul_of_nonneg_left (by linarith) hρa
      _ = μ ^ (2 * K) * ((ρ * μ) ^ (2 * N - 2 * K) * G (2 * K)) := by rw [hsplit]; ring
  exact le_of_mul_le_mul_left h1 (pow_pos hμ _)

/-- **The residual part from scalar bounds.** With `μ ≤ U`, `ρ ≤ ρ̄` and `G(2K) ≤ Γ`, the part of
`G(2N)` outside the leading line is at most `(ρ̄ U)^{2N-2K} Γ`. -/
theorem residual_le_of_bounds (h : OpenMomentWitness Z G μ c ρ s) (hμ : 0 < μ) {K N : ℕ}
    (hK : 1 ≤ K) (hKN : K ≤ N) {U ρ' Γ : ℝ} (hμU : μ ≤ U) (hρ : ρ ≤ ρ') (hΓ : G (2 * K) ≤ Γ) :
    G (2 * N) - c * μ ^ (2 * N) ≤ (ρ' * U) ^ (2 * N - 2 * K) * Γ := by
  have hG : 0 ≤ G (2 * K) :=
    (mul_nonneg h.c_nonneg (pow_nonneg hμ.le _)).trans (h.le_G K hK)
  have hρμ : (ρ * μ) ^ (2 * N - 2 * K) ≤ (ρ' * U) ^ (2 * N - 2 * K) :=
    pow_le_pow_left₀ (mul_nonneg h.rho_nonneg hμ.le)
      (mul_le_mul hρ hμU hμ.le (h.rho_nonneg.trans hρ)) _
  calc G (2 * N) - c * μ ^ (2 * N) ≤ (ρ * μ) ^ (2 * N - 2 * K) * G (2 * K) :=
        h.residual_le hμ hK hKN
    _ ≤ (ρ' * U) ^ (2 * N - 2 * K) * Γ :=
        mul_le_mul hρμ hΓ hG (pow_nonneg (mul_nonneg (h.rho_nonneg.trans hρ)
          (hμ.le.trans hμU)) _)

/-- **The starting boundary defect** (Lemma in:seedbounds, first inequality). With the residual
bounds of `residual_le_of_bounds`, a ground weight `g ≤ G(2N)` with `g > 0` and
`(ρ̄ U)^{2N-2K} Γ ≤ ε g`, the boundary defect at length `2N` is at most `ε`:
`(1 - ε) G(2N) ≤ c μ^{2N}`. -/
theorem one_sub_mul_le_of_ground (h : OpenMomentWitness Z G μ c ρ s) (hμ : 0 < μ) {K N : ℕ}
    (hK : 1 ≤ K) (hKN : K ≤ N) {U ρ' Γ g ε : ℝ} (hμU : μ ≤ U) (hρ : ρ ≤ ρ')
    (hΓ : G (2 * K) ≤ Γ) (hg : 0 < g) (hgG : g ≤ G (2 * N))
    (hε : (ρ' * U) ^ (2 * N - 2 * K) * Γ ≤ ε * g) :
    (1 - ε) * G (2 * N) ≤ c * μ ^ (2 * N) := by
  have hres := h.residual_le_of_bounds hμ hK hKN hμU hρ hΓ
  have hD : 0 ≤ G (2 * N) - c * μ ^ (2 * N) := by linarith [h.le_G N (hK.trans hKN)]
  have hε0 : 0 ≤ ε := by
    by_contra hneg
    push Not at hneg
    nlinarith
  nlinarith [mul_le_mul_of_nonneg_left hgG hε0]

end OpenMomentWitness

/-- **The leading part from a polynomial bound.** If `c ≥ 0`, `c P(μ)² ≤ C` and the ratio
certificate `RatioCheck P k₀ v` holds, then for `v ≤ μ ≤ U`, `N ≥ max k₀ 1`,
`c μ^{2N} ≤ C U^{2N} / P(U)²`: the ratio `x^{2N} / P(x)²` increases on `[v, U]`. -/
theorem mul_pow_le_of_ratioCheck {p : List ℚ} {k₀ : ℕ} {v : ℚ} (hp : RatioCheck p k₀ v)
    (hv : 0 ≤ v) {N : ℕ} (hN : 1 ≤ N) (hk₀ : k₀ ≤ N) {c μ U C : ℝ} (hc : 0 ≤ c)
    (hμv : (v : ℝ) ≤ μ) (hμU : μ ≤ U) (hC : c * eval p μ ^ 2 ≤ C) :
    c * μ ^ (2 * N) ≤ C * (U ^ (2 * N) / eval p U ^ 2) := by
  have hP := hp.eval_pos hμv
  have hr := pow_div_sq_le_of_ratioCheck hp hv hN hk₀ hμv hμU
  have hr0 : 0 ≤ μ ^ (2 * N) / eval p μ ^ 2 := by
    have : 0 ≤ μ := (by exact_mod_cast hv : (0 : ℝ) ≤ v).trans hμv
    positivity
  calc c * μ ^ (2 * N) = c * eval p μ ^ 2 * (μ ^ (2 * N) / eval p μ ^ 2) := by
        field_simp
    _ ≤ C * (U ^ (2 * N) / eval p U ^ 2) :=
        mul_le_mul hC hr hr0 ((mul_nonneg hc (sq_nonneg _)).trans hC)

namespace OpenMomentWitness

variable {Z G : ℕ → ℝ} {μ c ρ s : ℝ}

/-- **An upper bound on an even boundary moment**: the leading part from a polynomial bound
`c P(μ)² ≤ C` and the residual part from `residual_le_of_bounds`,
`G(2N) ≤ C U^{2N} / P(U)² + (ρ̄ U)^{2N-2K} Γ`. -/
theorem le_add_of_bounds (h : OpenMomentWitness Z G μ c ρ s) (hμ : 0 < μ) {K N : ℕ}
    (hK : 1 ≤ K) (hKN : K ≤ N) {p : List ℚ} {k₀ : ℕ} {v : ℚ} (hp : RatioCheck p k₀ v)
    (hv : 0 ≤ v) (hk₀ : k₀ ≤ N) {U ρ' Γ C : ℝ} (hμv : (v : ℝ) ≤ μ) (hμU : μ ≤ U)
    (hρ : ρ ≤ ρ') (hΓ : G (2 * K) ≤ Γ) (hC : c * eval p μ ^ 2 ≤ C) :
    G (2 * N) ≤ C * (U ^ (2 * N) / eval p U ^ 2) + (ρ' * U) ^ (2 * N - 2 * K) * Γ := by
  have h1 := h.residual_le_of_bounds hμ hK hKN hμU hρ hΓ
  have h2 := mul_pow_le_of_ratioCheck hp hv (hK.trans hKN) hk₀ h.c_nonneg hμv hμU hC
  linarith

/-- **The boundary defect from one witness.** If `Z 4 > 0`, `ρ < 1`, `G k > 0` and
`(1 - E) G(k) ≤ c μᵏ`, then `e_k ≤ E`: the witness modulus and weight are the leading root and
the leading overlap of the traces. -/
theorem boundaryDefect_le (h : OpenMomentWitness Z G μ c ρ s) (hZ4 : 0 < Z 4) (hρ : ρ < 1)
    {k : ℕ} (hG : 0 < G k) {E : ℝ} (he : (1 - E) * G k ≤ c * μ ^ k) :
    boundaryDefect Z G k ≤ E := by
  rw [boundaryDefect, h.leadingRoot_eq hZ4, h.leadingOverlap_eq hZ4 hρ, sub_le_comm,
    le_div_iff₀ hG]
  linarith

end OpenMomentWitness

/-! ### Thermal purity from one ground weight -/

/-- **One ground weight.** Let `yᵢ ≥ 0`, with a chosen `y_{i₀} ≥ g ≥ 0` and `∑ yᵢ ≤ g (1 + R)`,
`R ≥ 0`. Then the other weights sum to at most `R y_{i₀}`, their `l`th powers to at most
`R^l y_{i₀}^l` (a sum of nonnegative `l`th powers is at most the `l`th power of the sum), and
`(∑ yᵢ^l)² ≤ (1 + R^l)² ∑ yᵢ^{2l}`. Uniqueness of the ground weight is not needed. -/
theorem sum_pow_sq_le_of_ground {κ : Type*} (s : Finset κ) (y : κ → ℝ)
    (hy : ∀ i ∈ s, 0 ≤ y i) {i₀ : κ} (hi₀ : i₀ ∈ s) {g R : ℝ} (hgy : g ≤ y i₀) (hR : 0 ≤ R)
    (hsum : ∑ i ∈ s, y i ≤ g * (1 + R)) {l : ℕ} (hl : 1 ≤ l) :
    (∑ i ∈ s, y i ^ l) ^ 2 ≤ (1 + R ^ l) ^ 2 * ∑ i ∈ s, y i ^ (2 * l) := by
  classical
  have hy0 : 0 ≤ y i₀ := hy i₀ hi₀
  have hy' : ∀ i ∈ s.erase i₀, 0 ≤ y i := fun i hi ↦ hy i (mem_of_mem_erase hi)
  have hsplit (n : ℕ) : ∑ i ∈ s, y i ^ n = y i₀ ^ n + ∑ i ∈ s.erase i₀, y i ^ n :=
    (add_sum_erase s (fun i ↦ y i ^ n) hi₀).symm
  have hsplit1 : ∑ i ∈ s, y i = y i₀ + ∑ i ∈ s.erase i₀, y i :=
    (add_sum_erase s y hi₀).symm
  -- The other weights carry at most `R y_{i₀}`.
  have hr1 : ∑ i ∈ s.erase i₀, y i ≤ R * y i₀ := by
    rw [hsplit1] at hsum; nlinarith
  have hr0 : 0 ≤ ∑ i ∈ s.erase i₀, y i := sum_nonneg hy'
  have hrl : ∑ i ∈ s.erase i₀, y i ^ l ≤ R ^ l * y i₀ ^ l := by
    have a := sum_pow_pow_le _ y hy' le_rfl hl
    simp only [pow_one] at a
    calc ∑ i ∈ s.erase i₀, y i ^ l ≤ (∑ i ∈ s.erase i₀, y i) ^ l := a
      _ ≤ (R * y i₀) ^ l := pow_le_pow_left₀ hr0 hr1 l
      _ = R ^ l * y i₀ ^ l := mul_pow _ _ _
  have hZ : ∑ i ∈ s, y i ^ l ≤ y i₀ ^ l * (1 + R ^ l) := by rw [hsplit]; linarith
  have hZ0 : 0 ≤ ∑ i ∈ s, y i ^ l := sum_nonneg fun i hi ↦ pow_nonneg (hy i hi) l
  have h2 : (y i₀ ^ l) ^ 2 ≤ ∑ i ∈ s, y i ^ (2 * l) := by
    rw [hsplit, ← pow_mul, mul_comm l 2]
    linarith [sum_nonneg fun i hi ↦ pow_nonneg (hy' i hi) (2 * l)]
  calc (∑ i ∈ s, y i ^ l) ^ 2 ≤ (y i₀ ^ l * (1 + R ^ l)) ^ 2 := pow_le_pow_left₀ hZ0 hZ 2
    _ = (1 + R ^ l) ^ 2 * (y i₀ ^ l) ^ 2 := by ring
    _ ≤ (1 + R ^ l) ^ 2 * ∑ i ∈ s, y i ^ (2 * l) := by gcongr

section Ground

variable {ι : ℕ → Type*} [∀ n, Fintype (ι n)] [∀ n, Nonempty (ι n)]

/-- **The thermal purity from a ground weight** (Lemma in:seedbounds, second inequality). Let
some energy `E_{i₀} ≤ e₀` of the chain of length `n` (a variational bound), and let the partition
function at `β₀ ≥ 0` satisfy `Z_n(β₀) ≤ e^{-β₀ e₀} (1 + R)` with `R ≥ 0`. Then at every
`β_J = 2^J β₀` the thermal purity is `P_τ(n, β_J) ≥ (1 + R^{2^J})⁻²`. -/
theorem thermalPurity_ge_of_ground (E : (n : ℕ) → ι n → ℝ) {n : ℕ} {β₀ e₀ R : ℝ}
    (hβ₀ : 0 ≤ β₀) {i₀ : ι n} (hi₀ : E n i₀ ≤ e₀) (hR : 0 ≤ R)
    (hZ : partitionFunction E n β₀ ≤ exp (-(β₀ * e₀)) * (1 + R)) (J : ℕ) :
    1 / (1 + R ^ 2 ^ J) ^ 2 ≤ thermalPurity E n (2 ^ J * β₀) := by
  -- The Boltzmann factors at `2^J β₀` and `2^{J+1} β₀` are powers of `yᵢ = e^{-β₀ Eᵢ}`.
  have hpow (k : ℕ) (e : ℝ) : exp (-((k * β₀) * e)) = exp (-(β₀ * e)) ^ k := by
    rw [← exp_nat_mul]; ring_nf
  have e1 : partitionFunction E n (2 ^ J * β₀) = ∑ i, exp (-(β₀ * E n i)) ^ 2 ^ J := by
    simp only [partitionFunction]
    exact sum_congr rfl fun i _ ↦ by rw [← hpow]; push_cast; ring_nf
  have e2 : partitionFunction E n (2 * (2 ^ J * β₀)) =
      ∑ i, exp (-(β₀ * E n i)) ^ (2 * 2 ^ J) := by
    simp only [partitionFunction]
    exact sum_congr rfl fun i _ ↦ by rw [← hpow]; push_cast; ring_nf
  have hg : exp (-(β₀ * e₀)) ≤ exp (-(β₀ * E n i₀)) :=
    exp_le_exp.2 (by nlinarith)
  have key := sum_pow_sq_le_of_ground univ (fun i ↦ exp (-(β₀ * E n i)))
    (fun _ _ ↦ (exp_pos _).le) (mem_univ i₀) hg hR hZ (Nat.one_le_two_pow (n := J))
  have hpos := partitionFunction_pos E n (2 ^ J * β₀)
  rw [thermalPurity, div_le_div_iff₀ (by positivity) (pow_pos hpos 2), one_mul, e1, e2]
  linarith

end Ground

end CertifiedGaps
