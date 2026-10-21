/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.MomentBounds

/-!
# The boundary-field gap criterion for open chains

This file formalizes Appendices A–C of the thermal-MPO note: a uniform gap for a family of open
chains whose partition functions `G n = Tr e^{-β H_n}` share a transfer operator with a family
of rings `Z L = Tr e^{-β H_L^{ring}}` (`OpenMomentBounds`). For the spin-one Heisenberg chain,
`H_n` has equal end fields on odd chains and opposite ones on even chains, and the theorem is
the boundary-field gap of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain*.

## Setting

The ring spectra are `Er L`, the open-chain spectra `Eo n`. The scales are `βᵢ = β₀ 2^i` and
`Nᵢ = N₀ 2^i`, with ring lengths `2Nᵢ` and open-chain lengths `Mᵢ = 4Nᵢ`; for the note's seeds
(`β₀ = 49`, `N₀ = 60`) the note's scale `j` is `i + 2`. At each `βᵢ` we fix a witness
`(μᵢ, cᵢ, ρᵢ, sᵢ)` of the boundary moment inequalities. The three claims of the note are, in
multiplicative form,
* (R) ring purities `P_s(2Nᵢ, βᵢ), P_τ(Mᵢ, βᵢ) ≥ 1 - εᵢ` with `εᵢ = (C u₀)^{2^i} / C`;
* (E) the boundary defect `e_{2Kᵢ}(βᵢ) ≤ Eᵢ`, written `(1 - Eᵢ) G(2Kᵢ) ≤ cᵢ μᵢ^{2Kᵢ}`, with the
  seed's `K₀` and `Kᵢ = Nᵢ` for `i ≥ 1` (`boundaryHalf`);
* (T) the open-chain thermal purity `P_τ(Mᵢ, βᵢ) ≥ 1 - ηᵢ`.

## Main results

* `OpenMomentWitness.tail_le`, `mu_le_sq`, `mu_pow_ge`: claim (R) as transfer bounds.
* `move_one`, `move_two`, `oddChain_purity`: the overlap recurrence, the length doubling and
  the equal-field purity of Appendix A, all in multiplicative form; `evenChain_purity` is the
  purity at even lengths.
* `boundaryDefect_le_of_traces`: the boundary seed from three finite traces.
* `open_step`: one step of the boundary induction.
* `gapGT_of_openChain`: the boundary-field gap `Δ_n > log(1/(10x)) / (β₀ 2^m)` for every
  `n ≥ 4 N₀ 2^m`, odd or even, from general seeds, `m` base steps (`OpenStep`) and the
  hand-over `x`. `Open/Certificate.lean` checks the base steps in rational arithmetic.
-/

open Filter Finset Real

namespace CertifiedGaps

/-! ### Scalar functions -/

/-- The tail bound `d(ε) = ε / (2 - 3ε)` of the note, eq. (B1). -/
noncomputable def tailBound (ε : ℝ) : ℝ := ε / (2 - 3 * ε)

/-- The overlap update `F(q, e, t) = q/(1 - q) · (1/((1 - e)² t) - 1)` of the note, eq. (C2). -/
noncomputable def overlapUpdate (q e t : ℝ) : ℝ := q / (1 - q) * (1 / ((1 - e) ^ 2 * t) - 1)

/-- `d` is monotone on `[0, 2/3)`. -/
theorem tailBound_mono {ε ε' : ℝ} (h : ε ≤ ε') (h1 : ε' < 2 / 3) :
    tailBound ε ≤ tailBound ε' := by
  unfold tailBound
  rw [div_le_div_iff₀ (by linarith) (by linarith)]
  nlinarith

/-- **A purity close to one forces one dominant weight** (Appendix B). If the largest weight
`a/Z` of a spectrum with total `Z` satisfies `Y ≤ a Z` and `Y ≤ a² + (Z - a)²`, where `Y / Z²` is
the purity, and the purity is at least `1 - ε` with `ε ≤ 1/10`, then the remaining weights are at
most `d(ε)` times the largest. -/
theorem sub_le_tailBound_mul {a Z Y ε : ℝ} (ha : 0 < a) (haZ : a ≤ Z) (hY₁ : Y ≤ a * Z)
    (hY₂ : Y ≤ a ^ 2 + (Z - a) ^ 2) (hε : ε ≤ 1 / 10) (hP : (1 - ε) * Z ^ 2 ≤ Y) :
    Z - a ≤ tailBound ε * a := by
  have hZ : 0 < Z := ha.trans_le haZ
  set u := Z - a with hu
  have hu0 : 0 ≤ u := by linarith
  have hZu : Z = a + u := by linarith
  rw [hZu] at hY₁ hP
  -- `(1 - ε)(a + u) ≤ a`, so `u ≤ a`.
  have h1 : (1 - ε) * (a + u) ≤ a := by
    have := hP.trans hY₁
    nlinarith
  have hua : u ≤ a := by nlinarith
  -- `2 a u ≤ ε (a + u)²`, hence `(2 - 3ε) u ≤ ε a`.
  have h2 : 2 * a * u ≤ ε * (a + u) ^ 2 := by nlinarith
  have h3 : (2 - 3 * ε) * u ≤ ε * a := by
    have : ε * u ^ 2 ≤ ε * a * u := by
      rcases le_or_gt 0 ε with hε0 | hε0
      · nlinarith [mul_nonneg hε0 (mul_nonneg hu0 (sub_nonneg.2 hua))]
      · nlinarith
    nlinarith
  unfold tailBound
  rw [div_mul_eq_mul_div, le_div_iff₀ (by linarith)]
  linarith

/-- **The bulk factor on a shorter length.** If `R μ^{2(2K)} ≤ μ'^{2K}` and `μ' ≤ μ²`, then
`R μ^{2a} ≤ μ'^a` for every `a ≤ 2K`. -/
theorem bulk_pow_le {μ μ' R : ℝ} (hμ' : 0 < μ') (hr : μ' ≤ μ ^ 2) (hR0 : 0 ≤ R) {a K : ℕ}
    (ha : a ≤ 2 * K)
    (hR : R * μ ^ (2 * (2 * K)) ≤ μ' ^ (2 * K)) : R * (μ ^ a) ^ 2 ≤ μ' ^ a := by
  set b := 2 * K - a
  have hab : 2 * K = a + b := by omega
  have hμb : μ' ^ b ≤ (μ ^ b) ^ 2 := by
    rw [← pow_mul, mul_comm, pow_mul]; exact pow_le_pow_left₀ hμ'.le hr _
  have h4 : μ ^ (2 * (2 * K)) = (μ ^ a) ^ 2 * (μ ^ b) ^ 2 := by
    rw [hab, ← mul_pow, ← pow_add, ← pow_mul, mul_comm]
  rw [h4, hab, pow_add] at hR
  have hμb0 : 0 < μ' ^ b := pow_pos hμ' _
  have : R * (μ ^ a) ^ 2 * μ' ^ b ≤ μ' ^ a * μ' ^ b := by
    calc R * (μ ^ a) ^ 2 * μ' ^ b ≤ R * (μ ^ a) ^ 2 * (μ ^ b) ^ 2 := by gcongr
      _ = R * ((μ ^ a) ^ 2 * (μ ^ b) ^ 2) := by ring
      _ ≤ μ' ^ a * μ' ^ b := hR
  exact le_of_mul_le_mul_right this hμb0

/-! ### Transfer bounds from one witness -/

namespace OpenMomentWitness

variable {Z G : ℕ → ℝ} {μ c ρ s : ℝ}

/-- The largest weight of `σ_{2N}`: `P_s(2N) ≥ 1 - ε` gives `(1 - ε) Z(2N) ≤ μ^{2N}`. -/
theorem weight_ge (h : OpenMomentWitness Z G μ c ρ s) {N : ℕ} (hN : 1 ≤ N)
    (hZ : 0 < Z (2 * N)) {ε : ℝ} (hP : 1 - ε ≤ Z (2 * (2 * N)) / Z (2 * N) ^ 2) :
    (1 - ε) * Z (2 * N) ≤ μ ^ (2 * N) := by
  rw [le_div_iff₀ (by positivity)] at hP
  have := h.Z_two_mul_le N hN
  nlinarith

/-- **The residual ratio from a purity** (the note, eq. (B1)): `P_s(2N) ≥ 1 - ε` with
`ε ≤ 1/10` gives `ρ^{2N} ≤ d(ε)`. -/
theorem tail_le (h : OpenMomentWitness Z G μ c ρ s) (hμ : 0 < μ) {N : ℕ} (hN : 1 ≤ N)
    (hZ : 0 < Z (2 * N)) {ε : ℝ} (hε : ε ≤ 1 / 10)
    (hP : 1 - ε ≤ Z (2 * (2 * N)) / Z (2 * N) ^ 2) : ρ ^ (2 * N) ≤ tailBound ε := by
  rw [le_div_iff₀ (by positivity)] at hP
  have hμN : 0 < μ ^ (2 * N) := pow_pos hμ _
  have e : μ ^ (2 * (2 * N)) = (μ ^ (2 * N)) ^ 2 := by ring
  have key := sub_le_tailBound_mul hμN (h.pow_le_Z N hN) (h.Z_two_mul_le N hN)
    (by rw [← e]; exact h.Z_two_mul_le_sq N hN) hε hP
  have htail := h.tail N hN
  rw [mul_pow] at htail
  exact le_of_mul_le_mul_right (htail.trans key) hμN

/-- `ρ < 1` once `P_s(2N) ≥ 1 - ε` with `ε ≤ 1/10`. -/
theorem rho_lt_one (h : OpenMomentWitness Z G μ c ρ s) (hμ : 0 < μ) {N : ℕ} (hN : 1 ≤ N)
    (hZ : 0 < Z (2 * N)) {ε : ℝ} (hε : ε ≤ 1 / 10)
    (hP : 1 - ε ≤ Z (2 * (2 * N)) / Z (2 * N) ^ 2) : ρ < 1 := by
  have h1 := h.tail_le hμ hN hZ hε hP
  have hd : tailBound ε < 1 := by
    unfold tailBound; rw [div_lt_one (by linarith)]; linarith
  by_contra hρ
  replace hρ := not_lt.1 hρ
  linarith [one_le_pow₀ (n := 2 * N) hρ]

/-- **The leading eigenvalue does not grow faster than squaring** (Appendix B): if
`Z'(2N) ≤ Z(2N)²` at every length, as for the traces at `2β` and `β`, then `μ' ≤ μ²`. -/
theorem mu_le_sq {Z' G' : ℕ → ℝ} {μ' c' ρ' s' : ℝ} (h : OpenMomentWitness Z G μ c ρ s)
    (h' : OpenMomentWitness Z' G' μ' c' ρ' s') (hμ : 0 < μ)
    (hZZ : ∀ N, 1 ≤ N → Z' (2 * N) ≤ Z (2 * N) ^ 2) : μ' ≤ μ ^ 2 := by
  set K := Z 2 / μ ^ 2
  have hμ2 : 0 < μ ^ 2 := by positivity
  -- `Z(2N) ≤ μ^{2N} K`, so `(μ'/μ²)^{2N} ≤ K²` for every `N ≥ 1`.
  have hbound : ∀ N, 1 ≤ N → (μ' / μ ^ 2) ^ (2 * N) ≤ K ^ 2 := by
    intro N hN
    have h1 := h.Z_anti 1 N le_rfl hN
    have hZN : Z (2 * N) ≤ μ ^ (2 * N) * K := by
      simp only [K]; rw [mul_div_assoc', le_div_iff₀ hμ2]; simpa [mul_comm] using h1
    have hZ0 : 0 ≤ Z (2 * N) := (pow_nonneg h.mu_nonneg _).trans (h.pow_le_Z N hN)
    have h2 : μ' ^ (2 * N) ≤ (μ ^ (2 * N) * K) ^ 2 :=
      (h'.pow_le_Z N hN).trans ((hZZ N hN).trans (pow_le_pow_left₀ hZ0 hZN 2))
    rw [div_pow, div_le_iff₀ (by positivity)]
    calc μ' ^ (2 * N) ≤ (μ ^ (2 * N) * K) ^ 2 := h2
      _ = K ^ 2 * (μ ^ 2) ^ (2 * N) := by ring
  by_contra hlt
  replace hlt := not_le.1 hlt
  have hr : 1 < μ' / μ ^ 2 := (one_lt_div hμ2).2 hlt
  obtain ⟨N, hN⟩ := pow_unbounded_of_one_lt (K ^ 2) hr
  have := (pow_le_pow_right₀ hr.le (by omega : N ≤ 2 * (N + 1))).trans
    (hbound (N + 1) (by omega))
  linarith

/-- **The bulk ratio is close to one** (Appendix B): with `P_τ(2N, β) ≥ 1 - ε` (that is
`(1 - ε) Z(2N)² ≤ Z'(2N)`) and `P_s(2N, 2β) ≥ 1 - ε'`,
`(1 - ε)(1 - ε') μ^{2(2N)} ≤ μ'^{2N}`. -/
theorem mu_pow_ge {Z' G' : ℕ → ℝ} {μ' c' ρ' s' : ℝ} (h : OpenMomentWitness Z G μ c ρ s)
    (h' : OpenMomentWitness Z' G' μ' c' ρ' s') {N : ℕ} (hN : 1 ≤ N) (hZ' : 0 < Z' (2 * N))
    {ε ε' : ℝ} (hε : ε ≤ 1) (hε' : ε' ≤ 1) (hτ : (1 - ε) * Z (2 * N) ^ 2 ≤ Z' (2 * N))
    (hP' : 1 - ε' ≤ Z' (2 * (2 * N)) / Z' (2 * N) ^ 2) :
    (1 - ε) * (1 - ε') * μ ^ (2 * (2 * N)) ≤ μ' ^ (2 * N) := by
  have hw := h'.weight_ge hN hZ' hP'
  have hμN := h.pow_le_Z N hN
  have hsq : μ ^ (2 * (2 * N)) ≤ Z (2 * N) ^ 2 := by
    rw [mul_comm 2 (2 * N), pow_mul]
    exact pow_le_pow_left₀ (pow_nonneg h.mu_nonneg _) hμN 2
  have h1 : (1 - ε) * μ ^ (2 * (2 * N)) ≤ Z' (2 * N) :=
    (mul_le_mul_of_nonneg_left hsq (by linarith)).trans hτ
  calc (1 - ε) * (1 - ε') * μ ^ (2 * (2 * N)) = (1 - ε') * ((1 - ε) * μ ^ (2 * (2 * N))) := by
        ring
    _ ≤ (1 - ε') * Z' (2 * N) := mul_le_mul_of_nonneg_left h1 (by linarith)
    _ ≤ μ' ^ (2 * N) := by linarith

/-! ### The three moves of Appendix A

Throughout, the unprimed witness is at `β` and the primed one at `2β`, with `μ' ≤ μ²`. The
open-chain thermal purity `A(m) = G'(m)/G(m)²` enters through `G'(m) ≤ G(m)²` and lower bounds
`t G(m)² ≤ G'(m)`. -/

variable {Z' G' : ℕ → ℝ} {μ' c' ρ' s' : ℝ}

/-- **Move 1: the overlap recurrence**, eqs. (A18)–(A19). Let `k = 2K < m = 2N`. If
`e_k ≤ E < 1` at `β`, the thermal purity at length `m` is at least `t > 0`, and
`ρ'^{m-k} ≤ q < 1`, then `e'_m ≤ F(q, E, t)` at `2β`, in the form
`(1 - F) G'(m) ≤ c' μ'^m`. -/
theorem move_one (h : OpenMomentWitness Z G μ c ρ s) (h' : OpenMomentWitness Z' G' μ' c' ρ' s')
    (hμ' : 0 < μ') (hr : μ' ≤ μ ^ 2) {K N : ℕ} (hK : 1 ≤ K) (hKN : K ≤ N)
    (hAk : G' (2 * K) ≤ G (2 * K) ^ 2) {t : ℝ} (ht : 0 < t)
    (hAm : t * G (2 * N) ^ 2 ≤ G' (2 * N)) {E : ℝ} (hE : E < 1)
    (he : (1 - E) * G (2 * K) ≤ c * μ ^ (2 * K)) {q : ℝ} (hq0 : 0 ≤ q) (hq1 : q < 1)
    (hq : ρ' ^ (2 * N - 2 * K) ≤ q) :
    (1 - overlapUpdate q E t) * G' (2 * N) ≤ c' * μ' ^ (2 * N) := by
  set p := μ' ^ (2 * N - 2 * K)
  have hp0 : 0 < p := pow_pos hμ' _
  have hsplit : μ' ^ (2 * N) = μ' ^ (2 * K) * p := by rw [← pow_add]; congr 1; omega
  have ha : 0 < μ' ^ (2 * K) := pow_pos hμ' _
  have hDk' : 0 ≤ G' (2 * K) - c' * μ' ^ (2 * K) := by linarith [h'.le_G K hK]
  -- The decay at `2β`: `G'(m) - q p G'(k) ≤ (1 - q) c' μ'^m`.
  have key1 : G' (2 * N) - q * p * G' (2 * K) ≤ (1 - q) * (c' * μ' ^ (2 * N)) := by
    have hd := h'.G_decay K N hK hKN
    have hd' : μ' ^ (2 * K) * (G' (2 * N) - c' * μ' ^ (2 * N)) ≤
        q * μ' ^ (2 * N) * (G' (2 * K) - c' * μ' ^ (2 * K)) :=
      hd.trans (by gcongr)
    rw [hsplit] at hd' ⊢
    have : μ' ^ (2 * K) * (G' (2 * N) - q * p * G' (2 * K)) ≤
        μ' ^ (2 * K) * ((1 - q) * (c' * (μ' ^ (2 * K) * p))) := by nlinarith
    exact le_of_mul_le_mul_left this ha
  -- The ratio at `β`: `(1 - E)² t p G'(k) ≤ G'(m)`.
  have key2 : (1 - E) ^ 2 * t * (p * G' (2 * K)) ≤ G' (2 * N) := by
    have hGk : 0 ≤ G (2 * K) := (mul_nonneg h.c_nonneg (pow_nonneg h.mu_nonneg _)).trans
      (h.le_G K hK)
    have hμk : μ ^ (2 * N) = μ ^ (2 * K) * μ ^ (2 * N - 2 * K) := by
      rw [← pow_add]; congr 1; omega
    have h1 : (1 - E) * G (2 * K) * μ ^ (2 * N - 2 * K) ≤ G (2 * N) := by
      calc (1 - E) * G (2 * K) * μ ^ (2 * N - 2 * K)
          ≤ c * μ ^ (2 * K) * μ ^ (2 * N - 2 * K) := by
            gcongr; exact pow_nonneg h.mu_nonneg _
        _ = c * μ ^ (2 * N) := by rw [hμk]; ring
        _ ≤ G (2 * N) := h.le_G N (hK.trans hKN)
    have h10 : 0 ≤ (1 - E) * G (2 * K) * μ ^ (2 * N - 2 * K) := by
      have := pow_nonneg h.mu_nonneg (2 * N - 2 * K)
      have : 0 ≤ 1 - E := by linarith
      positivity
    have h2 : ((1 - E) * G (2 * K) * μ ^ (2 * N - 2 * K)) ^ 2 ≤ G (2 * N) ^ 2 :=
      pow_le_pow_left₀ h10 h1 2
    have hpμ : p ≤ (μ ^ (2 * N - 2 * K)) ^ 2 := by
      rw [← pow_mul, mul_comm, pow_mul]
      exact pow_le_pow_left₀ hμ'.le hr _
    have hGk' : 0 ≤ G' (2 * K) := by linarith [mul_nonneg h'.c_nonneg ha.le]
    calc (1 - E) ^ 2 * t * (p * G' (2 * K))
        ≤ (1 - E) ^ 2 * t * ((μ ^ (2 * N - 2 * K)) ^ 2 * G (2 * K) ^ 2) := by gcongr
      _ = t * ((1 - E) * G (2 * K) * μ ^ (2 * N - 2 * K)) ^ 2 := by ring
      _ ≤ t * G (2 * N) ^ 2 := by gcongr
      _ ≤ G' (2 * N) := hAm
  -- Combine: `(1 - F) G'(m) = (G'(m) - q G'(m) / w) / (1 - q)` with `w = (1 - E)² t`.
  have hw : 0 < (1 - E) ^ 2 * t := mul_pos (pow_pos (by linarith) 2) ht
  have h1q : 0 < 1 - q := by linarith
  have hF : (1 - overlapUpdate q E t) * G' (2 * N) =
      (G' (2 * N) - q * (G' (2 * N) / ((1 - E) ^ 2 * t))) / (1 - q) := by
    unfold overlapUpdate
    field_simp
    ring
  have hle : q * (p * G' (2 * K)) ≤ q * (G' (2 * N) / ((1 - E) ^ 2 * t)) := by
    gcongr
    rw [le_div_iff₀ hw]; linarith
  rw [hF, div_le_iff₀ h1q]
  linarith

/-- **Move 2: doubling the length**, eq. (A20). If `e'_{2M} ≤ E'` at `2β` and the bulk factor
satisfies `R μ^{2(2M)} ≤ μ'^{2M}`, then `A(4M) ≥ R (1 - E') A(2M)`, in the form
`R (1 - E') G'(2M) G(4M)² ≤ G'(4M) G(2M)²`. -/
theorem move_two (h : OpenMomentWitness Z G μ c ρ s) (h' : OpenMomentWitness Z' G' μ' c' ρ' s')
    (hμ : 0 < μ) {M : ℕ} (hM : 1 ≤ M) {E' : ℝ} (hE' : E' ≤ 1)
    (he' : (1 - E') * G' (2 * M) ≤ c' * μ' ^ (2 * M)) {R : ℝ} (hR0 : 0 ≤ R)
    (hR : R * μ ^ (2 * (2 * M)) ≤ μ' ^ (2 * M)) :
    R * (1 - E') * G' (2 * M) * G (2 * (2 * M)) ^ 2 ≤ G' (2 * (2 * M)) * G (2 * M) ^ 2 := by
  have hμM : 0 < μ ^ (2 * M) := pow_pos hμ _
  have hsq : μ ^ (2 * (2 * M)) = μ ^ (2 * M) * μ ^ (2 * M) := by rw [← pow_add]; ring_nf
  -- `G(4M) ≤ μ^{2M} G(2M)`: the normalized trace decreases, since `ρ ≤ 1`.
  have hdec : G (2 * (2 * M)) ≤ μ ^ (2 * M) * G (2 * M) := by
    have hd := h.G_decay M (2 * M) hM (by omega)
    have hρ : ρ ^ (2 * (2 * M) - 2 * M) ≤ 1 := pow_le_one₀ h.rho_nonneg h.rho_le_one
    have hD : 0 ≤ G (2 * M) - c * μ ^ (2 * M) := by linarith [h.le_G M hM]
    have : μ ^ (2 * M) * (G (2 * (2 * M)) - c * μ ^ (2 * (2 * M))) ≤
        μ ^ (2 * (2 * M)) * (G (2 * M) - c * μ ^ (2 * M)) := by
      refine hd.trans ?_
      have := mul_le_mul_of_nonneg_right hρ (mul_nonneg (pow_nonneg hμ.le (2 * (2 * M))) hD)
      linarith [this]
    rw [hsq] at this
    nlinarith
  have hG4 : 0 ≤ G (2 * (2 * M)) := (mul_nonneg h.c_nonneg (pow_nonneg hμ.le _)).trans
    (h.le_G (2 * M) (by omega))
  have hG4sq : G (2 * (2 * M)) ^ 2 ≤ μ ^ (2 * (2 * M)) * G (2 * M) ^ 2 := by
    rw [hsq]
    nlinarith [pow_le_pow_left₀ hG4 hdec 2]
  -- `G'(4M) ≥ c' μ'^{4M} ≥ (1 - E') G'(2M) μ'^{2M}`.
  have hG'4 : (1 - E') * G' (2 * M) * μ' ^ (2 * M) ≤ G' (2 * (2 * M)) := by
    have := h'.le_G (2 * M) (by omega)
    have hsq' : μ' ^ (2 * (2 * M)) = μ' ^ (2 * M) * μ' ^ (2 * M) := by
      rw [← pow_add]; ring_nf
    rw [hsq'] at this
    have hμ'M := pow_nonneg h'.mu_nonneg (2 * M)
    nlinarith [mul_le_mul_of_nonneg_right he' hμ'M]
  have hG'2 : 0 ≤ G' (2 * M) := by
    linarith [h'.le_G M hM, mul_nonneg h'.c_nonneg (pow_nonneg h'.mu_nonneg (2 * M))]
  have h1E : 0 ≤ 1 - E' := by linarith
  calc R * (1 - E') * G' (2 * M) * G (2 * (2 * M)) ^ 2
      ≤ R * (1 - E') * G' (2 * M) * (μ ^ (2 * (2 * M)) * G (2 * M) ^ 2) := by gcongr
    _ = (1 - E') * G' (2 * M) * (R * μ ^ (2 * (2 * M))) * G (2 * M) ^ 2 := by ring
    _ ≤ (1 - E') * G' (2 * M) * μ' ^ (2 * M) * G (2 * M) ^ 2 := by gcongr
    _ ≤ G' (2 * (2 * M)) * G (2 * M) ^ 2 := by gcongr

/-- **The equal-field purity at odd lengths**, eqs. (A23)–(A26). Let `M = 2K` and let `n` be
odd with `M ≤ n ≤ 2M`. If `e'_M ≤ E' < 1/2` at `2β`, `A(M) ≥ 1 - η` and
`R μ^{2M} ≤ μ'^M`, then `A(n) ≥ (1 - 2E')(1 - η) R`. The sign `s'` is forced to be `+1` by
`G'(n) > 0`. -/
theorem oddChain_purity (h : OpenMomentWitness Z G μ c ρ s)
    (h' : OpenMomentWitness Z' G' μ' c' ρ' s') (hμ' : 0 < μ') (hr : μ' ≤ μ ^ 2) {K N : ℕ}
    (hK : 1 ≤ K) (hKN : K ≤ N) (hN2 : 2 * N + 1 ≤ 2 * (2 * K)) (hGn' : 0 < G' (2 * N + 1))
    {E' : ℝ} (hE' : E' < 1 / 2) (he' : (1 - E') * G' (2 * K) ≤ c' * μ' ^ (2 * K)) {η : ℝ}
    (hη1 : η ≤ 1) (hη : (1 - η) * G (2 * K) ^ 2 ≤ G' (2 * K)) {R : ℝ} (hR0 : 0 ≤ R)
    (hR : R * μ ^ (2 * (2 * K)) ≤ μ' ^ (2 * K)) :
    (1 - 2 * E') * (1 - η) * R * G (2 * N + 1) ^ 2 ≤ G' (2 * N + 1) := by
  set a := 2 * N + 1 - 2 * K with ha
  have hμ := h.mu_nonneg
  have hsplit : ∀ x : ℝ, x ^ (2 * N + 1) = x ^ a * x ^ (2 * K) := fun x ↦ by
    rw [← pow_add]; congr 1; omega
  have hGk : 0 ≤ G (2 * K) := (mul_nonneg h.c_nonneg (pow_nonneg hμ _)).trans (h.le_G K hK)
  have hG'k : 0 ≤ G' (2 * K) :=
    (mul_nonneg h'.c_nonneg (pow_nonneg hμ'.le _)).trans (h'.le_G K hK)
  -- `|G(n)| ≤ μ^a G(2K)` at `β`.
  have habs : |G (2 * N + 1)| ≤ μ ^ a * G (2 * K) := by
    have hodd := h.G_odd K N hK hKN
    have hs : |s| = 1 := by
      have := h.s_sq; rw [← sq_abs] at this
      nlinarith [abs_nonneg s, sq_nonneg (|s| - 1), sq_nonneg (|s| + 1)]
    have hscn : |s * c * μ ^ (2 * N + 1)| = c * μ ^ (2 * N + 1) := by
      rw [abs_mul, abs_mul, hs, abs_of_nonneg h.c_nonneg, abs_of_nonneg (pow_nonneg hμ _),
        one_mul]
    calc |G (2 * N + 1)| ≤ |G (2 * N + 1) - s * c * μ ^ (2 * N + 1)| +
          |s * c * μ ^ (2 * N + 1)| := by
          simpa using abs_add_le (G (2 * N + 1) - s * c * μ ^ (2 * N + 1))
            (s * c * μ ^ (2 * N + 1))
      _ ≤ μ ^ a * (G (2 * K) - c * μ ^ (2 * K)) + c * μ ^ (2 * N + 1) := by
          rw [hscn]; exact add_le_add_left hodd _
      _ = μ ^ a * G (2 * K) := by rw [hsplit]; ring
  have hGn2 : G (2 * N + 1) ^ 2 ≤ (μ ^ a) ^ 2 * G (2 * K) ^ 2 := by
    rw [← mul_pow, ← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) habs 2
  -- At `2β` the sign is `+1`, and `G'(n) ≥ (1 - 2E') μ'^a G'(2K)`.
  have hodd' := h'.G_odd K N hK hKN
  rw [← ha] at hodd'
  have hc' : (1 - E') * G' (2 * K) ≤ c' * μ' ^ (2 * K) := he'
  have hμa' : 0 < μ' ^ a := pow_pos hμ' _
  have hs' : s' = 1 := by
    have hsq := h'.s_sq
    have : (s' - 1) * (s' + 1) = 0 := by linear_combination hsq
    rcases mul_eq_zero.1 this with h1 | h1
    · linarith
    · exfalso
      have hs' : s' = -1 := by linarith
      rw [hs', hsplit] at hodd'
      have := (abs_le.1 hodd').2
      nlinarith [mul_le_mul_of_nonneg_left hc' hμa'.le,
        mul_nonneg (mul_nonneg hμa'.le hG'k) (by linarith : (0 : ℝ) ≤ 1 - 2 * E')]
  rw [hs', one_mul, hsplit] at hodd'
  have hlow : (1 - 2 * E') * (μ' ^ a * G' (2 * K)) ≤ G' (2 * N + 1) := by
    have := (abs_le.1 hodd').1
    nlinarith [mul_le_mul_of_nonneg_left hc' hμa'.le]
  -- The bulk factor: `R μ^{2a} ≤ μ'^a`, from `R μ^{4K} ≤ μ'^{2K}`, `a ≤ 2K` and `μ' ≤ μ²`.
  have hbulk := bulk_pow_le hμ' hr hR0 (by omega : a ≤ 2 * K) hR
  have h12 : 0 ≤ 1 - 2 * E' := by linarith
  have h1η : 0 ≤ 1 - η := by linarith
  calc (1 - 2 * E') * (1 - η) * R * G (2 * N + 1) ^ 2
      ≤ (1 - 2 * E') * (1 - η) * R * ((μ ^ a) ^ 2 * G (2 * K) ^ 2) := by gcongr
    _ = (1 - 2 * E') * (R * (μ ^ a) ^ 2) * ((1 - η) * G (2 * K) ^ 2) := by ring
    _ ≤ (1 - 2 * E') * μ' ^ a * G' (2 * K) := by gcongr
    _ = (1 - 2 * E') * (μ' ^ a * G' (2 * K)) := by ring
    _ ≤ G' (2 * N + 1) := hlow

/-- **The open-chain purity at even lengths**, the even analogue of `oddChain_purity`. Let
`M = 2K` and let `n = 2N` with `M ≤ n ≤ 2M`. If `e'_M ≤ E' ≤ 1` at `2β`, `A(M) ≥ 1 - η` and
`R μ^{2M} ≤ μ'^M`, then `A(n) ≥ (1 - E')(1 - η) R`. At even lengths the leading weight has no
sign, so the factor is `1 - E'` rather than `1 - 2E'`. -/
theorem evenChain_purity (h : OpenMomentWitness Z G μ c ρ s)
    (h' : OpenMomentWitness Z' G' μ' c' ρ' s') (hμ' : 0 < μ') (hr : μ' ≤ μ ^ 2) {K N : ℕ}
    (hK : 1 ≤ K) (hKN : K ≤ N) (hN2 : N ≤ 2 * K) {E' : ℝ} (hE' : E' ≤ 1)
    (he' : (1 - E') * G' (2 * K) ≤ c' * μ' ^ (2 * K)) {η : ℝ} (hη1 : η ≤ 1)
    (hη : (1 - η) * G (2 * K) ^ 2 ≤ G' (2 * K)) {R : ℝ} (hR0 : 0 ≤ R)
    (hR : R * μ ^ (2 * (2 * K)) ≤ μ' ^ (2 * K)) :
    (1 - E') * (1 - η) * R * G (2 * N) ^ 2 ≤ G' (2 * N) := by
  set a := 2 * N - 2 * K with ha
  have hμ := h.mu_nonneg
  have hsplit : ∀ x : ℝ, x ^ (2 * N) = x ^ a * x ^ (2 * K) := fun x ↦ by
    rw [← pow_add]; congr 1; omega
  have hμpos : 0 < μ := by
    rcases hμ.lt_or_eq with h0 | h0
    · exact h0
    · rw [← h0] at hr; nlinarith
  have hGk : 0 ≤ G (2 * K) := (mul_nonneg h.c_nonneg (pow_nonneg hμ _)).trans (h.le_G K hK)
  have hG'k : 0 ≤ G' (2 * K) :=
    (mul_nonneg h'.c_nonneg (pow_nonneg hμ'.le _)).trans (h'.le_G K hK)
  -- `G(2N) ≤ μ^a G(2K)` at `β`: the residual part decays, since `ρ ≤ 1`.
  have hdec : G (2 * N) ≤ μ ^ a * G (2 * K) := by
    have hd := h.G_decay K N hK hKN
    have hρ : ρ ^ (2 * N - 2 * K) ≤ 1 := pow_le_one₀ h.rho_nonneg h.rho_le_one
    have hD : 0 ≤ G (2 * K) - c * μ ^ (2 * K) := by linarith [h.le_G K hK]
    have h1 : μ ^ (2 * K) * (G (2 * N) - c * μ ^ (2 * N)) ≤
        μ ^ (2 * K) * (μ ^ a * (G (2 * K) - c * μ ^ (2 * K))) := by
      refine hd.trans ?_
      rw [hsplit μ]
      have := mul_le_mul_of_nonneg_right hρ (mul_nonneg (mul_nonneg (pow_nonneg hμ a)
        (pow_nonneg hμ (2 * K))) hD)
      nlinarith [this]
    have h2 := le_of_mul_le_mul_left h1 (pow_pos hμpos _)
    rw [hsplit μ] at h2
    nlinarith [h.c_nonneg]
  have hGn2 : G (2 * N) ^ 2 ≤ (μ ^ a) ^ 2 * G (2 * K) ^ 2 := by
    have hGn : 0 ≤ G (2 * N) :=
      (mul_nonneg h.c_nonneg (pow_nonneg hμ _)).trans (h.le_G N (hK.trans hKN))
    rw [← mul_pow]; exact pow_le_pow_left₀ hGn hdec 2
  -- `G'(2N) ≥ c' μ'^{2N} ≥ (1 - E') μ'^a G'(2K)` at `2β`.
  have hlow : (1 - E') * (μ' ^ a * G' (2 * K)) ≤ G' (2 * N) := by
    have := h'.le_G N (hK.trans hKN)
    rw [hsplit] at this
    nlinarith [mul_le_mul_of_nonneg_left he' (pow_nonneg hμ'.le a)]
  have hbulk := bulk_pow_le hμ' hr hR0 (by omega : a ≤ 2 * K) hR
  have h1E : 0 ≤ 1 - E' := by linarith
  have h1η : 0 ≤ 1 - η := by linarith
  calc (1 - E') * (1 - η) * R * G (2 * N) ^ 2
      ≤ (1 - E') * (1 - η) * R * ((μ ^ a) ^ 2 * G (2 * K) ^ 2) := by gcongr
    _ = (1 - E') * (R * (μ ^ a) ^ 2) * ((1 - η) * G (2 * K) ^ 2) := by ring
    _ ≤ (1 - E') * μ' ^ a * G' (2 * K) := by gcongr
    _ = (1 - E') * (μ' ^ a * G' (2 * K)) := by ring
    _ ≤ G' (2 * N) := hlow

end OpenMomentWitness

/-! ### The boundary seed from finite traces -/

/-- **The boundary defect from finite traces.** The boundary defect `e_{2K}` involves the
leading modulus and weight of the transfer operator, but it is bounded by finitely many traces:
the decay inequality at `2K' = 2K - 2N₀ r` gives `e_{2K} ≤ (ρμ)^{2N₀ r} G(2K') / G(2K)`, and a ring
purity `P_s(2N₀) ≥ 1 - ε` gives `ρ^{2N₀} ≤ d(ε)` and `μ^{2N₀} ≤ Z(2N₀)`. So
`G(2K') (d(ε) Z(2N₀))^r ≤ E G(2K)` implies `e_{2K} ≤ E`. -/
theorem boundaryDefect_le_of_traces {Z G : ℕ → ℝ} (hb : OpenMomentBounds Z G)
    (hZ : ∀ L, 0 < Z L) (hG : ∀ n, 0 < G n) {N₀ K r : ℕ} (hN₀ : 1 ≤ N₀)
    (hK : N₀ * r + 1 ≤ K) {ε E : ℝ} (hε : ε ≤ 1 / 10)
    (hP : 1 - ε ≤ Z (2 * (2 * N₀)) / Z (2 * N₀) ^ 2)
    (hE : G (2 * (K - N₀ * r)) * (tailBound ε * Z (2 * N₀)) ^ r ≤ E * G (2 * K)) :
    boundaryDefect Z G (2 * K) ≤ E := by
  obtain ⟨μ, c, ρ, s, w⟩ := hb
  have hμ := w.mu_pos (hZ 4)
  have hρ1 := w.rho_lt_one hμ hN₀ (hZ _) hε hP
  have htail := w.tail_le hμ hN₀ (hZ _) hε hP
  have hμN := w.pow_le_Z N₀ hN₀
  rw [boundaryDefect, w.leadingRoot_eq (hZ 4), w.leadingOverlap_eq (hZ 4) hρ1]
  set K' := K - N₀ * r with hK'def
  have hK' : 1 ≤ K' := by omega
  have h2 : 2 * N₀ * r = 2 * (N₀ * r) := by ring
  have ha : 2 * K - 2 * K' = 2 * N₀ * r := by omega
  have hsplit : μ ^ (2 * K) = μ ^ (2 * K') * μ ^ (2 * N₀ * r) := by
    rw [← pow_add]; congr 1; omega
  have hρ0 := w.rho_nonneg
  -- The decay at `β`: `G(2K) - cμ^{2K} ≤ (ρμ)^{2N₀ r} G(2K')`.
  have hdec : G (2 * K) - c * μ ^ (2 * K) ≤ (ρ * μ) ^ (2 * N₀ * r) * G (2 * K') := by
    have hd := w.G_decay K' K hK' (by omega)
    rw [ha] at hd
    have hc0 : 0 ≤ c * μ ^ (2 * K') := mul_nonneg w.c_nonneg (pow_nonneg hμ.le _)
    have hρa : 0 ≤ ρ ^ (2 * N₀ * r) := pow_nonneg hρ0 _
    have h1 : μ ^ (2 * K') * (G (2 * K) - c * μ ^ (2 * K)) ≤
        μ ^ (2 * K') * ((ρ * μ) ^ (2 * N₀ * r) * G (2 * K')) := by
      refine hd.trans ?_
      rw [hsplit, mul_pow]
      have := mul_le_mul_of_nonneg_left (by linarith : G (2 * K') - c * μ ^ (2 * K') ≤
        G (2 * K')) (mul_nonneg hρa (mul_nonneg (pow_nonneg hμ.le (2 * K'))
          (pow_nonneg hμ.le (2 * N₀ * r))))
      linarith
    exact le_of_mul_le_mul_left h1 (pow_pos hμ _)
  -- `(ρμ)^{2N₀ r} ≤ (d(ε) Z(2N₀))^r`.
  have hρμ : (ρ * μ) ^ (2 * N₀ * r) ≤ (tailBound ε * Z (2 * N₀)) ^ r := by
    rw [pow_mul, mul_pow]
    exact pow_le_pow_left₀ (mul_nonneg (pow_nonneg hρ0 _) (pow_nonneg hμ.le _))
      (mul_le_mul htail hμN (pow_nonneg hμ.le _) ((pow_nonneg hρ0 _).trans htail)) r
  have hG2K := hG (2 * K)
  rw [sub_le_iff_le_add, ← sub_le_iff_le_add', le_div_iff₀ hG2K]
  nlinarith [mul_le_mul_of_nonneg_right hρμ (hG (2 * K')).le]

/-! ### One step of the boundary induction -/

section Step

variable {ιo : ℕ → Type*} [∀ n, Fintype (ιo n)] [∀ n, Nonempty (ιo n)]
  (Eo : (n : ℕ) → ιo n → ℝ)

/-- **One step of the boundary induction**, eqs. (C3)–(C5) of the note: Moves 1–3 in a row.
At `β` let the boundary defect at length `k = 2K` be at most `E` and the open-chain thermal
purity at length `m = 2N` at least `1 - η`. With the bulk factor `R`, the suppression `q` and
`E' ≥ F(q, E, 1 - η)`, `η' ≥ f(1 - (1 - η) R (1 - E'))`, at `2β` the boundary defect at length
`m` is at most `E'` and the thermal purity at length `2m` is at least `1 - η'`. -/
theorem open_step {β : ℝ} {Z₀ Z₁ : ℕ → ℝ} {μ₀ c₀ ρ₀ s₀ μ₁ c₁ ρ₁ s₁ : ℝ}
    (w₀ : OpenMomentWitness Z₀ (fun n ↦ partitionFunction Eo n β) μ₀ c₀ ρ₀ s₀)
    (w₁ : OpenMomentWitness Z₁ (fun n ↦ partitionFunction Eo n (2 * β)) μ₁ c₁ ρ₁ s₁)
    (hμ₀ : 0 < μ₀) (hμ₁ : 0 < μ₁) (hr : μ₁ ≤ μ₀ ^ 2) {K N : ℕ} (hK : 1 ≤ K) (hKN : K ≤ N)
    {R : ℝ} (hR0 : 0 ≤ R) (hR : R * μ₀ ^ (2 * (2 * N)) ≤ μ₁ ^ (2 * N))
    {q : ℝ} (hq0 : 0 ≤ q) (hq1 : q < 1) (hq : ρ₁ ^ (2 * N - 2 * K) ≤ q)
    {E η E' η' : ℝ} (hE : E < 1)
    (he : (1 - E) * partitionFunction Eo (2 * K) β ≤ c₀ * μ₀ ^ (2 * K))
    (hA : 1 - η ≤ thermalPurity Eo (2 * N) β) (hη : η < 1)
    (hE' : overlapUpdate q E (1 - η) ≤ E') (hE'1 : E' ≤ 1)
    (hprod : 0 < (1 - η) * R * (1 - E'))
    (hη' : squaringDefect (1 - (1 - η) * R * (1 - E')) ≤ η') :
    (1 - E') * partitionFunction Eo (2 * N) (2 * β) ≤ c₁ * μ₁ ^ (2 * N) ∧
      1 - η' ≤ thermalPurity Eo (2 * (2 * N)) (2 * β) := by
  have hG₀ : ∀ n, 0 < partitionFunction Eo n β := fun n ↦ partitionFunction_pos Eo n β
  have hG₁ : ∀ n, 0 < partitionFunction Eo n (2 * β) := fun n ↦ partitionFunction_pos Eo n _
  have hAk : partitionFunction Eo (2 * K) (2 * β) ≤ partitionFunction Eo (2 * K) β ^ 2 := by
    have := thermalPurity_le_one Eo (2 * K) β
    rwa [thermalPurity, div_le_one (pow_pos (hG₀ _) 2)] at this
  have hAm : (1 - η) * partitionFunction Eo (2 * N) β ^ 2 ≤ partitionFunction Eo (2 * N) (2 * β)
      := by
    rw [thermalPurity, le_div_iff₀ (pow_pos (hG₀ _) 2)] at hA; linarith
  -- Move 1.
  have h1 := w₀.move_one w₁ hμ₁ hr hK hKN hAk (by linarith) hAm hE he hq0 hq1 hq
  have hfirst : (1 - E') * partitionFunction Eo (2 * N) (2 * β) ≤ c₁ * μ₁ ^ (2 * N) :=
    le_trans (mul_le_mul_of_nonneg_right (by linarith) (hG₁ _).le) h1
  refine ⟨hfirst, ?_⟩
  -- Move 2: `A(4N, β) ≥ (1 - η) R (1 - E')`.
  have h2 := w₀.move_two w₁ hμ₀ (M := N) (hK.trans hKN) hE'1 hfirst hR0 hR
  have hA4 : 1 - (1 - (1 - η) * R * (1 - E')) ≤ thermalPurity Eo (2 * (2 * N)) β := by
    rw [sub_sub_cancel, thermalPurity, le_div_iff₀ (pow_pos (hG₀ _) 2)]
    have hRE : 0 ≤ R * (1 - E') := mul_nonneg hR0 (by linarith)
    have : (1 - η) * R * (1 - E') * partitionFunction Eo (2 * (2 * N)) β ^ 2 *
        partitionFunction Eo (2 * N) β ^ 2 ≤
          partitionFunction Eo (2 * (2 * N)) (2 * β) * partitionFunction Eo (2 * N) β ^ 2 := by
      calc _ = R * (1 - E') * ((1 - η) * partitionFunction Eo (2 * N) β ^ 2) *
            partitionFunction Eo (2 * (2 * N)) β ^ 2 := by ring
        _ ≤ R * (1 - E') * partitionFunction Eo (2 * N) (2 * β) *
            partitionFunction Eo (2 * (2 * N)) β ^ 2 :=
          mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hAm hRE) (sq_nonneg _)
        _ ≤ _ := h2
    exact le_of_mul_le_mul_right this (pow_pos (hG₀ _) 2)
  -- Move 3: squaring the Gibbs state.
  have h3 := one_sub_thermalPurity_two_mul_le Eo (2 * (2 * N)) β
  have hmono := squaringDefect_le_of_le (thermalPurity_le_one Eo _ β) hA4 (by linarith)
  linarith

end Step

/-! ### Scalar inequalities of the uniform phase (the note, eq. (C6)) -/

/-- `d(6x²) ≤ (2x)²` for `0 ≤ x ≤ 1/100`. -/
theorem uniform_tail {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 1 / 100) :
    tailBound (6 * x ^ 2) ≤ (2 * x) ^ 2 := by
  unfold tailBound
  rw [div_le_iff₀ (by nlinarith)]
  have hx2 : x ^ 2 ≤ 1 / 36 := by nlinarith
  nlinarith [mul_le_mul_of_nonneg_left hx2 (sq_nonneg x)]

/-- `F(2x, x, 1 - x) ≤ 10 x²` for `0 ≤ x ≤ 1/100`. -/
theorem uniform_overlap {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 1 / 100) :
    overlapUpdate (2 * x) x (1 - x) ≤ 10 * x ^ 2 := by
  set u := (1 - x) ^ 2 * (1 - x) with hu
  have hu0 : 97 / 100 ≤ u := by rw [hu]; nlinarith [mul_nonneg hx0 hx0]
  have hu1 : 1 - u ≤ 3 * x := by
    rw [hu]; nlinarith [mul_nonneg hx0 hx0, mul_nonneg (mul_nonneg hx0 hx0) hx0]
  have h2x : 0 < 1 - 2 * x := by linarith
  have e : overlapUpdate (2 * x) x (1 - x) = 2 * x * (1 - u) / ((1 - 2 * x) * u) := by
    unfold overlapUpdate
    rw [← hu]
    field_simp
  rw [e, div_le_iff₀ (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hu1 (by positivity : (0 : ℝ) ≤ 2 * x),
    mul_nonneg (mul_nonneg hx0 hx0) (by linarith : (0 : ℝ) ≤ u - 97 / 100)]

/-- `f(2x + 16x²) ≤ 10 x²` for `0 ≤ x ≤ 1/100`. -/
theorem uniform_squaring {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 1 / 100) :
    squaringDefect (2 * x + 16 * x ^ 2) ≤ 10 * x ^ 2 := by
  unfold squaringDefect
  have h1 : 9 / 10 ≤ 1 - (2 * x + 16 * x ^ 2) := by nlinarith
  rw [div_le_iff₀ (by nlinarith)]
  have h2 : (2 * x + 16 * x ^ 2) ^ 2 ≤ (54 / 25) ^ 2 * x ^ 2 := by
    have : 2 * x + 16 * x ^ 2 ≤ 54 / 25 * x := by nlinarith
    calc (2 * x + 16 * x ^ 2) ^ 2 ≤ (54 / 25 * x) ^ 2 := pow_le_pow_left₀ (by positivity) this 2
      _ = (54 / 25) ^ 2 * x ^ 2 := by ring
  nlinarith [pow_le_pow_left₀ (by norm_num) h1 2, sq_nonneg x]

/-- The four-factor product of eq. (C5) in the uniform phase: its defect is at most
`2x + 16x²`. -/
theorem uniform_product {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 1 / 100) :
    1 - (1 - x) * ((1 - x) * (1 - 6 * x ^ 2)) * (1 - 10 * x ^ 2) ≤ 2 * x + 16 * x ^ 2 := by
  nlinarith [mul_nonneg hx0 hx0, mul_nonneg (mul_nonneg hx0 hx0) hx0,
    mul_nonneg (mul_nonneg hx0 hx0) (mul_nonneg hx0 hx0),
    mul_nonneg (mul_nonneg (mul_nonneg hx0 hx0) (mul_nonneg hx0 hx0)) hx0,
    mul_nonneg (mul_nonneg (mul_nonneg hx0 hx0) (mul_nonneg hx0 hx0)) (mul_nonneg hx0 hx0)]

/-- The final equal-field purity: `(1 - 20x²)(1 - x)(1 - x)(1 - 6x²) ≥ 1 - 4x`. -/
theorem uniform_final {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 1 / 100) :
    1 - 4 * x ≤ (1 - 2 * (10 * x ^ 2)) * (1 - x) * ((1 - x) * (1 - 6 * x ^ 2)) := by
  nlinarith [mul_nonneg hx0 hx0, mul_nonneg (mul_nonneg hx0 hx0) hx0,
    mul_nonneg (mul_nonneg hx0 hx0) (mul_nonneg hx0 hx0),
    mul_nonneg (mul_nonneg (mul_nonneg hx0 hx0) (mul_nonneg hx0 hx0)) hx0]

/-- From `ρ^{2m} ≤ d ≤ q²` to `ρ^m ≤ q`. -/
theorem pow_le_of_pow_two_mul_le {ρ q d : ℝ} (hρ : 0 ≤ ρ) (hq : 0 ≤ q) {m : ℕ}
    (h : ρ ^ (2 * m) ≤ d) (hd : d ≤ q ^ 2) : ρ ^ m ≤ q := by
  refine (pow_le_pow_iff_left₀ (pow_nonneg hρ m) hq two_ne_zero).1 ?_
  rw [← pow_mul, mul_comm]
  exact h.trans hd

/-- From `ρ^b ≤ d` and `d^{a/g} ≤ q^{b/g}` with `g = gcd(a, b)` to `ρ^a ≤ q`: the exponents are
reduced by their common divisor before the rational comparison. -/
theorem pow_le_of_pow_le_gcd {ρ d q : ℝ} (hρ : 0 ≤ ρ) (hq : 0 ≤ q) {a b : ℕ} (hb : 0 < b)
    (h : ρ ^ b ≤ d) (hd : d ^ (a / a.gcd b) ≤ q ^ (b / a.gcd b)) : ρ ^ a ≤ q := by
  have hg : 0 < a.gcd b := Nat.gcd_pos_of_pos_right a hb
  have hbg : b / a.gcd b ≠ 0 := (Nat.div_pos (Nat.gcd_le_right a hb) hg).ne'
  refine (pow_le_pow_iff_left₀ (pow_nonneg hρ a) hq hbg).1 ?_
  have e : a * (b / a.gcd b) = b * (a / a.gcd b) := by
    rw [← Nat.mul_div_assoc a (Nat.gcd_dvd_right a b),
      ← Nat.mul_div_assoc b (Nat.gcd_dvd_left a b), mul_comm]
  calc (ρ ^ a) ^ (b / a.gcd b) = (ρ ^ b) ^ (a / a.gcd b) := by rw [← pow_mul, ← pow_mul, e]
    _ ≤ d ^ (a / a.gcd b) := pow_le_pow_left₀ (pow_nonneg hρ b) h _
    _ ≤ q ^ (b / a.gcd b) := hd

/-! ### The boundary-field gap -/

/-- The half-length `Kᵢ` of the boundary defect at scale `i`: the seed's `K₀` at `i = 0`, and the
ring length `N₀ 2^i` of that scale afterwards. -/
def boundaryHalf (N₀ K₀ i : ℕ) : ℕ := if i = 0 then K₀ else N₀ * 2 ^ i

/-- The scalar conditions of one base step of the boundary induction (the note, eqs.
(C3)–(C5)). The ring defects are `ε` and `ε'` at the two scales, the suppression `q` bounds
`ρ'^a` through `d(ε')^{a/g} ≤ q^{b/g}`, and `(E, η) ↦ (E', η')` follows the overlap update and the
squaring bound. -/
structure OpenStep (ε ε' q E η E' η' : ℝ) (a b : ℕ) : Prop where
  q_nonneg : 0 ≤ q
  q_lt_one : q < 1
  tail : tailBound ε' ^ (a / a.gcd b) ≤ q ^ (b / a.gcd b)
  E_lt_one : E < 1
  η_lt_one : η < 1
  overlap : overlapUpdate q E (1 - η) ≤ E'
  E'_le_one : E' ≤ 1
  prod_pos : 0 < (1 - η) * ((1 - ε) * (1 - ε')) * (1 - E')
  squaring : squaringDefect (1 - (1 - η) * ((1 - ε) * (1 - ε')) * (1 - E')) ≤ η'

section Main

variable {ιr ιo : ℕ → Type*} [∀ L, Fintype (ιr L)] [∀ L, Nonempty (ιr L)]
  [∀ n, Fintype (ιo n)] [∀ n, Nonempty (ιo n)] (Er : (L : ℕ) → ιr L → ℝ)
  (Eo : (n : ℕ) → ιo n → ℝ)

/-- **The boundary-field gap criterion** (Theorem 2 of the note, with general seeds). Let the ring
spectra `Er` satisfy the spatial moment inequalities and the ring and open-chain traces the
boundary moment inequalities at every `β > 0`. The scales are `βᵢ = β₀ 2^i` and
`Nᵢ = N₀ 2^i`, with ring lengths `2Nᵢ` and open-chain lengths `Mᵢ = 4Nᵢ`. Assume

* the ring seeds `P_s(2N₀, β₀), P_τ(4N₀, β₀) ≥ 1 - u₀`, with `G(u₀) ≤ C ≤ 6` and `C u₀ ≤ 1`, so
  that the ring defects are `εᵢ = (C u₀)^{2^i} / C`;
* the boundary seed `e_{2K₀}(β₀) ≤ E₀` and the open-chain seed `P_τ(M₀, β₀) ≥ 1 - η₀`;
* `m ≥ 1` base steps `OpenStep` from `(Eᵢ, ηᵢ)` to `(Eᵢ₊₁, ηᵢ₊₁)`, with boundary length
  `2Kᵢ` (`boundaryHalf`);
* the hand-over `Eₘ, ηₘ, εₘ ≤ x ≤ 1/100`.

Then every open chain of `n ≥ 4 N₀ 2^m` sites, odd or even, has a nondegenerate ground level
and gap `Δ_n > log(1/(10x)) / (β₀ 2^m)`. -/
theorem gapGT_of_openChain {N₀ K₀ m : ℕ} (hN₀ : 1 ≤ N₀) (hK₀ : 1 ≤ K₀) (hK₀N : K₀ ≤ 2 * N₀)
    (hm : 1 ≤ m) (hring : HasSpatialMomentBounds Er (2 * N₀))
    (hopen : ∀ β : ℝ, 0 < β → OpenMomentBounds (fun L ↦ partitionFunction Er L β)
      (fun n ↦ partitionFunction Eo n β))
    {β₀ u₀ C x : ℝ} (hβ₀ : 0 < β₀) (hu₀ : 0 < u₀) (hu₁ : u₀ ≤ 1 / 10)
    (hCoef : doublingCoefficient u₀ ≤ C) (hC : C ≤ 6) (hCu : C * u₀ ≤ 1) (E η q : ℕ → ℝ)
    (hstep : ∀ i < m, OpenStep (coefficientDefect C u₀ i) (coefficientDefect C u₀ (i + 1))
      (q i) (E i) (η i) (E (i + 1)) (η (i + 1))
      (2 * (2 * (N₀ * 2 ^ i)) - 2 * boundaryHalf N₀ K₀ i) (2 * (2 * (N₀ * 2 ^ i))))
    (hx0 : 0 < x) (hx : x ≤ 1 / 100) (hEm : E m ≤ x) (hηm : η m ≤ x)
    (hεm : coefficientDefect C u₀ m ≤ x)
    (hs : 1 - u₀ ≤ spatialPurity Er (2 * N₀) β₀)
    (ht : 1 - u₀ ≤ thermalPurity Er (2 * (2 * N₀)) β₀)
    (he : boundaryDefect (fun L ↦ partitionFunction Er L β₀) (fun n ↦ partitionFunction Eo n β₀)
      (2 * K₀) ≤ E 0)
    (hη : 1 - η 0 ≤ thermalPurity Eo (2 * (2 * N₀)) β₀) {n : ℕ} (hn : 4 * N₀ * 2 ^ m ≤ n) :
    GapGT (Eo n) (log (1 / (10 * x)) / (β₀ * 2 ^ m)) := by
  set β : ℕ → ℝ := fun i ↦ β₀ * 2 ^ i with hβdef
  have hβpos : ∀ i, 0 < β i := fun i ↦ by simp only [β]; positivity
  have hβs : ∀ i, β (i + 1) = 2 * β i := fun i ↦ by simp only [β]; ring
  have hβ0 : β 0 = β₀ := by simp [β]
  set Nn : ℕ → ℕ := fun i ↦ N₀ * 2 ^ i with hNndef
  have hNn1 : ∀ i, 1 ≤ Nn i := fun i ↦ by
    simp only [Nn]; exact Nat.one_le_iff_ne_zero.2 (by positivity)
  have hNns : ∀ i, Nn (i + 1) = 2 * Nn i := fun i ↦ by simp only [Nn]; ring
  have hNn0 : Nn 0 = N₀ := by simp [Nn]
  set ε : ℕ → ℝ := coefficientDefect C u₀ with hεdef
  have hCpos : 0 < C := (doublingCoefficient_pos (by linarith)).trans_le hCoef
  have hε0 : ∀ i, 0 ≤ ε i := fun i ↦ (coefficientDefect_mem hCpos hu₀ hCu i).1.le
  have hεu : ∀ i, ε i ≤ u₀ := fun i ↦ (coefficientDefect_mem hCpos hu₀ hCu i).2
  have hε10 : ∀ i, ε i ≤ 1 / 10 := fun i ↦ (hεu i).trans hu₁
  have hεs : ∀ i, ε (i + 1) = C * ε i ^ 2 := coefficientDefect_succ hCpos.ne' u₀
  -- The witnesses at every scale.
  choose μ c ρ s hw using fun i ↦ hopen (β i) (hβpos i)
  have hw' : ∀ i, OpenMomentWitness (fun L ↦ partitionFunction Er L (2 * β i))
      (fun n ↦ partitionFunction Eo n (2 * β i)) (μ (i + 1)) (c (i + 1)) (ρ (i + 1))
      (s (i + 1)) := fun i ↦ by
    have := hw (i + 1); rwa [hβs] at this
  have hμ : ∀ i, 0 < μ i := fun i ↦ (hw i).mu_pos (partitionFunction_pos Er 4 _)
  have hZ : ∀ L β', 0 < partitionFunction Er L β' := partitionFunction_pos Er
  have hG : ∀ m β', 0 < partitionFunction Eo m β' := partitionFunction_pos Eo
  /- The ring: the purities at every scale, the residual ratios `ρᵢ^{2Nᵢ} ≤ d(εᵢ)`, the growth
  bound `μᵢ₊₁ ≤ μᵢ²` and the bulk factor `Rᵢ = (1 - εᵢ)(1 - εᵢ₊₁)`. -/
  have hR := purities_ge_coefficientDefect Er (n₀ := 2 * N₀) ⟨N₀, by ring⟩ hβ₀ hring hu₀
    (by linarith) hCoef hCu hs ht
  have hRs : ∀ i, 1 - ε i ≤ partitionFunction Er (2 * (2 * Nn i)) (β i) /
      partitionFunction Er (2 * Nn i) (β i) ^ 2 := fun i ↦ by
    have := (hR i).1
    rwa [spatialPurity, show 2 * N₀ * 2 ^ i = 2 * Nn i by simp only [Nn]; ring] at this
  have hRt : ∀ i, (1 - ε i) * partitionFunction Er (2 * (2 * Nn i)) (β i) ^ 2 ≤
      partitionFunction Er (2 * (2 * Nn i)) (2 * β i) := fun i ↦ by
    have := (hR i).2
    rwa [thermalPurity, show 2 * N₀ * 2 ^ i = 2 * Nn i by simp only [Nn]; ring,
      le_div_iff₀ (pow_pos (hZ _ _) 2)] at this
  have htail : ∀ i, ρ i ^ (2 * Nn i) ≤ tailBound (ε i) := fun i ↦
    (hw i).tail_le (hμ i) (hNn1 i) (hZ _ _) (hε10 i) (hRs i)
  have hr : ∀ i, μ (i + 1) ≤ μ i ^ 2 := fun i ↦
    (hw i).mu_le_sq (hw' i) (hμ i) fun N _ ↦ by
      have := thermalPurity_le_one Er (2 * N) (β i)
      rwa [thermalPurity, div_le_one (pow_pos (hZ _ _) 2)] at this
  have hRbulk : ∀ i, (1 - ε i) * (1 - ε (i + 1)) * μ i ^ (2 * (2 * (2 * Nn i))) ≤
      μ (i + 1) ^ (2 * (2 * Nn i)) := fun i ↦ by
    have hP' := hRs (i + 1)
    rw [hNns, hβs] at hP'
    exact (hw i).mu_pow_ge (hw' i) (N := 2 * Nn i) (by have := hNn1 i; omega) (hZ _ _)
      ((hε10 i).trans (by norm_num)) ((hε10 _).trans (by norm_num)) (hRt i) hP'
  have hR0 : ∀ i, 0 ≤ (1 - ε i) * (1 - ε (i + 1)) := fun i ↦
    mul_nonneg (by linarith [hε10 i]) (by linarith [hε10 (i + 1)])
  -- The boundary lengths `2Kᵢ`.
  have hK1 : ∀ i, 1 ≤ boundaryHalf N₀ K₀ i := fun i ↦ by
    unfold boundaryHalf; split_ifs
    · exact hK₀
    · exact hNn1 i
  have hKN : ∀ i, boundaryHalf N₀ K₀ i ≤ 2 * Nn i := fun i ↦ by
    unfold boundaryHalf; split_ifs with h
    · subst h; rw [hNn0]; exact hK₀N
    · have := hNn1 i; simp only [Nn] at this ⊢; omega
  have hKs : ∀ i, i ≠ 0 → boundaryHalf N₀ K₀ i = Nn i := fun i hi ↦ by
    simp only [boundaryHalf, hi, ite_false, Nn]
  /- One step from scale `i` to `i + 1`: the boundary defect at length `2Kᵢ` and the open-chain
  thermal purity at length `4Nᵢ` at `βᵢ` give both at `βᵢ₊₁`, at lengths `2Nᵢ₊₁` and `4Nᵢ₊₁`. -/
  have step : ∀ i {q E η E' η' : ℝ}, 0 ≤ q → q < 1 →
      ρ (i + 1) ^ (2 * (2 * Nn i) - 2 * boundaryHalf N₀ K₀ i) ≤ q → E < 1 → η < 1 →
      overlapUpdate q E (1 - η) ≤ E' → E' ≤ 1 →
      0 < (1 - η) * ((1 - ε i) * (1 - ε (i + 1))) * (1 - E') →
      squaringDefect (1 - (1 - η) * ((1 - ε i) * (1 - ε (i + 1))) * (1 - E')) ≤ η' →
      (1 - E) * partitionFunction Eo (2 * boundaryHalf N₀ K₀ i) (β i) ≤
        c i * μ i ^ (2 * boundaryHalf N₀ K₀ i) →
      1 - η ≤ thermalPurity Eo (2 * (2 * Nn i)) (β i) →
      (1 - E') * partitionFunction Eo (2 * Nn (i + 1)) (β (i + 1)) ≤
          c (i + 1) * μ (i + 1) ^ (2 * Nn (i + 1)) ∧
        1 - η' ≤ thermalPurity Eo (2 * (2 * Nn (i + 1))) (β (i + 1)) := by
    intro i q E η E' η' hq0 hq1 hq hE hη hE' hE'1 hprod hη' he hA
    have h := open_step Eo (hw i) (hw' i) (hμ i) (hμ _) (hr i) (K := boundaryHalf N₀ K₀ i)
      (N := 2 * Nn i) (hK1 i) (hKN i) (hR0 i) (hRbulk i) hq0 hq1 hq hE he hA hη hE' hE'1 hprod
      hη'
    rw [hNns, hβs]
    exact h
  -- The base phase: the seeds and `m` certified steps.
  have hbase : ∀ i ≤ m, (1 - E i) * partitionFunction Eo (2 * boundaryHalf N₀ K₀ i) (β i) ≤
      c i * μ i ^ (2 * boundaryHalf N₀ K₀ i) ∧
        1 - η i ≤ thermalPurity Eo (2 * (2 * Nn i)) (β i) := by
    intro i hi
    induction i with
    | zero =>
      have w0 := hw 0
      rw [hβ0] at w0
      have hρ0 : ρ 0 < 1 := by
        have := hRs 0
        rw [hβ0] at this
        exact w0.rho_lt_one (hμ 0) (hNn1 0) (hZ _ _) (hε10 0) this
      rw [boundaryDefect, w0.leadingRoot_eq (hZ 4 β₀), w0.leadingOverlap_eq (hZ 4 β₀) hρ0] at he
      have hK0 : boundaryHalf N₀ K₀ 0 = K₀ := rfl
      rw [hK0, hNn0, hβ0]
      refine ⟨?_, hη⟩
      have : 1 - E 0 ≤ c 0 * μ 0 ^ (2 * K₀) / partitionFunction Eo (2 * K₀) β₀ := by linarith
      rw [le_div_iff₀ (hG _ _)] at this
      exact this
    | succ i ih =>
      obtain ⟨he, hA⟩ := ih (by omega)
      have hi' := hstep i (by omega)
      have hq : ρ (i + 1) ^ (2 * (2 * Nn i) - 2 * boundaryHalf N₀ K₀ i) ≤ q i := by
        have := htail (i + 1)
        rw [hNns] at this
        exact pow_le_of_pow_le_gcd (hw _).rho_nonneg hi'.q_nonneg
          (by have := hNn1 i; omega) this hi'.tail
      have h := step i hi'.q_nonneg hi'.q_lt_one hq hi'.E_lt_one hi'.η_lt_one hi'.overlap
        hi'.E'_le_one hi'.prod_pos hi'.squaring he hA
      rw [hKs (i + 1) (Nat.succ_ne_zero i)]
      exact h
  /- The uniform phase: with `x_j = (10x)^{2^j} / 10`, all three defects at scale `m + j` are at
  most `x_j`, and `x_{j+1} = 10 x_j²`. -/
  set xs : ℕ → ℝ := coefficientDefect 10 x with hxsdef
  have hxs0 : ∀ j, 0 < xs j := fun j ↦
    (coefficientDefect_mem (by norm_num) hx0 (by linarith) j).1
  have hxsle : ∀ j, xs j ≤ 1 / 100 := fun j ↦
    (coefficientDefect_mem (by norm_num) hx0 (by linarith) j).2.trans hx
  have hxss : ∀ j, xs (j + 1) = 10 * xs j ^ 2 := coefficientDefect_succ (by norm_num) x
  have hεsq : ∀ i {y : ℝ}, ε i ≤ y → ε (i + 1) ≤ 6 * y ^ 2 := fun i y hy ↦ by
    rw [hεs]
    have := pow_le_pow_left₀ (hε0 i) hy 2
    nlinarith [mul_le_mul_of_nonneg_right hC (sq_nonneg (ε i))]
  have hunif : ∀ j, ((1 - xs j) * partitionFunction Eo (2 * Nn (m + j)) (β (m + j)) ≤
      c (m + j) * μ (m + j) ^ (2 * Nn (m + j)) ∧
        1 - xs j ≤ thermalPurity Eo (2 * (2 * Nn (m + j))) (β (m + j))) ∧ ε (m + j) ≤ xs j := by
    intro j
    induction j with
    | zero =>
      have hx' : xs 0 = x := by simp [xs, coefficientDefect]
      obtain ⟨he, hA⟩ := hbase m le_rfl
      rw [hKs m (by omega)] at he
      simp only [add_zero, hx']
      exact ⟨⟨le_trans (mul_le_mul_of_nonneg_right (by linarith) (hG _ _).le) he,
        le_trans (by linarith) hA⟩, hεm⟩
    | succ j ih =>
      obtain ⟨⟨he, hA⟩, hεj⟩ := ih
      set y := xs j
      have hy0 := (hxs0 j).le
      have hy := hxsle j
      have hysq : y ^ 2 ≤ 1 / 10000 := by
        have := pow_le_pow_left₀ hy0 hy 2; norm_num at this ⊢; linarith
      have hysq0 : 0 ≤ y ^ 2 := sq_nonneg y
      have hεs' := hεsq (m + j) hεj
      set Ri := (1 - ε (m + j)) * (1 - ε (m + j + 1)) with hRi
      have hRlow : (1 - y) * (1 - 6 * y ^ 2) ≤ Ri :=
        mul_le_mul (by linarith) (by linarith) (by linarith) (by linarith [hε10 (m + j)])
      have hRle : Ri ≤ 1 := (mul_le_mul (by linarith [hε0 (m + j)])
        (by linarith [hε0 (m + j + 1)]) (by linarith [hε10 (m + j + 1)]) zero_le_one).trans_eq
        (mul_one 1)
      have h10 : 0 ≤ 1 - 10 * y ^ 2 := by linarith
      have hlow : (1 - y) * ((1 - y) * (1 - 6 * y ^ 2)) * (1 - 10 * y ^ 2) ≤
          (1 - y) * Ri * (1 - 10 * y ^ 2) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hRlow (by linarith)) h10
      have hpos : 0 < (1 - y) * ((1 - y) * (1 - 6 * y ^ 2)) * (1 - 10 * y ^ 2) :=
        mul_pos (mul_pos (by linarith) (mul_pos (by linarith) (by linarith))) (by linarith)
      have hle1 : (1 - y) * Ri * (1 - 10 * y ^ 2) ≤ 1 :=
        (mul_le_mul ((mul_le_mul (by linarith) hRle (hR0 _) zero_le_one).trans_eq (mul_one 1))
          (by linarith) h10 zero_le_one).trans_eq (mul_one 1)
      have hKj := hKs (m + j) (by omega)
      have hq : ρ (m + j + 1) ^ (2 * (2 * Nn (m + j)) - 2 * boundaryHalf N₀ K₀ (m + j)) ≤
          2 * y := by
        rw [hKj, show 2 * (2 * Nn (m + j)) - 2 * Nn (m + j) = 2 * Nn (m + j) by omega]
        have := htail (m + j + 1)
        rw [hNns] at this
        exact pow_le_of_pow_two_mul_le (hw _).rho_nonneg (by linarith) this
          ((tailBound_mono hεs' (by linarith)).trans (uniform_tail hy0 hy))
      rw [← hKj] at he
      have h := step (m + j) (q := 2 * y) (E := y) (η := y) (E' := 10 * y ^ 2)
        (η' := 10 * y ^ 2) (by linarith) (by linarith) hq (by linarith) (by linarith)
        (uniform_overlap hy0 hy) (by linarith) (hpos.trans_le hlow)
        (le_trans (squaringDefect_mono (by linarith)
          (by linarith [uniform_product hy0 hy]) (by linarith)) (uniform_squaring hy0 hy))
        he hA
      rw [hxss]
      exact ⟨h, hεs'.trans (by linarith)⟩
  -- The final scale: `Mᵢ ≤ n ≤ 2 Mᵢ` with `i = m + j`, at either parity of `n`.
  obtain ⟨j, hj1, hj2⟩ := exists_dyadic_window (n₀ := 2 * (2 * Nn m))
    (by have := hNn1 m; omega) ((le_of_eq (by simp only [Nn]; ring)).trans hn)
  have hM : 2 * (2 * Nn m) * 2 ^ j = 2 * (2 * Nn (m + j)) := by
    simp only [Nn]; rw [pow_add]; ring
  rw [hM] at hj1 hj2
  obtain ⟨⟨-, hA⟩, hεj⟩ := hunif j
  obtain ⟨⟨he', -⟩, -⟩ := hunif (j + 1)
  set y := xs j
  have hy0 := (hxs0 j).le
  have hy := hxsle j
  have hysq : y ^ 2 ≤ 1 / 10000 := by
    have := pow_le_pow_left₀ hy0 hy 2; norm_num at this ⊢; linarith
  have hysq0 : 0 ≤ y ^ 2 := sq_nonneg y
  rw [show m + (j + 1) = m + j + 1 by omega, hNns, hβs] at he'
  rw [thermalPurity, le_div_iff₀ (pow_pos (hG _ _) 2)] at hA
  -- `(1 - 2x')(1 - x) Rᵢ ≥ 1 - 4x`, with `x' = 10x²` the boundary defect at `2βᵢ`.
  have hlow : 1 - 4 * y ≤
      (1 - 2 * xs (j + 1)) * (1 - y) * ((1 - ε (m + j)) * (1 - ε (m + j + 1))) := by
    have hεs' := hεsq (m + j) hεj
    have hRlow : (1 - y) * (1 - 6 * y ^ 2) ≤ (1 - ε (m + j)) * (1 - ε (m + j + 1)) :=
      mul_le_mul (by linarith) (by linarith) (by linarith) (by linarith [hε10 (m + j)])
    rw [hxss]
    exact (uniform_final hy0 hy).trans
      (mul_le_mul_of_nonneg_left hRlow (mul_nonneg (by linarith) (by linarith)))
  have hK1 : 1 ≤ 2 * Nn (m + j) := by have := hNn1 (m + j); omega
  -- `P_τ(n, βᵢ) ≥ 1 - 4x` at both parities.
  have hP : 1 - 4 * y ≤ thermalPurity Eo n (β (m + j)) := by
    rw [thermalPurity, le_div_iff₀ (pow_pos (hG _ _) 2)]
    have hG2 := pow_pos (hG n (β (m + j))) 2
    obtain ⟨N, hN | hN⟩ := Nat.even_or_odd' n <;> subst hN
    · have hpur := (hw (m + j)).evenChain_purity (hw' (m + j)) (hμ _) (hr _)
        (K := 2 * Nn (m + j)) (N := N) hK1 (by omega) (by omega) (E' := xs (j + 1))
        (by rw [hxss]; linarith) he' (η := y) (by linarith) hA (hR0 _) (hRbulk _)
      have hle : (1 - 2 * xs (j + 1)) * (1 - y) * ((1 - ε (m + j)) * (1 - ε (m + j + 1))) ≤
          (1 - xs (j + 1)) * (1 - y) * ((1 - ε (m + j)) * (1 - ε (m + j + 1))) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right
          (by linarith [hxs0 (j + 1)]) (by linarith)) (hR0 _)
      exact (mul_le_mul_of_nonneg_right (hlow.trans hle) hG2.le).trans hpur
    · have hpur := (hw (m + j)).oddChain_purity (hw' (m + j)) (hμ _) (hr _)
        (K := 2 * Nn (m + j)) (N := N) hK1 (by omega) (by omega) (hG _ _)
        (E' := xs (j + 1)) (by rw [hxss]; linarith) he' (η := y) (by linarith) hA (hR0 _)
        (hRbulk _)
      exact (mul_le_mul_of_nonneg_right hlow hG2.le).trans hpur
  /- From the purity to the gap:
  `e^{β_{m+j} log(1/(10x)) / β_m} · 4x_j = (10x)^{-2^j} · 4 (10x)^{2^j} / 10`. -/
  refine gapGT_of_thermalPurity Eo n (hβpos _) (by linarith) hP ?_
  have h10x : 0 < 10 * x := by positivity
  have hexp : exp (β (m + j) * (log (1 / (10 * x)) / (β₀ * 2 ^ m))) =
      (1 / (10 * x)) ^ (2 ^ j) := by
    have : β (m + j) * (log (1 / (10 * x)) / (β₀ * 2 ^ m)) =
        ((2 ^ j : ℕ) : ℝ) * log (1 / (10 * x)) := by
      simp only [β]; push_cast; rw [pow_add]; field_simp
    rw [this, exp_nat_mul, exp_log (by positivity)]
  have hx10 : (1 / (10 * x)) ^ (2 ^ j) * y = 1 / 10 := by
    simp only [y, xs, coefficientDefect]
    rw [mul_div_assoc', ← mul_pow, one_div_mul_cancel h10x.ne', one_pow]
  rw [hexp, show (1 / (10 * x)) ^ 2 ^ j * (4 * y) = 4 * ((1 / (10 * x)) ^ 2 ^ j * y) by ring,
    hx10]
  linarith

end Main

end CertifiedGaps
