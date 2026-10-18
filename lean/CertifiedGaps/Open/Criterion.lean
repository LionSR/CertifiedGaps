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

The ring spectra are `Er L`, the open-chain spectra `Eo n`. The scales are
`βᵢ = 49 · 2^i`, `nᵢ = 120 · 2^i` (ring lengths) and `Mᵢ = 2 nᵢ` (open-chain lengths); the
note's scale `j` is `i + 2`. At each `βᵢ` we fix a witness `(μᵢ, cᵢ, ρᵢ, sᵢ)` of the boundary
moment inequalities. The three claims of the note are, in multiplicative form,
* (R) ring purities `P_s(nᵢ, βᵢ), P_τ(Mᵢ, βᵢ) ≥ 1 - εᵢ` with `ε₀ = 0.052`, `εᵢ₊₁ = 6 εᵢ²`;
* (E) the boundary defect `e_{kᵢ}(i) ≤ Eᵢ`, written `(1 - Eᵢ) G(kᵢ) ≤ cᵢ μᵢ^{kᵢ}`, with
  `k₀ = 142` and `kᵢ₊₁ = Mᵢ`;
* (T) the open-chain thermal purity `P_τ(Mᵢ, βᵢ) ≥ 1 - ηᵢ`.

## Main results

* `ring_doubling_six`: the ring doubling step with prefactor `6` for defects `≤ 0.052`.
* `OpenMomentWitness.tail_le`, `mu_le_sq`, `mu_pow_ge`: claim (R) as transfer bounds.
* `move_one`, `move_two`, `oddChain_purity`: the overlap recurrence, the length doubling and
  the equal-field purity of Appendix A, all in multiplicative form.
* `gapGT_of_openSeeds`: the boundary-field gap `Δ_n > (log 10)/392` for every odd `n ≥ 1921`,
  from the four inputs of the note's eq. (C1).
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

/-- The coefficient of the ring doubling at defect `ε ≤ 0.052`, the note's eq. (B3):
`G(ε) ≤ 6`, so `f(1 - (1 - ε)³) ≤ 6 ε²`. -/
theorem doublingCoefficient_le_six {ε : ℝ} (h1 : ε ≤ 13 / 250) : doublingCoefficient ε ≤ 6 :=
  (doublingCoefficient_mono h1 (by norm_num)).trans (by unfold doublingCoefficient; norm_num)

/-! ### Claim (R): the ring at the open-chain scales -/

section Ring

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)] (E : (L : ℕ) → ι L → ℝ)

/-- **Ring doubling with prefactor 6** (Appendix B). For a defect `ε ≤ 0.052` (above the
threshold `1/100` of `doubling`), `P_s(n, β), P_τ(2n, β) ≥ 1 - ε` give
`P_s(2n, 2β), P_τ(4n, 2β) ≥ 1 - 6ε²`. -/
theorem ring_doubling_six {n₀ : ℕ} (hZ : HasSpatialMomentBounds E n₀) {n : ℕ} (hn : n₀ ≤ n)
    (hne : Even n) {β : ℝ} (hβ : 0 < β) {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 13 / 250)
    (hs : 1 - ε ≤ spatialPurity E n β) (ht : 1 - ε ≤ thermalPurity E (2 * n) β) :
    1 - 6 * ε ^ 2 ≤ spatialPurity E (2 * n) (2 * β) ∧
      1 - 6 * ε ^ 2 ≤ thermalPurity E (4 * n) (2 * β) :=
  doubling_of_coefficient E hZ hn hne hβ h0 (by linarith) (doublingCoefficient_le_six h1)
    (by linarith) hs ht

end Ring

/-- The ring defects `ε₀ = 0.052`, `εᵢ₊₁ = 6 εᵢ²` (the note's `ε_{i+2}`). -/
noncomputable def ringDefect : ℕ → ℝ
  | 0 => 13 / 250
  | i + 1 => 6 * ringDefect i ^ 2

/-- `εᵢ ≥ 0`. -/
theorem ringDefect_nonneg (i : ℕ) : 0 ≤ ringDefect i := by
  cases i <;> simp only [ringDefect] <;> positivity

/-- `εᵢ ≤ 0.052`, since `6 εᵢ² ≤ εᵢ`. -/
theorem ringDefect_le (i : ℕ) : ringDefect i ≤ 13 / 250 := by
  induction i with
  | zero => simp [ringDefect]
  | succ i ih =>
    have := ringDefect_nonneg i
    simp only [ringDefect]
    nlinarith

section Ring

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)] (E : (L : ℕ) → ι L → ℝ)

/-- **Claim (R).** From `P_s(120, 49), P_τ(240, 49) ≥ 0.948`, induction on the scale gives
`P_s(120 · 2^i, 49 · 2^i), P_τ(240 · 2^i, 49 · 2^i) ≥ 1 - εᵢ`. -/
theorem ring_scales (hZ : HasSpatialMomentBounds E 2)
    (hs : 1 - 13 / 250 ≤ spatialPurity E 120 49) (ht : 1 - 13 / 250 ≤ thermalPurity E 240 49)
    (i : ℕ) :
    1 - ringDefect i ≤ spatialPurity E (120 * 2 ^ i) (49 * 2 ^ i) ∧
      1 - ringDefect i ≤ thermalPurity E (240 * 2 ^ i) (49 * 2 ^ i) := by
  induction i with
  | zero => simp only [ringDefect, pow_zero, mul_one]; exact ⟨hs, ht⟩
  | succ i ih =>
    have h := ring_doubling_six E hZ (n := 120 * 2 ^ i)
      (by have := Nat.one_le_two_pow (n := i); omega) ⟨60 * 2 ^ i, by ring⟩
      (by positivity : (0 : ℝ) < 49 * 2 ^ i) (ringDefect_nonneg i) (ringDefect_le i) ih.1
      (by simpa [show 2 * (120 * 2 ^ i) = 240 * 2 ^ i by ring] using ih.2)
    simp only [ringDefect]
    have e1 : 2 * (120 * 2 ^ i) = 120 * 2 ^ (i + 1) := by ring
    have e2 : 4 * (120 * 2 ^ i) = 240 * 2 ^ (i + 1) := by ring
    have e3 : (2 : ℝ) * (49 * 2 ^ i) = 49 * 2 ^ (i + 1) := by ring
    rw [e1, e2, e3] at h
    exact h

end Ring

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
  have hbulk : R * (μ ^ a) ^ 2 ≤ μ' ^ a := by
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
  have h12 : 0 ≤ 1 - 2 * E' := by linarith
  have h1η : 0 ≤ 1 - η := by linarith
  calc (1 - 2 * E') * (1 - η) * R * G (2 * N + 1) ^ 2
      ≤ (1 - 2 * E') * (1 - η) * R * ((μ ^ a) ^ 2 * G (2 * K) ^ 2) := by gcongr
    _ = (1 - 2 * E') * (R * (μ ^ a) ^ 2) * ((1 - η) * G (2 * K) ^ 2) := by ring
    _ ≤ (1 - 2 * E') * μ' ^ a * G' (2 * K) := by gcongr
    _ = (1 - 2 * E') * (μ' ^ a * G' (2 * K)) := by ring
    _ ≤ G' (2 * N + 1) := hlow

end OpenMomentWitness

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

/-! ### The boundary-field gap -/

section Main

variable {ιr ιo : ℕ → Type*} [∀ L, Fintype (ιr L)] [∀ L, Nonempty (ιr L)]
  [∀ n, Fintype (ιo n)] [∀ n, Nonempty (ιo n)] (Er : (L : ℕ) → ιr L → ℝ)
  (Eo : (n : ℕ) → ιo n → ℝ)

/-- **The boundary-field gap criterion** (Theorem 2 of the note). Let the ring spectra `Er`
satisfy the spatial moment inequalities and the ring and open-chain traces the boundary moment
inequalities at every `β > 0`. If

* `P_s(120, 49) ≥ 0.948` and `P_τ(240, 49) ≥ 0.948` for the ring,
* the boundary defect `e₁₄₂ ≤ 0.14` at `β = 49`, and
* the open-chain thermal purity `P_τ(240, 49) ≥ 0.79`,

then every odd open chain of `n ≥ 1921` sites has a nondegenerate ground level and gap
`Δ_n > (log 10)/392 ≈ 0.005874`. -/
theorem gapGT_of_openSeeds (hring : HasSpatialMomentBounds Er 2)
    (hopen : ∀ β : ℝ, 0 < β → OpenMomentBounds (fun L ↦ partitionFunction Er L β)
      (fun n ↦ partitionFunction Eo n β))
    (hs : 948 / 1000 ≤ spatialPurity Er 120 49) (ht : 948 / 1000 ≤ thermalPurity Er 240 49)
    (he : boundaryDefect (fun L ↦ partitionFunction Er L 49) (fun n ↦ partitionFunction Eo n 49)
      142 ≤ 14 / 100)
    (hη : 79 / 100 ≤ thermalPurity Eo 240 49) {n : ℕ} (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (Eo n) (log 10 / 392) := by
  /- The scales `βᵢ = 49 · 2^i` and `Nᵢ = 60 · 2^i`: ring length `nᵢ = 2Nᵢ`, open-chain length
  `Mᵢ = 2nᵢ`, boundary-defect length `kᵢ = nᵢ` for `i ≥ 1` (and `k₀ = 142`). -/
  set β : ℕ → ℝ := fun i ↦ 49 * 2 ^ i with hβdef
  have hβpos : ∀ i, 0 < β i := fun i ↦ by simp only [β]; positivity
  have hβs : ∀ i, β (i + 1) = 2 * β i := fun i ↦ by simp only [β]; ring
  have hβ0 : β 0 = 49 := by simp [β]
  set Nn : ℕ → ℕ := fun i ↦ 60 * 2 ^ i with hNndef
  have hNn1 : ∀ i, 1 ≤ Nn i := fun i ↦ by
    simp only [Nn]; have := Nat.one_le_two_pow (n := i); omega
  have hNns : ∀ i, Nn (i + 1) = 2 * Nn i := fun i ↦ by simp only [Nn]; ring
  set ε := ringDefect
  have hε0 : ∀ i, 0 ≤ ε i := ringDefect_nonneg
  have hε1 : ∀ i, ε i ≤ 13 / 250 := ringDefect_le
  -- The witnesses at every scale.
  choose μ c ρ s hw using fun i ↦ hopen (β i) (hβpos i)
  have hw' : ∀ i, OpenMomentWitness (fun L ↦ partitionFunction Er L (2 * β i))
      (fun n ↦ partitionFunction Eo n (2 * β i)) (μ (i + 1)) (c (i + 1)) (ρ (i + 1))
      (s (i + 1)) := fun i ↦ by
    have := hw (i + 1); rwa [hβs] at this
  have hμ : ∀ i, 0 < μ i := fun i ↦ (hw i).mu_pos (partitionFunction_pos Er 4 _)
  have hZ : ∀ L β', 0 < partitionFunction Er L β' := partitionFunction_pos Er
  have hG : ∀ m β', 0 < partitionFunction Eo m β' := partitionFunction_pos Eo
  /- Claim (R): the ring purities, the residual ratios `ρᵢ^{nᵢ} ≤ d(εᵢ)`, the growth bound
  `μᵢ₊₁ ≤ μᵢ²` and the bulk factor `R_i = (1 - εᵢ)(1 - εᵢ₊₁)`. -/
  have hR := ring_scales Er hring (by linarith) (by linarith)
  have hRs : ∀ i, 1 - ε i ≤ partitionFunction Er (2 * (2 * Nn i)) (β i) /
      partitionFunction Er (2 * Nn i) (β i) ^ 2 := fun i ↦ by
    have := (hR i).1
    rwa [spatialPurity, show 120 * 2 ^ i = 2 * Nn i by simp only [Nn]; ring] at this
  have hRt : ∀ i, (1 - ε i) * partitionFunction Er (2 * (2 * Nn i)) (β i) ^ 2 ≤
      partitionFunction Er (2 * (2 * Nn i)) (2 * β i) := fun i ↦ by
    have := (hR i).2
    rwa [thermalPurity, show 240 * 2 ^ i = 2 * (2 * Nn i) by simp only [Nn]; ring,
      le_div_iff₀ (pow_pos (hZ _ _) 2)] at this
  have hε10 : ∀ i, ε i ≤ 1 / 10 := fun i ↦ (hε1 i).trans (by norm_num)
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
      ((hε1 i).trans (by norm_num)) ((hε1 _).trans (by norm_num)) (hRt i) hP'
  have hR0 : ∀ i, 0 ≤ (1 - ε i) * (1 - ε (i + 1)) := fun i ↦
    mul_nonneg (by linarith [hε1 i]) (by linarith [hε1 (i + 1)])
  /- Claims (E) and (T) at scale `i ≥ 1`: `(1 - E) G(nᵢ) ≤ cᵢ μᵢ^{nᵢ}` and
  `P_τ^{open}(Mᵢ, βᵢ) ≥ 1 - η`. -/
  have step : ∀ i, 1 ≤ i → ∀ {q E η E' η' : ℝ}, 0 ≤ q → q < 1 →
      tailBound (ε (i + 1)) ≤ q ^ 2 → E < 1 → η < 1 → overlapUpdate q E (1 - η) ≤ E' →
      E' ≤ 1 → 0 < (1 - η) * ((1 - ε i) * (1 - ε (i + 1))) * (1 - E') →
      squaringDefect (1 - (1 - η) * ((1 - ε i) * (1 - ε (i + 1))) * (1 - E')) ≤ η' →
      (1 - E) * partitionFunction Eo (2 * Nn i) (β i) ≤ c i * μ i ^ (2 * Nn i) →
      1 - η ≤ thermalPurity Eo (2 * (2 * Nn i)) (β i) →
      (1 - E') * partitionFunction Eo (2 * Nn (i + 1)) (β (i + 1)) ≤
          c (i + 1) * μ (i + 1) ^ (2 * Nn (i + 1)) ∧
        1 - η' ≤ thermalPurity Eo (2 * (2 * Nn (i + 1))) (β (i + 1)) := by
    intro i _ q E η E' η' hq0 hq1 hqd hE hη hE' hE'1 hprod hη' he hA
    have hq : ρ (i + 1) ^ (2 * (2 * Nn i) - 2 * Nn i) ≤ q := by
      have := htail (i + 1)
      rw [hNns] at this
      rw [show 2 * (2 * Nn i) - 2 * Nn i = 2 * Nn i by omega]
      exact pow_le_of_pow_two_mul_le (hw (i + 1)).rho_nonneg hq0 this hqd
    have h := open_step Eo (hw i) (hw' i) (hμ i) (hμ _) (hr i) (K := Nn i) (N := 2 * Nn i)
      (hNn1 i) (by omega) (hR0 i) (hRbulk i) hq0 hq1 hq hE he hA hη hE' hE'1 hprod hη'
    rw [hNns, hβs]
    exact h
  -- The base case from the seeds: `e₁₄₂ ≤ 0.14` and `P_τ^{open}(240, 49) ≥ 0.79`.
  have w0 := hw 0
  rw [hβ0] at w0
  have hρ0 : ρ 0 < 1 := by
    have := hRs 0
    rw [hβ0] at this
    exact w0.rho_lt_one (hμ 0) (hNn1 0) (hZ _ _) (hε10 0) this
  have hS0e : (1 - 14 / 100) * partitionFunction Eo (2 * 71) (β 0) ≤ c 0 * μ 0 ^ (2 * 71) := by
    have hZ4 := hZ 4 49
    rw [boundaryDefect, w0.leadingRoot_eq hZ4, w0.leadingOverlap_eq hZ4 hρ0] at he
    have hG142 := hG 142 49
    rw [hβ0, show 2 * 71 = 142 from rfl]
    have : 86 / 100 ≤ c 0 * μ 0 ^ 142 / partitionFunction Eo 142 49 := by linarith
    rw [le_div_iff₀ hG142] at this
    linarith
  have hS0t : 1 - 21 / 100 ≤ thermalPurity Eo (2 * (2 * Nn 0)) (β 0) := by
    rw [show 2 * (2 * Nn 0) = 240 by simp [Nn], hβ0]; linarith
  -- Step `0 → 1`, with `k₀ = 142`, `q₀ = 0.142`.
  have hS1 : (1 - 12 / 100) * partitionFunction Eo (2 * Nn 1) (β 1) ≤ c 1 * μ 1 ^ (2 * Nn 1) ∧
      1 - 155 / 1000 ≤ thermalPurity Eo (2 * (2 * Nn 1)) (β 1) := by
    have hq : ρ 1 ^ (2 * (2 * Nn 0) - 2 * 71) ≤ 71 / 500 := by
      have h240 := htail 1
      rw [show 2 * Nn 1 = 240 by simp [Nn]] at h240
      rw [show 2 * (2 * Nn 0) - 2 * 71 = 98 by simp [Nn]]
      have hρ1 := (hw 1).rho_nonneg
      refine (pow_le_pow_iff_left₀ (pow_nonneg hρ1 _) (by norm_num) (by norm_num : 120 ≠ 0)).1 ?_
      calc (ρ 1 ^ 98) ^ 120 = (ρ 1 ^ 240) ^ 49 := by ring
        _ ≤ tailBound (ε 1) ^ 49 := pow_le_pow_left₀ (pow_nonneg hρ1 _) h240 49
        _ ≤ (71 / 500) ^ 120 := by simp only [ε, ringDefect, tailBound]; norm_num
    have h := open_step Eo (hw 0) (hw' 0) (hμ 0) (hμ _) (hr 0) (K := 71) (N := 2 * Nn 0)
      (by norm_num) (by simp [Nn]) (hR0 0) (hRbulk 0) (by norm_num) (by norm_num) hq
      (by norm_num) hS0e hS0t (by norm_num) (E' := 12 / 100) (η' := 155 / 1000)
      (by simp only [overlapUpdate]; norm_num) (by norm_num)
      (by simp only [ε, ringDefect]; norm_num)
      (by simp only [ε, ringDefect, squaringDefect]; norm_num)
    rw [← hβs, show 2 * (2 * Nn 0) = 2 * Nn 1 by simp [Nn]] at h
    exact h
  -- Steps `1 → 2 → 3`, with `q₁ = 0.0282`, `q₂ = 0.0028`.
  have hS2 := step 1 le_rfl (q := 141 / 5000) (E' := 16 / 1000) (η' := 26 / 1000) (by norm_num)
    (by norm_num) (by simp only [ε, ringDefect, tailBound]; norm_num) (by norm_num) (by norm_num)
    (by simp only [overlapUpdate]; norm_num) (by norm_num) (by simp only [ε, ringDefect]; norm_num)
    (by simp only [ε, ringDefect, squaringDefect]; norm_num) hS1.1 hS1.2
  have hS3 := step 2 (by norm_num) (q := 7 / 2500) (E' := 1 / 1000) (η' := 1 / 1000)
    (by norm_num) (by norm_num) (by simp only [ε, ringDefect, tailBound]; norm_num)
    (by norm_num) (by norm_num) (by simp only [overlapUpdate]; norm_num) (by norm_num)
    (by simp only [ε, ringDefect]; norm_num)
    (by simp only [ε, ringDefect, squaringDefect]; norm_num) hS2.1 hS2.2
  /- The uniform phase: with `x_j = 10^{-(2^j + 1)}`, all three defects at scale `3 + j` are at
  most `x_j`, and `x_{j+1} = 10 x_j²`. -/
  set x : ℕ → ℝ := fun j ↦ (1 / 10) ^ (2 ^ j + 1) with hxdef
  have hx0 : ∀ j, 0 < x j := fun j ↦ by simp only [x]; positivity
  have hxle : ∀ j, x j ≤ 1 / 100 := fun j ↦ by
    simp only [x]
    calc ((1 : ℝ) / 10) ^ (2 ^ j + 1) ≤ (1 / 10) ^ 2 :=
          pow_le_pow_of_le_one (by norm_num) (by norm_num)
            (by have := Nat.one_le_two_pow (n := j); omega)
      _ = 1 / 100 := by norm_num
  have hxs : ∀ j, x (j + 1) = 10 * x j ^ 2 := fun j ↦ by
    simp only [x]; rw [pow_succ 2 j, ← pow_mul]; ring
  have hunif : ∀ j, ((1 - x j) * partitionFunction Eo (2 * Nn (3 + j)) (β (3 + j)) ≤
      c (3 + j) * μ (3 + j) ^ (2 * Nn (3 + j)) ∧
        1 - x j ≤ thermalPurity Eo (2 * (2 * Nn (3 + j))) (β (3 + j))) ∧ ε (3 + j) ≤ x j := by
    intro j
    induction j with
    | zero =>
      have hx : x 0 = 1 / 100 := by simp [x]; norm_num
      refine ⟨⟨le_trans (mul_le_mul_of_nonneg_right (by rw [hx]; norm_num) (hG _ _).le)
        hS3.1, le_trans (by rw [hx]; norm_num) hS3.2⟩, ?_⟩
      rw [hx]; simp only [ε, ringDefect]; norm_num
    | succ j ih =>
      obtain ⟨⟨he, hA⟩, hεj⟩ := ih
      set y := x j
      have hy0 := (hx0 j).le
      have hy := hxle j
      have hysq : y ^ 2 ≤ 1 / 10000 := by
        have := pow_le_pow_left₀ hy0 hy 2; norm_num at this ⊢; linarith
      have hysq0 : 0 ≤ y ^ 2 := sq_nonneg y
      have hεs : ε (3 + j + 1) ≤ 6 * y ^ 2 := by
        change 6 * ε (3 + j) ^ 2 ≤ 6 * y ^ 2
        have := pow_le_pow_left₀ (hε0 _) hεj 2
        linarith
      set Ri := (1 - ε (3 + j)) * (1 - ε (3 + j + 1)) with hRi
      have hRlow : (1 - y) * (1 - 6 * y ^ 2) ≤ Ri :=
        mul_le_mul (by linarith) (by linarith) (by linarith) (by linarith [hε1 (3 + j)])
      have hRle : Ri ≤ 1 := (mul_le_mul (by linarith [hε0 (3 + j)]) (by linarith [hε0 (3 + j + 1)])
        (by linarith [hε1 (3 + j + 1)]) zero_le_one).trans_eq (mul_one 1)
      have h10 : 0 ≤ 1 - 10 * y ^ 2 := by linarith
      have hlow : (1 - y) * ((1 - y) * (1 - 6 * y ^ 2)) * (1 - 10 * y ^ 2) ≤
          (1 - y) * Ri * (1 - 10 * y ^ 2) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hRlow (by linarith)) h10
      have hpos : 0 < (1 - y) * ((1 - y) * (1 - 6 * y ^ 2)) * (1 - 10 * y ^ 2) :=
        mul_pos (mul_pos (by linarith) (mul_pos (by linarith) (by linarith))) (by linarith)
      have hle1 : (1 - y) * Ri * (1 - 10 * y ^ 2) ≤ 1 :=
        (mul_le_mul ((mul_le_mul (by linarith) hRle (hR0 _) zero_le_one).trans_eq (mul_one 1))
          (by linarith) h10 zero_le_one).trans_eq (mul_one 1)
      have h := step (3 + j) (by omega) (q := 2 * y) (E := y) (η := y) (E' := 10 * y ^ 2)
        (η' := 10 * y ^ 2) (by linarith) (by linarith)
        ((tailBound_mono hεs (by linarith)).trans (uniform_tail hy0 hy)) (by linarith)
        (by linarith) (uniform_overlap hy0 hy) (by linarith) (hpos.trans_le hlow)
        (le_trans (squaringDefect_mono (by linarith)
          (by linarith [uniform_product hy0 hy]) (by linarith)) (uniform_squaring hy0 hy))
        he hA
      rw [hxs]
      exact ⟨h, hεs.trans (by linarith)⟩
  -- The final scale: `Mᵢ ≤ n ≤ 2 Mᵢ` with `i = 3 + j`.
  obtain ⟨j, hj1, hj2⟩ := exists_dyadic_window (n₀ := 1920) (by norm_num) (by omega : 1920 ≤ n)
  obtain ⟨N, rfl⟩ := hno
  have hM : 2 * (2 * Nn (3 + j)) = 1920 * 2 ^ j := by simp only [Nn]; rw [pow_add]; ring
  obtain ⟨⟨-, hA⟩, hεj⟩ := hunif j
  obtain ⟨⟨he', -⟩, -⟩ := hunif (j + 1)
  set y := x j
  have hy0 := (hx0 j).le
  have hy := hxle j
  have hysq : y ^ 2 ≤ 1 / 10000 := by
    have := pow_le_pow_left₀ hy0 hy 2; norm_num at this ⊢; linarith
  have hysq0 : 0 ≤ y ^ 2 := sq_nonneg y
  rw [show 3 + (j + 1) = 3 + j + 1 by omega, hNns, hβs] at he'
  rw [thermalPurity, le_div_iff₀ (pow_pos (hG _ _) 2)] at hA
  have hpur := (hw (3 + j)).oddChain_purity (hw' (3 + j)) (hμ _) (hr _) (K := 2 * Nn (3 + j))
    (N := N) (by have := hNn1 (3 + j); omega) (by omega) (by omega) (hG _ _)
    (E' := x (j + 1)) (by rw [hxs]; linarith) he' (η := y) (by linarith) hA (hR0 _)
    (hRbulk _)
  -- `P_τ(n, βᵢ) ≥ (1 - 2x')(1 - x) Rᵢ ≥ 1 - 4x`.
  have hP : 1 - 4 * y ≤ thermalPurity Eo (2 * N + 1) (β (3 + j)) := by
    rw [thermalPurity, le_div_iff₀ (pow_pos (hG _ _) 2)]
    have hεs : ε (3 + j + 1) ≤ 6 * y ^ 2 := by
      change 6 * ε (3 + j) ^ 2 ≤ 6 * y ^ 2
      have := pow_le_pow_left₀ (hε0 _) hεj 2
      linarith
    have hRlow : (1 - y) * (1 - 6 * y ^ 2) ≤ (1 - ε (3 + j)) * (1 - ε (3 + j + 1)) :=
      mul_le_mul (by linarith) (by linarith) (by linarith) (by linarith [hε1 (3 + j)])
    have hfin := uniform_final hy0 hy
    have hlow : (1 - 2 * (10 * y ^ 2)) * (1 - y) * ((1 - y) * (1 - 6 * y ^ 2)) ≤
        (1 - 2 * x (j + 1)) * (1 - y) * ((1 - ε (3 + j)) * (1 - ε (3 + j + 1))) := by
      rw [hxs]
      exact mul_le_mul_of_nonneg_left hRlow (mul_nonneg (by linarith) (by linarith))
    have hG2 := pow_pos (hG (2 * N + 1) (β (3 + j))) 2
    calc (1 - 4 * y) * partitionFunction Eo (2 * N + 1) (β (3 + j)) ^ 2
        ≤ (1 - 2 * x (j + 1)) * (1 - y) * ((1 - ε (3 + j)) * (1 - ε (3 + j + 1))) *
            partitionFunction Eo (2 * N + 1) (β (3 + j)) ^ 2 :=
          mul_le_mul_of_nonneg_right (hfin.trans hlow) hG2.le
      _ ≤ partitionFunction Eo (2 * N + 1) (2 * β (3 + j)) := hpur
  -- From the purity to the gap: `e^{βᵢ (log 10)/392} · 4x = 10^{2^j} · 4 · 10^{-(2^j + 1)}`.
  refine gapGT_of_thermalPurity Eo (2 * N + 1) (hβpos _) (by linarith) hP ?_
  have hexp : exp (β (3 + j) * (log 10 / 392)) = 10 ^ (2 ^ j) := by
    have : β (3 + j) * (log 10 / 392) = ((2 ^ j : ℕ) : ℝ) * log 10 := by
      simp only [β]; push_cast; rw [pow_add]; ring
    rw [this, exp_nat_mul, exp_log (by norm_num)]
  have hx10 : (10 : ℝ) ^ (2 ^ j) * y = 1 / 10 := by
    simp only [y, x]; rw [pow_succ, ← mul_assoc, ← mul_pow]; norm_num
  rw [hexp, show (10 : ℝ) ^ 2 ^ j * (4 * y) = 4 * (10 ^ 2 ^ j * y) by ring, hx10]
  linarith

end Main

end CertifiedGaps
