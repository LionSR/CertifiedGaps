/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.SectorTransferC
import CertifiedGaps.Ring.EnergyBound

/-!
# Uniform energy bounds for the spin-one Heisenberg ring

The Lemma "Energy bounds" of the spin-one Haldane-gap paper: for every rotation `R` of the
closing bond and every ring length `L ≥ 3`,

  `-(3L/2) ≤ H_L(R) ≤ L`,

where `H_L(R) = ∑_{j < L-1} Sⱼ · Sⱼ₊₁ + ∑_{αβ} R_{αβ} S^α_{L-1} S^β_0` (`rotTwistedHamiltonian`).
The untwisted ring (`hamiltonian`), the ring twisted by the `π` rotation about `z`
(`twistedHamiltonian`) and the ring twisted by `C = C_z(2π/3)` (`cTwistedHamiltonian`) are
special cases.

**Upper bound.** In the Cartesian basis the exchange is `S · S = SWAP - |Φ⟩⟨Φ|` with
`Φ = ∑ₛ |s s⟩`, an integer matrix (`exchange_eq_map`) satisfying
`(h - 1)(h + 1)(h + 2) = 0` (`exchangeInt_poly`, checked by `decide`), so `h ≤ 1`. The seam bond
is the conjugate of `h` by `1 ⊗ R` (`conj_exchange_right`), hence also `≤ 1`, and summing the `L`
bonds gives `H_L(R) ≤ L`.

**Lower bound.** For three distinct sites, `T = h₁₂ + h₂₃` has no diagonal entries, and with the
weights `q(aaa) = 4`, `q(aab) = q(abb) = 3`, `q(aba) = q(abc) = 2` of the three letters every
weighted absolute row sum is at most `3q` (`triad_rowSum`, checked by `decide`). The weighted
Schur test gives `|T| ≤ 3` (`abs_qform_triad`). A pair of bonds containing the seam is the
conjugate of an untwisted pair by a product of one-site rotations: rotate the sites after the
seam (`abs_qform_pair_conj`). Summing the `L` pairs counts every bond twice, so
`H_L(R) ≥ -(3/2) L`.

## Main results

* `qform_exchange_le`, `qform_rotExchange_le`: `S · S ≤ 1` and the same for the twisted bond.
* `abs_qform_triad`: `|h₁₂ + h₂₃| ≤ 3` on three distinct sites of a ring.
* `le_qform_rotTwistedHamiltonian`, `qform_rotTwistedHamiltonian_le`: the energy bounds for
  every proper rotation `R` (`Rᵀ R = 1`, `PreservesCross R`).
* The corollaries for `hamiltonian`, `twistedHamiltonian` and `cTwistedHamiltonian`.

## Comparison with the paper

The paper states the lemma for `n ≥ 4`; the proof only uses that the three sites of a pair of
consecutive bonds are distinct, so the lower bound is proved for `L ≥ 3` and the upper bound for
`L ≥ 2`. The paper's Schur test bounds the spectral radius of `D⁻¹ T D`; here the same weights
bound the quadratic form directly (`Spectral.abs_qform_le_of_schur`).
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

open Spectral SpinOne

/-! ### The exchange as an integer matrix -/

/-- The exchange `S · S` in the Cartesian basis, `SWAP - |Φ⟩⟨Φ|` with `Φ = ∑ₛ |s s⟩`:
the entry `⟨s t| S · S |s' t'⟩` is `[s = t', t = s'] - [s = t, s' = t']`. -/
def exchangeInt : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℤ :=
  of fun p p' ↦ (if p.1 = p'.2 ∧ p.2 = p'.1 then 1 else 0) -
    (if p.1 = p.2 ∧ p'.1 = p'.2 then 1 else 0)

/-- The contraction `∑_α ε_{α s s'} ε_{α t t'} = δ_{st} δ_{s't'} - δ_{st'} δ_{s't}`. -/
theorem sum_leviCivita_mul (s t s' t' : Fin 3) :
    ∑ α, leviCivita α s s' * leviCivita α t t' =
      (if s = t ∧ s' = t' then 1 else 0) - (if s = t' ∧ t = s' then 1 else 0) := by
  revert s t s' t'
  decide

/-- The Cartesian exchange is the integer matrix `exchangeInt`. -/
theorem exchange_eq_map : exchange = exchangeInt.map (↑) := by
  ext ⟨s, t⟩ ⟨s', t'⟩
  have hk := congrArg (fun z : ℤ ↦ (z : ℂ)) (sum_leviCivita_mul s t s' t')
  simp only [Fin.sum_univ_three, Int.cast_add, Int.cast_mul, Int.cast_sub] at hk
  simp only [exchange, Matrix.sum_apply, kroneckerMap_apply, S, of_apply, map_apply, exchangeInt,
    Fin.sum_univ_three, Int.cast_sub]
  linear_combination -hk + ((leviCivita 0 s s' : ℂ) * leviCivita 0 t t' +
    (leviCivita 1 s s' : ℂ) * leviCivita 1 t t' + (leviCivita 2 s s' : ℂ) * leviCivita 2 t t') *
    Complex.I_sq

/-- The entries of the exchange have absolute values `0` or `1`. -/
theorem norm_exchange (p p' : Fin 3 × Fin 3) :
    ‖exchange p p'‖ = ((|exchangeInt p p'| : ℤ) : ℝ) := by
  rw [exchange_eq_map, map_apply, Complex.norm_intCast, Int.cast_abs]

/-- **The spectrum of `S · S` is `{1, -1, -2}`**: `(h - 1)(h + 1)(h + 2) = 0`, checked on the
integer matrix. -/
theorem exchangeInt_poly : (exchangeInt - 1) * ((exchangeInt + 1) * (exchangeInt + 2)) = 0 := by
  have h : ∀ p p', ((exchangeInt - 1) * ((exchangeInt + 1) * (exchangeInt + 2))) p p' = 0 := by
    decide +kernel
  ext p p'
  exact h p p'

/-- The polynomial identity of the exchange in the form of
`Spectral.qform_le_of_prod_eq_zero`. -/
theorem exchange_poly :
    ([1, -1, -2].map fun r : ℝ ↦ exchange - (r : ℂ) • (1 : Matrix (Fin 3 × Fin 3) _ ℂ)).prod
      = 0 := by
  have h := congrArg (Int.castRingHom ℂ).mapMatrix exchangeInt_poly
  rw [map_mul, map_mul, map_sub, map_add, map_add, map_one, map_ofNat, map_zero,
    RingHom.mapMatrix_apply] at h
  rw [exchange_eq_map]
  convert h using 1
  simp only [List.map_cons, List.map_nil, List.prod_cons, List.prod_nil, Matrix.mul_one]
  congr 2 <;> simp [sub_eq_add_neg, two_smul, one_add_one_eq_two]

/-- **The exchange is at most one**: `S · S ≤ 1`. -/
theorem qform_exchange_le (y : Fin 3 × Fin 3 → ℂ) : qform exchange y ≤ 1 * nsq y :=
  qform_le_of_prod_eq_zero exchange_isHermitian exchange_poly
    (by intro r hr; simp only [List.mem_cons, List.not_mem_nil, or_false] at hr
        rcases hr with rfl | rfl | rfl <;> norm_num) y

/-! ### The twisted bond -/

/-- **The twisted bond is a rotated exchange.** Rotating the second site by `R` turns `S · S`
into `∑_{αβ} R_{αβ} S^α ⊗ S^β`. -/
theorem conj_exchange_right {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : PreservesCross R) :
    (1 ⊗ₖ R.map (↑))ᴴ * exchange * (1 ⊗ₖ R.map (↑)) = rotExchange R := by
  simp only [exchange, Finset.mul_sum, Finset.sum_mul, rotExchange]
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [conjTranspose_kronecker, conjTranspose_one, ← mul_kronecker_mul, ← mul_kronecker_mul,
    Matrix.one_mul, Matrix.mul_one, conj_S_of_preservesCross hR]
  simp only [Fin.sum_univ_three, kronecker_add, kronecker_smul]

/-- `1 ⊗ R` is unitary for a real orthogonal `R`. -/
theorem one_kronecker_map_unitary {R : Matrix (Fin 3) (Fin 3) ℝ} (hRR : Rᵀ * R = 1) :
    (1 ⊗ₖ R.map (↑) : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ)ᴴ *
      (1 ⊗ₖ R.map (↑) : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ)ᴴᴴ = 1 := by
  rw [conjTranspose_conjTranspose, conjTranspose_kronecker, conjTranspose_one,
    ← mul_kronecker_mul, Matrix.one_mul, map_conjTranspose_mul_map hRR, one_kronecker_one]

/-- **The twisted bond is at most one** for a proper rotation `R`. -/
theorem qform_rotExchange_le {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : PreservesCross R)
    (hRR : Rᵀ * R = 1) (y : Fin 3 × Fin 3 → ℂ) : qform (rotExchange R) y ≤ 1 * nsq y := by
  have := qform_conj_le (one_kronecker_map_unitary hRR) qform_exchange_le y
  rwa [conjTranspose_conjTranspose, conj_exchange_right hR] at this

/-! ### Two consecutive bonds -/

/-- The Schur weight of a word of three letters: `4` for `aaa`, `3` for `aab` and `abb`, `2` for
`aba` and `abc`. -/
def triadWeight (a b c : Fin 3) : ℤ :=
  if a = b ∧ b = c then 4 else if a = b ∨ b = c then 3 else 2

/-- **The local Schur test.** For every word `s₁ s₂ s₃`, the absolute row sum of
`h₁₂ + h₂₃` weighted by `triadWeight` is at most three times the weight of the word. -/
theorem triad_rowSum (s₁ s₂ s₃ : Fin 3) :
    ∑ p : Fin 3 × Fin 3, |exchangeInt (s₁, s₂) p| * triadWeight p.1 p.2 s₃ +
      ∑ p : Fin 3 × Fin 3, |exchangeInt (s₂, s₃) p| * triadWeight s₁ p.1 p.2 ≤
        3 * triadWeight s₁ s₂ s₃ := by
  revert s₁ s₂ s₃
  decide

/-- **Two bonds of three distinct sites**: `|h₁₂ + h₂₃| ≤ 3` on a ring, by the weighted Schur
test with weights depending only on the letters at the three sites. -/
theorem abs_qform_triad {L : ℕ} {i₁ i₂ i₃ : Fin L} (h12 : i₁ ≠ i₂) (h23 : i₂ ≠ i₃)
    (h13 : i₁ ≠ i₃) (x : (Fin L → Fin 3) → ℂ) :
    |qform (bondOp i₁ i₂ exchange + bondOp i₂ i₃ exchange) x| ≤ 3 * nsq x := by
  have hH₁ := bondOp_isHermitian i₁ i₂ exchange_isHermitian
  have hH₂ := bondOp_isHermitian i₂ i₃ exchange_isHermitian
  refine abs_qform_le_of_schur
    (K := fun σ τ ↦ ‖bondOp i₁ i₂ exchange σ τ‖ + ‖bondOp i₂ i₃ exchange σ τ‖)
    (q := fun σ ↦ (triadWeight (σ i₁) (σ i₂) (σ i₃) : ℝ)) (fun σ τ ↦ norm_add_le _ _)
    (fun σ τ ↦ ?_) (fun σ ↦ ?_) (fun σ ↦ ?_) x
  · rw [← hH₁.apply σ τ, norm_star, ← hH₂.apply σ τ, norm_star]
  · simp only [triadWeight]
    split_ifs <;> norm_num
  · simp only [add_mul, Finset.sum_add_distrib]
    rw [sum_norm_bondOp_mul h12, sum_norm_bondOp_mul h23]
    simp only [Function.update_self, Function.update_of_ne h12, Function.update_of_ne h13,
      Function.update_of_ne h23, Function.update_of_ne h23.symm, Function.update_of_ne h13.symm,
      norm_exchange]
    exact_mod_cast triad_rowSum (σ i₁) (σ i₂) (σ i₃)

/-- **Rotated pairs of bonds.** Conjugating `h₁₂ + h₂₃` by a product `⨂ᵢ aᵢ` of one-site unitaries
keeps the bound `|·| ≤ 3`; the two bonds become `(a₁ ⊗ a₂) h (a₁ ⊗ a₂)ᴴ` and
`(a₂ ⊗ a₃) h (a₂ ⊗ a₃)ᴴ`. -/
theorem abs_qform_pair_conj {L : ℕ} (a : Fin L → Matrix (Fin 3) (Fin 3) ℂ)
    (ha : ∀ i, a i * (a i)ᴴ = 1) {i₁ i₂ i₃ : Fin L} (h12 : i₁ ≠ i₂) (h23 : i₂ ≠ i₃)
    (h13 : i₁ ≠ i₃) {B₁ B₂ : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ}
    (hB₁ : (a i₁ ⊗ₖ a i₂) * exchange * (a i₁ ⊗ₖ a i₂)ᴴ = B₁)
    (hB₂ : (a i₂ ⊗ₖ a i₃) * exchange * (a i₂ ⊗ₖ a i₃)ᴴ = B₂) (x : (Fin L → Fin 3) → ℂ) :
    |qform (bondOp i₁ i₂ B₁ + bondOp i₂ i₃ B₂) x| ≤ 3 * nsq x := by
  have hconj : piTensor a * (bondOp i₁ i₂ exchange + bondOp i₂ i₃ exchange) * (piTensor a)ᴴ =
      bondOp i₁ i₂ B₁ + bondOp i₂ i₃ B₂ := by
    rw [conjTranspose_piTensor, Matrix.mul_add, Matrix.add_mul, conj_bondOp a ha h12,
      conj_bondOp a ha h23, hB₁, hB₂]
  have hM : piTensor a * (piTensor a)ᴴ = 1 := by
    rw [conjTranspose_piTensor, piTensor_mul]
    rw [show (fun i ↦ a i * (a i)ᴴ) = fun _ ↦ 1 from funext ha, piTensor_one]
  rw [← hconj]
  exact abs_qform_conj_le hM (abs_qform_triad h12 h23 h13) x

/-- The three sites `j, j + 1, j + 2` of a ring of `L ≥ 3` sites are distinct. -/
theorem ne_finRotate_finRotate {L : ℕ} (hL : 3 ≤ L) (j : Fin L) :
    j ≠ finRotate L (finRotate L j) := by
  intro h
  have h1 := val_finRotate j
  have h2 := val_finRotate (finRotate L j)
  rw [← h] at h2
  have := j.isLt
  split_ifs at h1 h2 <;> omega

/-- **Pairs of bonds of the twisted ring.** For a proper rotation `R` and `L ≥ 3`, every pair of
consecutive bonds of `H_L(R)` satisfies `-3 ≤ B_{j,j+1} + B_{j+1,j+2}`. A pair containing the
seam is rotated into an untwisted pair by `R` on the sites after the seam. -/
theorem le_qform_pair {L : ℕ} (hL : 3 ≤ L) {R : Matrix (Fin 3) (Fin 3) ℝ}
    (hR : PreservesCross R) (hRR : Rᵀ * R = 1) (j : Fin L) (x : (Fin L → Fin 3) → ℂ) :
    -3 * nsq x ≤ qform (twistedBond exchange (rotExchange R) j +
      twistedBond exchange (rotExchange R) (finRotate L j)) x := by
  set m := finRotate L j with hm
  set m' := finRotate L m with hm'
  have hjm : j ≠ m := ne_finRotate (by omega) j
  have hmm' : m ≠ m' := ne_finRotate (by omega) m
  have hjm' : j ≠ m' := ne_finRotate_finRotate hL j
  set V : Matrix (Fin 3) (Fin 3) ℂ := R.map (↑) with hV_def
  have hV : Vᴴ * V = 1 := map_conjTranspose_mul_map hRR
  have hVh : Vᴴ * Vᴴᴴ = 1 := by rw [conjTranspose_conjTranspose, hV]
  -- The three local conjugations: none, the second site, and both sites.
  have F1 : (1 ⊗ₖ Vᴴ) * exchange * (1 ⊗ₖ Vᴴ)ᴴ = rotExchange R := by
    rw [← conj_exchange_right hR, conjTranspose_kronecker, conjTranspose_conjTranspose,
      conjTranspose_one, conjTranspose_kronecker, conjTranspose_one]
  have F2 : (Vᴴ ⊗ₖ Vᴴ) * exchange * (Vᴴ ⊗ₖ Vᴴ)ᴴ = exchange := by
    rw [conjTranspose_kronecker, conjTranspose_conjTranspose, ← conjTranspose_kronecker,
      ← rotExchange_one, conj_rotExchange hR, Matrix.mul_one, hRR]
  have F3 : (1 ⊗ₖ 1 : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ) * exchange * (1 ⊗ₖ 1)ᴴ =
      exchange := by
    rw [one_kronecker_one, conjTranspose_one, Matrix.one_mul, Matrix.mul_one]
  have hlow : ∀ y : ℝ, |y| ≤ 3 * nsq x → -3 * nsq x ≤ y := fun y hy ↦ by
    linarith [neg_abs_le y]
  simp only [twistedBond]
  by_cases hj : j.val + 1 = L
  · -- The first bond is the seam: rotate the sites `m = 0` and `m' = 1`.
    have hmL : ¬ (m.val + 1 = L) := by
      rw [hm, val_finRotate, ite_eq_left hj]; omega
    rw [ite_eq_left hj, ite_eq_right hmL]
    refine hlow _ (abs_qform_pair_conj (fun i ↦ if i = j then 1 else Vᴴ) (fun i ↦ ?_) hjm hmm'
      hjm' ?_ ?_ x)
    · split_ifs
      · rw [conjTranspose_one, Matrix.mul_one]
      · exact hVh
    · simp only [ite_eq_right hjm.symm]; exact F1
    · simp only [ite_eq_right hjm.symm, ite_eq_right hjm'.symm]; exact F2
  · by_cases hmL : m.val + 1 = L
    · -- The second bond is the seam: rotate the site `m' = 0`.
      rw [ite_eq_right hj, ite_eq_left hmL]
      refine hlow _ (abs_qform_pair_conj (fun i ↦ if i = m' then Vᴴ else 1) (fun i ↦ ?_) hjm
        hmm' hjm' ?_ ?_ x)
      · split_ifs
        · exact hVh
        · rw [conjTranspose_one, Matrix.mul_one]
      · simp only [ite_eq_right hjm', ite_eq_right hmm']; exact F3
      · simp only [ite_eq_right hmm']; exact F1
    · rw [ite_eq_right hj, ite_eq_right hmL]
      exact hlow _ (abs_qform_triad hjm hmm' hjm' x)

/-! ### The energy bounds -/

/-- **Lower energy bound** (Lemma "Energy bounds", lower half). For every proper rotation `R`
and every ring of `L ≥ 3` sites, `-(3L/2) ≤ H_L(R)`. -/
theorem le_qform_rotTwistedHamiltonian {L : ℕ} (hL : 3 ≤ L) {R : Matrix (Fin 3) (Fin 3) ℝ}
    (hR : PreservesCross R) (hRR : Rᵀ * R = 1) (x : (Fin L → Fin 3) → ℂ) :
    -(3 * L / 2) * nsq x ≤ qform (rotTwistedHamiltonian R L) x := by
  have := le_qform_twistedRingHamiltonian_of_pairs (a := -3)
    (fun j y ↦ le_qform_pair hL hR hRR j y) x
  rw [rotTwistedHamiltonian, show -(3 * (L : ℝ) / 2) = L * -3 / 2 by ring]
  exact this

/-- **Upper energy bound** (Lemma "Energy bounds", upper half). For every proper rotation `R`
and every ring of `L ≥ 2` sites, `H_L(R) ≤ L`. -/
theorem qform_rotTwistedHamiltonian_le {L : ℕ} (hL : 2 ≤ L) {R : Matrix (Fin 3) (Fin 3) ℝ}
    (hR : PreservesCross R) (hRR : Rᵀ * R = 1) (x : (Fin L → Fin 3) → ℂ) :
    qform (rotTwistedHamiltonian R L) x ≤ L * nsq x := by
  have := qform_twistedRingHamiltonian_le hL exchange_isHermitian (rotExchange_isHermitian R)
    qform_exchange_le (qform_rotExchange_le hR hRR) x
  rwa [mul_one] at this

/-! ### The untwisted, `P`-twisted and `C`-twisted rings -/

/-- The identity preserves the cross product. -/
theorem preservesCross_one : PreservesCross (1 : Matrix (Fin 3) (Fin 3) ℝ) := by
  intro α b c
  fin_cases α <;> fin_cases b <;> fin_cases c <;> simp [leviCivita, one_apply]

/-- The untwisted ring is the ring twisted by the identity. -/
theorem hamiltonian_eq_rotTwisted (L : ℕ) : hamiltonian L = rotTwistedHamiltonian 1 L := by
  rw [rotTwistedHamiltonian, rotExchange_one, twistedRingHamiltonian_self]
  rfl

/-- **Energy bounds of the periodic chain**: `-(3L/2) ≤ H_L` for `L ≥ 3`. -/
theorem le_qform_hamiltonian {L : ℕ} (hL : 3 ≤ L) (x : (Fin L → Fin 3) → ℂ) :
    -(3 * L / 2) * nsq x ≤ qform (hamiltonian L) x := by
  rw [hamiltonian_eq_rotTwisted]
  exact le_qform_rotTwistedHamiltonian hL preservesCross_one (by simp) x

/-- **Energy bounds of the periodic chain**: `H_L ≤ L` for `L ≥ 2`. -/
theorem qform_hamiltonian_le {L : ℕ} (hL : 2 ≤ L) (x : (Fin L → Fin 3) → ℂ) :
    qform (hamiltonian L) x ≤ L * nsq x := by
  rw [hamiltonian_eq_rotTwisted]
  exact qform_rotTwistedHamiltonian_le hL preservesCross_one (by simp) x

/-- The `π` rotation about `z` as a rotation matrix, `diag(-1, -1, 1)`. -/
def piRotZ : Matrix (Fin 3) (Fin 3) ℝ := !![-1, 0, 0; 0, -1, 0; 0, 0, 1]

theorem piRotZ_transpose_mul_self : piRotZᵀ * piRotZ = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [piRotZ, mul_apply, Fin.sum_univ_three]

theorem piRotZ_preservesCross : PreservesCross piRotZ := by
  intro α c e
  fin_cases α <;> fin_cases c <;> fin_cases e <;> simp [piRotZ, Fin.sum_univ_three, leviCivita]

/-- The `π`-twisted ring is the ring twisted by `piRotZ`. -/
theorem twistedHamiltonian_eq_rotTwisted (L : ℕ) :
    twistedHamiltonian L = rotTwistedHamiltonian piRotZ L := by
  have h : twistedExchange = rotExchange piRotZ := by
    simp [twistedExchange, rotExchange, piRotZ, twistSign, Fin.sum_univ_three]
  rw [rotTwistedHamiltonian, ← h]
  rfl

/-- **Energy bounds of the `P`-twisted ring**: `-(3L/2) ≤ H_L(P)` for `L ≥ 3`. -/
theorem le_qform_twistedHamiltonian {L : ℕ} (hL : 3 ≤ L) (x : (Fin L → Fin 3) → ℂ) :
    -(3 * L / 2) * nsq x ≤ qform (twistedHamiltonian L) x := by
  rw [twistedHamiltonian_eq_rotTwisted]
  exact le_qform_rotTwistedHamiltonian hL piRotZ_preservesCross piRotZ_transpose_mul_self x

/-- **Energy bounds of the `P`-twisted ring**: `H_L(P) ≤ L` for `L ≥ 2`. -/
theorem qform_twistedHamiltonian_le {L : ℕ} (hL : 2 ≤ L) (x : (Fin L → Fin 3) → ℂ) :
    qform (twistedHamiltonian L) x ≤ L * nsq x := by
  rw [twistedHamiltonian_eq_rotTwisted]
  exact qform_rotTwistedHamiltonian_le hL piRotZ_preservesCross piRotZ_transpose_mul_self x

-- One simp set serves all entries; each entry uses only some of the root identities.
set_option linter.unusedSimpArgs false in
/-- `cRot` is orthogonal. -/
theorem cRot_transpose_mul_self : cRotᵀ * cRot = 1 := by
  have hb : √3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  generalize √3 = b at hb ⊢
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [cRot, mul_apply, Fin.sum_univ_three, one_apply] <;> ring_nf <;> simp [hb] <;> ring_nf

set_option linter.unusedSimpArgs false in
/-- `cRot` preserves the cross product: it is a proper rotation. -/
theorem cRot_preservesCross : PreservesCross cRot := by
  have hb : √3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  generalize √3 = b at hb ⊢
  intro α c e
  fin_cases α <;> fin_cases c <;> fin_cases e <;>
    simp [cRot, Fin.sum_univ_three, leviCivita] <;> ring_nf <;> simp [hb] <;> ring_nf

/-- **Energy bounds of the `C`-twisted ring**: `-(3L/2) ≤ H_L(C)` for `L ≥ 3`. -/
theorem le_qform_cTwistedHamiltonian {L : ℕ} (hL : 3 ≤ L) (x : (Fin L → Fin 3) → ℂ) :
    -(3 * L / 2) * nsq x ≤ qform (cTwistedHamiltonian L) x :=
  le_qform_rotTwistedHamiltonian hL cRot_preservesCross cRot_transpose_mul_self x

/-- **Energy bounds of the `C`-twisted ring**: `H_L(C) ≤ L` for `L ≥ 2`. -/
theorem qform_cTwistedHamiltonian_le {L : ℕ} (hL : 2 ≤ L) (x : (Fin L → Fin 3) → ℂ) :
    qform (cTwistedHamiltonian L) x ≤ L * nsq x :=
  qform_rotTwistedHamiltonian_le hL cRot_preservesCross cRot_transpose_mul_self x

end CertifiedGaps.SpinOneHeisenberg
