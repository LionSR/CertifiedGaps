/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryData
import CertifiedGaps.Heisenberg.BoundaryTraces
import CertifiedGaps.Heisenberg.EnergyBound
import CertifiedGaps.Heisenberg.SzBridge
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Trace.WeightedFrobenius

/-!
# The fifteen thermal traces from the integer recurrence

This file proves `BoundaryTraceTable` (`Heisenberg/BoundaryTraces.lean`) from the hypothesis
`BoundaryTotals` (`Heisenberg/BoundaryData.lean`) that the integer recurrence of the thermal
appendix of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026)
reproduces the printed integers `k_{n,g}` (eq. (appth:integer-table)). The chain of the proof
follows the subsections "An integer matrix", "The exact recurrence", "Error of the matrix
approximation" and "Rational intervals from the printed table" of the appendix.

1. **The `S^z` basis** (`conj_hamiltonian`, `trace_eq_gibbsTrace_Z`): `Z_{n,g}` is the Gibbs
   trace of `H^z_{n,g} + a n`, the Hamiltonian in the `S^z` basis.
2. **The contraction** (`energy_bounds`, `abs_qform_contraction`): `0 ≤ H_{n,g} + C_n ≤ 16` for
   `2 ≤ n ≤ 6`, so `A = 1 - (H + C_n)/8` satisfies `|⟨x, A x⟩| ≤ ‖x‖²`. Consecutive bonds are
   paired (`abs_qform_triad`, `≥ -3`), a left-over bond is `≥ -2`, every bond is `≤ 1`, and each
   boundary field has norm `3/5`.
3. **The integer matrix** (`similarity_apply`, `matrixD_similar`): the banded matrix `band n g`
   of the recurrence, relabelled by the word numbers `∑ⱼ σⱼ 3ʲ`, is `80 J A J⁻¹` with
   `J_σ = √(1 + p_σ)`. The entries are matched hop by hop (`bondOp_exchangeZ_apply`,
   `siteOp_spinZ_two_apply`, `siteOp_spinZ_zero_apply`, `band_entry`).
4. **The error** (`trace_bounds_of_run`): `Trace.gibbsTrace_enclosure_weighted` with
   `Q = 2⁵⁰`, mean `49`, `K = 192`, the weighted Horner bound and the Poisson tail
   `η(49, 192) ≤ 3⁴⁹/2¹⁹²` give `|√W - ‖E‖_F| ≤ 2·10⁻¹⁰`, where `W = N/(2Q²)` is computed from
   the packed recurrence (`weighted_sum_eq`, through `Trace.Packed.run_numerator`).
5. **The table** (`boundaryTraceTable_of_totals`): `k ≤ 10¹⁶ W < k + 1` and the exponential
   enclosure `expLo`, `expHi` (Taylor polynomials of degree `40`) give the rational intervals;
   the fifteen rational comparisons are one `decide +kernel` (`entryCheck_all`).

## Main results

* `boundaryTraceTable_of_totals : BoundaryTotals → BoundaryTraceTable`.
* `trace_bounds_of_run`: the enclosure of one trace from one successful packed run.

## Differences from the paper

* The rounding error is bounded in the weighted norm `‖J⁻¹ ·‖`, in which `D/80` is a
  contraction, column by column (`Trace/WeightedHorner.lean`); the paper bounds it after
  undoing the similarity with `‖J⁻¹‖ ‖J‖ ≤ √2`. Both give the paper's `2·10⁻¹⁰`: here
  `27 (3/2 · 193 · 28 + 2⁵⁰ η) ≤ 2·10⁻¹⁰ · 2⁵⁰`.
* The Poisson tail is bounded by `η ≤ e^m/2^K ≤ 3⁴⁹/2¹⁹²` (`Trace.poissonError_le_three_pow_div`)
  instead of the paper's `2R`, `R < 3.152·10⁻³³` (eq. (appth:tail)); both are negligible.
* The exponential `e^{b₀(C_n - a n)}` is enclosed with the Taylor polynomial of degree `40` and
  the tail `y⁴⁰/40!` (`Trace.exp_mem_Icc_taylor`) instead of degree `80` with the tail
  `2y⁸¹/81!` (eq. (appth:scalar-exp)).
* The lower bound uses `max(0, √W - ε)² ≥ W - 2ε√W`, which does not subtract the paper's
  `4·10⁻²⁰`; the upper bound is the paper's.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryThermal

open Spectral SpinOne BoundaryTraces Trace.Packed

/-! ### Energy bounds of the open chains -/

/-- The spin-one matrices satisfy `(Sᵅ)³ = Sᵅ`: their eigenvalues are `-1, 0, 1`. -/
theorem S_cube (α : Fin 3) : S α * S α * S α = S α := by
  rw [S_eq]
  fin_cases α <;> ext i j <;> fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_three]

/-- `Sᵅ ⊗ 1` is annihilated by `X (X - 1)(X + 1)`. -/
theorem kron_one_poly (α : Fin 3) :
    ([1, 0, -1].map fun r : ℝ ↦
      S α ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℂ) - (r : ℂ) • (1 : Matrix (Fin 3 × Fin 3) _ ℂ)).prod
      = 0 := by
  have h : (S α ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℂ)) * (S α ⊗ₖ 1) * (S α ⊗ₖ 1) = S α ⊗ₖ 1 := by
    rw [← mul_kronecker_mul, ← mul_kronecker_mul, S_cube, Matrix.mul_one, Matrix.mul_one]
  have e : ∀ B : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ,
      (B - ((1 : ℝ) : ℂ) • 1) * ((B - ((0 : ℝ) : ℂ) • 1) * ((B - ((-1 : ℝ) : ℂ) • 1) * 1)) =
        B * B * B - B := fun B ↦ by
    push_cast; simp only [one_smul, zero_smul, sub_zero, neg_smul, sub_neg_eq_add]; noncomm_ring
  simp only [List.map_cons, List.map_nil, List.prod_cons, List.prod_nil]
  rw [e, h, sub_self]

/-- A one-site operator is a two-site operator with the identity on a second site. -/
theorem siteOp_eq_bondOp {L : ℕ} {j k : Fin L} (hjk : j ≠ k) (A : Matrix (Fin 3) (Fin 3) ℂ) :
    siteOp j A = bondOp j k (A ⊗ₖ 1) := by
  rw [bondOp_kronecker hjk, siteOp]
  congr 1; funext i; split_ifs <;> rfl

/-- The roots `1, 0, -1` lie in `[-1, 1]`. -/
private theorem roots_mem {r : ℝ} (hr : r ∈ [(1 : ℝ), 0, -1]) : -1 ≤ r ∧ r ≤ 1 := by
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hr
  rcases hr with rfl | rfl | rfl <;> norm_num

/-- **One-site fields**: `|⟨x, c Sᵅⱼ x⟩| ≤ |c| ‖x‖²` on a chain of at least two sites. -/
theorem abs_qform_siteOp {L : ℕ} {j k : Fin L} (hjk : j ≠ k) (c : ℝ) (α : Fin 3)
    (x : (Fin L → Fin 3) → ℂ) : |qform (siteOp j ((c : ℂ) • S α)) x| ≤ |c| * nsq x := by
  have hB : (S α ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℂ)).IsHermitian := by
    change (S α ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℂ))ᴴ = _
    rw [conjTranspose_kronecker, (S_isHermitian α).eq, conjTranspose_one]
  rw [siteOp_eq_bondOp hjk, smul_kronecker, bondOp_smul,
    show qform ((c : ℂ) • bondOp j k (S α ⊗ₖ 1)) x = c * qform (bondOp j k (S α ⊗ₖ 1)) x from
      qform_smul_left c _ x, abs_mul]
  refine mul_le_mul_of_nonneg_left (abs_le.2 ⟨?_, ?_⟩) (abs_nonneg _)
  · have := le_qform_bondOp hjk hB (a := -1) (fun y ↦ le_qform_of_prod_eq_zero hB
      (kron_one_poly α) (fun r hr ↦ (roots_mem hr).1) y) x
    linarith
  · have := qform_bondOp_le hjk hB (b := 1) (fun y ↦ qform_le_of_prod_eq_zero hB
      (kron_one_poly α) (fun r hr ↦ (roots_mem hr).2) y) x
    linarith

theorem qform_zero_left' {L : ℕ} (x : (Fin L → Fin 3) → ℂ) :
    qform (0 : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ) x = 0 := by simp [qform]

/-- Bond `j → j + 1` of a chain of `L` sites (zero if `j + 1 ≥ L`). -/
noncomputable def chainBond (L j : ℕ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  if h : j + 1 < L then bondOp (⟨j, by omega⟩ : Fin L) ⟨j + 1, h⟩ exchange else 0

/-- A single bond is at least `-2`. -/
theorem le_qform_chainBond (L j : ℕ) (x : (Fin L → Fin 3) → ℂ) :
    -2 * nsq x ≤ qform (chainBond L j) x := by
  unfold chainBond
  split_ifs with h
  · refine le_qform_bondOp (by simp [Fin.ext_iff]) exchange_isHermitian (fun y ↦
      le_qform_of_prod_eq_zero exchange_isHermitian exchange_poly ?_ y) x
    intro r hr; simp only [List.mem_cons, List.not_mem_nil, or_false] at hr
    rcases hr with rfl | rfl | rfl <;> norm_num
  · have := nsq_nonneg x
    simp only [qform_zero_left']
    linarith

/-- A single bond is at most `1`. -/
theorem qform_chainBond_le (L j : ℕ) (x : (Fin L → Fin 3) → ℂ) :
    qform (chainBond L j) x ≤ 1 * nsq x := by
  unfold chainBond
  split_ifs with h
  · exact qform_bondOp_le (by simp [Fin.ext_iff]) exchange_isHermitian qform_exchange_le x
  · simp only [qform_zero_left', one_mul]; exact nsq_nonneg x

/-- **The bonds of an open chain**: the first `m` bonds are at least
`-(3⌊m/2⌋ + 2 (m mod 2))`, pairing consecutive bonds (`abs_qform_triad`). -/
theorem le_qform_chainBonds {L : ℕ} : ∀ m, m + 1 ≤ L → ∀ x : (Fin L → Fin 3) → ℂ,
    -((3 * (m / 2) + 2 * (m % 2) : ℕ) : ℝ) * nsq x ≤
      qform (∑ j ∈ Finset.range m, chainBond L j) x
  | 0, _, x => by simp [qform_zero_left']
  | 1, _, x => by simpa using le_qform_chainBond L 0 x
  | m + 2, h, x => by
    rw [Finset.sum_range_succ, Finset.sum_range_succ, add_assoc, qform_add_left]
    have ih := le_qform_chainBonds m (by omega) x
    have hp : -3 * nsq x ≤ qform (chainBond L m + chainBond L (m + 1)) x := by
      simp only [chainBond, show m + 1 < L by omega, show m + 1 + 1 < L by omega,
        ↓reduceDIte]
      have := abs_qform_triad (i₁ := (⟨m, by omega⟩ : Fin L)) (i₂ := ⟨m + 1, by omega⟩)
        (i₃ := ⟨m + 1 + 1, by omega⟩) (by simp [Fin.ext_iff]) (by simp [Fin.ext_iff])
        (by simp only [ne_eq, Fin.mk.injEq]; omega) x
      linarith [neg_abs_le (qform (bondOp (⟨m, by omega⟩ : Fin L) ⟨m + 1, by omega⟩ exchange +
        bondOp (⟨m + 1, by omega⟩ : Fin L) ⟨m + 1 + 1, by omega⟩ exchange) x)]
    have e1 : (m + 2) / 2 = m / 2 + 1 := by omega
    have e2 : (m + 2) % 2 = m % 2 := by omega
    rw [e1, e2]
    push_cast at ih ⊢
    linarith

/-- The first `m` bonds are at most `m`. -/
theorem qform_chainBonds_le {L : ℕ} (m : ℕ) (x : (Fin L → Fin 3) → ℂ) :
    qform (∑ j ∈ Finset.range m, chainBond L j) x ≤ m * nsq x := by
  have := qform_sum_le (Finset.range m) (b := fun _ ↦ 1) (fun j _ y ↦ qform_chainBond_le L j y) x
  simpa using this

/-- The bonds of `H_{n,g}` are the chain bonds. -/
theorem sum_bondOp_eq_chainBonds (n : ℕ) :
    ∑ j : Fin n, bondOp j.castSucc j.succ exchange = ∑ j ∈ Finset.range n, chainBond (n + 1) j := by
  rw [← Fin.sum_univ_eq_sum_range]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [chainBond, show (j : ℕ) + 1 < n + 1 by omega, ↓reduceDIte]
  rfl

/-- The right field `(3/5) R_g` is a real multiple `c Sᵅ` with `|c| ≤ 3/5`. -/
theorem rightField_eq {g : ℤ} (hg : g ∈ fields) : ∃ (c : ℝ) (α : Fin 3), |c| ≤ 3 / 5 ∧
    ((3 / 5 : ℝ) : ℂ) • rightField g = (c : ℂ) • S α := by
  simp only [fields, List.mem_cons, List.not_mem_nil, or_false] at hg
  rcases hg with rfl | rfl | rfl
  · refine ⟨3 / 5, 2, by norm_num [abs_of_pos], ?_⟩
    simp [rightField]
  · refine ⟨-(3 / 5), 2, by norm_num [abs_of_pos], ?_⟩
    simp [rightField]
  · exact ⟨3 / 5, 0, by norm_num [abs_of_pos], by simp [rightField]⟩

/-- **Energy bounds of the boundary chains** (the paper's `0 ≤ H_{n,g} + C_n ≤ (n - 1 + 6/5 +
C_n)`): `-C_{n+1} ≤ H_{n+1,g} ≤ n + 6/5` for `n ≥ 1` bonds, with
`10 C_{n+1} = c10 (n + 1)`. -/
theorem energy_bounds {g : ℤ} (hg : g ∈ fields) {n : ℕ} (hn : 1 ≤ n)
    (x : (Fin (n + 1) → Fin 3) → ℂ) :
    -((c10 (n + 1) : ℝ) / 10) * nsq x ≤ qform (BoundaryTraces.hamiltonian g (n + 1)) x ∧
      qform (BoundaryTraces.hamiltonian g (n + 1)) x ≤ (n + 6 / 5) * nsq x := by
  obtain ⟨c, α, hc, hR⟩ := rightField_eq hg
  have h0l : (0 : Fin (n + 1)) ≠ Fin.last n := by simp [Fin.ext_iff]; omega
  have hz : -((3 / 5 : ℝ) : ℂ) • S 2 = (((-(3 / 5)) : ℝ) : ℂ) • S 2 := by push_cast; rfl
  have h1 := abs_qform_siteOp h0l (-(3 / 5)) 2 x
  have h2 := abs_qform_siteOp h0l.symm c α x
  have hb1 := le_qform_chainBonds (L := n + 1) n le_rfl x
  have hb2 := qform_chainBonds_le (L := n + 1) n x
  have hnx := nsq_nonneg x
  simp only [BoundaryTraces.hamiltonian, sum_bondOp_eq_chainBonds, qform_add_left, hz, hR]
  rw [abs_le] at h1 h2
  have hc' : |c| * nsq x ≤ 3 / 5 * nsq x := mul_le_mul_of_nonneg_right hc hnx
  have h35 : |(-(3 / 5) : ℝ)| = 3 / 5 := by norm_num [abs_of_pos]
  rw [h35] at h1
  have hC : ((c10 (n + 1) : ℕ) : ℝ) / 10 = ((3 * (n / 2) + 2 * (n % 2) : ℕ) : ℝ) + 6 / 5 := by
    simp only [c10, Nat.add_sub_cancel]; push_cast; ring
  rw [hC]
  constructor
  · nlinarith [h1.1, h2.1]
  · nlinarith [h1.2, h2.2]

/-! ### The `S^z` basis -/

/-- The right field `R_g` in the `S^z` basis. -/
noncomputable def rightFieldZ (g : ℤ) : Matrix (Fin 3) (Fin 3) ℂ :=
  if g = 0 then spinZ 0 else (g : ℂ) • spinZ 2

/-- The boundary Hamiltonian `H_{n,g}` in the `S^z` basis. -/
noncomputable def hamiltonianZ (g : ℤ) : (n : ℕ) → Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℂ
  | 0 => 0
  | n + 1 => ∑ j : Fin n, bondOp j.castSucc j.succ (exchangeZ.map (↑)) +
      siteOp 0 (-((3 / 5 : ℝ) : ℂ) • spinZ 2) +
        siteOp (Fin.last n) (((3 / 5 : ℝ) : ℂ) • rightFieldZ g)

/-- Changing the basis of a one-site operator with the same unitary on every site. -/
theorem conj_siteOp_const {L : ℕ} (V : Matrix (Fin 3) (Fin 3) ℂ) (hV : Vᴴ * V = 1) (j : Fin L)
    (A : Matrix (Fin 3) (Fin 3) ℂ) :
    (piTensor fun _ : Fin L ↦ V)ᴴ * siteOp j A * piTensor (fun _ : Fin L ↦ V) =
      siteOp j (Vᴴ * A * V) := by
  have h := conj_siteOp (fun _ : Fin L ↦ Vᴴ) (fun _ ↦ by rw [conjTranspose_conjTranspose, hV]) j A
  simp only [conjTranspose_conjTranspose] at h
  rw [conjTranspose_piTensor, h]

/-- A scalar multiple of a spin matrix in the `S^z` basis. -/
theorem conj_smul_S (c : ℂ) (α : Fin 3) : szBasisᴴ * (c • S α) * szBasis = c • spinZ α := by
  rw [Matrix.mul_smul, Matrix.smul_mul, conj_S]

/-- **The boundary Hamiltonians in the `S^z` basis.** -/
theorem conj_hamiltonian (g : ℤ) (n : ℕ) :
    (szFrame n)ᴴ * BoundaryTraces.hamiltonian g n * szFrame n = hamiltonianZ g n := by
  cases n with
  | zero => simp [BoundaryTraces.hamiltonian, hamiltonianZ]
  | succ n =>
    simp only [BoundaryTraces.hamiltonian, hamiltonianZ, Matrix.mul_add, Matrix.add_mul,
      Finset.mul_sum, Finset.sum_mul, szFrame]
    congr 1
    · congr 1
      · refine Finset.sum_congr rfl fun j _ ↦ ?_
        rw [conj_bondOp_const szBasis szBasis_conjTranspose_mul_self
          (by simp [Fin.ext_iff, Fin.val_succ]), conj_exchange]
      · rw [conj_siteOp_const szBasis szBasis_conjTranspose_mul_self, conj_smul_S]
    · rw [conj_siteOp_const szBasis szBasis_conjTranspose_mul_self]
      congr 1
      unfold rightField rightFieldZ
      split_ifs
      · exact conj_smul_S _ 0
      · rw [smul_smul, conj_smul_S, smul_smul]

theorem szFrame_unitary (n : ℕ) : (szFrame n)ᴴ * szFrame n = 1 := szFrame_conjTranspose_mul_self n

/-- The trace `Z_{n,g}` in the `S^z` basis. -/
theorem trace_eq_gibbsTrace_Z (n : ℕ) (g : ℤ) :
    trace n g = gibbsTrace (hamiltonianZ g n + ((shift * n : ℝ) : ℂ) • 1) (49 / 4) := by
  unfold BoundaryTraces.trace
  rw [← gibbsTrace_conj (szFrame n) _ (szFrame_unitary n), Matrix.mul_add, Matrix.add_mul,
    conj_hamiltonian, Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_one, szFrame_unitary]

/-- The boundary Hamiltonians are Hermitian. -/
theorem hamiltonian_isHermitian {g : ℤ} (hg : g ∈ fields) (n : ℕ) :
    (BoundaryTraces.hamiltonian g n).IsHermitian := by
  cases n with
  | zero => exact isHermitian_zero
  | succ n =>
    obtain ⟨c, α, -, hR⟩ := rightField_eq hg
    have hS : ∀ (c : ℝ) (α : Fin 3), ((c : ℂ) • S α).IsHermitian := fun c α ↦ by
      change ((c : ℂ) • S α)ᴴ = _
      rw [conjTranspose_smul, (S_isHermitian α).eq, Complex.star_def, Complex.conj_ofReal]
    have hz : -((3 / 5 : ℝ) : ℂ) • S 2 = (((-(3 / 5)) : ℝ) : ℂ) • S 2 := by push_cast; rfl
    simp only [BoundaryTraces.hamiltonian, hz, hR]
    exact ((isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian _ _ exchange_isHermitian).add
      (siteOp_isHermitian _ (hS _ _))).add (siteOp_isHermitian _ (hS _ _))

/-- The boundary Hamiltonians in the `S^z` basis are Hermitian. -/
theorem hamiltonianZ_isHermitian {g : ℤ} (hg : g ∈ fields) (n : ℕ) :
    (hamiltonianZ g n).IsHermitian := by
  rw [← conj_hamiltonian]
  exact isHermitian_conjTranspose_mul_mul _ (hamiltonian_isHermitian hg n)

/-- `C_n = c10 n / 10`. -/
noncomputable def cShift (n : ℕ) : ℝ := (c10 n : ℝ) / 10

/-- The contraction `A = 1 - (H_{n,g} + C_n)/8` in the `S^z` basis (eq. (appth:C)). -/
noncomputable def contraction (g : ℤ) (n : ℕ) : Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℂ :=
  1 - ((1 / 8 : ℝ) : ℂ) • (hamiltonianZ g n + ((cShift n : ℝ) : ℂ) • 1)

theorem contraction_isHermitian {g : ℤ} (hg : g ∈ fields) (n : ℕ) :
    (contraction g n).IsHermitian := by
  have hc : ∀ c : ℝ, ((c : ℂ) • (1 : Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℂ)).IsHermitian :=
    fun c ↦ isHermitian_one.smul (Complex.conj_ofReal c)
  exact isHermitian_one.sub (((hamiltonianZ_isHermitian hg n).add (hc _)).smul
    (Complex.conj_ofReal _))

/-- **`A` is a contraction** for `2 ≤ n ≤ 6`: `0 ≤ H_{n,g} + C_n ≤ 16` (the paper's upper
bounds `5.4, 7.4, 10.4, 12.4, 15.4`). -/
theorem abs_qform_contraction {g : ℤ} (hg : g ∈ fields) {n : ℕ} (hn : n ∈ lengths)
    (x : (Fin n → Fin 3) → ℂ) : |qform (contraction g n) x| ≤ nsq x := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by simp [lengths] at hn; omega⟩
  have hm1 : 1 ≤ m := by simp [lengths] at hn; omega
  have hm5 : m ≤ 5 := by simp [lengths] at hn; omega
  -- `qform H^Z x = qform H (U x)` and `‖U x‖ = ‖x‖`.
  have hq : qform (hamiltonianZ g (m + 1)) x =
      qform (BoundaryTraces.hamiltonian g (m + 1)) (szFrame (m + 1) *ᵥ x) := by
    rw [qform_mulVec, conj_hamiltonian]
  have hn' : nsq (szFrame (m + 1) *ᵥ x) = nsq x := nsq_mulVec_of_unitary (szFrame_unitary _) x
  obtain ⟨hlo, hhi⟩ := energy_bounds hg hm1 (szFrame (m + 1) *ᵥ x)
  rw [hn', ← hq] at hlo hhi
  have hA : qform (contraction g (m + 1)) x = nsq x - 1 / 8 *
      (qform (hamiltonianZ g (m + 1)) x + cShift (m + 1) * nsq x) := by
    rw [contraction, qform_sub_left, qform_one,
      show qform (((1 / 8 : ℝ) : ℂ) • (hamiltonianZ g (m + 1) + ((cShift (m + 1) : ℝ) : ℂ) • 1)) x
        = 1 / 8 * qform (hamiltonianZ g (m + 1) + ((cShift (m + 1) : ℝ) : ℂ) • 1) x from
        qform_smul_left _ _ x, qform_add_left,
      show qform (((cShift (m + 1) : ℝ) : ℂ) • (1 : Matrix (Fin (m + 1) → Fin 3) _ ℂ)) x =
        cShift (m + 1) * qform 1 x from qform_smul_left _ _ x, qform_one]
  have hC : cShift (m + 1) + (m + 6 / 5) ≤ 16 := by
    unfold cShift c10
    interval_cases m <;> norm_num
  have hnx := nsq_nonneg x
  rw [hA, abs_le]
  unfold cShift at hC ⊢
  constructor <;> nlinarith

/-! ### Word numbers -/

/-- The word number `e σ = ∑ᵢ σᵢ 3ⁱ` (site `0` least significant). -/
abbrev wordNum (n : ℕ) : (Fin n → Fin 3) ≃ Fin (3 ^ n) := finFunctionFinEquiv

theorem wordNum_val {n : ℕ} (σ : Fin n → Fin 3) :
    ((wordNum n σ : ℕ) : ℤ) = ∑ i, ((σ i : ℕ) : ℤ) * 3 ^ (i : ℕ) := by
  rw [wordNum, finFunctionFinEquiv_apply]; push_cast; rfl

/-- The digits of a word number are the letters. -/
theorem digit_wordNum {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) :
    digit (wordNum n σ) j = σ j := by
  have := finFunctionFinEquiv_symm_apply_val (m := 3) (wordNum n σ) j
  rw [wordNum, Equiv.symm_apply_apply] at this
  exact this.symm

/-- Changing two letters changes a weighted letter sum by the weighted differences. -/
theorem sum_update₂ {n : ℕ} (w : Fin n → ℤ) (σ : Fin n → Fin 3) {j k : Fin n} (hjk : j ≠ k)
    (a b : Fin 3) :
    ∑ i, ((Function.update (Function.update σ j a) k b i : ℕ) : ℤ) * w i =
      ∑ i, ((σ i : ℕ) : ℤ) * w i + ((a : ℤ) - σ j) * w j + ((b : ℤ) - σ k) * w k := by
  have h : ∀ i, ((Function.update (Function.update σ j a) k b i : ℕ) : ℤ) * w i =
      ((σ i : ℕ) : ℤ) * w i + (if i = j then ((a : ℤ) - σ j) * w j else 0) +
        (if i = k then ((b : ℤ) - σ k) * w k else 0) := by
    intro i
    by_cases hk : i = k
    · subst hk; simp [Ne.symm hjk]; ring
    · by_cases hj : i = j
      · subst hj; simp [hk]; ring
      · simp [hk, hj]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

/-- Changing one letter changes a weighted letter sum by the weighted difference. -/
theorem sum_update {n : ℕ} (w : Fin n → ℤ) (σ : Fin n → Fin 3) (j : Fin n) (a : Fin 3) :
    ∑ i, ((Function.update σ j a i : ℕ) : ℤ) * w i =
      ∑ i, ((σ i : ℕ) : ℤ) * w i + ((a : ℤ) - σ j) * w j := by
  have h : ∀ i, ((Function.update σ j a i : ℕ) : ℤ) * w i =
      ((σ i : ℕ) : ℤ) * w i + (if i = j then ((a : ℤ) - σ j) * w j else 0) := by
    intro i
    by_cases hj : i = j
    · subst hj; simp; ring
    · simp [hj]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

/-- Changing two letters changes the word number by the differences of the place values. -/
theorem wordNum_update₂ {n : ℕ} (σ : Fin n → Fin 3) {j k : Fin n} (hjk : j ≠ k) (a b : Fin 3) :
    ((wordNum n (Function.update (Function.update σ j a) k b) : ℕ) : ℤ) =
      wordNum n σ + ((a : ℤ) - σ j) * 3 ^ (j : ℕ) + ((b : ℤ) - σ k) * 3 ^ (k : ℕ) := by
  rw [wordNum_val, wordNum_val, sum_update₂ (fun i ↦ 3 ^ (i : ℕ)) σ hjk]

/-- Changing one letter changes the word number by the difference of the place value. -/
theorem wordNum_update {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) (a : Fin 3) :
    ((wordNum n (Function.update σ j a) : ℕ) : ℤ) =
      wordNum n σ + ((a : ℤ) - σ j) * 3 ^ (j : ℕ) := by
  rw [wordNum_val, wordNum_val, sum_update (fun i ↦ 3 ^ (i : ℕ)) σ j]

/-- The parity `p_σ = (∑ᵢ σᵢ + n) mod 2` of the word, equal to `(∑ᵢ wᵢ) mod 2`. -/
def parity {n : ℕ} (σ : Fin n → Fin 3) : ℕ := (∑ i, (σ i : ℕ) + n) % 2

theorem list_sum_range {M : Type*} [AddCommMonoid M] (f : ℕ → M) (m : ℕ) :
    ((List.range m).map f).sum = ∑ i ∈ Finset.range m, f i := by
  induction m with
  | zero => simp
  | succ m ih => rw [List.sum_range_succ, ih, Finset.sum_range_succ]

theorem par_wordNum {n : ℕ} (σ : Fin n → Fin 3) : par n (wordNum n σ) = parity σ := by
  rw [par, parity, list_sum_range, ← Fin.sum_univ_eq_sum_range]
  simp only [digit_wordNum]

theorem letter_wordNum {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) :
    letter (wordNum n σ) j = mz (σ j) := by
  rw [letter, digit_wordNum, mz]

theorem parity_lt_two {n : ℕ} (σ : Fin n → Fin 3) : parity σ < 2 := Nat.mod_lt _ (by norm_num)

/-- Two letters changing in opposite directions keep the parity. -/
theorem parity_update₂ {n : ℕ} (σ : Fin n → Fin 3) {j k : Fin n} (hjk : j ≠ k) (a b : Fin 3)
    (hab : (a : ℕ) + b = σ j + σ k) :
    parity (Function.update (Function.update σ j a) k b) = parity σ := by
  have h := sum_update₂ (fun _ ↦ 1) σ hjk a b
  simp only [mul_one] at h
  have h' : (∑ i, (Function.update (Function.update σ j a) k b i : ℕ) : ℤ) = ∑ i, (σ i : ℕ) := by
    push_cast; rw [h]; have : ((a : ℕ) : ℤ) + b = σ j + σ k := by exact_mod_cast hab
    linarith
  unfold parity
  rw [show ∑ i, (Function.update (Function.update σ j a) k b i : ℕ) = ∑ i, (σ i : ℕ) by
    exact_mod_cast h']

/-- One letter changing by one flips the parity. -/
theorem parity_update {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) (a : Fin 3)
    (ha : (a : ℕ) = σ j + 1 ∨ (a : ℕ) + 1 = σ j) :
    parity (Function.update σ j a) = 1 - parity σ := by
  have h := sum_update (fun _ ↦ 1) σ j a
  simp only [mul_one] at h
  unfold parity
  have hs : (∑ i, (Function.update σ j a i : ℕ) : ℤ) = ∑ i, (σ i : ℕ) + (a - σ j) := by
    push_cast; rw [h]
  rcases ha with ha | ha
  · have : ∑ i, (Function.update σ j a i : ℕ) = ∑ i, (σ i : ℕ) + 1 := by
      have : ((a : ℕ) : ℤ) - σ j = 1 := by omega
      rw [this] at hs; exact_mod_cast hs
    rw [this]; omega
  · have h1 : 1 ≤ ∑ i, (σ i : ℕ) := le_trans (by omega)
      (Finset.single_le_sum (f := fun i ↦ (σ i : ℕ)) (fun _ _ ↦ Nat.zero_le _)
        (Finset.mem_univ j))
    have : ∑ i, (Function.update σ j a i : ℕ) + 1 = ∑ i, (σ i : ℕ) := by
      have : ((a : ℕ) : ℤ) - σ j = -1 := by omega
      rw [this] at hs; have : ((∑ i, (Function.update σ j a i : ℕ) : ℕ) : ℤ) + 1 =
        ∑ i, (σ i : ℕ) := by push_cast at hs ⊢; linarith
      exact_mod_cast this
    omega

/-! ### Hops between words -/

/-- The bond hop `σⱼ ↦ σⱼ - 1`, `σₖ ↦ σₖ + 1`, allowed when `σⱼ ≥ 1` and `σₖ ≤ 1`. -/
def hopUp {n : ℕ} (j k : Fin n) (σ : Fin n → Fin 3) : Option (Fin n → Fin 3) :=
  if h : 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1 then
    some (Function.update (Function.update σ j ⟨σ j - 1, by have := (σ j).isLt; omega⟩) k
      ⟨σ k + 1, by omega⟩)
  else none

/-- The end hop `σₗ ↦ σₗ + 1`, allowed when `σₗ ≤ 1`. -/
def endUp' {n : ℕ} (l : Fin n) (σ : Fin n → Fin 3) : Option (Fin n → Fin 3) :=
  if h : (σ l : ℕ) ≤ 1 then some (Function.update σ l ⟨σ l + 1, by omega⟩) else none

/-- The end hop `σₗ ↦ σₗ - 1`, allowed when `σₗ ≥ 1`. -/
def endDown' {n : ℕ} (l : Fin n) (σ : Fin n → Fin 3) : Option (Fin n → Fin 3) :=
  if h : 1 ≤ (σ l : ℕ) then
    some (Function.update σ l ⟨σ l - 1, by have := (σ l).isLt; omega⟩) else none

theorem hopUp_eq_some {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (σ τ : Fin n → Fin 3) :
    hopUp j k σ = some τ ↔ (∀ i, i ≠ j → i ≠ k → σ i = τ i) ∧
      (σ j : ℕ) = τ j + 1 ∧ (σ k : ℕ) + 1 = τ k := by
  unfold hopUp
  constructor
  · intro h
    split_ifs at h with hc
    cases h
    refine ⟨fun i hj hk ↦ by simp [Function.update_of_ne hk, Function.update_of_ne hj], ?_, ?_⟩
    · rw [Function.update_of_ne hjk, Function.update_self]; dsimp only; omega
    · simp
  · rintro ⟨hag, h1, h2⟩
    have hτk := (τ k).isLt
    rw [dite_eq_left ⟨by omega, by omega⟩]
    congr 1
    funext i
    by_cases hk : i = k
    · subst hk; rw [Function.update_self]; exact Fin.ext (by dsimp only; omega)
    · by_cases hj : i = j
      · subst hj; rw [Function.update_of_ne hk, Function.update_self]
        exact Fin.ext (by dsimp only; omega)
      · simp [Function.update_of_ne hk, Function.update_of_ne hj, hag i hj hk]

theorem endUp'_eq_some {n : ℕ} (l : Fin n) (σ τ : Fin n → Fin 3) :
    endUp' l σ = some τ ↔ (∀ i, i ≠ l → σ i = τ i) ∧ (σ l : ℕ) + 1 = τ l := by
  unfold endUp'
  constructor
  · intro h
    split_ifs at h with hc
    cases h
    exact ⟨fun i hi ↦ by simp [Function.update_of_ne hi], by simp⟩
  · rintro ⟨hag, h1⟩
    have := (τ l).isLt
    rw [dite_eq_left (by omega)]
    congr 1
    funext i
    by_cases hi : i = l
    · subst hi; rw [Function.update_self]; exact Fin.ext (by dsimp only; omega)
    · simp [Function.update_of_ne hi, hag i hi]

theorem endDown'_eq_some {n : ℕ} (l : Fin n) (σ τ : Fin n → Fin 3) :
    endDown' l σ = some τ ↔ (∀ i, i ≠ l → σ i = τ i) ∧ (σ l : ℕ) = τ l + 1 := by
  unfold endDown'
  constructor
  · intro h
    split_ifs at h with hc
    cases h
    exact ⟨fun i hi ↦ by simp [Function.update_of_ne hi], by simp; omega⟩
  · rintro ⟨hag, h1⟩
    rw [dite_eq_left (by omega)]
    congr 1
    funext i
    by_cases hi : i = l
    · subst hi; rw [Function.update_self]; exact Fin.ext (by dsimp only; omega)
    · simp [Function.update_of_ne hi, hag i hi]

/-- **The integer bond in the `S^z` basis, entry by entry**: the diagonal `mz mz` and one unit
entry for each of the two bond hops. -/
theorem bondOp_exchangeZ_apply {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (σ τ : Fin n → Fin 3) :
    bondOp j k exchangeZ σ τ = (if σ = τ then mz (σ j) * mz (σ k) else 0) +
      (if hopUp j k σ = some τ then 1 else 0) + (if hopUp k j σ = some τ then 1 else 0) := by
  simp only [bondOp, of_apply, exchangeZ]
  by_cases hag : ∀ i, i ≠ j → i ≠ k → σ i = τ i
  · rw [ite_eq_left hag]
    have e1 : (σ j, σ k) = (τ j, τ k) ↔ σ = τ := by
      constructor
      · intro h
        simp only [Prod.mk.injEq] at h
        funext i
        by_cases hj : i = j
        · subst hj; exact h.1
        · by_cases hk : i = k
          · subst hk; exact h.2
          · exact hag i hj hk
      · rintro rfl; rfl
    have e2 : ((σ j : ℕ) = τ j + 1 ∧ (σ k : ℕ) + 1 = τ k) ↔ hopUp j k σ = some τ := by
      rw [hopUp_eq_some hjk]; exact ⟨fun h ↦ ⟨hag, h⟩, fun h ↦ h.2⟩
    have e3 : ((σ j : ℕ) + 1 = τ j ∧ (σ k : ℕ) = τ k + 1) ↔ hopUp k j σ = some τ := by
      rw [hopUp_eq_some (Ne.symm hjk)]
      exact ⟨fun h ↦ ⟨fun i hk hj ↦ hag i hj hk, h.2, h.1⟩, fun h ↦ ⟨h.2.2, h.2.1⟩⟩
    simp only [e1, e2, e3]
  · rw [ite_eq_right hag]
    have h1 : σ ≠ τ := fun h ↦ hag fun i _ _ ↦ by rw [h]
    have h2 : hopUp j k σ ≠ some τ := fun h ↦ hag ((hopUp_eq_some hjk σ τ).1 h).1
    have h3 : hopUp k j σ ≠ some τ := fun h ↦
      hag fun i hj hk ↦ ((hopUp_eq_some (Ne.symm hjk) σ τ).1 h).1 i hk hj
    simp [h1, h2, h3]

/-- The field `c Sᶻ` on one site, entry by entry. -/
theorem siteOp_spinZ_two_apply {n : ℕ} (j : Fin n) (c : ℂ) (σ τ : Fin n → Fin 3) :
    siteOp j (c • spinZ 2) σ τ = if σ = τ then c * mz (σ j) else 0 := by
  rw [siteOp_apply]
  by_cases h : σ = τ
  · subst h
    simp only [implies_true, ite_true, Matrix.smul_apply, smul_eq_mul]
    congr 1
    generalize σ j = a
    fin_cases a <;> simp [spinZ, mz]
  · rw [ite_eq_right h]
    split_ifs with hag
    · have : σ j ≠ τ j := fun e ↦ h (funext fun i ↦ by
        by_cases hi : i = j
        · subst hi; exact e
        · exact hag i hi)
      simp only [Matrix.smul_apply, smul_eq_mul]
      have h0 : spinZ 2 (σ j) (τ j) = 0 := by
        revert this; generalize σ j = a; generalize τ j = b
        intro hab; fin_cases a <;> fin_cases b <;> simp_all [spinZ]
      rw [h0, mul_zero]
    · rfl

/-- The field `c Sˣ` on one site, entry by entry: the two end hops with `c/√2`. -/
theorem siteOp_spinZ_zero_apply {n : ℕ} (l : Fin n) (c : ℂ) (σ τ : Fin n → Fin 3) :
    siteOp l (c • spinZ 0) σ τ = c * invSqrtTwo *
      ((if endUp' l σ = some τ then 1 else 0) + (if endDown' l σ = some τ then 1 else 0)) := by
  rw [siteOp_apply]
  by_cases hag : ∀ i, i ≠ l → σ i = τ i
  · rw [ite_eq_left hag]
    have e1 : endUp' l σ = some τ ↔ (σ l : ℕ) + 1 = τ l := by
      rw [endUp'_eq_some]; exact ⟨fun h ↦ h.2, fun h ↦ ⟨hag, h⟩⟩
    have e2 : endDown' l σ = some τ ↔ (σ l : ℕ) = τ l + 1 := by
      rw [endDown'_eq_some]; exact ⟨fun h ↦ h.2, fun h ↦ ⟨hag, h⟩⟩
    simp only [e1, e2, Matrix.smul_apply, smul_eq_mul]
    rw [mul_assoc]
    congr 1
    generalize σ l = a; generalize τ l = b
    fin_cases a <;> fin_cases b <;> simp [spinZ]
  · rw [ite_eq_right hag]
    have h1 : endUp' l σ ≠ some τ := fun h ↦ hag ((endUp'_eq_some l σ τ).1 h).1
    have h2 : endDown' l σ ≠ some τ := fun h ↦ hag ((endDown'_eq_some l σ τ).1 h).1
    simp [h1, h2]

/-! ### The band entries in terms of hops -/

theorem getD_map_range {α : Type*} (f : ℕ → α) (d k : ℕ) (a : α) (hk : k < d) :
    ((List.range d).map f).getD k a = f k := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_range hk]; rfl

theorem le_wordNum {n : ℕ} (σ : Fin n → Fin 3) (j : Fin n) (hj : 1 ≤ (σ j : ℕ)) :
    3 ^ (j : ℕ) ≤ (wordNum n σ : ℕ) := by
  have h := wordNum_val σ
  have h1 : ((σ j : ℕ) : ℤ) * 3 ^ (j : ℕ) ≤ ∑ i, ((σ i : ℕ) : ℤ) * 3 ^ (i : ℕ) :=
    Finset.single_le_sum (f := fun i ↦ ((σ i : ℕ) : ℤ) * 3 ^ (i : ℕ))
      (fun i _ ↦ by positivity) (Finset.mem_univ j)
  have h2 : (3 : ℤ) ^ (j : ℕ) ≤ ((σ j : ℕ) : ℤ) * 3 ^ (j : ℕ) :=
    le_mul_of_one_le_left (by positivity) (by exact_mod_cast hj)
  have : ((3 ^ (j : ℕ) : ℕ) : ℤ) ≤ ((wordNum n σ : ℕ) : ℤ) := by push_cast; linarith
  exact_mod_cast this

theorem wordNum_eq_iff {n : ℕ} {σ τ : Fin n → Fin 3} :
    ((wordNum n σ : ℕ) : ℤ) = wordNum n τ ↔ σ = τ := by
  constructor
  · intro h; exact (wordNum n).injective (Fin.ext (by exact_mod_cast h))
  · rintro rfl; rfl

/-- The up band of a bond in terms of the hop. -/
theorem bondUp_term {n : ℕ} {j k : Fin n} (hjk : (k : ℕ) = j + 1) (σ τ : Fin n → Fin 3)
    (c : ℕ) : (if (wordNum n τ : ℕ) = wordNum n σ + 2 * 3 ^ (j : ℕ) then
      (if 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1 then c else 0) else 0) =
      if hopUp j k σ = some τ then c else 0 := by
  have hne : j ≠ k := by intro h; rw [h] at hjk; omega
  by_cases hc : 1 ≤ (σ j : ℕ) ∧ (σ k : ℕ) ≤ 1
  · have hh : hopUp j k σ = some _ := dite_eq_left hc
    rw [hh, ite_eq_left hc]
    set τ' := Function.update (Function.update σ j ⟨σ j - 1, by have := (σ j).isLt; omega⟩) k
      ⟨σ k + 1, by omega⟩
    have he : ((wordNum n τ' : ℕ) : ℤ) = wordNum n σ + 2 * 3 ^ (j : ℕ) := by
      rw [wordNum_update₂ σ hne, hjk]
      push_cast [Nat.cast_sub hc.1]
      ring
    have : (wordNum n τ : ℕ) = wordNum n σ + 2 * 3 ^ (j : ℕ) ↔ some τ' = some τ := by
      rw [Option.some_inj, ← wordNum_eq_iff (σ := τ') (τ := τ), he, ← Nat.cast_inj (R := ℤ)]
      push_cast
      constructor <;> intro h <;> linarith
    simp only [this]
  · have hh : hopUp j k σ = none := dite_eq_right hc
    rw [hh, ite_eq_right hc]
    simp

/-- The down band of a bond in terms of the hop. -/
theorem bondDown_term {n : ℕ} {j k : Fin n} (hjk : (k : ℕ) = j + 1) (σ τ : Fin n → Fin 3)
    (c : ℕ) : (if (wordNum n τ : ℕ) = wordNum n σ - 2 * 3 ^ (j : ℕ) then
      (if (σ j : ℕ) ≤ 1 ∧ 1 ≤ (σ k : ℕ) then c else 0) else 0) =
      if hopUp k j σ = some τ then c else 0 := by
  have hne : k ≠ j := by intro h; rw [h] at hjk; omega
  by_cases hc : 1 ≤ (σ k : ℕ) ∧ (σ j : ℕ) ≤ 1
  · have hh : hopUp k j σ = some _ := dite_eq_left hc
    rw [hh, ite_eq_left (show (σ j : ℕ) ≤ 1 ∧ 1 ≤ (σ k : ℕ) from ⟨hc.2, hc.1⟩)]
    set τ' := Function.update (Function.update σ k ⟨σ k - 1, by have := (σ k).isLt; omega⟩) j
      ⟨σ j + 1, by omega⟩
    have he : ((wordNum n τ' : ℕ) : ℤ) = wordNum n σ - 2 * 3 ^ (j : ℕ) := by
      rw [wordNum_update₂ σ hne, hjk]
      push_cast [Nat.cast_sub hc.1]
      ring
    have hle : 2 * 3 ^ (j : ℕ) ≤ (wordNum n σ : ℕ) := by
      have := le_wordNum σ k hc.1
      rw [hjk, pow_succ] at this
      omega
    have : (wordNum n τ : ℕ) = wordNum n σ - 2 * 3 ^ (j : ℕ) ↔ some τ' = some τ := by
      rw [Option.some_inj, ← wordNum_eq_iff (σ := τ') (τ := τ), he, ← Nat.cast_inj (R := ℤ),
        Nat.cast_sub hle]
      push_cast
      constructor <;> intro h <;> linarith
    simp only [this]
  · have hh : hopUp k j σ = none := dite_eq_right hc
    rw [hh, ite_eq_right (show ¬((σ j : ℕ) ≤ 1 ∧ 1 ≤ (σ k : ℕ)) from fun h ↦ hc ⟨h.2, h.1⟩)]
    simp

/-- The up end band in terms of the hop. -/
theorem endUp_term {n : ℕ} (l : Fin n) (σ τ : Fin n → Fin 3) (c : ℕ) :
    (if (wordNum n τ : ℕ) = wordNum n σ + 3 ^ (l : ℕ) then
      (if (σ l : ℕ) ≤ 1 then c else 0) else 0) = if endUp' l σ = some τ then c else 0 := by
  by_cases hc : (σ l : ℕ) ≤ 1
  · have hh : endUp' l σ = some _ := dite_eq_left hc
    rw [hh, ite_eq_left hc]
    set τ' := Function.update σ l ⟨σ l + 1, by omega⟩
    have he : ((wordNum n τ' : ℕ) : ℤ) = wordNum n σ + 3 ^ (l : ℕ) := by
      rw [wordNum_update]; push_cast; ring
    have : (wordNum n τ : ℕ) = wordNum n σ + 3 ^ (l : ℕ) ↔ some τ' = some τ := by
      rw [Option.some_inj, ← wordNum_eq_iff (σ := τ') (τ := τ), he, ← Nat.cast_inj (R := ℤ)]
      push_cast
      constructor <;> intro h <;> linarith
    simp only [this]
  · have hh : endUp' l σ = none := dite_eq_right hc
    rw [hh, ite_eq_right hc]
    simp

/-- The down end band in terms of the hop. -/
theorem endDown_term {n : ℕ} (l : Fin n) (σ τ : Fin n → Fin 3) (c : ℕ) :
    (if (wordNum n τ : ℕ) = wordNum n σ - 3 ^ (l : ℕ) then
      (if 1 ≤ (σ l : ℕ) then c else 0) else 0) = if endDown' l σ = some τ then c else 0 := by
  by_cases hc : 1 ≤ (σ l : ℕ)
  · have hh : endDown' l σ = some _ := dite_eq_left hc
    rw [hh, ite_eq_left hc]
    set τ' := Function.update σ l ⟨σ l - 1, by have := (σ l).isLt; omega⟩
    have he : ((wordNum n τ' : ℕ) : ℤ) = wordNum n σ - 3 ^ (l : ℕ) := by
      rw [wordNum_update]; push_cast [Nat.cast_sub hc]; ring
    have hle := le_wordNum σ l hc
    have : (wordNum n τ : ℕ) = wordNum n σ - 3 ^ (l : ℕ) ↔ some τ' = some τ := by
      rw [Option.some_inj, ← wordNum_eq_iff (σ := τ') (τ := τ), he, ← Nat.cast_inj (R := ℤ),
        Nat.cast_sub hle]
      push_cast
      constructor <;> intro h <;> linarith
    simp only [this]
  · have hh : endDown' l σ = none := dite_eq_right hc
    rw [hh, ite_eq_right hc]
    simp

theorem sum_map_flatMap {α β M : Type*} [AddCommMonoid M] (L : List α) (f : α → List β)
    (g : β → M) : ((L.flatMap f).map g).sum = (L.map fun a ↦ ((f a).map g).sum).sum := by
  induction L with
  | nil => simp
  | cons a L ih => simp [List.flatMap_cons, List.map_append, List.sum_append, ih]

/-- **The entries of the integer matrix in terms of hops.** -/
theorem band_entry (g : ℤ) (m : ℕ) (σ τ : Fin (m + 1) → Fin 3) :
    (band (m + 1) g).matrix (3 ^ (m + 1)) (wordNum (m + 1) σ) (wordNum (m + 1) τ) =
      (if σ = τ then eta (m + 1) g (wordNum (m + 1) σ) else 0) -
      ∑ j : Fin m, ((if hopUp j.castSucc j.succ σ = some τ then (10 : ℤ) else 0) +
        (if hopUp j.succ j.castSucc σ = some τ then (10 : ℤ) else 0)) -
      (if g = 0 then (if endUp' (Fin.last m) σ = some τ then ((3 * (1 + parity σ) : ℕ) : ℤ)
          else 0) + (if endDown' (Fin.last m) σ = some τ then ((3 * (1 + parity σ) : ℕ) : ℤ)
          else 0) else 0) := by
  have hk := (wordNum (m + 1) σ).isLt
  simp only [Band.matrix, band, BoundaryThermal.hops, (wordNum (m + 1)).injective.eq_iff,
    List.map_append, List.sum_append, sum_map_flatMap]
  rw [getD_map_range _ _ _ _ hk, getD_map_range _ _ _ _ hk, sub_sub]
  congr 1
  · split_ifs <;> simp
  have hm : m + 1 - 1 = m := by omega
  rw [hm, list_sum_range, ← Fin.sum_univ_eq_sum_range]
  congr 1
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero, bondUp,
      bondDown, Hop.nb, ite_true, Bool.false_eq_true, ite_false]
    rw [getD_map_range _ _ _ _ hk, getD_map_range _ _ _ _ hk]
    have d1 : digit (wordNum (m + 1) σ) j = σ j.castSucc := digit_wordNum σ j.castSucc
    have d2 : digit (wordNum (m + 1) σ) (j + 1) = σ j.succ := digit_wordNum σ j.succ
    rw [d1, d2]
    have h1 := congrArg (fun x : ℕ ↦ (x : ℤ))
      (bondUp_term (j := j.castSucc) (k := j.succ) (by simp) σ τ 10)
    have h2 := congrArg (fun x : ℕ ↦ (x : ℤ))
      (bondDown_term (j := j.castSucc) (k := j.succ) (by simp) σ τ 10)
    simp only [Fin.val_castSucc, Nat.cast_ite, Nat.cast_zero, Nat.cast_ofNat] at h1 h2 ⊢
    rw [h1, h2]
  · by_cases hg : g = 0
    · simp only [hg, ↓reduceIte, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil,
        add_zero, endUp, endDown, Hop.nb, Bool.false_eq_true, hm]
      rw [getD_map_range _ _ _ _ hk, getD_map_range _ _ _ _ hk]
      have d1 : digit (wordNum (m + 1) σ) m = σ (Fin.last m) := digit_wordNum σ (Fin.last m)
      have hp : par (m + 1) (wordNum (m + 1) σ) = parity σ := par_wordNum σ
      rw [d1, hp]
      have h1 := congrArg (fun x : ℕ ↦ (x : ℤ))
        (endUp_term (Fin.last m) σ τ (3 * (1 + parity σ)))
      have h2 := congrArg (fun x : ℕ ↦ (x : ℤ))
        (endDown_term (Fin.last m) σ τ (3 * (1 + parity σ)))
      simp only [Fin.val_last, Nat.cast_ite, Nat.cast_zero] at h1 h2 ⊢
      rw [h1, h2]
    · simp [hg]

/-! ### The diagonal similarity -/

/-- The diagonal similarity `J_σ = √(1 + p_σ)`. -/
noncomputable def jv {n : ℕ} (σ : Fin n → Fin 3) : ℝ := Real.sqrt (1 + parity σ)

theorem jv_pos {n : ℕ} (σ : Fin n → Fin 3) : 0 < jv σ := Real.sqrt_pos.2 (by positivity)

theorem parity_of_hopUp {n : ℕ} {j k : Fin n} (hjk : j ≠ k) {σ τ : Fin n → Fin 3}
    (h : hopUp j k σ = some τ) : parity σ = parity τ := by
  unfold hopUp at h
  split_ifs at h with hc
  cases h
  rw [parity_update₂ σ hjk]
  have := (σ j).isLt; simp only; omega

theorem parity_of_endUp' {n : ℕ} {l : Fin n} {σ τ : Fin n → Fin 3}
    (h : endUp' l σ = some τ) : parity τ = 1 - parity σ := by
  unfold endUp' at h
  split_ifs at h with hc
  cases h
  exact parity_update σ l _ (Or.inl rfl)

theorem parity_of_endDown' {n : ℕ} {l : Fin n} {σ τ : Fin n → Fin 3}
    (h : endDown' l σ = some τ) : parity τ = 1 - parity σ := by
  unfold endDown' at h
  split_ifs at h with hc
  cases h
  exact parity_update σ l _ (Or.inr (by simp; omega))

/-- Conjugating an entry that couples words of equal parity does not change it. -/
theorem jv_conj_ite {n : ℕ} {σ τ : Fin n → Fin 3} {c : Prop} [Decidable c]
    (h : c → parity σ = parity τ) (a : ℂ) :
    (jv σ : ℂ) * (if c then a else 0) * ((jv τ : ℂ))⁻¹ = if c then a else 0 := by
  split_ifs with hc
  · have hne : ((jv τ : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (jv_pos τ).ne'
    have heq : jv σ = jv τ := by rw [jv, jv, h hc]
    rw [heq]; field_simp
  · simp

theorem jv_flip_real {n : ℕ} {σ τ : Fin n → Fin 3} (hτ : parity τ = 1 - parity σ) :
    jv σ * (Real.sqrt 2)⁻¹ * (jv τ)⁻¹ = (1 + parity σ : ℕ) / 2 := by
  have hp := parity_lt_two σ
  have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
  have h2' : Real.sqrt 2 ≠ 0 := (Real.sqrt_pos.2 (by norm_num)).ne'
  rw [jv, jv, hτ]
  interval_cases hps : parity σ
  · norm_num
    field_simp
    linarith
  · norm_num

/-- Conjugating an end-hop entry `a/√2` between words of opposite parity gives
`(1 + p_σ) a / 2`. -/
theorem jv_conj_flip {n : ℕ} {σ τ : Fin n → Fin 3} {c : Prop} [Decidable c]
    (h : c → parity τ = 1 - parity σ) (a : ℂ) :
    (jv σ : ℂ) * (invSqrtTwo * if c then a else 0) * ((jv τ : ℂ))⁻¹ =
      if c then ((1 + parity σ : ℕ) : ℂ) / 2 * a else 0 := by
  split_ifs with hc
  · have e : (jv σ : ℂ) * (invSqrtTwo * a) * ((jv τ : ℂ))⁻¹ =
        ((jv σ * (Real.sqrt 2)⁻¹ * (jv τ)⁻¹ : ℝ) : ℂ) * a := by
      rw [invSqrtTwo]; push_cast; ring
    rw [e, jv_flip_real (h hc)]
    push_cast; ring
  · simp

/-- **The entries of `H_{n,g}` in the `S^z` basis.** -/
theorem hamiltonianZ_apply (g : ℤ) (m : ℕ) (σ τ : Fin (m + 1) → Fin 3) :
    hamiltonianZ g (m + 1) σ τ =
      ∑ j : Fin m, ((bondOp j.castSucc j.succ exchangeZ σ τ : ℤ) : ℂ) +
        (if σ = τ then -((3 / 5 : ℝ) : ℂ) * mz (σ 0) else 0) +
        (if g = 0 then ((3 / 5 : ℝ) : ℂ) * invSqrtTwo *
            ((if endUp' (Fin.last m) σ = some τ then 1 else 0) +
              (if endDown' (Fin.last m) σ = some τ then 1 else 0))
          else if σ = τ then ((3 / 5 : ℝ) : ℂ) * g * mz (σ (Fin.last m)) else 0) := by
  simp only [hamiltonianZ, Matrix.add_apply, Matrix.sum_apply]
  congr 1
  · congr 1
    · refine Finset.sum_congr rfl fun j _ ↦ ?_
      rw [← bondOp_map _ _ _ _ Int.cast_zero]; rfl
    · rw [siteOp_spinZ_two_apply]
  · unfold rightFieldZ
    by_cases hg : g = 0
    · simp only [hg, ↓reduceIte]
      rw [siteOp_spinZ_zero_apply]
    · simp only [hg, ↓reduceIte]
      rw [smul_smul, siteOp_spinZ_two_apply, mul_assoc]

theorem hopUp_ne_self {n : ℕ} {j k : Fin n} (hjk : j ≠ k) (σ : Fin n → Fin 3) :
    hopUp j k σ ≠ some σ := fun h ↦ by have := ((hopUp_eq_some hjk σ σ).1 h).2.1; omega

theorem endUp'_ne_self {n : ℕ} (l : Fin n) (σ : Fin n → Fin 3) : endUp' l σ ≠ some σ :=
  fun h ↦ by have := ((endUp'_eq_some l σ σ).1 h).2; omega

theorem endDown'_ne_self {n : ℕ} (l : Fin n) (σ : Fin n → Fin 3) : endDown' l σ ≠ some σ :=
  fun h ↦ by have := ((endDown'_eq_some l σ σ).1 h).2; omega

/-- The diagonal entry `η` in terms of the letters of the word. -/
theorem eta_wordNum (g : ℤ) (m : ℕ) (σ : Fin (m + 1) → Fin 3) :
    eta (m + 1) g (wordNum (m + 1) σ) = 80 - (c10 (m + 1) : ℤ) -
      10 * ∑ j : Fin m, mz (σ j.castSucc) * mz (σ j.succ) + 6 * mz (σ 0) -
        (if g = 0 then 0 else 6 * g * mz (σ (Fin.last m))) := by
  rw [eta, Nat.add_sub_cancel, list_sum_range, ← Fin.sum_univ_eq_sum_range]
  have h0 : letter (wordNum (m + 1) σ) 0 = mz (σ 0) := letter_wordNum σ 0
  have hl : letter (wordNum (m + 1) σ) m = mz (σ (Fin.last m)) := letter_wordNum σ (Fin.last m)
  have hj : ∀ j : Fin m, letter (wordNum (m + 1) σ) j * letter (wordNum (m + 1) σ) (j + 1) =
      mz (σ j.castSucc) * mz (σ j.succ) := fun j ↦ by
    rw [show (j : ℕ) = (j.castSucc : Fin (m + 1)) from rfl, letter_wordNum,
      show (j.castSucc : ℕ) + 1 = (j.succ : Fin (m + 1)) by simp, letter_wordNum]
  simp only [h0, hl, hj]

/-- **The integer matrix is `D = 80 J A J⁻¹`**, entry by entry (eq. (appth:diagonal) and the
hop entries `-10`, `-3(1 + p_w)`). -/
theorem similarity_apply (g : ℤ) (m : ℕ) (σ τ : Fin (m + 1) → Fin 3) :
    (((band (m + 1) g).matrix (3 ^ (m + 1)) (wordNum (m + 1) σ) (wordNum (m + 1) τ) : ℤ) : ℂ) =
      80 * ((jv σ : ℂ) * contraction g (m + 1) σ τ * ((jv τ : ℂ))⁻¹) := by
  rw [band_entry, contraction]
  simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.add_apply, one_apply, smul_eq_mul,
    hamiltonianZ_apply]
  have hjv : ((jv τ : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (jv_pos τ).ne'
  by_cases hστ : σ = τ
  · subst hστ
    have hb : ∀ j : Fin m, ((bondOp j.castSucc j.succ exchangeZ σ σ : ℤ) : ℂ) =
        (mz (σ j.castSucc) * mz (σ j.succ) : ℤ) := fun j ↦ by
      have hne : j.castSucc ≠ j.succ := by simp [Fin.ext_iff]
      rw [bondOp_exchangeZ_apply hne, ite_eq_left rfl, ite_eq_right (hopUp_ne_self hne σ),
        ite_eq_right (hopUp_ne_self hne.symm σ)]
      simp
    have hu : ∀ j : Fin m, hopUp j.castSucc j.succ σ ≠ some σ := fun j ↦
      hopUp_ne_self (by simp [Fin.ext_iff]) σ
    have hd : ∀ j : Fin m, hopUp j.succ j.castSucc σ ≠ some σ := fun j ↦
      hopUp_ne_self (by simp [Fin.ext_iff]) σ
    simp only [hb, hu, hd, endUp'_ne_self, endDown'_ne_self, ite_true, ite_false, add_zero,
      Finset.sum_const_zero, sub_zero, ite_self, eta_wordNum]
    field_simp
    unfold cShift
    push_cast
    split_ifs <;> ring
  · have hb : ∀ j : Fin m, ((bondOp j.castSucc j.succ exchangeZ σ τ : ℤ) : ℂ) =
        (if hopUp j.castSucc j.succ σ = some τ then 1 else 0) +
          (if hopUp j.succ j.castSucc σ = some τ then 1 else 0) := fun j ↦ by
      have hne : j.castSucc ≠ j.succ := by simp [Fin.ext_iff]
      rw [bondOp_exchangeZ_apply hne, ite_eq_right hστ]
      push_cast; ring
    simp only [hb, hστ, ite_false, mul_zero, add_zero, zero_sub, mul_neg, neg_mul]
    have hcu : ∀ j : Fin m, (jv σ : ℂ) * (if hopUp j.castSucc j.succ σ = some τ then 1 else 0) *
        ((jv τ : ℂ))⁻¹ = if hopUp j.castSucc j.succ σ = some τ then 1 else 0 := fun j ↦
      jv_conj_ite (parity_of_hopUp (by simp [Fin.ext_iff])) 1
    have hcd : ∀ j : Fin m, (jv σ : ℂ) * (if hopUp j.succ j.castSucc σ = some τ then 1 else 0) *
        ((jv τ : ℂ))⁻¹ = if hopUp j.succ j.castSucc σ = some τ then 1 else 0 := fun j ↦
      jv_conj_ite (parity_of_hopUp (by simp [Fin.ext_iff])) 1
    set r : ℂ := (jv σ : ℂ) * ((jv τ : ℂ))⁻¹ with hr
    have hconj : ∀ X : ℂ, (jv σ : ℂ) * X * ((jv τ : ℂ))⁻¹ = r * X := fun X ↦ by rw [hr]; ring
    have hS : r * ∑ j : Fin m, ((if hopUp j.castSucc j.succ σ = some τ then (1 : ℂ) else 0) +
        (if hopUp j.succ j.castSucc σ = some τ then 1 else 0)) =
        ∑ j : Fin m, ((if hopUp j.castSucc j.succ σ = some τ then (1 : ℂ) else 0) +
          (if hopUp j.succ j.castSucc σ = some τ then 1 else 0)) := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun j _ ↦ ?_
      rw [mul_add, ← hconj, ← hconj, hcu, hcd]
    have h10 : ∑ j : Fin m, ((if hopUp j.castSucc j.succ σ = some τ then (10 : ℂ) else 0) +
        (if hopUp j.succ j.castSucc σ = some τ then 10 else 0)) =
        10 * ∑ j : Fin m, ((if hopUp j.castSucc j.succ σ = some τ then (1 : ℂ) else 0) +
          (if hopUp j.succ j.castSucc σ = some τ then 1 else 0)) := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun j _ ↦ ?_
      split_ifs <;> ring
    push_cast
    rw [h10, hconj]
    by_cases hg : g = 0
    · simp only [hg, ↓reduceIte]
      have hfu := jv_conj_flip (σ := σ) (τ := τ) (c := endUp' (Fin.last m) σ = some τ)
        parity_of_endUp' 1
      have hfd := jv_conj_flip (σ := σ) (τ := τ) (c := endDown' (Fin.last m) σ = some τ)
        parity_of_endDown' 1
      rw [hconj] at hfu hfd
      push_cast at hfu hfd
      split_ifs at hfu hfd ⊢ <;> linear_combination (10 : ℂ) * hS + 6 * hfu + 6 * hfd
    · simp only [hg, ↓reduceIte]
      linear_combination (10 : ℂ) * hS

/-! ### Rational enclosures of the exponential -/

/-- The Taylor polynomial `∑_{j ≤ N} yʲ/j!`. -/
def taylorQ (y : ℚ) (N : ℕ) : ℚ := ((List.range (N + 1)).map fun j ↦ y ^ j / (j.factorial : ℚ)).sum

/-- A rational lower bound of `e^y` for `|y| ≤ 20`. -/
def expLo (y : ℚ) : ℚ :=
  if 0 ≤ y then taylorQ y 40 else 1 / (taylorQ (-y) 40 + (-y) ^ 40 / (Nat.factorial 40 : ℚ))

/-- A rational upper bound of `e^y` for `|y| ≤ 20`. -/
def expHi (y : ℚ) : ℚ :=
  if 0 ≤ y then taylorQ y 40 + y ^ 40 / (Nat.factorial 40 : ℚ) else 1 / taylorQ (-y) 40

theorem taylorQ_cast (y : ℚ) (N : ℕ) :
    (taylorQ y N : ℝ) = ∑ j ∈ Finset.range (N + 1), (y : ℝ) ^ j / (j.factorial : ℝ) := by
  rw [taylorQ, list_sum_range]; push_cast; rfl

theorem exp_mem_taylor {y : ℚ} (hy0 : 0 ≤ y) (hy : y ≤ 20) :
    (taylorQ y 40 : ℝ) ≤ Real.exp y ∧
      Real.exp y ≤ (taylorQ y 40 : ℝ) + (y : ℝ) ^ 40 / (Nat.factorial 40 : ℝ) := by
  rw [taylorQ_cast]
  exact Trace.exp_mem_Icc_taylor (x := (y : ℝ)) (by exact_mod_cast hy0) (N := 40)
    (by have : (y : ℝ) ≤ 20 := by exact_mod_cast hy
        push_cast; linarith)

/-- **The rational enclosure of the exponential.** -/
theorem exp_mem_expLo_expHi {y : ℚ} (hy : |y| ≤ 20) :
    (expLo y : ℝ) ≤ Real.exp y ∧ Real.exp y ≤ expHi y := by
  rw [abs_le] at hy
  unfold expLo expHi
  split_ifs with h
  · have := exp_mem_taylor h hy.2
    push_cast; exact this
  · have hn : 0 ≤ -y := by linarith
    obtain ⟨h1, h2⟩ := exp_mem_taylor hn (by linarith)
    have hpos : (0 : ℝ) < taylorQ (-y) 40 := by
      rw [taylorQ_cast]
      refine Finset.sum_pos' (fun j _ ↦ by positivity) ⟨0, by simp, by simp⟩
    have he : Real.exp (y : ℝ) = (Real.exp (-(y : ℝ)))⁻¹ := by
      rw [Real.exp_neg, inv_inv]
    push_cast at h1 h2 ⊢
    rw [he]
    constructor
    · rw [one_div]; exact inv_anti₀ (Real.exp_pos _) h2
    · rw [one_div]; exact inv_anti₀ hpos h1

/-! ### From the integer to the interval -/

/-- **The enclosure of `Z` from a square-root error bound.** If `e^y max(0, √W - ε)² ≤ Z ≤
e^y (√W + ε)²`, `L ≤ W ≤ U`, `√W ≤ r` and `s ≤ e^y ≤ s'`, then
`s (L - 2εr) ≤ Z ≤ s' (U + 2εr + ε²)` (eq. (appth:partition-interval)). -/
theorem interval_of_bounds {Z W ey ε L U r s s' : ℝ} (hW : 0 ≤ W) (hL : L ≤ W) (hU : W ≤ U)
    (hr : Real.sqrt W ≤ r) (hε : 0 ≤ ε) (hs : s ≤ ey) (hs' : ey ≤ s') (hs0 : 0 ≤ s)
    (hL' : 0 ≤ L - 2 * ε * r) (hlo : ey * max 0 (Real.sqrt W - ε) ^ 2 ≤ Z)
    (hhi : Z ≤ ey * (Real.sqrt W + ε) ^ 2) :
    s * (L - 2 * ε * r) ≤ Z ∧ Z ≤ s' * (U + 2 * ε * r + ε ^ 2) := by
  have hsq := Real.sq_sqrt hW
  have hrt := Real.sqrt_nonneg W
  have hey : 0 ≤ ey := hs0.trans hs
  constructor
  · -- `max(0, √W - ε)² ≥ W - 2ε√W ≥ L - 2εr`.
    have h1 : W - 2 * ε * Real.sqrt W ≤ max 0 (Real.sqrt W - ε) ^ 2 := by
      rcases le_total (Real.sqrt W - ε) 0 with h | h
      · rw [max_eq_left h]; nlinarith
      · rw [max_eq_right h]; nlinarith
    have h2 : L - 2 * ε * r ≤ W - 2 * ε * Real.sqrt W := by nlinarith
    calc s * (L - 2 * ε * r) ≤ ey * (L - 2 * ε * r) := mul_le_mul_of_nonneg_right hs hL'
      _ ≤ ey * max 0 (Real.sqrt W - ε) ^ 2 :=
          mul_le_mul_of_nonneg_left (h2.trans h1) hey
      _ ≤ Z := hlo
  · have h1 : (Real.sqrt W + ε) ^ 2 ≤ U + 2 * ε * r + ε ^ 2 := by nlinarith
    have h0 : 0 ≤ U + 2 * ε * r + ε ^ 2 := by nlinarith
    calc Z ≤ ey * (Real.sqrt W + ε) ^ 2 := hhi
      _ ≤ ey * (U + 2 * ε * r + ε ^ 2) := mul_le_mul_of_nonneg_left h1 hey
      _ ≤ s' * (U + 2 * ε * r + ε ^ 2) := mul_le_mul_of_nonneg_right hs' h0

/-! ### The integer matrix on words -/

/-- The integer matrix `D = 80 J A J⁻¹` on words, the banded matrix relabelled by the word
numbers. -/
def matrixD (n : ℕ) (g : ℤ) : Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℤ :=
  ((band n g).matrix (3 ^ n)).submatrix (wordNum n) (wordNum n)

/-- **`D = 80 J A J⁻¹`.** -/
theorem matrixD_similar (g : ℤ) (m : ℕ) :
    (matrixD (m + 1) g).map ((↑) : ℤ → ℂ) = ((80 : ℕ) : ℂ) • (diagonal (fun σ ↦ (jv σ : ℂ)) *
      contraction g (m + 1) * diagonal (fun σ ↦ ((jv σ)⁻¹ : ℂ))) := by
  ext σ τ
  simp only [matrixD, Matrix.map_apply, submatrix_apply, Matrix.smul_apply, diagonal_mul,
    mul_diagonal, smul_eq_mul]
  rw [similarity_apply]
  push_cast
  ring

/-- The Horner coefficient `tᵢ = ⌊Q sᵢ / S_K⌋`. -/
def coeffFun (i : ℕ) : ℤ := ((2 ^ 50 * poissonNum i / poissonDen : ℕ) : ℤ)

theorem coeffs_eq : coeffs.map ((↑) : ℕ → ℤ) = (List.range (192 + 1)).map coeffFun := by
  rw [coeffs, List.map_map]; rfl

theorem poissonDen_pos : 0 < poissonDen := by
  unfold poissonDen
  rw [list_sum_range]
  exact Finset.sum_pos' (fun i _ ↦ Nat.zero_le _) ⟨0, by simp, Nat.mul_pos (by positivity)
    (Nat.div_pos (Nat.factorial_le (Nat.zero_le _)) (Nat.factorial_pos 0))⟩

/-- `tᵢ` is within `1` of `Q pᵢ`, `pᵢ` the Poisson weight of mean `49` conditioned to
`i ≤ 192`. -/
theorem abs_coeffFun_sub_le {k : ℕ} (hk : k ≤ 192) :
    |(coeffFun k : ℝ) - 2 ^ 50 * (Trace.poissonWeight 49 192 k : ℝ)| ≤ 1 := by
  -- `pₖ = sₖ / S_K` with the natural numerators `sₖ = 49ᵏ (192!/k!)`.
  have hnum : ∀ l ≤ 192, ((poissonNum l : ℕ) : ℚ) = 49 ^ l * (Nat.factorial 192 : ℚ) /
      (Nat.factorial l : ℚ) := fun l hl ↦ by
    rw [poissonNum, Nat.cast_mul, Nat.cast_div (Nat.factorial_dvd_factorial hl)
      (by positivity)]
    push_cast; ring
  have hden : ((poissonDen : ℕ) : ℚ) = ∑ l ∈ Finset.range (192 + 1),
      49 ^ l * (Nat.factorial 192 : ℚ) / (Nat.factorial l : ℚ) := by
    rw [poissonDen, list_sum_range, Nat.cast_sum]
    exact Finset.sum_congr rfl fun l hl ↦ hnum l (by simp at hl; omega)
  have hp : (Trace.poissonWeight 49 192 k : ℝ) = ((poissonNum k : ℕ) : ℝ) / (poissonDen : ℕ) := by
    rw [Trace.poissonWeight_eq, ← hden, ← hnum k hk]; push_cast; rfl
  have hd : (0 : ℝ) < (poissonDen : ℕ) := by exact_mod_cast poissonDen_pos
  rw [hp, coeffFun, Int.cast_natCast]
  set a := 2 ^ 50 * poissonNum k
  have h1 : ((a / poissonDen : ℕ) : ℝ) ≤ (a : ℝ) / (poissonDen : ℕ) := by
    rw [le_div_iff₀ hd]; exact_mod_cast Nat.div_mul_le_self a poissonDen
  have h2 : (a : ℝ) / (poissonDen : ℕ) < ((a / poissonDen : ℕ) : ℝ) + 1 := by
    rw [div_lt_iff₀ hd]
    have := Nat.lt_div_mul_add (a := a) poissonDen_pos
    have : (a : ℝ) < ((a / poissonDen : ℕ) : ℝ) * (poissonDen : ℕ) + (poissonDen : ℕ) := by
      exact_mod_cast this
    linarith
  have ha : (2 : ℝ) ^ 50 * ((poissonNum k : ℕ) : ℝ) / (poissonDen : ℕ) =
      (a : ℝ) / (poissonDen : ℕ) := by simp only [a, Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat]
  rw [← mul_div_assoc, ha, abs_le]
  constructor <;> linarith

/-! ### The enclosure of one trace -/

/-- The half-step Gibbs factor is the Poisson target `exp (49 (A - 1))`:
`exp[-(b₀/2)(H + C_n)] = e^{-49} exp(49 A)`. -/
theorem gibbsHalf_eq (g : ℤ) (n : ℕ) :
    Trace.gibbsHalf (hamiltonianZ g n + ((shift * n : ℝ) : ℂ) • 1) (49 / 4)
        (cShift n - shift * n) = NormedSpace.exp (((49 : ℚ) : ℂ) • (contraction g n - 1)) := by
  rw [Trace.gibbsHalf_add_smul_one, show shift * n + (cShift n - shift * n) = cShift n by ring,
    Trace.gibbsHalf, contraction, sub_sub_cancel_left, smul_neg, smul_smul]
  congr 1
  rw [← neg_smul]
  congr 1
  push_cast; ring

theorem jv_le (n : ℕ) (σ : Fin n → Fin 3) : jv σ ≤ 3 / 2 := by
  rw [jv, Real.sqrt_le_left (by norm_num)]
  have := parity_lt_two σ
  have : (parity σ : ℝ) ≤ 1 := by exact_mod_cast (by omega : parity σ ≤ 1)
  linarith

theorem one_le_jv (n : ℕ) (σ : Fin n → Fin 3) : 1 ≤ jv σ := by
  rw [jv, Real.le_sqrt (by norm_num) (by positivity)]
  have : (0 : ℝ) ≤ parity σ := Nat.cast_nonneg _
  linarith

theorem inv_jv_le_one (n : ℕ) (σ : Fin n → Fin 3) : (jv σ)⁻¹ ≤ 1 :=
  inv_le_one_of_one_le₀ (one_le_jv n σ)

theorem sqrt_card_le {n : ℕ} (hn : n ∈ lengths) :
    Real.sqrt (Fintype.card (Fin n → Fin 3)) ≤ 27 := by
  rw [Fintype.card_fun, Fintype.card_fin, Fintype.card_fin, Real.sqrt_le_left (by norm_num)]
  have : 3 ^ n ≤ 729 := by
    simp only [lengths, List.mem_cons, List.not_mem_nil, or_false] at hn
    rcases hn with rfl | rfl | rfl | rfl | rfl <;> norm_num
  exact_mod_cast this

/-- The weights `(J_r / J_w)² = (2 - p_w)(1 + p_r)/2`. -/
theorem jv_div_sq {n : ℕ} (r w : Fin n → Fin 3) :
    (jv r / jv w) ^ 2 = ((2 - parity w : ℕ) : ℝ) * ((1 + parity r : ℕ) : ℝ) / 2 := by
  rw [div_pow, jv, jv, Real.sq_sqrt (by positivity), Real.sq_sqrt (by positivity)]
  have := parity_lt_two w
  interval_cases parity w <;> push_cast <;> ring

/-- The exponent `b₀ (C_n - a n)` as a rational number. -/
def expo (n : ℕ) : ℚ := 49 / 4 * ((c10 n : ℚ) / 10 - 700741 / 500000 * n)

theorem expo_cast (n : ℕ) : ((expo n : ℚ) : ℝ) = 49 / 4 * (cShift n - shift * n) := by
  rw [expo, cShift, shift]; push_cast; ring

/-- The column recursion of the banded matrix computes the Horner columns of `D` on words. -/
theorem hornerFloor_matrixD {n : ℕ} {g : ℤ} (hD : (band n g).InRange (3 ^ n))
    (r w : Fin n → Fin 3) :
    Trace.hornerFloor (matrixD n g) 80 coeffFun r 192 w =
      specCol (band n g) 80 (wordNum n r) (coeffs.map (↑)) (wordNum n w) := by
  unfold matrixD
  rw [hornerFloor_submatrix, hornerFloor_eq_specCol hD, coeffs_eq]

/-- A sum over words is a sum over word numbers. -/
theorem sum_words {n : ℕ} {M : Type*} [AddCommMonoid M] (f : ℕ → M) :
    ∑ σ : Fin n → Fin 3, f (wordNum n σ) = ∑ k ∈ Finset.range (3 ^ n), f k := by
  rw [(wordNum n).sum_comp (fun k ↦ f k), Fin.sum_univ_eq_sum_range]

/-- **`2 W = N`**: the weighted Frobenius norm of the computed columns is half the integer
numerator of eq. (appth:quotient). -/
theorem weighted_sum_eq {n : ℕ} {g : ℤ} (hchk : (band n g).check params
      ((band n g).cst params (3 ^ n)) (3 ^ n) = true) {ys : List ℕ}
    (hrun : runBoundary n g = some ys) :
    ∑ r, ∑ w, (jv r / jv w) ^ 2 * (Trace.hornerFloor (matrixD n g) 80 coeffFun r 192 w : ℝ) ^ 2
      = (numerator params (rowWeights n) (colWeights n) ys : ℝ) / 2 := by
  obtain ⟨hD, hR⟩ := Band.check_spec hchk
  have hP : params.Ok := by decide +kernel
  have hcs : ∀ t ∈ coeffs, t ≤ params.B := by decide +kernel
  have hN := run_numerator hP hD hR hcs hrun (u := rowWeights n) (v := colWeights n)
    (by simp [rowWeights]) (by simp [colWeights])
  simp only [zero_add] at hN
  have hN' : (numerator params (rowWeights n) (colWeights n) ys : ℝ) =
      ((numerator params (rowWeights n) (colWeights n) ys : ℤ) : ℝ) := by push_cast; rfl
  rw [hN', hN]
  simp only [hornerFloor_matrixD hD, jv_div_sq, ← par_wordNum]
  set F : ℕ → ℕ → ℝ := fun a b ↦ ((2 - par n b : ℕ) : ℝ) * ((1 + par n a : ℕ) : ℝ) / 2 *
    (specCol (band n g) 80 a (coeffs.map (↑)) b : ℝ) ^ 2 with hF
  have hL : ∑ r : Fin n → Fin 3, ∑ w : Fin n → Fin 3, F (wordNum n r) (wordNum n w) =
      ∑ a ∈ Finset.range (3 ^ n), ∑ b ∈ Finset.range (3 ^ n), F a b := by
    rw [← sum_words (fun a ↦ ∑ b ∈ Finset.range (3 ^ n), F a b)]
    exact Finset.sum_congr rfl fun r _ ↦ sum_words (fun b ↦ F (wordNum n r) b)
  rw [hL, Finset.sum_comm]
  push_cast
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl fun k hk ↦ ?_
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl fun r hr ↦ ?_
  have hk' := Finset.mem_range.1 hk
  have hr' := Finset.mem_range.1 hr
  simp only [hF, rowWeights, colWeights, getD_map_range _ _ _ _ hk', getD_map_range _ _ _ _ hr']
  have h2 : par n k ≤ 2 := (Nat.mod_lt _ (by norm_num)).le
  push_cast [h2]
  change _ = _ * _ * (specCol (band n g) 80 r (coeffs.map (↑)) k : ℝ) ^ 2 / 2
  ring

/-- **One trace from one successful run.** With `W = N / (2 Q²)` the paper's `W_{n,g}` and
`ε = 2·10⁻¹⁰` (Lemma appth:matrix-error),
`e^{b₀(C_n - a n)} max(0, √W - ε)² ≤ Z_{n,g} ≤ e^{b₀(C_n - a n)} (√W + ε)²`
(eq. (appth:partition-error)). -/
theorem trace_bounds_of_run {n : ℕ} (hn : n ∈ lengths) {g : ℤ} (hg : g ∈ fields)
    (hchk : (band n g).check params ((band n g).cst params (3 ^ n)) (3 ^ n) = true)
    {ys : List ℕ} (hrun : runBoundary n g = some ys) :
    Real.exp (expo n) * max 0 (Real.sqrt ((numerator params (rowWeights n) (colWeights n) ys : ℝ)
        / (2 * (2 ^ 50) ^ 2)) - 2 / 10 ^ 10) ^ 2 ≤ trace n g ∧
      trace n g ≤ Real.exp (expo n) * (Real.sqrt ((numerator params (rowWeights n)
        (colWeights n) ys : ℝ) / (2 * (2 ^ 50) ^ 2)) + 2 / 10 ^ 10) ^ 2 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by simp [lengths] at hn; omega⟩
  have hH : (hamiltonianZ g (m + 1) + ((shift * (m + 1 : ℕ) : ℝ) : ℂ) • 1).IsHermitian :=
    (hamiltonianZ_isHermitian hg _).add (isHermitian_one.smul (Complex.conj_ofReal _))
  have hε : Real.sqrt (Fintype.card (Fin (m + 1) → Fin 3)) *
      (3 / 2 * ((192 + 1) * (27 + 1)) + 2 ^ 50 * Trace.poissonError ((49 : ℚ) : ℝ) 192) ≤
        2 / 10 ^ 10 * 2 ^ 50 := by
    have h1 := sqrt_card_le hn
    have h2 := Trace.poissonError_le_three_pow_div 49 192
    have h3 := Trace.poissonError_nonneg (m := ((49 : ℚ) : ℝ)) (by norm_num) 192
    push_cast at h2 h3 ⊢
    have h4 : (2 : ℝ) ^ 50 * Trace.poissonError 49 192 ≤ 2 ^ 50 * (3 ^ 49 / 2 ^ 192) :=
      mul_le_mul_of_nonneg_left h2 (by positivity)
    have h5 : (0 : ℝ) ≤ 3 / 2 * ((192 + 1) * (27 + 1)) + 2 ^ 50 * Trace.poissonError 49 192 := by
      positivity
    calc Real.sqrt (Fintype.card (Fin (m + 1) → Fin 3)) *
          (3 / 2 * ((192 + 1) * (27 + 1)) + 2 ^ 50 * Trace.poissonError 49 192)
        ≤ 27 * (3 / 2 * ((192 + 1) * (27 + 1)) + 2 ^ 50 * (3 ^ 49 / 2 ^ 192)) := by
          gcongr
      _ ≤ 2 / 10 ^ 10 * 2 ^ 50 := by norm_num
  have hδ : Real.sqrt (∑ σ : Fin (m + 1) → Fin 3, (jv σ)⁻¹ ^ 2) ≤ 27 := by
    refine le_trans (Real.sqrt_le_sqrt ?_) (sqrt_card_le hn)
    calc ∑ σ : Fin (m + 1) → Fin 3, (jv σ)⁻¹ ^ 2 ≤ ∑ _σ : Fin (m + 1) → Fin 3, (1 : ℝ) :=
          Finset.sum_le_sum fun σ _ ↦ by
            have := inv_jv_le_one _ σ
            have : 0 ≤ (jv σ)⁻¹ := inv_nonneg.2 (jv_pos σ).le
            nlinarith
      _ = Fintype.card (Fin (m + 1) → Fin 3) := by simp
  have henc := Trace.gibbsTrace_enclosure_weighted hH (49 / 4)
    (cShift (m + 1) - shift * (m + 1 : ℕ)) (contraction_isHermitian hg _)
    (abs_qform_contraction hg hn) (m := 49) (by norm_num)
    (gibbsHalf_eq g (m + 1)) jv_pos (jmax := 3 / 2) (by norm_num) (jv_le _)
    (D := matrixD (m + 1) g) (s := 80) (by norm_num) (matrixD_similar g m) 192
    (Q := 2 ^ 50) (by positivity) coeffFun (fun k hk ↦ abs_coeffFun_sub_le hk) hδ
    (by norm_num) (inv_jv_le_one _) hε
  rw [weighted_sum_eq hchk hrun] at henc
  rw [trace_eq_gibbsTrace_Z, expo_cast]
  set N : ℝ := (numerator params (rowWeights (m + 1)) (colWeights (m + 1)) ys : ℝ)
  have hN : 0 ≤ N := Nat.cast_nonneg _
  -- `√(N/2)/Q = √(N/(2Q²))` and `ε / Q = 2·10⁻¹⁰`.
  have hQ : (0 : ℝ) < 2 ^ 50 := by positivity
  have hs : Real.sqrt (N / 2) = 2 ^ 50 * Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) := by
    have : Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) = Real.sqrt (N / 2) / 2 ^ 50 := by
      rw [show N / (2 * (2 ^ 50) ^ 2) = (N / 2) / (2 ^ 50) ^ 2 by ring,
        Real.sqrt_div (by positivity), Real.sqrt_sq hQ.le]
    rw [this]; field_simp
  have hcast : ((shift * (m + 1 : ℕ) : ℝ)) = shift * ((m + 1 : ℕ) : ℝ) := rfl
  rw [hs] at henc
  obtain ⟨h1, h2⟩ := henc
  constructor
  · refine le_trans (le_of_eq ?_) h1
    rw [show 2 ^ 50 * Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) - 2 / 10 ^ 10 * 2 ^ 50 =
      2 ^ 50 * (Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) - 2 / 10 ^ 10) by ring]
    rw [show max 0 ((2 : ℝ) ^ 50 * (Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) - 2 / 10 ^ 10)) =
      2 ^ 50 * max 0 (Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) - 2 / 10 ^ 10) by
        rw [mul_max_of_nonneg _ _ hQ.le, mul_zero]]
    field_simp
  · refine le_trans h2 (le_of_eq ?_)
    rw [show 2 ^ 50 * Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) + 2 / 10 ^ 10 * 2 ^ 50 =
      2 ^ 50 * (Real.sqrt (N / (2 * (2 ^ 50) ^ 2)) + 2 / 10 ^ 10) by ring]
    field_simp

/-! ### The table -/

/-- The rational checks of one table entry: from the printed integer `k = k_{n,g}`, with
`L = k/10¹⁶`, `U = (k + 1)/10¹⁶`, `r = (⌊√(k+1)⌋ + 1)/10⁸`, `ε = 2·10⁻¹⁰` (eq.
(appth:interval-data)), the enclosure `traceLo ≤ e̲ (L - 2εr)` and `ē (U + 2εr + ε²) ≤ traceHi`
with `e̲ ≤ e^{b₀(C_n - a n)} ≤ ē`. -/
def entryCheck (n : ℕ) (g : ℤ) : Bool :=
  let k := paperK n g
  let L : ℚ := k / 10 ^ 16
  let U : ℚ := (k + 1) / 10 ^ 16
  let r : ℚ := ((Nat.sqrt (k + 1) + 1 : ℕ) : ℚ) / 10 ^ 8
  let ε : ℚ := 2 / 10 ^ 10
  decide (|expo n| ≤ 20) && decide (0 ≤ L - 2 * ε * r) && decide (0 ≤ expLo (expo n)) &&
    decide (traceLo n g ≤ expLo (expo n) * (L - 2 * ε * r)) &&
    decide (expHi (expo n) * (U + 2 * ε * r + ε ^ 2) ≤ traceHi n g)

/-- All fifteen rational checks hold. -/
theorem entryCheck_all : ∀ n ∈ lengths, ∀ g ∈ fields, entryCheck n g = true := by
  decide +kernel

/-- All fifteen banded matrices satisfy the band conditions. -/
theorem bandCheck_all : ∀ n ∈ lengths, ∀ g ∈ fields,
    (band n g).check params ((band n g).cst params (3 ^ n)) (3 ^ n) = true := by
  decide +kernel

/-- **The fifteen thermal traces from the integer table.** If the packed integer recurrence
reproduces the paper's printed integers `k_{n,g}` (`BoundaryTotals`), the traces satisfy the
rational enclosures of `BoundaryTraceTable`. -/
theorem boundaryTraceTable_of_totals (h : BoundaryTotals) : BoundaryTraceTable := by
  intro n hn g hg
  obtain ⟨ys, hrun, hk⟩ := Option.map_eq_some_iff.1 (h n hn g hg)
  obtain ⟨hlo, hhi⟩ := trace_bounds_of_run hn hg (bandCheck_all n hn g hg) hrun
  have he := entryCheck_all n hn g hg
  simp only [entryCheck, Bool.and_eq_true, decide_eq_true_eq] at he
  obtain ⟨⟨⟨⟨hy, hL0⟩, hs0⟩, hclo⟩, hchi⟩ := he
  set k := paperK n g
  set N := numerator params (rowWeights n) (colWeights n) ys
  -- `k/10¹⁶ ≤ W < (k + 1)/10¹⁶` from the integer quotient.
  have hD : (0 : ℝ) < 2 * (2 ^ 50) ^ 2 := by positivity
  have hkN : k * (2 * (2 ^ 50) ^ 2) ≤ 10 ^ 16 * N ∧ 10 ^ 16 * N < (k + 1) * (2 * (2 ^ 50) ^ 2) := by
    rw [← hk]
    refine ⟨Nat.div_mul_le_self _ _, ?_⟩
    have := Nat.lt_div_mul_add (a := 10 ^ 16 * N) (b := 2 * (2 ^ 50) ^ 2) (by positivity)
    rw [add_mul, one_mul]; exact this
  set W : ℝ := (N : ℝ) / (2 * (2 ^ 50) ^ 2)
  have hW0 : 0 ≤ W := by positivity
  have hWL : ((k / 10 ^ 16 : ℚ) : ℝ) ≤ W := by
    have : (k : ℝ) * (2 * (2 ^ 50) ^ 2) ≤ 10 ^ 16 * N := by exact_mod_cast hkN.1
    push_cast
    rw [div_le_div_iff₀ (by positivity) hD]
    linarith
  have hWU : W ≤ (((k + 1) / 10 ^ 16 : ℚ) : ℝ) := by
    have : (10 : ℝ) ^ 16 * N < (k + 1) * (2 * (2 ^ 50) ^ 2) := by exact_mod_cast hkN.2
    push_cast
    rw [div_le_div_iff₀ hD (by positivity)]
    linarith
  -- `√W ≤ (⌊√(k+1)⌋ + 1)/10⁸`.
  have hWr : Real.sqrt W ≤ ((((Nat.sqrt (k + 1) + 1 : ℕ) : ℚ) / 10 ^ 8 : ℚ) : ℝ) := by
    have hs : (k + 1 : ℝ) ≤ ((Nat.sqrt (k + 1) + 1 : ℕ) : ℝ) ^ 2 := by
      exact_mod_cast (Nat.lt_succ_sqrt' (k + 1)).le
    rw [Real.sqrt_le_left (by positivity)]
    refine hWU.trans ?_
    push_cast at hs ⊢
    rw [div_pow, div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith
  have hexp := exp_mem_expLo_expHi hy
  have c0 := (Rat.cast_le (K := ℝ)).2 hs0
  have c1 := (Rat.cast_le (K := ℝ)).2 hL0
  have c2 := (Rat.cast_le (K := ℝ)).2 hclo
  have c3 := (Rat.cast_le (K := ℝ)).2 hchi
  push_cast at c0 c1 c2 c3 hWL hWU hWr hexp
  have hint := interval_of_bounds (ε := 2 / 10 ^ 10) hW0 hWL hWU hWr (by norm_num) hexp.1 hexp.2
    c0 c1 hlo hhi
  exact ⟨c2.trans hint.1, hint.2.trans c3⟩

end BoundaryThermal

end CertifiedGaps.SpinOneHeisenberg
