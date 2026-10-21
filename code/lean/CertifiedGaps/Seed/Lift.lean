/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Hamiltonian
import CertifiedGaps.Seed.CappedMass
import CertifiedGaps.Seed.Limits

/-!
# Temperature lifting of Boltzmann sums

A trace table is computed at one inverse temperature `β₀`, but a seed may be needed at a larger
`β = (l/j) β₀`. For a sum of Boltzmann factors `Z(β) = ∑ᵢ e^{-β Eᵢ}`, write `yᵢ = e^{-β₀ Eᵢ / j}`:
then `Z(β₀) = ∑ yᵢ^j` and `Z(β) = ∑ yᵢ^l`, and the power-mean inequality

  `(∑ yᵢ^l)^j ≤ (∑ yᵢ^j)^l`  for `yᵢ ≥ 0`, `1 ≤ j ≤ l`

lifts every upper bound `Z(β₀) ≤ m` to `Z(β)^j ≤ m^l` (`BoltzmannRep.lift`). This is the step
`Z(b) ≤ Z(b₀)^{5/2}` of the initialization at length sixty in OpenAI's *The periodic spin-one
Haldane gap* (Proposition 5.4, with `j = 2`, `l = 5`), stated with integer powers only.

The same inequality, applied to all factors but one, gives a thermal purity: if one factor
exceeds `1` at `β₀` and the sum is at most `q` with `(q - 1)^l < R^j`, then the factors other than
the largest carry less than `R` times its weight at `β`, so `Z(2β) / Z(β)² ≥ (1 + R)⁻²`
(`BoltzmannRep.sq_le_of_dominant`).

Twisted partition functions `Re Tr e^{-β H(g)}` with a Hermitian twisted Hamiltonian `H(g)` are
Boltzmann sums as well (`boltzmannRep_gibbsTrace`), so the lifting applies to every twist.

## Main results

* `sum_pow_pow_le`: `(∑ yᵢ^l)^j ≤ (∑ yᵢ^j)^l`.
* `BoltzmannRep`: `f β = ∑ᵢ e^{-β Eᵢ}` for finitely many real `Eᵢ`; `BoltzmannRep.shift`,
  `BoltzmannRep.lift`.
* `sum_pow_sq_le_of_dominant`, `boltzmann_sq_le_of_dominant`: the thermal purity from one
  dominant factor.
* `shiftedTrace_sq_le_iff`: such a purity bound is invariant under an energy shift.
-/

open Finset Real

namespace CertifiedGaps

/-! ### The power-mean inequality with integer exponents -/

/-- **Power means.** Nonnegative weights satisfy `(∑ yᵢ^l)^j ≤ (∑ yᵢ^j)^l` for `1 ≤ j ≤ l`: with
`M = max yᵢ` and `S = ∑ yᵢ^j ≥ M^j`, `∑ yᵢ^l ≤ M^{l-j} S`. -/
theorem sum_pow_pow_le {κ : Type*} (s : Finset κ) (y : κ → ℝ) (hy : ∀ i ∈ s, 0 ≤ y i)
    {j l : ℕ} (hj : 1 ≤ j) (hjl : j ≤ l) :
    (∑ i ∈ s, y i ^ l) ^ j ≤ (∑ i ∈ s, y i ^ j) ^ l := by
  classical
  rcases s.eq_empty_or_nonempty with rfl | hs
  · simp [zero_pow (by omega : j ≠ 0), zero_pow (by omega : l ≠ 0)]
  obtain ⟨i₀, hi₀, hmax⟩ := s.exists_max_image y hs
  set S := ∑ i ∈ s, y i ^ j
  have hM0 : 0 ≤ y i₀ := hy i₀ hi₀
  have hMS : y i₀ ^ j ≤ S :=
    single_le_sum (f := fun i ↦ y i ^ j) (fun i hi ↦ pow_nonneg (hy i hi) j) hi₀
  have h1 : ∑ i ∈ s, y i ^ l ≤ y i₀ ^ (l - j) * S := by
    rw [mul_sum]
    refine sum_le_sum fun i hi ↦ ?_
    rw [show y i ^ l = y i ^ (l - j) * y i ^ j by rw [← pow_add, Nat.sub_add_cancel hjl]]
    exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (hy i hi) (hmax i hi) _)
      (pow_nonneg (hy i hi) j)
  have hl0 : 0 ≤ ∑ i ∈ s, y i ^ l := sum_nonneg fun i hi ↦ pow_nonneg (hy i hi) l
  calc (∑ i ∈ s, y i ^ l) ^ j ≤ (y i₀ ^ (l - j) * S) ^ j := pow_le_pow_left₀ hl0 h1 j
    _ = (y i₀ ^ j) ^ (l - j) * S ^ j := by ring
    _ ≤ S ^ (l - j) * S ^ j := by
        gcongr; exact pow_nonneg (sum_nonneg fun i hi ↦ pow_nonneg (hy i hi) j) j
    _ = S ^ l := by rw [← pow_add, Nat.sub_add_cancel hjl]

/-! ### Boltzmann sums -/

/-- `f` is a sum of Boltzmann factors: `f β = ∑ᵢ e^{-β Eᵢ}` for finitely many real energies
`Eᵢ`, such as the eigenvalues of a Hermitian Hamiltonian. -/
def BoltzmannRep (f : ℝ → ℝ) : Prop :=
  ∃ (ι : Type) (_ : Fintype ι) (E : ι → ℝ), ∀ β, f β = ∑ i, exp (-(β * E i))

/-- The Gibbs trace of a Hermitian matrix is a Boltzmann sum over its eigenvalues. -/
theorem boltzmannRep_gibbsTrace {n : Type} [Fintype n] [DecidableEq n] {H : Matrix n n ℂ}
    (hH : H.IsHermitian) : BoltzmannRep (gibbsTrace H) :=
  ⟨n, inferInstance, hH.eigenvalues, gibbsTrace_eq_sum hH⟩

namespace BoltzmannRep

variable {f : ℝ → ℝ}

/-- Boltzmann sums are nonnegative. -/
theorem nonneg (h : BoltzmannRep f) (β : ℝ) : 0 ≤ f β := by
  obtain ⟨ι, _, E, hE⟩ := h
  rw [hE]
  exact sum_nonneg fun _ _ ↦ (exp_pos _).le

/-- The shifted Boltzmann sum: the energies `Eᵢ + κ n`. -/
theorem shiftedTrace_eq {ι : Type*} [Fintype ι] (E : ι → ℝ) (β κ : ℝ) (n : ℕ) :
    shiftedTrace β κ n (∑ i, exp (-(β * E i))) = ∑ i, exp (-(β * (E i + κ * n))) := by
  rw [shiftedTrace, mul_sum]
  refine sum_congr rfl fun i _ ↦ ?_
  rw [← exp_add]
  ring_nf

/-- An energy shift by `κ` per bond keeps a Boltzmann sum a Boltzmann sum. -/
theorem shift (h : BoltzmannRep f) (κ : ℝ) (n : ℕ) :
    BoltzmannRep fun β ↦ shiftedTrace β κ n (f β) := by
  obtain ⟨ι, _, E, hE⟩ := h
  exact ⟨ι, inferInstance, fun i ↦ E i + κ * n, fun β ↦ by
    change shiftedTrace β κ n (f β) = _
    rw [hE, shiftedTrace_eq]⟩

/-- The Boltzmann factors at `β` and `(l/j) β` as powers of `yᵢ = e^{-β Eᵢ / j}`. -/
theorem exp_eq_pow {j : ℕ} (hj : 1 ≤ j) (l : ℕ) (β e : ℝ) :
    exp (-((l / j * β) * e)) = exp (-(β * e) / j) ^ l := by
  rw [← exp_nat_mul]
  have : (j : ℝ) ≠ 0 := by exact_mod_cast (show j ≠ 0 by omega)
  congr 1
  field_simp

/-- **Temperature lifting.** A Boltzmann sum satisfies `f((l/j) β)^j ≤ f(β)^l` for
`1 ≤ j ≤ l`. -/
theorem lift (h : BoltzmannRep f) {j l : ℕ} (hj : 1 ≤ j) (hjl : j ≤ l) (β : ℝ) :
    f (l / j * β) ^ j ≤ f β ^ l := by
  obtain ⟨ι, _, E, hE⟩ := h
  have e₁ : f (l / j * β) = ∑ i, exp (-(β * E i) / j) ^ l := by
    rw [hE]; exact sum_congr rfl fun i _ ↦ exp_eq_pow hj l β (E i)
  have e₂ : f β = ∑ i, exp (-(β * E i) / j) ^ j := by
    rw [hE]
    refine sum_congr rfl fun i _ ↦ ?_
    have := exp_eq_pow hj j β (E i)
    have hj' : (j : ℝ) / j = 1 := div_self (by exact_mod_cast (show j ≠ 0 by omega))
    rwa [hj', one_mul] at this
  rw [e₁, e₂]
  exact sum_pow_pow_le _ _ (fun _ _ ↦ (exp_pos _).le) hj hjl

/-- An upper bound `f β ≤ m` lifts to `f ((l/j) β) < D` when `m^l < D^j` and `D ≥ 0`. -/
theorem lt_of_lift (h : BoltzmannRep f) {j l : ℕ} (hj : 1 ≤ j) (hjl : j ≤ l) {β m D : ℝ}
    (hm : f β ≤ m) (hD0 : 0 ≤ D) (hmD : m ^ l < D ^ j) : f (l / j * β) < D := by
  by_contra hc
  push Not at hc
  have h1 : D ^ j ≤ f (l / j * β) ^ j := pow_le_pow_left₀ hD0 hc j
  have h2 : f β ^ l ≤ m ^ l := pow_le_pow_left₀ (h.nonneg β) hm l
  linarith [h.lift hj hjl β]

end BoltzmannRep

/-! ### The thermal purity from one dominant factor -/

/-- **One dominant weight.** Let `yᵢ ≥ 0` with `y_{i₀}^j > 1`, `∑ yᵢ^j ≤ q` and
`(q - 1)^l < R^j`, `R ≥ 0`, `1 ≤ j ≤ l`. Then the weights other than `y_{i₀}` carry
`∑_{i ≠ i₀} yᵢ^l < R ≤ R y_{i₀}^l`, so `(∑ yᵢ^l)² ≤ (1 + R)² ∑ yᵢ^{2l}`. -/
theorem sum_pow_sq_le_of_dominant {κ : Type*} (s : Finset κ) (y : κ → ℝ)
    (hy : ∀ i ∈ s, 0 ≤ y i) {j l : ℕ} (hj : 1 ≤ j) (hjl : j ≤ l) {i₀ : κ} (hi₀ : i₀ ∈ s)
    (h1 : 1 < y i₀ ^ j) {q R : ℝ} (hR : 0 ≤ R) (hq : ∑ i ∈ s, y i ^ j ≤ q)
    (hqR : (q - 1) ^ l < R ^ j) :
    (∑ i ∈ s, y i ^ l) ^ 2 ≤ (1 + R) ^ 2 * ∑ i ∈ s, y i ^ (2 * l) := by
  classical
  have hy0 : 0 ≤ y i₀ := hy i₀ hi₀
  have hy1 : 1 ≤ y i₀ := by
    by_contra hc
    push Not at hc
    have := pow_le_one₀ hy0 hc.le (n := j)
    linarith
  have hyl : 1 ≤ y i₀ ^ l := one_le_pow₀ hy1
  have hy' : ∀ i ∈ s.erase i₀, 0 ≤ y i := fun i hi ↦ hy i (mem_of_mem_erase hi)
  -- The weights other than `y i₀`.
  set rj := ∑ i ∈ s.erase i₀, y i ^ j
  set rl := ∑ i ∈ s.erase i₀, y i ^ l
  have hsplit (n : ℕ) : ∑ i ∈ s, y i ^ n = y i₀ ^ n + ∑ i ∈ s.erase i₀, y i ^ n :=
    (add_sum_erase s (fun i ↦ y i ^ n) hi₀).symm
  have hrj0 : 0 ≤ rj := sum_nonneg fun i hi ↦ pow_nonneg (hy' i hi) j
  have hrl0 : 0 ≤ rl := sum_nonneg fun i hi ↦ pow_nonneg (hy' i hi) l
  have hrjq : rj ≤ q - 1 := by rw [hsplit] at hq; linarith
  have hrl : rl < R := by
    by_contra hc
    push Not at hc
    have a := sum_pow_pow_le _ y hy' hj hjl
    have b := pow_le_pow_left₀ hrj0 hrjq l
    have c := pow_le_pow_left₀ hR hc j
    linarith
  have hZ : ∑ i ∈ s, y i ^ l ≤ y i₀ ^ l * (1 + R) := by
    rw [hsplit]; nlinarith
  have hZ0 : 0 ≤ ∑ i ∈ s, y i ^ l := sum_nonneg fun i hi ↦ pow_nonneg (hy i hi) l
  have h2 : (y i₀ ^ l) ^ 2 ≤ ∑ i ∈ s, y i ^ (2 * l) := by
    rw [hsplit, ← pow_mul, mul_comm l 2]
    linarith [sum_nonneg fun i hi ↦ pow_nonneg (hy' i hi) (2 * l)]
  calc (∑ i ∈ s, y i ^ l) ^ 2 ≤ (y i₀ ^ l * (1 + R)) ^ 2 := pow_le_pow_left₀ hZ0 hZ 2
    _ = (1 + R) ^ 2 * (y i₀ ^ l) ^ 2 := by ring
    _ ≤ (1 + R) ^ 2 * ∑ i ∈ s, y i ^ (2 * l) := by gcongr

/-- **The thermal purity from one dominant Boltzmann factor.** If `e^{-β₀ E_{i₀}} > 1`,
`∑ᵢ e^{-β₀ Eᵢ} ≤ q` and `(q - 1)^l < R^j` with `R ≥ 0`, `1 ≤ j ≤ l`, then at `β = (l/j) β₀`,
`Z(β)² ≤ (1 + R)² Z(2β)`. -/
theorem boltzmann_sq_le_of_dominant {ι : Type*} [Fintype ι] (E : ι → ℝ) {j l : ℕ} (hj : 1 ≤ j)
    (hjl : j ≤ l) {β₀ q R : ℝ} {i₀ : ι} (h1 : 1 < exp (-(β₀ * E i₀))) (hR : 0 ≤ R)
    (hq : ∑ i, exp (-(β₀ * E i)) ≤ q) (hqR : (q - 1) ^ l < R ^ j) :
    (∑ i, exp (-((l / j * β₀) * E i))) ^ 2 ≤
      (1 + R) ^ 2 * ∑ i, exp (-((2 * (l / j * β₀)) * E i)) := by
  classical
  have hjj : (j : ℝ) / j = 1 := div_self (by exact_mod_cast (show j ≠ 0 by omega))
  have ej (e : ℝ) : exp (-(β₀ * e)) = exp (-(β₀ * e) / j) ^ j := by
    have := BoltzmannRep.exp_eq_pow hj j β₀ e
    rwa [hjj, one_mul] at this
  have e2 (e : ℝ) : exp (-((2 * (l / j * β₀)) * e)) = exp (-(β₀ * e) / j) ^ (2 * l) := by
    have := BoltzmannRep.exp_eq_pow hj (2 * l) β₀ e
    rw [← this]; push_cast; ring_nf
  simp only [BoltzmannRep.exp_eq_pow hj l, e2]
  simp only [ej] at h1 hq
  exact sum_pow_sq_le_of_dominant _ (fun i ↦ exp (-(β₀ * E i) / j)) (fun _ _ ↦ (exp_pos _).le)
    hj hjl (mem_univ i₀) h1 hR hq hqR

/-- One Boltzmann factor of the shifted energies `Eᵢ + κ n` exceeds `1` when `E_{i₀} < -κ n`, so
the shifted Boltzmann sum exceeds `1` at every `β > 0`. -/
theorem one_lt_shiftedTrace_of_lt {ι : Type*} [Fintype ι] (E : ι → ℝ) {β κ : ℝ} {n : ℕ}
    (hβ : 0 < β) {i₀ : ι} (h : E i₀ < -(κ * n)) :
    1 < shiftedTrace β κ n (∑ i, exp (-(β * E i))) := by
  rw [BoltzmannRep.shiftedTrace_eq]
  calc (1 : ℝ) < exp (-(β * (E i₀ + κ * n))) := by
        rw [← exp_zero]; exact exp_lt_exp.2 (by nlinarith)
    _ ≤ ∑ i, exp (-(β * (E i + κ * n))) :=
        single_le_sum (f := fun i ↦ exp (-(β * (E i + κ * n)))) (fun _ _ ↦ (exp_pos _).le)
          (mem_univ i₀)

/-- A thermal-purity bound `Z(β)² ≤ c Z(2β)` holds for the shifted traces `e^{-βκn} Z` exactly
when it holds for the traces: the shift factor at `2β` is the square of the one at `β`. -/
theorem shiftedTrace_sq_le_iff (β κ : ℝ) (n : ℕ) (z₁ z₂ c : ℝ) :
    shiftedTrace β κ n z₁ ^ 2 ≤ c * shiftedTrace (2 * β) κ n z₂ ↔ z₁ ^ 2 ≤ c * z₂ := by
  have e : exp (-(2 * β * κ * n)) = exp (-(β * κ * n)) ^ 2 := by
    rw [← exp_nat_mul]; push_cast; ring_nf
  rw [shiftedTrace, shiftedTrace, e, mul_pow, mul_left_comm]
  exact mul_le_mul_iff_of_pos_left (by positivity)

end CertifiedGaps
