/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Spectral
import Mathlib.Analysis.SpecialFunctions.Exponential

/-!
# Partition functions from a low spectrum

A spectrum `e : ι → ℝ` has a *low spectrum* `(a₀, b₀, a₁, ρ)` if one eigenvalue lies in
`[a₀, b₀]`, a second one is at least `a₁`, and all others are at least `ρ` (`LowSpectrum`). This
is what a residual bound on two trial vectors (`exists_near_of_residual`) and counting by
comparison (`card_eigenvalues_lt_le`, `le_re_quadratic_of_weights`) give for any Hermitian matrix
(`lowSpectrum_of_near`). This file turns a low spectrum into bounds on Boltzmann sums
`Z(β) = ∑ᵢ e^{-β eᵢ}` and on the three ratios a seed needs, for any spectra:

* `thermal_ratio_ge`: `Z(2β) / Z(β)² ≥ (1 - 4βr) / (1 + t)²`;
* `spatial_ratio_ge`: `Z₂(β) / Z₁(β)² ≥ (1 - β(θ₂ + r - 2(θ₁ - r))) / (1 + t)²`;
* `trace_le`: `Z₁ (τ Z₂)^m ≤ E Z₃` from upper bounds on `Z₁, Z₂` and a lower bound on `Z₃`;

where `t` bounds the excited-state sum `e^{-β(a₁ - a₀)} + N e^{-β(ρ - a₀)}`. Nothing here depends
on the model.
-/

open Matrix Real

namespace SeedKernel

/-- One eigenvalue in `[a₀, b₀]`, a second one at least `a₁`, all others at least `ρ`. -/
def LowSpectrum {ι : Type*} (e : ι → ℝ) (a₀ b₀ a₁ ρ : ℝ) : Prop :=
  ∃ i₀ i₁, i₀ ≠ i₁ ∧ a₀ ≤ e i₀ ∧ e i₀ ≤ b₀ ∧ a₁ ≤ e i₁ ∧ ∀ i, i ≠ i₀ → i ≠ i₁ → ρ ≤ e i

theorem sum_exp_ge {ι : Type*} [Fintype ι] {e : ι → ℝ} {a₀ b₀ a₁ ρ β : ℝ} (hβ : 0 ≤ β)
    (h : LowSpectrum e a₀ b₀ a₁ ρ) : exp (-(β * b₀)) ≤ ∑ i, exp (-(β * e i)) := by
  obtain ⟨i₀, -, -, -, hb, -, -⟩ := h
  calc exp (-(β * b₀)) ≤ exp (-(β * e i₀)) := exp_le_exp.2 (by nlinarith)
    _ ≤ ∑ i, exp (-(β * e i)) :=
      Finset.single_le_sum (f := fun i ↦ exp (-(β * e i))) (fun i _ ↦ (exp_pos _).le)
        (Finset.mem_univ i₀)

theorem sum_exp_le {ι : Type*} [Fintype ι] {e : ι → ℝ} {a₀ b₀ a₁ ρ β : ℝ}
    (hβ : 0 ≤ β) (h : LowSpectrum e a₀ b₀ a₁ ρ) :
    ∑ i, exp (-(β * e i)) ≤
      exp (-(β * a₀)) + exp (-(β * a₁)) + Fintype.card ι * exp (-(β * ρ)) := by
  classical
  obtain ⟨i₀, i₁, hne, ha, -, ha₁, hρ⟩ := h
  have hi₁ : i₁ ∈ Finset.univ.erase i₀ := Finset.mem_erase.2 ⟨hne.symm, Finset.mem_univ _⟩
  rw [← Finset.add_sum_erase _ _ (Finset.mem_univ i₀), ← Finset.add_sum_erase _ _ hi₁,
    ← add_assoc]
  refine add_le_add (add_le_add (exp_le_exp.2 (by nlinarith)) (exp_le_exp.2 (by nlinarith))) ?_
  calc ∑ i ∈ (Finset.univ.erase i₀).erase i₁, exp (-(β * e i))
      ≤ ∑ i ∈ (Finset.univ.erase i₀).erase i₁, exp (-(β * ρ)) :=
        Finset.sum_le_sum fun i hi ↦ by
          have hi' := Finset.mem_erase.1 hi
          have hi'' := Finset.mem_erase.1 hi'.2
          exact exp_le_exp.2 (by nlinarith [hρ i hi''.1 hi'.1])
    _ ≤ Fintype.card ι * exp (-(β * ρ)) := by
        rw [Finset.sum_const, nsmul_eq_mul]
        gcongr
        exact_mod_cast (Finset.card_le_univ _)

/-- The ratio step: from `Z₂ ≥ A² e^{-x}`, `Z₁ ≤ A (1 + T)`, `T ≤ t` to
`Z₂ / Z₁² ≥ (1-x)/(1+t)²`. -/
theorem ratio_ge {Z₁ Z₂ A x T t : ℝ} (hA : 0 < A) (hZ₁ : 0 < Z₁) (h₁ : Z₁ ≤ A * (1 + T))
    (h₂ : A ^ 2 * exp (-x) ≤ Z₂) (hT : 0 ≤ T) (hTt : T ≤ t) :
    (1 - x) / (1 + t) ^ 2 ≤ Z₂ / Z₁ ^ 2 := by
  have hx : 1 - x ≤ exp (-x) := by linarith [add_one_le_exp (-x)]
  have hZ : Z₁ ^ 2 ≤ A ^ 2 * (1 + t) ^ 2 := by
    rw [← mul_pow]; gcongr; linarith [mul_le_mul_of_nonneg_left hTt hA.le]
  have ht : 0 < 1 + t := by linarith
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  by_cases hx0 : 1 - x ≤ 0
  · nlinarith [sq_nonneg Z₁, exp_pos (-x), sq_nonneg A, sq_nonneg (1 + t)]
  · calc (1 - x) * Z₁ ^ 2 ≤ (1 - x) * (A ^ 2 * (1 + t) ^ 2) := by gcongr; linarith
      _ ≤ exp (-x) * A ^ 2 * (1 + t) ^ 2 := by nlinarith [sq_nonneg A, sq_nonneg (1 + t)]
      _ ≤ Z₂ * (1 + t) ^ 2 := by nlinarith [sq_nonneg (1 + t)]

/-! ### A low spectrum from two eigenvalues and a count -/

/-- **A residual bound with a rational radius.** If `‖x‖² ≥ N > 0` and `‖A x - θ x‖² ≤ r² N`,
some eigenvalue of `A` lies within `r` of `θ`. -/
theorem exists_near_of_residual {n : Type*} [Fintype n] [DecidableEq n] [Nonempty n]
    {A : Matrix n n ℂ} (hA : A.IsHermitian) (x : n → ℂ) {θ r N : ℝ} (hr : 0 ≤ r) (hN : 0 < N)
    (hx : N ≤ normSqSum x) (hres : normSqSum (A *ᵥ x - (θ : ℂ) • x) ≤ r ^ 2 * N) :
    ∃ i, |hA.eigenvalues i - θ| ≤ r := by
  obtain ⟨i, hi⟩ := exists_eigenvalue_near hA x θ
  refine ⟨i, abs_le_of_sq_le_sq' ?_ hr |> fun h ↦ abs_le.2 h⟩
  have h1 : (hA.eigenvalues i - θ) ^ 2 * N ≤ r ^ 2 * N :=
    (mul_le_mul_of_nonneg_left hx (sq_nonneg _)).trans (hi.trans hres)
  exact le_of_mul_le_mul_right h1 hN

/-- **A low spectrum from two eigenvalues and a count.** Two eigenvalues within `r` of
`θ₀ + 2r < θ₁`, and at most two eigenvalues below `ρ > θ₁ + r`. -/
theorem lowSpectrum_of_near {ι : Type*} [Fintype ι] {e : ι → ℝ}
    {θ₀ θ₁ r ρ : ℝ} (h₀ : ∃ i, |e i - θ₀| ≤ r) (h₁ : ∃ i, |e i - θ₁| ≤ r) (h01 : θ₀ + 2 * r < θ₁)
    (h1ρ : θ₁ + r < ρ) (hcard : (Finset.univ.filter fun i ↦ e i < ρ).card ≤ 2) :
    LowSpectrum e (θ₀ - r) (θ₀ + r) (θ₁ - r) ρ := by
  classical
  obtain ⟨i₀, h₀⟩ := h₀
  obtain ⟨i₁, h₁⟩ := h₁
  have h₀' := abs_le.1 h₀
  have h₁' := abs_le.1 h₁
  have hne : i₀ ≠ i₁ := by rintro rfl; linarith
  refine ⟨i₀, i₁, hne, by linarith, by linarith, by linarith, fun i hi₀ hi₁ ↦ ?_⟩
  by_contra hlt
  push Not at hlt
  have hsub : ({i₀, i₁, i} : Finset _) ⊆ Finset.univ.filter fun j ↦ e j < ρ := by
    intro j hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    rcases hj with rfl | rfl | rfl <;> linarith
  have h3 : ({i₀, i₁, i} : Finset _).card = 3 := by
    rw [Finset.card_insert_of_notMem (by simp [hne, Ne.symm hi₀]),
      Finset.card_pair (Ne.symm hi₁)]
  have := (Finset.card_le_card hsub).trans hcard
  omega

/-- The lowest eigenvalue of a low spectrum `(θ - r, θ + r, a₁, ρ)` lies within `r` of `θ`. -/
theorem near_of_lowSpectrum {ι : Type*} {e : ι → ℝ} {θ r a₁ ρ : ℝ}
    (h : LowSpectrum e (θ - r) (θ + r) a₁ ρ) : ∃ i, |e i - θ| ≤ r := by
  obtain ⟨i₀, -, -, h₀, h₁, -, -⟩ := h
  exact ⟨i₀, abs_le.2 ⟨by linarith, by linarith⟩⟩

/-! ### Boltzmann sums and their ratios -/

variable {ι κ μ : Type*} [Fintype ι] [Fintype κ] [Fintype μ]

/-- `Z(β) ≤ e^{-β a₀} (1 + e^{-β (a₁ - a₀)} + N e^{-β (ρ - a₀)})` with `N = |ι|`. -/
theorem sum_exp_le_mul {e : ι → ℝ} {a₀ b₀ a₁ ρ β : ℝ} (hβ : 0 ≤ β)
    (h : LowSpectrum e a₀ b₀ a₁ ρ) :
    ∑ i, exp (-(β * e i)) ≤ exp (-(β * a₀)) *
      (1 + (exp (-(β * (a₁ - a₀))) + Fintype.card ι * exp (-(β * (ρ - a₀))))) := by
  refine (sum_exp_le hβ h).trans (le_of_eq ?_)
  simp only [mul_add, mul_one, ← exp_add, mul_left_comm (exp (-(β * a₀))) (Fintype.card ι : ℝ)]
  ring_nf

/-- An eigenvalue within `r` of `θ` gives `Z(β) ≥ e^{-β (θ + r)}`. -/
theorem exp_le_sum_of_near {e : ι → ℝ} {θ r β : ℝ} (hβ : 0 ≤ β) (h : ∃ i, |e i - θ| ≤ r) :
    exp (-(β * (θ + r))) ≤ ∑ i, exp (-(β * e i)) := by
  obtain ⟨i, hi⟩ := h
  calc exp (-(β * (θ + r))) ≤ exp (-(β * e i)) :=
        exp_le_exp.2 (by nlinarith [(abs_le.1 hi).2])
    _ ≤ _ := Finset.single_le_sum (f := fun i ↦ exp (-(β * e i)))
        (fun i _ ↦ (exp_pos _).le) (Finset.mem_univ i)

/-- The excited-state sum `e^{-β (a₁ - a₀)} + N e^{-β (ρ - a₀)}` of a low spectrum. -/
noncomputable def excited (N : ℕ) (β a₀ a₁ ρ : ℝ) : ℝ :=
  exp (-(β * (a₁ - a₀))) + N * exp (-(β * (ρ - a₀)))

theorem excited_nonneg (N : ℕ) (β a₀ a₁ ρ : ℝ) : 0 ≤ excited N β a₀ a₁ ρ := by
  unfold excited; positivity

theorem sum_exp_pos [Nonempty ι] (e : ι → ℝ) (β : ℝ) : 0 < ∑ i, exp (-(β * e i)) :=
  Finset.sum_pos (fun _ _ ↦ exp_pos _) Finset.univ_nonempty

/-- **The thermal ratio.** For a low spectrum `(θ - r, θ + r, a₁, ρ)` and
`excited N β (θ - r) a₁ ρ ≤ t`, `Z(2β) / Z(β)² ≥ (1 - 4βr) / (1 + t)²`. -/
theorem thermal_ratio_ge [Nonempty ι] {e : ι → ℝ} {θ r a₁ ρ β t : ℝ}
    (hβ : 0 ≤ β) (h : LowSpectrum e (θ - r) (θ + r) a₁ ρ)
    (ht : excited (Fintype.card ι) β (θ - r) a₁ ρ ≤ t) :
    (1 - 4 * β * r) / (1 + t) ^ 2 ≤
      (∑ i, exp (-(2 * β * e i))) / (∑ i, exp (-(β * e i))) ^ 2 := by
  have hu := sum_exp_le_mul hβ h
  have hl : exp (-(2 * β * (θ + r))) ≤ ∑ i, exp (-(2 * β * e i)) := by
    obtain ⟨i₀, -, -, h₀, h₁, -, -⟩ := h
    have := exp_le_sum_of_near (e := e) (θ := θ) (r := r) (β := 2 * β) (by positivity)
      ⟨i₀, abs_le.2 ⟨by linarith, by linarith⟩⟩
    simpa only [mul_assoc] using this
  refine ratio_ge (exp_pos _) (sum_exp_pos e β) hu (le_trans (le_of_eq ?_) hl)
    (excited_nonneg _ _ _ _ _) ht
  rw [← exp_nat_mul, ← exp_add]
  congr 1
  push_cast
  ring

/-- **The spatial ratio.** For a low spectrum `(θ₁ - r, θ₁ + r, a₁, ρ)` of `e₁`, an eigenvalue of
`e₂` within `r` of `θ₂`, and `excited N β (θ₁ - r) a₁ ρ ≤ t`,
`Z₂(β) / Z₁(β)² ≥ (1 - β (θ₂ + r - 2 (θ₁ - r))) / (1 + t)²`. -/
theorem spatial_ratio_ge [Nonempty ι] {e₁ : ι → ℝ} {e₂ : κ → ℝ}
    {θ₁ θ₂ r a₁ ρ β t : ℝ} (hβ : 0 ≤ β) (h₁ : LowSpectrum e₁ (θ₁ - r) (θ₁ + r) a₁ ρ)
    (h₂ : ∃ i, |e₂ i - θ₂| ≤ r) (ht : excited (Fintype.card ι) β (θ₁ - r) a₁ ρ ≤ t) :
    (1 - β * (θ₂ + r - 2 * (θ₁ - r))) / (1 + t) ^ 2 ≤
      (∑ i, exp (-(β * e₂ i))) / (∑ i, exp (-(β * e₁ i))) ^ 2 := by
  refine ratio_ge (exp_pos _) (sum_exp_pos e₁ β) (sum_exp_le_mul hβ h₁)
    (le_trans (le_of_eq ?_) (exp_le_sum_of_near hβ h₂)) (excited_nonneg _ _ _ _ _) ht
  rw [← exp_nat_mul, ← exp_add]
  congr 1
  push_cast
  ring

/-- **The trace bound.** From low spectra of `e₁` and `e₂` and an eigenvalue of `e₃` within `r`
of `θ₃`: if `e^{β (θ₃ + r - (θ₁ - r) - m (θ₂ - r))} (1 + s₁) ((1 + s₂) τ)^m ≤ E`, where `s₁, s₂`
bound the excited-state sums, then `Z₁ (τ Z₂)^m ≤ E Z₃`. -/
theorem trace_le {e₁ : ι → ℝ} {e₂ : κ → ℝ} {e₃ : μ → ℝ}
    {θ₁ θ₂ θ₃ r a₁ a₂ ρ₁ ρ₂ β τ s₁ s₂ E : ℝ} (m : ℕ) (hβ : 0 ≤ β)
    (h₁ : LowSpectrum e₁ (θ₁ - r) (θ₁ + r) a₁ ρ₁) (h₂ : LowSpectrum e₂ (θ₂ - r) (θ₂ + r) a₂ ρ₂)
    (h₃ : ∃ i, |e₃ i - θ₃| ≤ r) (hτ : 0 ≤ τ)
    (hs₁ : excited (Fintype.card ι) β (θ₁ - r) a₁ ρ₁ ≤ s₁)
    (hs₂ : excited (Fintype.card κ) β (θ₂ - r) a₂ ρ₂ ≤ s₂)
    (hE : exp (β * (θ₃ + r - (θ₁ - r) - m * (θ₂ - r))) * (1 + s₁) * ((1 + s₂) * τ) ^ m ≤ E) :
    (∑ i, exp (-(β * e₁ i))) * (τ * ∑ i, exp (-(β * e₂ i))) ^ m ≤
      E * ∑ i, exp (-(β * e₃ i)) := by
  have hZ₁ := sum_exp_le_mul hβ h₁
  have hZ₂ := sum_exp_le_mul hβ h₂
  have hZ₃ := exp_le_sum_of_near hβ h₃
  have hs₁0 : 0 ≤ s₁ := (excited_nonneg _ _ _ _ _).trans hs₁
  have hs₂0 : 0 ≤ s₂ := (excited_nonneg _ _ _ _ _).trans hs₂
  set A₁ := exp (-(β * (θ₁ - r)))
  set A₂ := exp (-(β * (θ₂ - r)))
  set A₃ := exp (-(β * (θ₃ + r)))
  have hsplit : A₁ * A₂ ^ m = exp (β * (θ₃ + r - (θ₁ - r) - m * (θ₂ - r))) * A₃ := by
    simp only [A₁, A₂, A₃, ← exp_nat_mul, ← exp_add]
    congr 1
    ring
  have hZ₂0 : 0 ≤ ∑ i, exp (-(β * e₂ i)) := Finset.sum_nonneg fun i _ ↦ (exp_pos _).le
  calc (∑ i, exp (-(β * e₁ i))) * (τ * ∑ i, exp (-(β * e₂ i))) ^ m
      ≤ (A₁ * (1 + s₁)) * (τ * (A₂ * (1 + s₂))) ^ m := by
        gcongr
        · exact hZ₁.trans (by gcongr; exact hs₁)
        · exact hZ₂.trans (by gcongr; exact hs₂)
    _ = (A₁ * A₂ ^ m) * ((1 + s₁) * ((1 + s₂) * τ) ^ m) := by
        rw [show τ * (A₂ * (1 + s₂)) = A₂ * ((1 + s₂) * τ) by ring, mul_pow]; ring
    _ = (exp (β * (θ₃ + r - (θ₁ - r) - m * (θ₂ - r))) * (1 + s₁) * ((1 + s₂) * τ) ^ m) * A₃ := by
        rw [hsplit]; ring
    _ ≤ E * ∑ i, exp (-(β * e₃ i)) :=
        mul_le_mul hE hZ₃ (exp_pos _).le (le_trans (by positivity) hE)

end SeedKernel
