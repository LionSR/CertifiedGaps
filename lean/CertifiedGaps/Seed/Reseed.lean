/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Certificate
import CertifiedGaps.Seed.Lift

/-!
# Reseeding a table seed at a doubled temperature

A table seed (`TableCertificate`) at `(n₀, β₀)` starts a doubling chain. Some uses need a second
seed at a larger scale, with smaller defects than the chain itself provides there. This file
derives such a seed at `2β₀` from three facts: the table at `β₀`, one row of the chain at `2β₀`,
and trial states. It is the step "Reseeding at `(120, b₂)`" of Proposition "Concentration at the
doubling scales" in OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026).

The argument has four steps.

1. **Trace bounds at `β₀`.** Sector data that obey the table have each `|x|` below its cap and
   each sector's `q`th moment below its mass (`TableCertificate.abs_lt_cap_of_table`,
   `TableCertificate.sum_pow_le_mass_of_table`). The capped-mass bound of every sector, the
   dominant one included, bounds the trace at every length `q r` (`SectorSeed.traceBound`,
   `TableCertificate.shiftedTrace_le`). For the dominant eigenvalue itself,
   `SectorSeed.exists_window` gives `ℓ < |x i₀| < cap₀` and the residual moments
   `restBound r` at every `r ≥ 1`.
2. **Physical squaring.** A Boltzmann sum has `Z(2β) ≤ Z(β)²`, also after an energy shift
   (`BoltzmannRep.shiftedTrace_two_mul_le_sq`). This bounds the traces at `2β₀`.
3. **The spatial seed at `2β₀`** (`ReseedData.spatial`). Write `Z_n` for the traces at `2β₀`,
   approximated by power sums `∑ xᵢⁿ` of real numbers (`PowerSumApprox`). A spatial purity
   `Z_{2m}/Z_m² ≥ 1 - p` from the chain and a trial state `Z_m > 1` give one eigenvalue `x₀` with
   `x₀^m > 1 - c` (`exists_gt_of_purity`, the paper's `λ₁⁷² > .968`). A bound `Z_k ≤ B` from
   step 2 then bounds the relative tail `∑_{i ≠ 0} (xᵢ/x₀)^k ≤ T`, so `|xᵢ| ≤ ρ |x₀|` with
   `T ≤ ρ^k`, and the tail at a length `n ≥ k` is at most `ρ^{n-k} T`. Hence
   `P_s(n, 2β₀) ≥ (1 + ρ^{n-k} T)⁻²`.
4. **The thermal seed at `2β₀`** (`thermalPurity_ge_of_dominant`). One Boltzmann factor at
   length `N` exceeds `1` (a trial state) and the shifted trace at `β₀` is at most `q`; if
   `(q - 1)² < R`, then `P_τ(N, 2β₀) ≥ (1 + R)⁻²`.

One coupled update of the doubling recurrence (`purities_ge_chainCertificate`) then gives the
seed at `(2n, 4β₀)`.

## Differences from the paper

* The trace bound of step 1 bounds the dominant sector by its cap and mass, `U^{qr} +
  (m_J - U^q)^r`, instead of the paper's `U^k (1 + D(k))`, which divides the residual moments by
  `ℓ^k`. The bound is sharper and needs no dominance.
* The tail length `k` is a multiple of `q` (`k = 36` for the spin-one seed instead of the paper's
  `k = 32`), so that every power is an integer power. The moment order is raised from `k` to `n`
  through the ratio bound `ρ` instead of the power `(∑ aᵢ^k)^{n/k}`.
* Step 3 never names the transfer operator: the eigenvalue `x₀` lives on the approximating power
  sums, and the seed passes to the traces as the approximation error tends to zero.

## Main results

* `TableCertificate.shiftedTrace_le`: the trace bound at every length `q r`.
* `SectorSeed.exists_window`: the dominant eigenvalue window and the residual moments.
* `exists_gt_of_purity`: a dominant weight from a purity.
* `ReseedData.exists_dominant_of_sums`: the dominant eigenvalue, the tail and the ratio bound.
* `ReseedData.spatial`: the spatial seed at the doubled temperature.
* `thermalPurity_ge_of_dominant`: the thermal seed at the doubled temperature.
-/

open Finset Real

namespace CertifiedGaps

open PolynomialFilter

/-! ### Caps, masses and trace bounds from a table -/

namespace SectorSeed

variable {m : ℕ} (S : SectorSeed m)

/-- The capped-mass bound `∑_s μ_s · sectorBound s r` on the trace at length `q r`. The dominant
sector is bounded like the others, by its cap and its mass. -/
def traceBound (r : ℕ) : ℚ := ∑ s, S.mult s * S.sectorBound s r

/-- Every sector, the dominant one included, admits its capped-mass bound: a nonnegative mass
with one full bar under its cap. -/
def TraceValid : Prop := ∀ s, 0 ≤ S.mass s ∧ OneBar S.q (S.cap s) (S.mass s)

instance : Decidable S.TraceValid := by unfold TraceValid; infer_instance

variable {S} {ι : Type*} {sec : Fin (m + 1) → Finset ι} {x : ι → ℝ}

/-- **The trace bound of sector data** at length `q r` from the caps and the masses. -/
theorem sectorSum_le_traceBound (hq : Even S.q) (hT : S.TraceValid)
    (hcap : ∀ s v, S.cap s = some v → ∀ i ∈ sec s, |x i| < v)
    (hmass : ∀ s, ∑ i ∈ sec s, x i ^ S.q ≤ S.mass s) {r : ℕ} (hr : 1 ≤ r) :
    S.sectorSum sec x (S.q * r) ≤ S.traceBound r := by
  simp only [sectorSum, traceBound, Rat.cast_sum, Rat.cast_mul, Rat.cast_natCast]
  exact sum_le_sum fun s _ ↦ mul_le_mul_of_nonneg_left
    (sum_le_sectorBound hq (by exact_mod_cast (hT s).1) (hT s).2 (hcap s) (hmass s) hr)
    (Nat.cast_nonneg _)

/-- The trace of sector data without one eigenvalue `i₀` of the dominant sector. -/
def restSum [DecidableEq ι] (S : SectorSeed m) (sec : Fin (m + 1) → Finset ι) (x : ι → ℝ)
    (i₀ : ι) (n : ℕ) : ℝ :=
  ∑ i ∈ (sec 0).erase i₀, x i ^ n + ∑ i : Fin m, (S.mult i.succ : ℝ) * ∑ j ∈ sec i.succ, x j ^ n

/-- The trace splits into the eigenvalue `i₀` and the rest. -/
theorem sectorSum_eq_add_restSum [DecidableEq ι] (hm0 : S.mult 0 = 1) {i₀ : ι}
    (hi₀ : i₀ ∈ sec 0) (n : ℕ) : S.sectorSum sec x n = x i₀ ^ n + S.restSum sec x i₀ n := by
  simp only [sectorSum, Fin.sum_univ_succ, hm0, Nat.cast_one, one_mul, restSum,
    ← add_sum_erase (sec 0) (fun i ↦ x i ^ n) hi₀]
  ring

/-- **The residual moments.** Without an eigenvalue `|x i₀| > ℓ` of the dominant sector, the
moment of order `q r` is at most `restBound r`. -/
theorem restSum_le_restBound [DecidableEq ι] (hS : S.Valid)
    (hcap : ∀ s v, S.cap s = some v → ∀ i ∈ sec s, |x i| < v)
    (hmass : ∀ s, ∑ i ∈ sec s, x i ^ S.q ≤ S.mass s) {i₀ : ι} (hi₀ : i₀ ∈ sec 0)
    (hdom : (S.ell : ℝ) < |x i₀|) {r : ℕ} (hr : 1 ≤ r) :
    S.restSum sec x i₀ (S.q * r) ≤ S.restBound r := by
  obtain ⟨hq, hq0, -, hell, -, -, -, -, -, hN, -⟩ := hS
  -- The dominant sector without `i₀` carries at most `mass₀ - ℓ^q`.
  have hJ : ∑ i ∈ (sec 0).erase i₀, x i ^ (S.q * r) ≤ ((S.mass 0 : ℝ) - (S.ell : ℝ) ^ S.q) ^ r := by
    simp only [pow_mul]
    have hsp : x i₀ ^ S.q + ∑ i ∈ (sec 0).erase i₀, x i ^ S.q = ∑ i ∈ sec 0, x i ^ S.q :=
      add_sum_erase (sec 0) (fun i ↦ x i ^ S.q) hi₀
    have h0 : 0 ≤ ∑ i ∈ (sec 0).erase i₀, x i ^ S.q := sum_nonneg fun i _ ↦ hq.pow_nonneg _
    have hℓ : (S.ell : ℝ) ^ S.q ≤ x i₀ ^ S.q := by
      rw [← hq.pow_abs (x i₀)]
      exact pow_le_pow_left₀ (by exact_mod_cast hell.le) hdom.le _
    refine (sum_pow_le_pow_sum _ (fun i ↦ x i ^ S.q) (fun i _ ↦ hq.pow_nonneg _) hr).trans
      (pow_le_pow_left₀ h0 ?_ r)
    linarith [hmass 0]
  have hrest : ∑ i : Fin m, (S.mult i.succ : ℝ) * ∑ j ∈ sec i.succ, x j ^ (S.q * r) ≤
      ∑ i : Fin m, (S.mult i.succ : ℝ) * S.sectorBound i.succ r :=
    sum_le_sum fun i _ ↦ mul_le_mul_of_nonneg_left
      (sum_le_sectorBound hq (by exact_mod_cast (hN i).1) (hN i).2 (hcap _) (hmass _) hr)
      (Nat.cast_nonneg _)
  simp only [restSum, restBound, Rat.cast_add, Rat.cast_pow, Rat.cast_sub, Rat.cast_sum,
    Rat.cast_mul, Rat.cast_natCast]
  linarith

/-- **The dominant window.** Under `S.Valid`, the caps, the masses and `Z_{2n₀} > 1`, one
eigenvalue of the dominant sector has `ℓ < |x i₀| < cap₀`, and the rest of the spectrum has
moments of order `q r` at most `restBound r` for every `r ≥ 1`. -/
theorem exists_window [DecidableEq ι] (hS : S.Valid)
    (hcap : ∀ s v, S.cap s = some v → ∀ i ∈ sec s, |x i| < v)
    (hmass : ∀ s, ∑ i ∈ sec s, x i ^ S.q ≤ S.mass s)
    (h2 : 1 < S.sectorSum sec x (S.q * (2 * S.a))) :
    ∃ i₀ ∈ sec 0, (S.ell : ℝ) < |x i₀| ∧ |x i₀| < S.cap₀ ∧
      ∀ r, 1 ≤ r → S.restSum sec x i₀ (S.q * r) ≤ S.restBound r := by
  obtain ⟨i₀, hi₀, hdom⟩ := exists_dominant hS hcap hmass h2
  obtain ⟨v₀, hv₀⟩ := Option.isSome_iff_exists.1 hS.2.2.2.2.2.1
  have hcapv : S.cap₀ = v₀ := by simp [cap₀, hv₀]
  exact ⟨i₀, hi₀, hdom, by rw [hcapv]; exact hcap 0 v₀ hv₀ i₀ hi₀,
    fun r hr ↦ restSum_le_restBound hS hcap hmass hi₀ hdom hr⟩

end SectorSeed

namespace TableCertificate

variable {C : TableCertificate} {ι : Type*} {sec : Fin (C.scheme.m + 1) → Finset ι}
  {x : ι → ℝ}

/-- **The caps.** Sector data that obey the table of a valid certificate have `|x i| < v` in
every sector with a cap `v`: the filters exclude larger eigenvalues. -/
theorem abs_lt_cap_of_table (hC : C.Valid) (h : C.table.Holds (C.scheme.trace sec x)) :
    ∀ s v, C.seed.cap s = some v → ∀ i ∈ sec s, |x i| < v := by
  obtain ⟨hS, -, -, hfilter, -⟩ := hC
  have hM := TraceBounds.momentBounds hS h
  intro s v hv i hi
  have hf := hfilter s
  unfold FilterCheck at hf
  rw [hv] at hf
  rcases hfs : C.filter s with _ | ⟨a, qf⟩
  · rw [hfs] at hf; exact hf.elim
  · rw [hfs] at hf
    obtain ⟨B, hB, hcapB⟩ := CapCheckOpt.exists hf
    exact abs_lt_of_filter_bounds (sec s) x (eval_filterPoly_nonneg _ _) (hM s) hB hcapB hi

/-- **The masses.** Sector data that obey the table of a valid certificate have `q`th sector
moments at most the masses. -/
theorem sum_pow_le_mass_of_table (hC : C.Valid) (h : C.table.Holds (C.scheme.trace sec x)) :
    ∀ s, ∑ i ∈ sec s, x i ^ C.seed.q ≤ C.seed.mass s := by
  obtain ⟨hS, -, -, -, hmass, -⟩ := hC
  have hM := TraceBounds.momentBounds hS h
  intro s
  obtain ⟨b, hb, hbm⟩ := (hmass s).exists
  exact ((hM s C.seed.q b).2 hb).trans (by exact_mod_cast hbm)

/-- The untwisted trace of a valid certificate is the seed's sector sum. -/
theorem trace_zero_eq_sectorSum (hC : C.Valid) (n : ℕ) :
    C.scheme.trace sec x 0 n = C.seed.sectorSum sec x n := by
  simp only [SectorScheme.trace, SectorSeed.sectorSum]
  refine sum_congr rfl fun s _ ↦ ?_
  rw [← hC.2.2.1 s]; push_cast; ring

/-- **The trace bound of sector data** that obey the table, at every length `q r`, `r ≥ 1`. -/
theorem trace_le_of_table (hC : C.Valid) (hT : C.seed.TraceValid)
    (h : C.table.Holds (C.scheme.trace sec x)) {r : ℕ} (hr : 1 ≤ r) :
    C.scheme.trace sec x 0 (C.seed.q * r) ≤ C.seed.traceBound r := by
  rw [trace_zero_eq_sectorSum hC]
  exact SectorSeed.sectorSum_le_traceBound hC.2.1.1 hT (abs_lt_cap_of_table hC h)
    (sum_pow_le_mass_of_table hC h) hr

/-- **The dominant window of sector data** that obey the table, with `Z_{2n₀} > 1` (a trial
state): one eigenvalue of the dominant sector has `ℓ < |x i₀| < cap₀`, and the rest of the
spectrum has moments of order `q r` at most `restBound r`. -/
theorem exists_window_of_table [DecidableEq ι] (hC : C.Valid)
    (h : C.table.Holds (C.scheme.trace sec x)) (h2 : 1 < C.scheme.trace sec x 0 (2 * C.n₀)) :
    ∃ i₀ ∈ sec 0, (C.seed.ell : ℝ) < |x i₀| ∧ |x i₀| < C.seed.cap₀ ∧
      ∀ r, 1 ≤ r → C.seed.restSum sec x i₀ (C.seed.q * r) ≤ C.seed.restBound r := by
  refine SectorSeed.exists_window hC.2.1 (abs_lt_cap_of_table hC h)
    (sum_pow_le_mass_of_table hC h) ?_
  rw [two_mul_n₀, ← trace_zero_eq_sectorSum hC]
  exact h2

/-- **The trace bound at every length `q r`.** For a valid certificate whose sectors all admit
their capped-mass bounds, traces with the sector representation on `ns ⊇ table` that obey the
table strictly have `e^{-βκ q r} Z(q r, β) ≤ traceBound r` for every `q r ∈ ns`, `r ≥ 1`. -/
theorem shiftedTrace_le (hC : C.Valid) (hT : C.seed.TraceValid)
    {Z : Fin (C.scheme.t + 1) → ℕ → ℝ → ℝ} {ns : Finset ℕ}
    (happ : SectorApprox C.scheme (fun g n ↦ Z g n C.beta) ns)
    (hns : ∀ n ∈ C.table.lengths, n ∈ ns)
    (htab : C.table.HoldsStrict (fun g n ↦ shiftedTrace C.beta C.kappa n (Z g n C.beta)))
    {r : ℕ} (hr : 1 ≤ r) (hmem : C.seed.q * r ∈ ns) :
    shiftedTrace C.beta C.kappa (C.seed.q * r) (Z 0 (C.seed.q * r) C.beta) ≤
      C.seed.traceBound r := by
  set s := exp (-((C.beta : ℝ) * C.kappa))
  have hs0 : 0 < s := exp_pos _
  have hsh (n : ℕ) (z : ℝ) : shiftedTrace C.beta C.kappa n z = s ^ n * z :=
    shiftedTrace_eq_pow _ _ n z
  have htab' : C.table.HoldsStrict (fun g n ↦ s ^ n * Z g n C.beta) := by
    simpa only [hsh] using htab
  obtain ⟨μ, hμ, hmargin⟩ := C.table.exists_margin C.table.lengths htab'
  rw [hsh]
  refine le_of_forall_small one_pos hμ fun δ hδ hδμ ↦ ?_
  obtain ⟨ι', _, sec', x', hx⟩ := happ.scale hs0 δ hδ
  have htable : C.table.Holds (C.scheme.trace sec' x') := by
    intro n hn g b
    have a := abs_le.1 (hx n (hns n hn) g)
    have hm := hmargin n hn g b
    exact ⟨fun hb ↦ by linarith [hm.1 hb], fun hb ↦ by linarith [hm.2 hb]⟩
  have h1 := trace_le_of_table hC hT htable hr
  have c := abs_le.1 (hx _ hmem 0)
  linarith

end TableCertificate

/-! ### Physical squaring -/

/-- **Physical squaring of shifted traces.** A Boltzmann sum has `Z(2β) ≤ Z(β)²`, and the shift
factor at `2β` is the square of the one at `β`: `e^{-2βκn} Z(2β) ≤ (e^{-βκn} Z(β))²`. -/
theorem BoltzmannRep.shiftedTrace_two_mul_le_sq {f : ℝ → ℝ} (h : BoltzmannRep f) (β κ : ℝ)
    (n : ℕ) : shiftedTrace (2 * β) κ n (f (2 * β)) ≤ shiftedTrace β κ n (f β) ^ 2 := by
  have := (h.shift κ n).lift (j := 1) (l := 2) le_rfl (by norm_num) β
  simpa using this

/-! ### Power-sum approximations -/

/-- The values `Z n` are approximated by power sums: for every `δ > 0` there are finitely many
real numbers `x i` with `|∑ᵢ xᵢⁿ - Z n| ≤ δ` at every length `n ∈ ns`. The transfer eigenvalues
of a reflection-positive ring give this at every `β ≥ 0` (`RPRing.exists_transfer_approx`). -/
def PowerSumApprox (Z : ℕ → ℝ) (ns : Finset ℕ) : Prop :=
  ∀ δ : ℝ, 0 < δ → ∃ (D : Type) (_ : Fintype D) (x : D → ℝ), ∀ n ∈ ns, |∑ i, x i ^ n - Z n| ≤ δ

/-- Rescaling the values by `sⁿ` with `s > 0` (an energy shift) rescales the power sums. -/
theorem PowerSumApprox.scale {Z : ℕ → ℝ} {ns : Finset ℕ} (h : PowerSumApprox Z ns) {s : ℝ}
    (hs0 : 0 < s) : PowerSumApprox (fun n ↦ s ^ n * Z n) ns := by
  intro δ hδ
  set M := 1 + ∑ n ∈ ns, s ^ n
  have hM : 0 < M := by positivity
  have hsM {n : ℕ} (hn : n ∈ ns) : s ^ n ≤ M := by
    have := single_le_sum (fun m _ ↦ (pow_pos hs0 m).le) hn
    linarith
  obtain ⟨D, _, x, hx⟩ := h (δ / M) (div_pos hδ hM)
  refine ⟨D, inferInstance, fun i ↦ s * x i, fun n hn ↦ ?_⟩
  have e : ∑ i, (s * x i) ^ n - s ^ n * Z n = s ^ n * (∑ i, x i ^ n - Z n) := by
    simp only [mul_pow, ← mul_sum, mul_sub]
  rw [e, abs_mul, abs_of_nonneg (pow_nonneg hs0.le n)]
  calc s ^ n * |∑ i, x i ^ n - Z n| ≤ M * (δ / M) :=
        mul_le_mul (hsM hn) (hx n hn) (abs_nonneg _) hM.le
    _ = δ := by field_simp

/-! ### A dominant weight from a purity -/

/-- **A dominant weight from a purity** (Lemma "Squaring and concentration" of the boundary-field
paper). Nonnegative weights `vᵢ` with total `A > 0` and purity `∑ vᵢ² ≥ (1 - u) A²`, where
`u < 2c(1 - c)` and `0 ≤ c ≤ 1/2`, have one weight `v_{i₀} > (1 - c) A`. The largest weight `v₀`
is at least `(1 - u) A > A/2`, and the rest `y = A - v₀` has `2 y (A - y) ≤ u A²`. -/
theorem exists_gt_of_purity {ι : Type*} (s : Finset ι) (v : ι → ℝ) (hv : ∀ i ∈ s, 0 ≤ v i)
    {u c : ℝ} (hc0 : 0 ≤ c) (hc : c ≤ 1 / 2) (hu : u < 2 * c * (1 - c))
    (hA : 0 < ∑ i ∈ s, v i) (hP : (1 - u) * (∑ i ∈ s, v i) ^ 2 ≤ ∑ i ∈ s, v i ^ 2) :
    ∃ i₀ ∈ s, (1 - c) * ∑ i ∈ s, v i < v i₀ := by
  classical
  have hs : s.Nonempty := nonempty_of_sum_ne_zero hA.ne'
  obtain ⟨i₀, hi₀, hmax⟩ := s.exists_max_image v hs
  refine ⟨i₀, hi₀, ?_⟩
  set A := ∑ i ∈ s, v i
  set y := ∑ i ∈ s.erase i₀, v i
  have hv' : ∀ i ∈ s.erase i₀, 0 ≤ v i := fun i hi ↦ hv i (mem_of_mem_erase hi)
  have hy0 : 0 ≤ y := sum_nonneg hv'
  have hAy : A = v i₀ + y := (add_sum_erase s v hi₀).symm
  -- `∑ v² ≤ v₀ A` and `∑ v² ≤ v₀² + y²`.
  have h2max : ∑ i ∈ s, v i ^ 2 ≤ v i₀ * A := by
    rw [mul_sum]
    exact sum_le_sum fun i hi ↦ by
      rw [sq]; exact mul_le_mul_of_nonneg_right (hmax i hi) (hv i hi)
  have h2sum : ∑ i ∈ s, v i ^ 2 ≤ v i₀ ^ 2 + y ^ 2 := by
    rw [← add_sum_erase s (fun i ↦ v i ^ 2) hi₀]
    have := sum_pow_le_pow_sum (s.erase i₀) v hv' (n := 2) (by norm_num)
    linarith
  -- `2c(1 - c) ≤ 1/2`, so `u < 1/2` and `v₀ ≥ (1 - u) A > A/2`.
  have hu2 : u < 1 / 2 := by nlinarith [sq_nonneg (1 - 2 * c)]
  have hv₀ : (1 - u) * A ≤ v i₀ := by
    have : (1 - u) * A * A ≤ v i₀ * A := by nlinarith
    exact le_of_mul_le_mul_right this hA
  by_contra hno
  push Not at hno
  -- Then `cA ≤ y < A/2`, so `y (A - y) ≥ c(1 - c) A²`, against `2 y (A - y) ≤ u A²`.
  have hyc : c * A ≤ y := by linarith
  have hyA : y < A / 2 := by nlinarith
  have hyy : c * (1 - c) * A ^ 2 ≤ y * (A - y) := by
    nlinarith [mul_nonneg (sub_nonneg.2 hyc) (by linarith : 0 ≤ (1 - c) * A - y)]
  have hP' : 2 * y * (A - y) ≤ u * A ^ 2 := by
    have : v i₀ = A - y := by linarith
    rw [this] at h2sum
    nlinarith
  have : 2 * c * (1 - c) * A ^ 2 ≤ u * A ^ 2 := by linarith
  have hA2 : 0 < A ^ 2 := by positivity
  nlinarith

/-! ### The spatial seed at the doubled temperature -/

/-- The rational data of the spatial reseeding step at one inverse temperature: a chain row
`1 - P_s(m) ≤ p` with `m = a k`, a trace bound `Z_k ≤ B`, and the target length `n ≥ k`. -/
structure ReseedData where
  /-- The tail length `k`. -/
  k : ℕ
  /-- The chain length is `m = a k`. -/
  a : ℕ
  /-- The target length `n ≥ k` of the spatial seed. -/
  n : ℕ
  /-- The spatial defect `1 - P_s(m) ≤ p` of the traces. -/
  p : ℚ
  /-- A spatial defect `p' > p` that the approximating power sums satisfy. -/
  p' : ℚ
  /-- The dominance: one weight `x₀^m / Z_m > 1 - c`. -/
  c : ℚ
  /-- The lower bound `x₀^k > lk`, with `lk^a ≤ 1 - c`. -/
  lk : ℚ
  /-- The upper bound `Z_k ≤ B`. -/
  B : ℚ
  /-- The bound on the relative tail `∑_{i ≠ 0} (xᵢ/x₀)^k`. -/
  T : ℚ
  /-- The ratio bound `|xᵢ| ≤ ρ |x₀|`, with `T ≤ ρ^k`. -/
  rho : ℚ
  /-- The spatial defect `1 - P_s(n)` to certify. -/
  spatialDefect : ℚ

namespace ReseedData

variable (R : ReseedData)

/-- The chain length `m = a k`. -/
def m : ℕ := R.a * R.k

/-- The lengths of the power-sum approximation: `k`, `m`, `2m`, `n`, `2n`. -/
def lengths : Finset ℕ := {R.k, R.m, 2 * R.m, R.n, 2 * R.n}

/-- The rational conditions of the reseeding step. -/
def Valid : Prop :=
  Even R.k ∧ 0 < R.k ∧ 1 ≤ R.a ∧ Even R.n ∧ R.k ≤ R.n ∧ 0 ≤ R.p ∧ R.p < R.p' ∧ 0 ≤ R.c ∧
    R.c ≤ 1 / 2 ∧ R.p' < 2 * R.c * (1 - R.c) ∧ 0 < R.lk ∧ R.lk ^ R.a ≤ 1 - R.c ∧ 0 ≤ R.T ∧
    R.B < (1 + R.T) * R.lk ∧ 0 ≤ R.rho ∧ R.T ≤ R.rho ^ R.k ∧ R.spatialDefect ≤ 1 ∧
    (1 - R.spatialDefect) * (1 + R.rho ^ (R.n - R.k) * R.T) ^ 2 ≤ 1

instance : Decidable R.Valid := by unfold Valid; infer_instance

variable {R}

/-- **The dominant eigenvalue and the tail for power sums.** Real numbers `xᵢ` with
`∑ xᵢ^m ≥ 1`, `(1 - p') (∑ xᵢ^m)² ≤ ∑ xᵢ^{2m}` and `∑ xᵢ^k ≤ (1 + T) lk` have one `x₀ = x i₀`
with `x₀^m > 1 - c` and `x₀^k > lk` (the paper's `λ₁⁷² > .968`), relative tail
`∑_{i ≠ i₀} xᵢ^k ≤ T x₀^k` at length `k`, and ratios `|xᵢ| ≤ ρ |x₀|` for `i ≠ i₀` (the paper's
`ρ₁`). -/
theorem exists_dominant_of_sums (hR : R.Valid) {D : Type*} [Fintype D] [DecidableEq D]
    (x : D → ℝ) (hm : 1 ≤ ∑ i, x i ^ R.m)
    (hpur : (1 - (R.p' : ℝ)) * (∑ i, x i ^ R.m) ^ 2 ≤ ∑ i, x i ^ (2 * R.m))
    (hk : ∑ i, x i ^ R.k ≤ (1 + (R.T : ℝ)) * R.lk) :
    ∃ i₀, 1 - (R.c : ℝ) < x i₀ ^ R.m ∧ (R.lk : ℝ) < x i₀ ^ R.k ∧
      ∑ i ∈ univ.erase i₀, x i ^ R.k ≤ R.T * x i₀ ^ R.k ∧
      ∀ i ∈ univ.erase i₀, |x i| ≤ R.rho * |x i₀| := by
  obtain ⟨hke, hk0, -, -, -, -, -, hc0, hc, hp'c, hlk, hlka, hT, -, hρ, hTρ, -⟩ := hR
  have hme : Even R.m := by rw [m]; exact hke.mul_left _
  have hc' : (R.c : ℝ) ≤ 1 / 2 := by
    have := (Rat.cast_le (K := ℝ)).2 hc
    simpa using this
  -- Step 1: one weight `x₀^m > (1 - c) ∑ x^m ≥ 1 - c`.
  have hv : ∀ i ∈ (univ : Finset D), 0 ≤ x i ^ R.m := fun i _ ↦ hme.pow_nonneg _
  have hsq : ∀ i, (x i ^ R.m) ^ 2 = x i ^ (2 * R.m) := fun i ↦ by rw [← pow_mul, mul_comm]
  obtain ⟨i₀, -, hi₀⟩ := exists_gt_of_purity univ (fun i ↦ x i ^ R.m) hv
    (u := (R.p' : ℝ)) (c := (R.c : ℝ)) (by exact_mod_cast hc0)
    hc' (by exact_mod_cast hp'c) (by linarith)
    (by simp only [hsq]; exact hpur)
  have hc1 : (0 : ℝ) ≤ 1 - R.c := by linarith
  have hx0m : (1 - (R.c : ℝ)) < x i₀ ^ R.m := by nlinarith
  -- Step 2: `x₀^k > lk`.
  set P := x i₀ ^ R.k with hPdef
  have hP0 : 0 ≤ P := hke.pow_nonneg _
  have hlk' : (0 : ℝ) < R.lk := by exact_mod_cast hlk
  have hPlk : (R.lk : ℝ) < P := by
    have h1 : (R.lk : ℝ) ^ R.a ≤ 1 - R.c := by exact_mod_cast hlka
    have h2 : x i₀ ^ R.m = P ^ R.a := by rw [m, mul_comm, pow_mul]
    exact lt_of_pow_lt_pow_left₀ R.a hP0 (by linarith)
  -- Step 3: the tail at length `k` is at most `T x₀^k`, so `|xᵢ| ≤ ρ |x₀|`.
  have hsplit (n : ℕ) : ∑ i, x i ^ n = x i₀ ^ n + ∑ i ∈ univ.erase i₀, x i ^ n :=
    (add_sum_erase univ (fun i ↦ x i ^ n) (mem_univ i₀)).symm
  have hT' : (0 : ℝ) ≤ R.T := by exact_mod_cast hT
  have htailk : ∑ i ∈ univ.erase i₀, x i ^ R.k ≤ R.T * P := by
    rw [hsplit] at hk
    nlinarith
  have hρ' : (0 : ℝ) ≤ R.rho := by exact_mod_cast hρ
  have hratio : ∀ i ∈ univ.erase i₀, |x i| ≤ R.rho * |x i₀| := by
    intro i hi
    have h1 : x i ^ R.k ≤ R.T * P :=
      (single_le_sum (fun j _ ↦ hke.pow_nonneg (x j)) hi).trans htailk
    have h2 : (R.T : ℝ) ≤ (R.rho : ℝ) ^ R.k := by exact_mod_cast hTρ
    have h3 : |x i| ^ R.k ≤ (R.rho * |x i₀|) ^ R.k := by
      rw [hke.pow_abs, mul_pow, hke.pow_abs]
      calc x i ^ R.k ≤ R.T * P := h1
        _ ≤ (R.rho : ℝ) ^ R.k * P := mul_le_mul_of_nonneg_right h2 hP0
    exact le_of_pow_le_pow_left₀ hk0.ne' (by positivity) h3
  exact ⟨i₀, hx0m, hPlk, htailk, hratio⟩

/-- **The reseeding step for power sums.** Real numbers `xᵢ` with `∑ xᵢ^m ≥ 1`,
`(1 - p') (∑ xᵢ^m)² ≤ ∑ xᵢ^{2m}` and `∑ xᵢ^k ≤ (1 + T) lk` satisfy
`(1 - p_s) (∑ xᵢⁿ)² ≤ ∑ xᵢ^{2n}`. -/
theorem spatial_of_sums (hR : R.Valid) {D : Type*} [Fintype D] (x : D → ℝ)
    (hm : 1 ≤ ∑ i, x i ^ R.m)
    (hpur : (1 - (R.p' : ℝ)) * (∑ i, x i ^ R.m) ^ 2 ≤ ∑ i, x i ^ (2 * R.m))
    (hk : ∑ i, x i ^ R.k ≤ (1 + (R.T : ℝ)) * R.lk) :
    (1 - (R.spatialDefect : ℝ)) * (∑ i, x i ^ R.n) ^ 2 ≤ ∑ i, x i ^ (2 * R.n) := by
  classical
  obtain ⟨i₀, -, -, htailk, hratio⟩ := exists_dominant_of_sums hR x hm hpur hk
  obtain ⟨hke, -, -, hne, hkn, -, -, -, -, -, -, -, hT, -, hρ, -, hps1, hps⟩ := hR
  set P := x i₀ ^ R.k with hPdef
  have hsplit (n : ℕ) : ∑ i, x i ^ n = x i₀ ^ n + ∑ i ∈ univ.erase i₀, x i ^ n :=
    (add_sum_erase univ (fun i ↦ x i ^ n) (mem_univ i₀)).symm
  have hρ' : (0 : ℝ) ≤ R.rho := by exact_mod_cast hρ
  -- Step 4: the tail at length `n` is at most `ρ^{n-k} T x₀ⁿ`.
  have hnk : Even (R.n - R.k) := (Nat.even_sub hkn).2 (iff_of_true hne hke)
  set Q := (R.rho : ℝ) ^ (R.n - R.k) * |x i₀| ^ (R.n - R.k)
  have hQ0 : 0 ≤ Q := by positivity
  have htailn : ∑ i ∈ univ.erase i₀, x i ^ R.n ≤ Q * (R.T * P) := by
    calc ∑ i ∈ univ.erase i₀, x i ^ R.n
        ≤ ∑ i ∈ univ.erase i₀, Q * x i ^ R.k := by
          refine sum_le_sum fun i hi ↦ ?_
          rw [show R.n = (R.n - R.k) + R.k by omega, pow_add]
          refine mul_le_mul_of_nonneg_right ?_ (hke.pow_nonneg _)
          rw [← hnk.pow_abs]
          calc |x i| ^ (R.n - R.k) ≤ (R.rho * |x i₀|) ^ (R.n - R.k) :=
                pow_le_pow_left₀ (abs_nonneg _) (hratio i hi) _
            _ = Q := by rw [mul_pow]
      _ = Q * ∑ i ∈ univ.erase i₀, x i ^ R.k := by rw [mul_sum]
      _ ≤ Q * (R.T * P) := mul_le_mul_of_nonneg_left htailk hQ0
  have hxn : x i₀ ^ R.n = |x i₀| ^ (R.n - R.k) * P := by
    rw [hPdef, ← hke.pow_abs, ← pow_add, Nat.sub_add_cancel hkn, hne.pow_abs]
  -- Step 5: `Z_n ≤ x₀ⁿ (1 + ρ^{n-k} T)` and `Z_{2n} ≥ x₀^{2n}`.
  set X := x i₀ ^ R.n
  have hX0 : 0 ≤ X := hne.pow_nonneg _
  have hZn : ∑ i, x i ^ R.n ≤ X * (1 + (R.rho : ℝ) ^ (R.n - R.k) * R.T) := by
    rw [hsplit]
    have : Q * (R.T * P) = X * ((R.rho : ℝ) ^ (R.n - R.k) * R.T) := by
      simp only [Q, X, hxn]; ring
    linarith
  have hZn0 : 0 ≤ ∑ i, x i ^ R.n := sum_nonneg fun i _ ↦ hne.pow_nonneg _
  have hZ2n : X ^ 2 ≤ ∑ i, x i ^ (2 * R.n) := by
    rw [hsplit, show X ^ 2 = x i₀ ^ (2 * R.n) by rw [← pow_mul, mul_comm]]
    have : 0 ≤ ∑ i ∈ univ.erase i₀, x i ^ (2 * R.n) :=
      sum_nonneg fun i _ ↦ (even_two_mul R.n).pow_nonneg _
    linarith
  have hps' : (1 - (R.spatialDefect : ℝ)) * (1 + (R.rho : ℝ) ^ (R.n - R.k) * R.T) ^ 2 ≤ 1 := by
    exact_mod_cast hps
  have hps0 : 0 ≤ 1 - (R.spatialDefect : ℝ) := by
    have : (R.spatialDefect : ℝ) ≤ 1 := by exact_mod_cast hps1
    linarith
  calc (1 - (R.spatialDefect : ℝ)) * (∑ i, x i ^ R.n) ^ 2
      ≤ (1 - R.spatialDefect) * (X * (1 + (R.rho : ℝ) ^ (R.n - R.k) * R.T)) ^ 2 := by gcongr
    _ = X ^ 2 * ((1 - R.spatialDefect) * (1 + (R.rho : ℝ) ^ (R.n - R.k) * R.T) ^ 2) := by ring
    _ ≤ X ^ 2 * 1 := by gcongr
    _ ≤ ∑ i, x i ^ (2 * R.n) := by linarith

/-- **The spatial seed at the doubled temperature.** Let `Z` be approximated by power sums at
the lengths `k, m, 2m, n, 2n`. If `Z_m > 1` (a trial state), `(1 - p) Z_m² ≤ Z_{2m}` (a chain
row), `Z_k ≤ B` (physical squaring of a table bound) and `Z_n > 0`, then
`(1 - p_s) Z_n² ≤ Z_{2n}`. -/
theorem spatial (hR : R.Valid) {Z : ℕ → ℝ} (happ : PowerSumApprox Z R.lengths)
    (hm : 1 < Z R.m) (hpur : (1 - (R.p : ℝ)) * Z R.m ^ 2 ≤ Z (2 * R.m)) (hk : Z R.k ≤ R.B)
    (hn : 0 < Z R.n) : (1 - (R.spatialDefect : ℝ)) * Z R.n ^ 2 ≤ Z (2 * R.n) := by
  have hR' := hR
  obtain ⟨-, -, -, -, -, hp0, hpp, -, -, hp'c, -, -, -, hB, -, -, hps1, -⟩ := hR'
  set Zm := Z R.m
  set ε : ℝ := (R.p' - R.p : ℝ) * Zm ^ 2 / (2 * Zm + 2)
  have hpp' : (0 : ℝ) < R.p' - R.p := by
    have : (R.p : ℝ) < R.p' := by exact_mod_cast hpp
    linarith
  have hε : 0 < ε := by positivity
  set ε' : ℝ := (1 + R.T) * R.lk - R.B
  have hε' : 0 < ε' := by
    have : ((R.B : ℚ) : ℝ) < (((1 + R.T) * R.lk : ℚ) : ℝ) := by exact_mod_cast hB
    push_cast at this
    simp only [ε']
    linarith
  set δ₀ := min (min (Zm - 1) ε) (min ε' 1)
  have hδ₀ : 0 < δ₀ := lt_min (lt_min (by linarith) hε) (lt_min hε' one_pos)
  have hps1' : (R.spatialDefect : ℝ) ≤ 1 := by exact_mod_cast hps1
  refine (seed_limit (B := Z (2 * R.n)) hps1' hn (lt_min hδ₀ hn) fun δ hδ hδ₀' ↦ ?_).1
  have hδ1 : δ < Zm - 1 := hδ₀'.trans_le ((min_le_left _ _).trans ((min_le_left _ _).trans
    (min_le_left _ _)))
  have hδ2 : δ < ε := hδ₀'.trans_le ((min_le_left _ _).trans ((min_le_left _ _).trans
    (min_le_right _ _)))
  have hδ3 : δ < ε' := hδ₀'.trans_le ((min_le_left _ _).trans ((min_le_right _ _).trans
    (min_le_left _ _)))
  have hδ4 : δ < 1 := hδ₀'.trans_le ((min_le_left _ _).trans ((min_le_right _ _).trans
    (min_le_right _ _)))
  have hδ5 : δ < Z R.n := hδ₀'.trans_le (min_le_right _ _)
  obtain ⟨D, _, x, hx⟩ := happ δ hδ
  have c (j : ℕ) (hj : j ∈ R.lengths) := abs_le.1 (hx j hj)
  have ck := c R.k (by simp [lengths])
  have cm := c R.m (by simp [lengths])
  have c2m := c (2 * R.m) (by simp [lengths])
  have cn := c R.n (by simp [lengths])
  have c2n := c (2 * R.n) (by simp [lengths])
  -- The approximating power sums satisfy the hypotheses of `spatial_of_sums`.
  have hm' : 1 ≤ ∑ i, x i ^ R.m := by linarith
  have hp'1 : 0 ≤ 1 - (R.p' : ℝ) := by
    have : (R.p' : ℝ) < 2 * R.c * (1 - R.c) := by exact_mod_cast hp'c
    nlinarith [sq_nonneg (1 - 2 * (R.c : ℝ))]
  have hpur' : (1 - (R.p' : ℝ)) * (∑ i, x i ^ R.m) ^ 2 ≤ ∑ i, x i ^ (2 * R.m) := by
    have hZm : 1 < Zm := hm
    have h1 : (∑ i, x i ^ R.m) ^ 2 ≤ (Zm + δ) ^ 2 :=
      pow_le_pow_left₀ (by linarith) (by linarith) 2
    have h2 : (Zm + δ) ^ 2 ≤ Zm ^ 2 + (2 * Zm + 1) * δ := by nlinarith
    have h3 : (2 * Zm + 2) * δ ≤ (R.p' - R.p : ℝ) * Zm ^ 2 := by
      have : (2 * Zm + 2) * ε = (R.p' - R.p : ℝ) * Zm ^ 2 := by
        simp only [ε]; field_simp
      nlinarith
    have h4 : (1 - (R.p' : ℝ)) * (Zm ^ 2 + (2 * Zm + 1) * δ) ≤
        (1 - R.p) * Zm ^ 2 - (R.p' - R.p) * Zm ^ 2 + (2 * Zm + 1) * δ := by
      have hp'0 : (0 : ℝ) ≤ R.p' := by
        have : (0 : ℝ) ≤ R.p := by exact_mod_cast hp0
        linarith
      nlinarith [mul_nonneg hp'0 (mul_nonneg (by linarith : (0 : ℝ) ≤ 2 * Zm + 1) hδ.le)]
    calc (1 - (R.p' : ℝ)) * (∑ i, x i ^ R.m) ^ 2 ≤ (1 - R.p') * (Zm + δ) ^ 2 :=
          mul_le_mul_of_nonneg_left h1 hp'1
      _ ≤ (1 - R.p') * (Zm ^ 2 + (2 * Zm + 1) * δ) := mul_le_mul_of_nonneg_left h2 hp'1
      _ ≤ (1 - R.p) * Zm ^ 2 - (R.p' - R.p) * Zm ^ 2 + (2 * Zm + 1) * δ := h4
      _ ≤ ∑ i, x i ^ (2 * R.m) := by linarith
  have hk' : ∑ i, x i ^ R.k ≤ (1 + (R.T : ℝ)) * R.lk := by
    simp only [ε'] at hδ3
    linarith
  have key := spatial_of_sums hR x hm' hpur' hk'
  refine ⟨?_, by linarith⟩
  have hps0 : 0 ≤ 1 - (R.spatialDefect : ℝ) := by linarith
  calc (1 - (R.spatialDefect : ℝ)) * (Z R.n - δ) ^ 2
      ≤ (1 - R.spatialDefect) * (∑ i, x i ^ R.n) ^ 2 :=
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by linarith) (by linarith) 2) hps0
    _ ≤ ∑ i, x i ^ (2 * R.n) := key
    _ ≤ Z (2 * R.n) + δ := by linarith

end ReseedData

/-! ### The thermal seed at the doubled temperature -/

/-- **The thermal seed at the doubled temperature from one dominant Boltzmann factor.** If one
energy lies below `-κN` (a trial state), the shifted Boltzmann sum at `β > 0` is at most `q`,
`(q - 1)² < R` with `R ≥ 0`, and `(1 - p)(1 + R)² ≤ 1`, then the thermal purity at `2β` is at
least `1 - p`: `Z(4β) / Z(2β)² ≥ 1 - p`. -/
theorem thermalPurity_ge_of_dominant {ι : Type*} [Fintype ι] (E : ι → ℝ) {β κ q R p : ℝ}
    {N : ℕ} (hβ : 0 < β) (htrial : ∃ i, E i < -(κ * N))
    (hq : shiftedTrace β κ N (∑ i, exp (-(β * E i))) ≤ q) (hR : 0 ≤ R) (hqR : (q - 1) ^ 2 < R)
    (hp : (1 - p) * (1 + R) ^ 2 ≤ 1) (hp1 : p ≤ 1) :
    1 - p ≤ (∑ i, exp (-(2 * (2 * β) * E i))) / (∑ i, exp (-(2 * β * E i))) ^ 2 := by
  obtain ⟨i₀, hi₀⟩ := htrial
  have h1 : 1 < exp (-(β * (E i₀ + κ * N))) := by
    rw [← exp_zero]; exact exp_lt_exp.2 (by nlinarith)
  rw [BoltzmannRep.shiftedTrace_eq] at hq
  have hT := boltzmann_sq_le_of_dominant (fun i ↦ E i + κ * N) (j := 1) (l := 2) le_rfl
    (by norm_num) h1 hR hq (by simpa using hqR)
  have e : ((2 : ℕ) : ℝ) / ((1 : ℕ) : ℝ) * β = 2 * β := by norm_num
  rw [e, ← BoltzmannRep.shiftedTrace_eq, ← BoltzmannRep.shiftedTrace_eq,
    shiftedTrace_sq_le_iff] at hT
  have hZ : 0 < ∑ i, exp (-(2 * β * E i)) :=
    sum_pos (fun _ _ ↦ exp_pos _) ⟨i₀, mem_univ _⟩
  have hZ2 : 0 ≤ ∑ i, exp (-(2 * (2 * β) * E i)) := sum_nonneg fun _ _ ↦ (exp_pos _).le
  rw [le_div_iff₀ (by positivity)]
  have hp0 : 0 ≤ 1 - p := by linarith
  calc (1 - p) * (∑ i, exp (-(2 * β * E i))) ^ 2
      ≤ (1 - p) * ((1 + R) ^ 2 * ∑ i, exp (-(2 * (2 * β) * E i))) :=
        mul_le_mul_of_nonneg_left hT hp0
    _ = ((1 - p) * (1 + R) ^ 2) * ∑ i, exp (-(2 * (2 * β) * E i)) := by ring
    _ ≤ 1 * ∑ i, exp (-(2 * (2 * β) * E i)) := mul_le_mul_of_nonneg_right hp hZ2
    _ = _ := one_mul _

end CertifiedGaps
