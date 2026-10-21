/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.PolynomialFilter
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.Matrix.Mul
import Mathlib.Data.Rat.BigOperators
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Tactic.FieldSimp

/-!
# Sector schemes and trace tables

A symmetry of the transfer operator splits its spectrum into sectors, and twisting the ring by a
symmetry element weights each sector by a character. A `SectorScheme` records this linear
algebra: sectors `s : Fin (m + 1)` (sector `0` carries the dominant eigenvalue), twists
`g : Fin (t + 1)` (twist `0` is the untwisted ring), and the integer matrix `coeff` with

  `Z_g(n) = ∑_s coeff g s · M_s(n)`,

where `M_s(n)` is the `n`th power sum of the eigenvalues in one copy of sector `s`. Row `0` gives
the multiplicities. The rational matrix `inv` is a left inverse, `M_s = ∑_g inv s g · Z_g`, so a
table of twisted traces bounds every sector moment. The schemes used so far:

* `plain`: one sector, no twist (`Z = M₀`);
* `d2`: `D₂` with a cyclic axis permutation, `Z = J + 3N`, `Z(P) = J - N`;
* `ion`: three twists `1, P, C` (`π` and `2π/3` rotations), `Z = I + 2O + 3N`,
  `Z(P) = I + 2O - N`, `Z(C) = I - O`.

A `TraceBounds t` is a table of optional lower and upper bounds on the `t + 1` traces at finitely
many lengths; a missing side is never used (`termBound`). The sector bounds `loS`, `hiS` follow
by interval arithmetic through `inv`.

`SectorApprox S Z ns` is the spectral input: the traces `Z g n` are, to any accuracy, the
combinations `∑_s coeff g s · M_s(n)` of power sums of finitely many real numbers in the sectors.
-/

open Finset Matrix

namespace CertifiedGaps

open PolynomialFilter

/-! ### Interval bounds on linear combinations -/

/-- The upper bound on `∑ₖ w k · z k` from optional bounds `lo k ≤ z k ≤ hi k`; `none` if a
needed side is missing. -/
def linBound {k : ℕ} (w : Fin k → ℚ) (lo hi : Fin k → Option ℚ) : Option ℚ :=
  (List.finRange k).foldr
    (fun g acc ↦ (termBound (w g) (lo g) (hi g)).bind fun b ↦ acc.map (b + ·)) (some 0)

theorem linBound_spec_list {k : ℕ} (w : Fin k → ℚ) (lo hi : Fin k → Option ℚ) (z : Fin k → ℝ)
    (hlo : ∀ g b, lo g = some b → (b : ℝ) ≤ z g) (hhi : ∀ g b, hi g = some b → z g ≤ b) :
    ∀ (l : List (Fin k)) {B : ℚ},
      l.foldr (fun g acc ↦ (termBound (w g) (lo g) (hi g)).bind fun b ↦ acc.map (b + ·))
          (some 0) = some B → (l.map fun g ↦ (w g : ℝ) * z g).sum ≤ B
  | [], B, h => by simp only [List.foldr_nil, Option.some.injEq] at h; simp [← h]
  | g :: l, B, h => by
    simp only [List.foldr_cons] at h
    obtain ⟨b, hb, h'⟩ := Option.bind_eq_some_iff.1 h
    obtain ⟨r, hr, rfl⟩ := Option.map_eq_some_iff.1 h'
    have h1 := termBound_spec hb (hlo g) (hhi g)
    have h2 := linBound_spec_list w lo hi z hlo hhi l hr
    simp only [List.map_cons, List.sum_cons]
    push_cast
    linarith

/-- The combination is at most `linBound` when every value obeys its bounds. -/
theorem linBound_spec {k : ℕ} {w : Fin k → ℚ} {lo hi : Fin k → Option ℚ} {z : Fin k → ℝ}
    {B : ℚ} (h : linBound w lo hi = some B) (hlo : ∀ g b, lo g = some b → (b : ℝ) ≤ z g)
    (hhi : ∀ g b, hi g = some b → z g ≤ b) : ∑ g, (w g : ℝ) * z g ≤ B := by
  rw [Fin.sum_univ_def]
  exact linBound_spec_list w lo hi z hlo hhi _ h

/-! ### Sector schemes -/

/-- A sector scheme: `m + 1` sectors (sector `0` dominant), `t + 1` twists (twist `0`
untwisted), the integer matrix `coeff` of `Z_g = ∑_s coeff g s · M_s`, and a rational left
inverse `inv`. -/
structure SectorScheme where
  /-- The number of nontrivial sectors. -/
  m : ℕ
  /-- The number of nontrivial twists. -/
  t : ℕ
  /-- The sector content of each twisted trace. -/
  coeff : Fin (t + 1) → Fin (m + 1) → ℤ
  /-- A left inverse of `coeff`. -/
  inv : Fin (m + 1) → Fin (t + 1) → ℚ

namespace SectorScheme

variable (S : SectorScheme)

/-- The rational conditions: `inv` is a left inverse of `coeff`, the dominant sector appears
once in the untwisted trace, and every multiplicity is nonnegative. -/
def Valid : Prop :=
  (∀ s s', ∑ g, S.inv s g * S.coeff g s' = if s = s' then 1 else 0) ∧
    S.coeff 0 0 = 1 ∧ ∀ s, 0 ≤ S.coeff 0 s

instance : Decidable S.Valid := by unfold Valid; infer_instance

/-- One sector and no twist. -/
def plain : SectorScheme := ⟨0, 0, ![![1]], ![![1]]⟩

/-- `D₂` sectors from the untwisted and one `π`-twisted trace: `Z = J + 3N`, `Z(P) = J - N`. -/
def d2 : SectorScheme :=
  ⟨1, 1, ![![1, 3], ![1, -1]], ![![1 / 4, 3 / 4], ![1 / 4, -1 / 4]]⟩

/-- Sectors `I, O, N` from the twists `1, P, C`: `Z = I + 2O + 3N`, `Z(P) = I + 2O - N`,
`Z(C) = I - O`. -/
def ion : SectorScheme :=
  ⟨2, 2, ![![1, 2, 3], ![1, 2, -1], ![1, -1, 0]],
    ![![1 / 12, 3 / 12, 8 / 12], ![1 / 12, 3 / 12, -4 / 12], ![1 / 4, -1 / 4, 0]]⟩

theorem plain_valid : plain.Valid := by decide +kernel
theorem d2_valid : d2.Valid := by decide +kernel
theorem ion_valid : ion.Valid := by decide +kernel

/-- The multiplicity of a sector in the untwisted trace. -/
def mult (s : Fin (S.m + 1)) : ℤ := S.coeff 0 s

/-- The traces of sector data: `∑_s coeff g s · ∑_{i ∈ sec s} x i ^ n`. -/
def trace {ι : Type*} (sec : Fin (S.m + 1) → Finset ι) (x : ι → ℝ) (g : Fin (S.t + 1))
    (n : ℕ) : ℝ :=
  ∑ s, (S.coeff g s : ℝ) * ∑ i ∈ sec s, x i ^ n

/-- The untwisted trace of the scheme `plain` is the power sum of its one sector. -/
theorem plain_trace {ι : Type*} (sec : Fin (plain.m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) :
    plain.trace sec x 0 n = ∑ i ∈ sec 0, x i ^ n := by
  change ∑ s : Fin 1, ((![![(1 : ℤ)]] 0 s : ℤ) : ℝ) * ∑ i ∈ sec s, x i ^ n = _
  simp; rfl

/-- The untwisted trace of the scheme `d2`: `J + 3N`. -/
theorem d2_trace_zero {ι : Type*} (sec : Fin (d2.m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) :
    d2.trace sec x 0 n = ∑ i ∈ sec 0, x i ^ n + 3 * ∑ i ∈ sec 1, x i ^ n := by
  change ∑ s : Fin 2, ((![![(1 : ℤ), 3], ![1, -1]] 0 s : ℤ) : ℝ) * ∑ i ∈ sec s, x i ^ n = _
  simp [Fin.sum_univ_two]; rfl

/-- The twisted trace of the scheme `d2`: `J - N`. -/
theorem d2_trace_one {ι : Type*} (sec : Fin (d2.m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) :
    d2.trace sec x 1 n = ∑ i ∈ sec 0, x i ^ n - ∑ i ∈ sec 1, x i ^ n := by
  change ∑ s : Fin 2, ((![![(1 : ℤ), 3], ![1, -1]] 1 s : ℤ) : ℝ) * ∑ i ∈ sec s, x i ^ n = _
  simp [Fin.sum_univ_two]; rfl

/-- The untwisted trace of the scheme `ion`: `I + 2O + 3N`. -/
theorem ion_trace_zero {ι : Type*} (sec : Fin (ion.m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) :
    ion.trace sec x 0 n =
      ∑ i ∈ sec 0, x i ^ n + 2 * ∑ i ∈ sec 1, x i ^ n + 3 * ∑ i ∈ sec 2, x i ^ n := by
  change ∑ s : Fin 3, ((![![(1 : ℤ), 2, 3], ![1, 2, -1], ![1, -1, 0]] 0 s : ℤ) : ℝ) *
    ∑ i ∈ sec s, x i ^ n = _
  simp [Fin.sum_univ_three]; rfl

/-- The `P`-twisted trace of the scheme `ion`: `I + 2O - N`. -/
theorem ion_trace_one {ι : Type*} (sec : Fin (ion.m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) :
    ion.trace sec x 1 n =
      ∑ i ∈ sec 0, x i ^ n + 2 * ∑ i ∈ sec 1, x i ^ n - ∑ i ∈ sec 2, x i ^ n := by
  change ∑ s : Fin 3, ((![![(1 : ℤ), 2, 3], ![1, 2, -1], ![1, -1, 0]] 1 s : ℤ) : ℝ) *
    ∑ i ∈ sec s, x i ^ n = _
  simp [Fin.sum_univ_three]; ring_nf; rfl

/-- The `C`-twisted trace of the scheme `ion`: `I - O`. -/
theorem ion_trace_two {ι : Type*} (sec : Fin (ion.m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) :
    ion.trace sec x 2 n = ∑ i ∈ sec 0, x i ^ n - ∑ i ∈ sec 1, x i ^ n := by
  change ∑ s : Fin 3, ((![![(1 : ℤ), 2, 3], ![1, 2, -1], ![1, -1, 0]] 2 s : ℤ) : ℝ) *
    ∑ i ∈ sec s, x i ^ n = _
  simp [Fin.sum_univ_three]; ring_nf; rfl


/-- `ion_trace_zero` with the sectors indexed by `Fin 3`, the form that `rw` matches for families
written over `Fin 3`. -/
theorem ion_trace_zero_fin3 {ι : Type*} (sec : Fin 3 → Finset ι) (x : ι → ℝ) (n : ℕ) :
    ion.trace sec x (0 : Fin 3) n =
      ∑ i ∈ sec 0, x i ^ n + 2 * ∑ i ∈ sec 1, x i ^ n + 3 * ∑ i ∈ sec 2, x i ^ n :=
  ion_trace_zero sec x n

/-- `ion_trace_one` with the sectors indexed by `Fin 3`. -/
theorem ion_trace_one_fin3 {ι : Type*} (sec : Fin 3 → Finset ι) (x : ι → ℝ) (n : ℕ) :
    ion.trace sec x (1 : Fin 3) n =
      ∑ i ∈ sec 0, x i ^ n + 2 * ∑ i ∈ sec 1, x i ^ n - ∑ i ∈ sec 2, x i ^ n :=
  ion_trace_one sec x n

/-- `ion_trace_two` with the sectors indexed by `Fin 3`. -/
theorem ion_trace_two_fin3 {ι : Type*} (sec : Fin 3 → Finset ι) (x : ι → ℝ) (n : ℕ) :
    ion.trace sec x (2 : Fin 3) n = ∑ i ∈ sec 0, x i ^ n - ∑ i ∈ sec 1, x i ^ n :=
  ion_trace_two sec x n

/-- The inverse recovers the sector moments from the traces. -/
theorem sum_inv_mul_trace (hS : S.Valid) {ι : Type*} (sec : Fin (S.m + 1) → Finset ι)
    (x : ι → ℝ) (s : Fin (S.m + 1)) (n : ℕ) :
    ∑ g, (S.inv s g : ℝ) * S.trace sec x g n = ∑ i ∈ sec s, x i ^ n := by
  have key (s' : Fin (S.m + 1)) :
      ∑ g, (S.inv s g : ℝ) * (S.coeff g s' : ℝ) = if s = s' then 1 else 0 := by
    have h := hS.1 s s'
    have h' := congrArg (Rat.cast : ℚ → ℝ) h
    rw [Rat.cast_sum] at h'
    push_cast at h'
    rw [h']
    split_ifs <;> simp
  set M : Fin (S.m + 1) → ℝ := fun s' ↦ ∑ i ∈ sec s', x i ^ n with hM
  have ht (g : Fin (S.t + 1)) : S.trace sec x g n = ∑ s', (S.coeff g s' : ℝ) * M s' := rfl
  calc ∑ g, (S.inv s g : ℝ) * S.trace sec x g n
      = ∑ g, ∑ s', (S.inv s g : ℝ) * ((S.coeff g s' : ℝ) * M s') := by simp only [ht, mul_sum]
    _ = ∑ s', ∑ g, (S.inv s g : ℝ) * ((S.coeff g s' : ℝ) * M s') := sum_comm
    _ = ∑ s', (∑ g, (S.inv s g : ℝ) * (S.coeff g s' : ℝ)) * M s' := by
        simp only [sum_mul, mul_assoc]
    _ = M s := by simp [key]

end SectorScheme

/-! ### Approximation by sector power sums -/

/-- The traces `Z g n` are approximated by sector power sums: for every `δ > 0` there are finitely
many real numbers `x i`, split into sectors, with `|S.trace sec x g n - Z g n| ≤ δ` for every twist
`g` and length `n ∈ ns`. -/
def SectorApprox (S : SectorScheme) (Z : Fin (S.t + 1) → ℕ → ℝ) (ns : Finset ℕ) : Prop :=
  ∀ δ : ℝ, 0 < δ → ∃ (ι : Type) (_ : Fintype ι) (sec : Fin (S.m + 1) → Finset ι) (x : ι → ℝ),
    ∀ n ∈ ns, ∀ g, |S.trace sec x g n - Z g n| ≤ δ

namespace SectorApprox

variable {S : SectorScheme} {Z : Fin (S.t + 1) → ℕ → ℝ} {ns : Finset ℕ}

/-- An approximation at the lengths `ns` is one at every subset of them. -/
theorem mono (h : SectorApprox S Z ns) {ns' : Finset ℕ} (hsub : ns' ⊆ ns) :
    SectorApprox S Z ns' := fun δ hδ ↦ by
  obtain ⟨ι, i, sec, x, hx⟩ := h δ hδ
  exact ⟨ι, i, sec, x, fun n hn ↦ hx n (hsub hn)⟩

/-- Rescaling the traces by `sⁿ` with `s > 0` (an energy shift by `κ` per bond, `s = e^{-βκ}`)
rescales the sector data by `s`: the error `δ` is traded for `δ / M` with `M ≥ sⁿ` on `ns`. -/
theorem scale (h : SectorApprox S Z ns) {s : ℝ} (hs0 : 0 < s) :
    SectorApprox S (fun g n ↦ s ^ n * Z g n) ns := by
  intro δ hδ
  set M := 1 + ∑ n ∈ ns, s ^ n
  have hM : 0 < M := by positivity
  have hsM {n : ℕ} (hn : n ∈ ns) : s ^ n ≤ M := by
    have := single_le_sum (fun m _ ↦ (pow_pos hs0 m).le) hn
    linarith
  obtain ⟨ι, _, sec, x, hx⟩ := h (δ / M) (div_pos hδ hM)
  refine ⟨ι, inferInstance, sec, fun i ↦ s * x i, fun n hn g ↦ ?_⟩
  have e : S.trace sec (fun i ↦ s * x i) g n - s ^ n * Z g n =
      s ^ n * (S.trace sec x g n - Z g n) := by
    simp only [SectorScheme.trace, mul_pow, ← mul_sum, mul_sub]
    rw [mul_sum]
    congr 1
    refine sum_congr rfl fun _ _ ↦ by ring
  rw [e, abs_mul, abs_of_nonneg (pow_nonneg hs0.le n)]
  calc s ^ n * |S.trace sec x g n - Z g n| ≤ M * (δ / M) :=
        mul_le_mul (hsM hn) (hx n hn g) (abs_nonneg _) hM.le
    _ = δ := by field_simp

end SectorApprox

/-! ### Trace tables -/

/-- A table of optional lower and upper bounds on the `t + 1` traces at the lengths `lengths`. -/
structure TraceBounds (t : ℕ) where
  /-- The lengths of the table. -/
  lengths : List ℕ
  /-- Lower bounds on the traces. -/
  lo : Fin (t + 1) → ℕ → Option ℚ
  /-- Upper bounds on the traces. -/
  hi : Fin (t + 1) → ℕ → Option ℚ

namespace TraceBounds

variable {t : ℕ} (T : TraceBounds t)

/-- The symmetric table `c g n ± ε` at the lengths `ns`. -/
def centered (c : Fin (t + 1) → ℕ → ℚ) (ns : List ℕ) (ε : ℚ) : TraceBounds t where
  lengths := ns
  lo g n := some (c g n - ε)
  hi g n := some (c g n + ε)

/-- Values `Z g n` obey the table at every length. -/
def Holds (Z : Fin (t + 1) → ℕ → ℝ) : Prop :=
  ∀ n ∈ T.lengths, ∀ g (b : ℚ), (T.lo g n = some b → (b : ℝ) ≤ Z g n) ∧
    (T.hi g n = some b → Z g n ≤ b)

/-- Values `Z g n` obey the table strictly at every length. -/
def HoldsStrict (Z : Fin (t + 1) → ℕ → ℝ) : Prop :=
  ∀ n ∈ T.lengths, ∀ g (b : ℚ), (T.lo g n = some b → (b : ℝ) < Z g n) ∧
    (T.hi g n = some b → Z g n < b)

/-- A symmetric table within `ε` holds strictly for values within `εt < ε` of the centers. -/
theorem holdsStrict_centered {c : Fin (t + 1) → ℕ → ℚ} {ns : List ℕ} {ε εt : ℚ} (hε : εt < ε)
    {Z : Fin (t + 1) → ℕ → ℝ} (h : ∀ n ∈ ns, ∀ g, |Z g n - c g n| ≤ εt) :
    (centered c ns ε).HoldsStrict Z := by
  intro n hn g b
  have a := abs_le.1 (h n hn g)
  have hε' : (εt : ℝ) < ε := by exact_mod_cast hε
  simp only [centered, Option.some.injEq]
  refine ⟨?_, ?_⟩ <;> rintro rfl <;> push_cast <;> linarith

/-- The table bounds at one length hold with margin `δ`. -/
def HoldsWith (Z : Fin (t + 1) → ℕ → ℝ) (δ : ℝ) (n : ℕ) : Prop :=
  ∀ g (b : ℚ), (T.lo g n = some b → (b : ℝ) + δ ≤ Z g n) ∧ (T.hi g n = some b → Z g n + δ ≤ b)

/-- A strict optional bound has a positive margin. -/
theorem exists_margin_opt (o : Option ℚ) (z : ℝ) (sgn : ℝ) 
    (h : ∀ b, o = some b → sgn * b < sgn * z) :
    ∃ δ > 0, ∀ b, o = some b → sgn * b + δ ≤ sgn * z := by
  cases o with
  | none => exact ⟨1, one_pos, fun _ hb ↦ by cases hb⟩
  | some b =>
    have := h b rfl
    exact ⟨sgn * z - sgn * b, by linarith, fun c hc ↦ by cases hc; linarith⟩

/-- A table that holds strictly holds with a positive margin at finitely many lengths. -/
theorem exists_margin {Z : Fin (t + 1) → ℕ → ℝ} :
    ∀ l : List ℕ, (∀ n ∈ l, ∀ g (b : ℚ), (T.lo g n = some b → (b : ℝ) < Z g n) ∧
        (T.hi g n = some b → Z g n < b)) → ∃ δ > 0, ∀ n ∈ l, T.HoldsWith Z δ n
  | [], _ => ⟨1, one_pos, fun _ hn ↦ by cases hn⟩
  | n :: l, h => by
    obtain ⟨δ, hδ, hl⟩ := exists_margin l fun k hk ↦ h k (List.mem_cons_of_mem n hk)
    have hn := h n List.mem_cons_self
    -- One margin per twist and side, then the minimum over the finitely many twists.
    have hg : ∀ g : Fin (t + 1), ∃ ε > 0, ∀ b : ℚ,
        (T.lo g n = some b → (b : ℝ) + ε ≤ Z g n) ∧ (T.hi g n = some b → Z g n + ε ≤ b) := by
      intro g
      obtain ⟨ε₁, h₁, e₁⟩ := exists_margin_opt (T.lo g n) (Z g n) 1
        fun b hb ↦ by simpa using (hn g b).1 hb
      obtain ⟨ε₂, h₂, e₂⟩ := exists_margin_opt (T.hi g n) (Z g n) (-1)
        fun b hb ↦ by have := (hn g b).2 hb; linarith
      refine ⟨min ε₁ ε₂, lt_min h₁ h₂, fun b ↦ ⟨fun hb ↦ ?_, fun hb ↦ ?_⟩⟩
      · have := e₁ b hb; have := min_le_left ε₁ ε₂; linarith
      · have := e₂ b hb; have := min_le_right ε₁ ε₂; linarith
    choose ε hε hεb using hg
    set ε₀ := univ.inf' univ_nonempty ε
    have hε₀ : 0 < ε₀ := (lt_inf'_iff _).2 fun g _ ↦ hε g
    have hε₀g (g) : ε₀ ≤ ε g := inf'_le _ (mem_univ g)
    refine ⟨min δ ε₀, lt_min hδ hε₀, fun k hk ↦ ?_⟩
    rcases List.mem_cons.1 hk with rfl | hk
    · intro g b
      have := min_le_right δ ε₀
      have := hε₀g g
      exact ⟨fun hb ↦ by linarith [(hεb g b).1 hb], fun hb ↦ by linarith [(hεb g b).2 hb]⟩
    · intro g b
      have := min_le_left δ ε₀
      exact ⟨fun hb ↦ by linarith [(hl k hk g b).1 hb], fun hb ↦ by linarith [(hl k hk g b).2 hb]⟩

/-- The upper bound on the moment of sector `s` at length `n` from the table. -/
def hiS (S : SectorScheme) (T : TraceBounds S.t) (s : Fin (S.m + 1)) (n : ℕ) : Option ℚ :=
  if n ∈ T.lengths then linBound (S.inv s) (fun g ↦ T.lo g n) (fun g ↦ T.hi g n) else none

/-- The lower bound on the moment of sector `s` at length `n` from the table. -/
def loS (S : SectorScheme) (T : TraceBounds S.t) (s : Fin (S.m + 1)) (n : ℕ) : Option ℚ :=
  if n ∈ T.lengths then
    (linBound (fun g ↦ -S.inv s g) (fun g ↦ T.lo g n) (fun g ↦ T.hi g n)).map Neg.neg
  else none

end TraceBounds

/-- **Sector moments from the table.** If the traces of sector data obey a table, each sector's
power sums obey the bounds `loS`, `hiS`. -/
theorem TraceBounds.momentBounds {S : SectorScheme} (hS : S.Valid) {T : TraceBounds S.t}
    {ι : Type*} {sec : Fin (S.m + 1) → Finset ι} {x : ι → ℝ} (h : T.Holds (S.trace sec x))
    (s : Fin (S.m + 1)) :
    MomentBounds (fun n ↦ ∑ i ∈ sec s, x i ^ n) (T.loS S s) (T.hiS S s) := by
  intro n b
  simp only [TraceBounds.loS, TraceBounds.hiS]
  split_ifs with hn
  · have hlo : ∀ g (c : ℚ), T.lo g n = some c → (c : ℝ) ≤ S.trace sec x g n :=
      fun g c hc ↦ (h n hn g c).1 hc
    have hhi : ∀ g (c : ℚ), T.hi g n = some c → S.trace sec x g n ≤ c :=
      fun g c hc ↦ (h n hn g c).2 hc
    have e := S.sum_inv_mul_trace hS sec x s n
    constructor
    · intro hb
      obtain ⟨B, hB, rfl⟩ := Option.map_eq_some_iff.1 hb
      have := linBound_spec hB hlo hhi
      simp only [Rat.cast_neg, neg_mul, sum_neg_distrib] at this
      push_cast
      linarith
    · intro hb
      have := linBound_spec hb hlo hhi
      linarith
  · simp

end CertifiedGaps
