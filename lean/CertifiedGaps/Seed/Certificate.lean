/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Dominant
import CertifiedGaps.Seed.Limits
import CertifiedGaps.Seed.Scheme

/-!
# Seed certificates from trace tables

A `TableCertificate` is the finite rational data of a seed: a `SectorScheme` (which twisted
traces are tabulated and how they split into sectors), a table of bounds on those traces
(`TraceBounds`), an optional polynomial filter `X^(2a) q²` per sector, a `SectorSeed`, and the
physical parameters `β`, `κ` and the thermal defect. `TableCertificate.Valid` collects every
rational condition; it is decidable, so `decide` checks a concrete certificate.

For a model whose twisted partition functions `Z g n β` have the sector representation
(`SectorApprox`) at `β`, whose shifted traces `e^{-βκn} Z g n β` obey the table strictly, and
whose trial state gives `e^{-βκ·2n₀} Z₀(2n₀, β) > 1` at `β` and `2β`,
`TableCertificate.purities` gives the seed purities

  `1 - p ≤ P_s(n₀, β)`  and  `1 - q ≤ P_τ(2n₀, β)`.

The analytic input in `D₂` form, `D2Approx Z ZP ns` (provided by `SingleSite.D2Symmetry`), is the
`SectorApprox` of the scheme `SectorScheme.d2` (`D2Approx.toSectorApprox`).

## Main results

* `TableCertificate.seed_of_table`: the seed for sector data that obey the table.
* `TableCertificate.purities`: the seed purities of the traces.
-/

open Finset

namespace CertifiedGaps

open PolynomialFilter

/-! ### The `D₂` form of the spectral input -/

/-- The traces `Z n` (untwisted) and `ZP n` (twisted by one `π` rotation) are approximated by
`D₂` sector power sums: for every `δ > 0` there are real numbers `x i` in two finite index sets
with `|J_n + 3 N_n - Z n| ≤ δ` and `|J_n - N_n - ZP n| ≤ δ` at every length `n ∈ ns`. -/
def D2Approx (Z ZP : ℕ → ℝ) (ns : Finset ℕ) : Prop :=
  ∀ δ : ℝ, 0 < δ → ∃ (ι : Type) (_ : Fintype ι) (SJ SN : Finset ι) (x : ι → ℝ), ∀ n ∈ ns,
    |∑ i ∈ SJ, x i ^ n + 3 * ∑ i ∈ SN, x i ^ n - Z n| ≤ δ ∧
      |∑ i ∈ SJ, x i ^ n - ∑ i ∈ SN, x i ^ n - ZP n| ≤ δ

/-- `D2Approx` is the sector representation of the scheme `SectorScheme.d2`. -/
theorem D2Approx.toSectorApprox {Z ZP : ℕ → ℝ} {ns : Finset ℕ} (h : D2Approx Z ZP ns) :
    SectorApprox SectorScheme.d2 ![Z, ZP] ns := by
  intro δ hδ
  obtain ⟨ι, _, SJ, SN, x, hx⟩ := h δ hδ
  refine ⟨ι, inferInstance, ![SJ, SN], x, fun n hn g ↦ ?_⟩
  obtain ⟨h1, h2⟩ := hx n hn
  fin_cases g
  · simpa [SectorScheme.d2_trace_zero] using h1
  · show |SectorScheme.d2.trace ![SJ, SN] x 1 n - ![Z, ZP] 1 n| ≤ δ
    simpa [SectorScheme.d2_trace_one] using h2

/-! ### Certificates -/

/-- An optional upper bound that exists and is at most `b`. -/
def OptLE : Option ℚ → ℚ → Prop
  | some a, b => a ≤ b
  | none, _ => False

instance : (o : Option ℚ) → (b : ℚ) → Decidable (OptLE o b)
  | some a, b => inferInstanceAs (Decidable (a ≤ b))
  | none, _ => isFalse id

theorem OptLE.exists {o : Option ℚ} {b : ℚ} (h : OptLE o b) : ∃ a, o = some a ∧ a ≤ b := by
  cases o with
  | none => exact h.elim
  | some a => exact ⟨a, rfl, h⟩

/-- The rational data and parameters of a seed from a trace table. -/
structure TableCertificate where
  /-- The tabulated twists and their sector content. -/
  scheme : SectorScheme
  /-- The table of bounds on the shifted twisted traces. -/
  table : TraceBounds scheme.t
  /-- The filter `X^(2a) q²` of each sector, as `(a, q)`, if the sector is capped. -/
  filter : Fin (scheme.m + 1) → Option (ℕ × List ℚ)
  /-- The dominant-sector seed data. -/
  seed : SectorSeed scheme.m
  /-- The inverse temperature `β` of the table. -/
  beta : ℚ
  /-- The energy shift `κ` per bond of the table. -/
  kappa : ℚ
  /-- The thermal defect `q` to certify. -/
  thermalDefect : ℚ

namespace TableCertificate

variable (C : TableCertificate)

/-- The seed length `n₀ = q a`. -/
def n₀ : ℕ := C.seed.q * C.seed.a

/-- The lengths at which the sector representation is needed: the table, `n₀` and `2 n₀`. -/
def lengths : Finset ℕ := C.table.lengths.toFinset ∪ {C.n₀, 2 * C.n₀}

/-- The filter of a sector, with the table's bounds on its moments, caps it at the seed's cap;
a sector has a filter exactly when it has a cap. -/
def FilterCheck (s : Fin (C.scheme.m + 1)) : Prop :=
  match C.filter s, C.seed.cap s with
  | some f, some v => CapCheckOpt (filterPoly f.1 f.2) v
      (pairingBound (filterPoly f.1 f.2) (C.table.loS C.scheme s) (C.table.hiS C.scheme s))
  | none, none => True
  | _, _ => False

instance (s : Fin (C.scheme.m + 1)) : Decidable (C.FilterCheck s) := by
  unfold FilterCheck; split <;> infer_instance

/-- The rational conditions of a certificate. -/
def Valid : Prop :=
  C.scheme.Valid ∧ C.seed.Valid ∧ (∀ s, (C.seed.mult s : ℤ) = C.scheme.coeff 0 s) ∧
    (∀ s, C.FilterCheck s) ∧ (∀ s, OptLE (C.table.hiS C.scheme s C.seed.q) (C.seed.mass s)) ∧
    (∀ n ∈ C.table.lengths, 2 ≤ n) ∧ 2 ≤ C.n₀ ∧ 0 < C.beta ∧ C.thermalDefect ≤ 1 ∧
    (1 - C.thermalDefect) * C.seed.thermalBound ^ 2 ≤ 1

instance : Decidable C.Valid := by unfold Valid; infer_instance

theorem two_mul_n₀ : C.seed.q * (2 * C.seed.a) = 2 * C.n₀ := by rw [n₀]; ring

variable {C}

/-- Every length of the sector representation is at least `2`. -/
theorem two_le_of_mem_lengths (hC : C.Valid) {n : ℕ} (hn : n ∈ C.lengths) : 2 ≤ n := by
  obtain ⟨-, -, -, -, -, hlen, hn₀, -⟩ := hC
  rcases mem_union.1 hn with h | h
  · exact hlen n (List.mem_toFinset.1 h)
  · simp only [mem_insert, mem_singleton] at h
    rcases h with rfl | rfl <;> omega

/-- **The seed for sector data that obey the table.** -/
theorem seed_of_table (hC : C.Valid) {ι : Type*} {sec : Fin (C.scheme.m + 1) → Finset ι}
    {x : ι → ℝ} (h : C.table.Holds (C.scheme.trace sec x))
    (h2 : 1 < C.scheme.trace sec x 0 (2 * C.n₀)) :
    (1 - C.seed.spatialDefect : ℝ) * C.scheme.trace sec x 0 C.n₀ ^ 2 ≤
        C.scheme.trace sec x 0 (2 * C.n₀) ∧
      C.scheme.trace sec x 0 (2 * C.n₀) ≤ C.seed.thermalBound := by
  obtain ⟨hS, hseed, hmult, hfilter, hmass, -⟩ := hC
  have hM := TraceBounds.momentBounds hS h
  -- The untwisted trace is the sector sum with the seed's multiplicities.
  have htr (n : ℕ) : C.scheme.trace sec x 0 n = C.seed.sectorSum sec x n := by
    simp only [SectorScheme.trace, SectorSeed.sectorSum]
    refine sum_congr rfl fun s _ ↦ ?_
    rw [← hmult s]; push_cast; ring
  -- Caps from the filters.
  have hcap : ∀ s v, C.seed.cap s = some v → ∀ i ∈ sec s, |x i| < v := by
    intro s v hv i hi
    have hf := hfilter s
    unfold FilterCheck at hf
    rw [hv] at hf
    rcases hfs : C.filter s with _ | ⟨a, qf⟩
    · rw [hfs] at hf; exact hf.elim
    · rw [hfs] at hf
      obtain ⟨B, hB, hcapB⟩ := CapCheckOpt.exists hf
      exact abs_lt_of_filter_bounds (sec s) x (eval_filterPoly_nonneg _ _) (hM s) hB hcapB hi
  -- Masses from the table.
  have hmass' : ∀ s, ∑ i ∈ sec s, x i ^ C.seed.q ≤ C.seed.mass s := by
    intro s
    obtain ⟨b, hb, hbm⟩ := (hmass s).exists
    exact ((hM s C.seed.q b).2 hb).trans (by exact_mod_cast hbm)
  have := SectorSeed.seed hseed hcap hmass' (by rw [two_mul_n₀, ← htr]; exact h2)
  rw [two_mul_n₀] at this
  simp only [← htr] at this
  exact this

/-- **The seed purities from a certificate.** For a model whose twisted partition functions
`Z g n β` have the sector representation at `β`, whose shifted traces obey the table strictly,
with `Z₀(n₀, β) > 0` and the trial-state bounds at `β` and `2β`:
`1 - p ≤ Z₀(2n₀, β) / Z₀(n₀, β)²` and `1 - q ≤ Z₀(2n₀, 2β) / Z₀(2n₀, β)²`. -/
theorem purities (hC : C.Valid) {Z : Fin (C.scheme.t + 1) → ℕ → ℝ → ℝ}
    (happ : SectorApprox C.scheme (fun g n ↦ Z g n C.beta) C.lengths)
    (htab : C.table.HoldsStrict (fun g n ↦ shiftedTrace C.beta C.kappa n (Z g n C.beta)))
    (hpos : 0 < Z 0 C.n₀ C.beta)
    (htrial : 1 < shiftedTrace C.beta C.kappa (2 * C.n₀) (Z 0 (2 * C.n₀) C.beta))
    (htrial₂ : 1 < shiftedTrace (2 * C.beta) C.kappa (2 * C.n₀) (Z 0 (2 * C.n₀) (2 * C.beta))) :
    1 - (C.seed.spatialDefect : ℝ) ≤ Z 0 (2 * C.n₀) C.beta / Z 0 C.n₀ C.beta ^ 2 ∧
      1 - (C.thermalDefect : ℝ) ≤ Z 0 (2 * C.n₀) (2 * C.beta) / Z 0 (2 * C.n₀) C.beta ^ 2 := by
  have hC' := hC
  obtain ⟨-, hseed, -, -, -, -, -, -, hq1, hq⟩ := hC'
  set β : ℝ := (C.beta : ℝ)
  set κ : ℝ := (C.kappa : ℝ)
  set s := Real.exp (-(β * κ))
  have hs0 : 0 < s := Real.exp_pos _
  have hsh (n : ℕ) (z : ℝ) : shiftedTrace β κ n z = s ^ n * z := shiftedTrace_eq_pow β κ n z
  set Zs : Fin (C.scheme.t + 1) → ℕ → ℝ := fun g n ↦ s ^ n * Z g n β
  have hp1 : (C.seed.spatialDefect : ℝ) ≤ 1 := by exact_mod_cast hseed.2.2.2.2.1
  have htab' : C.table.HoldsStrict Zs := by simpa only [hsh] using htab
  obtain ⟨μ, hμ, hmargin⟩ := C.table.exists_margin C.table.lengths htab'
  have h₁ : C.n₀ ∈ C.lengths := mem_union_right _ (by simp)
  have h₂ : 2 * C.n₀ ∈ C.lengths := mem_union_right _ (by simp)
  have hone : 1 < Zs 0 (2 * C.n₀) := by simpa only [hsh] using htrial
  have hpos' : 0 < Zs 0 C.n₀ := mul_pos (pow_pos hs0 _) hpos
  -- For each small `δ`, the approximating sector data obey the table and the seed.
  have happrox : ∀ δ, 0 < δ → δ < min μ (min (Zs 0 (2 * C.n₀) - 1) (Zs 0 C.n₀)) →
      (1 - (C.seed.spatialDefect : ℝ)) * (Zs 0 C.n₀ - δ) ^ 2 ≤ Zs 0 (2 * C.n₀) + δ ∧
        Zs 0 (2 * C.n₀) ≤ C.seed.thermalBound + δ := by
    intro δ hδ hδ₀
    have hδ1 : δ < μ := lt_of_lt_of_le hδ₀ (min_le_left _ _)
    have hδ2 : δ < Zs 0 (2 * C.n₀) - 1 :=
      lt_of_lt_of_le hδ₀ ((min_le_right _ _).trans (min_le_left _ _))
    have hδ3 : δ < Zs 0 C.n₀ := lt_of_lt_of_le hδ₀ ((min_le_right _ _).trans (min_le_right _ _))
    obtain ⟨ι, _, sec, x, hx⟩ := happ.scale hs0 δ hδ
    have htable : C.table.Holds (C.scheme.trace sec x) := by
      intro n hn g b
      have a := abs_le.1 (hx n (mem_union_left _ (List.mem_toFinset.2 hn)) g)
      have hm := hmargin n hn g b
      exact ⟨fun hb ↦ by linarith [hm.1 hb], fun hb ↦ by linarith [hm.2 hb]⟩
    have c₁ := abs_le.1 (hx _ h₁ 0)
    have c₂ := abs_le.1 (hx _ h₂ 0)
    obtain ⟨hp, hq'⟩ := seed_of_table hC htable (by linarith)
    refine ⟨?_, by linarith⟩
    have hp0 : 0 ≤ 1 - (C.seed.spatialDefect : ℝ) := by linarith
    calc (1 - (C.seed.spatialDefect : ℝ)) * (Zs 0 C.n₀ - δ) ^ 2
        ≤ (1 - C.seed.spatialDefect) * C.scheme.trace sec x 0 C.n₀ ^ 2 :=
          mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by linarith) (by linarith) 2) hp0
      _ ≤ _ := hp
      _ ≤ Zs 0 (2 * C.n₀) + δ := by linarith
  obtain ⟨hp, hth⟩ := seed_limit hp1 hpos' (lt_min hμ (lt_min (by linarith) hpos')) happrox
  have hq' : (1 - (C.thermalDefect : ℝ)) * (C.seed.thermalBound : ℝ) ^ 2 ≤ 1 := by
    exact_mod_cast hq
  have hq1' : (C.thermalDefect : ℝ) ≤ 1 := by exact_mod_cast hq1
  exact purities_of_shifted hs0 hpos hp hth (by rw [← hsh]; exact htrial)
    (by rw [← shiftedTrace_two_mul]; exact htrial₂) hq' hq1'

end TableCertificate

end CertifiedGaps
