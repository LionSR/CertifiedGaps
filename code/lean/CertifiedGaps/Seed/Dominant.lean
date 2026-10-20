/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.CappedMass
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.Rat.Cast.Lemmas
import Mathlib.Data.Rat.Cast.Order
import Mathlib.Data.Rat.BigOperators

/-!
# Seeds from a dominant sector

Let `x i` be real numbers (transfer eigenvalues) split into sectors `sec s`, `s : Fin (m + 1)`,
where sector `s` appears `μ_s` times in the trace, `Zₙ = ∑_s μ_s ∑_{i ∈ sec s} x iⁿ`, and the
dominant sector `0` appears once. Suppose the `q`th moment of each sector is at most `mass s`,
each `|x i|` in sector `0` is below `cap 0`, and `Z_{2n₀} > 1` at `n₀ = q a` (a trial state).
Then one eigenvalue of sector `0` dominates and the seed defects at `n₀` are small (Appendix B of
the certificate note, eqs. (c-capped)–(c-seed)):

  `(1 - p) Z_{n₀}² ≤ Z_{2n₀}`  and  `Z_{2n₀} ≤ cap₀^{2n₀} + R(2a)`.

A nontrivial sector contributes `μ_s · sectorBound s r` to `R(r)`, the bound on its `q r`th moment:
`mass^r` without a cap (superadditivity), and the capped-mass bound `hʳ + (mass - h)ʳ` with
`h = cap^q` when a filter caps it and its mass exceeds `h`. The rational inequalities are
bundled in `SectorSeed.Valid`, which `decide` checks.

## Main results

* `SectorSeed.sum_le_sectorBound`: the moment bound of a nontrivial sector.
* `SectorSeed.exists_dominant`: some `|x i₀| > ℓ` in sector `0`.
* `SectorSeed.seed`: the two inequalities above.
-/

open Finset

namespace CertifiedGaps

/-- The rational data of a dominant-sector seed with `m` nontrivial sectors. -/
structure SectorSeed (m : ℕ) where
  /-- The moment order of the masses (`12` for the spin-one seed). -/
  q : ℕ
  /-- The seed length is `n₀ = q a`. -/
  a : ℕ
  /-- The multiplicity of each sector in the trace (`1` for sector `0`). -/
  mult : Fin (m + 1) → ℕ
  /-- The cap on `|x|` in each sector, if a filter provides one (required for sector `0`). -/
  cap : Fin (m + 1) → Option ℚ
  /-- The bound on the `q`th moment of each sector. -/
  mass : Fin (m + 1) → ℚ
  /-- The dominance threshold `ℓ`. -/
  ell : ℚ
  /-- The spatial defect `p` to certify. -/
  spatialDefect : ℚ

namespace SectorSeed

variable {m : ℕ} (S : SectorSeed m)

/-- The cap of the dominant sector. -/
def cap₀ : ℚ := (S.cap 0).getD 0

/-- The bound on the `q r`th moment of a nontrivial sector `s`. -/
def sectorBound (s : Fin (m + 1)) (r : ℕ) : ℚ :=
  match S.cap s with
  | none => S.mass s ^ r
  | some v => if S.mass s ≤ v ^ S.q then S.mass s ^ r else (v ^ S.q) ^ r + (S.mass s - v ^ S.q) ^ r

/-- The bound `R(r)` on the `q r`th moment of all eigenvalues but the dominant one. -/
def restBound (r : ℕ) : ℚ :=
  ∑ i : Fin m, S.mult i.succ * S.sectorBound i.succ r + (S.mass 0 - S.ell ^ S.q) ^ r

/-- The bound `cap₀^{2n₀} + R(2a)` on `Z_{2n₀}`. -/
def thermalBound : ℚ := (S.cap₀ ^ S.q) ^ (2 * S.a) + S.restBound (2 * S.a)

/-- A capped mass admits one full bar: `mass ≤ 2 cap^q` whenever there is a cap. -/
def OneBar (q : ℕ) : Option ℚ → ℚ → Prop
  | none, _ => True
  | some v, M => M ≤ 2 * v ^ q

instance (q : ℕ) : (o : Option ℚ) → (M : ℚ) → Decidable (OneBar q o M)
  | none, _ => isTrue trivial
  | some v, M => inferInstanceAs (Decidable (M ≤ 2 * v ^ q))

/-- The rational conditions: the shape, the dominant sector (a cap, multiplicity `1`, one full
bar at `ℓ`), the nontrivial masses, dominance at `ℓ` (`Z_{2n₀} < 1` otherwise) and the spatial
defect. -/
def Valid : Prop :=
  Even S.q ∧ 0 < S.q ∧ 1 ≤ S.a ∧ 0 < S.ell ∧ S.spatialDefect ≤ 1 ∧ (S.cap 0).isSome ∧
    S.mult 0 = 1 ∧ S.ell ^ S.q ≤ S.mass 0 ∧ S.mass 0 ≤ 2 * S.ell ^ S.q ∧
    (∀ i : Fin m, 0 ≤ S.mass i.succ ∧ OneBar S.q (S.cap i.succ) (S.mass i.succ)) ∧
    (S.ell ^ S.q) ^ (2 * S.a) + (S.mass 0 - S.ell ^ S.q) ^ (2 * S.a) +
        ∑ i : Fin m, S.mult i.succ * S.sectorBound i.succ (2 * S.a) < 1 ∧
    (1 - S.spatialDefect) * (1 + S.restBound S.a / (S.ell ^ S.q) ^ S.a) ^ 2 ≤ 1

instance : Decidable S.Valid := by unfold Valid; infer_instance

variable {S} {ι : Type*} {sec : Fin (m + 1) → Finset ι} {x : ι → ℝ}

/-- The trace `Zₙ = ∑_s μ_s ∑_{i ∈ sec s} x iⁿ` of sector data. -/
def sectorSum (S : SectorSeed m) (sec : Fin (m + 1) → Finset ι) (x : ι → ℝ) (n : ℕ) : ℝ :=
  ∑ s, (S.mult s : ℝ) * ∑ i ∈ sec s, x i ^ n

/-- A cap on `|y|` caps `y^q` for even `q`. -/
theorem pow_le_of_abs_le {q : ℕ} (hq : Even q) {y v : ℝ} (h : |y| ≤ v) : y ^ q ≤ v ^ q := by
  rw [← hq.pow_abs]
  exact pow_le_pow_left₀ (abs_nonneg _) h _

/-- **The moment bound of a sector**: `∑_{sec s} x^{q r} ≤ sectorBound s r` from the mass and,
if present, the cap. -/
theorem sum_le_sectorBound (hq : Even S.q) {s : Fin (m + 1)} (hmass0 : 0 ≤ S.mass s)
    (hbar : OneBar S.q (S.cap s) (S.mass s))
    (hcap : ∀ v, S.cap s = some v → ∀ i ∈ sec s, |x i| < v)
    (hmass : ∑ i ∈ sec s, x i ^ S.q ≤ S.mass s) {r : ℕ} (hr : 1 ≤ r) :
    ∑ i ∈ sec s, x i ^ (S.q * r) ≤ S.sectorBound s r := by
  simp only [pow_mul]
  have ht0 : ∀ i ∈ sec s, 0 ≤ x i ^ S.q := fun i _ ↦ hq.pow_nonneg _
  have hsum0 : 0 ≤ ∑ i ∈ sec s, x i ^ S.q := sum_nonneg ht0
  -- Without a usable cap: superadditivity and the mass.
  have hplain : ∑ i ∈ sec s, (x i ^ S.q) ^ r ≤ (S.mass s : ℝ) ^ r :=
    (sum_pow_le_pow_sum _ _ ht0 hr).trans (pow_le_pow_left₀ hsum0 hmass r)
  rcases hc : S.cap s with _ | v
  · simp only [sectorBound, hc]; push_cast; exact hplain
  · simp only [sectorBound, hc]
    rw [hc] at hbar
    split_ifs with hle
    · push_cast; exact hplain
    · have h1 : ((v : ℝ) ^ S.q) ≤ S.mass s := by exact_mod_cast (not_le.1 hle).le
      have h2 : (S.mass s : ℝ) ≤ 2 * (v : ℝ) ^ S.q := by exact_mod_cast hbar
      push_cast
      exact sum_pow_le_capped _ _ ht0
        (fun i hi ↦ pow_le_of_abs_le hq (hcap v hc i hi).le) hmass h1 h2 hr

/-- **The dominant eigenvalue.** Some `x i₀` in sector `0` has `|x i₀| > ℓ`, since otherwise the
capped-mass bounds give `Z_{2n₀} < 1`. -/
theorem exists_dominant (hS : S.Valid)
    (hcap : ∀ s v, S.cap s = some v → ∀ i ∈ sec s, |x i| < v)
    (hmass : ∀ s, ∑ i ∈ sec s, x i ^ S.q ≤ S.mass s)
    (h2 : 1 < S.sectorSum sec x (S.q * (2 * S.a))) :
    ∃ i₀ ∈ sec 0, (S.ell : ℝ) < |x i₀| := by
  by_contra hno
  push Not at hno
  have hV := hS
  obtain ⟨hq, -, ha, -, -, -, hm0, hJ1, hJ2, hN, hdom, -⟩ := hV
  have hJ : ∑ i ∈ sec 0, x i ^ (S.q * (2 * S.a)) ≤
      ((S.ell : ℝ) ^ S.q) ^ (2 * S.a) + (S.mass 0 - (S.ell : ℝ) ^ S.q) ^ (2 * S.a) := by
    rw [sum_congr rfl fun i _ ↦ pow_mul (x i) S.q (2 * S.a)]
    exact sum_pow_le_capped _ (fun i ↦ x i ^ S.q) (fun i _ ↦ hq.pow_nonneg _)
      (fun i hi ↦ pow_le_of_abs_le hq (hno i hi)) (hmass 0) (by exact_mod_cast hJ1)
      (by exact_mod_cast hJ2) (by omega)
  have hrest : ∑ i : Fin m, (S.mult i.succ : ℝ) * ∑ j ∈ sec i.succ, x j ^ (S.q * (2 * S.a)) ≤
      ∑ i : Fin m, (S.mult i.succ : ℝ) * S.sectorBound i.succ (2 * S.a) :=
    sum_le_sum fun i _ ↦ mul_le_mul_of_nonneg_left
      (sum_le_sectorBound hq (by exact_mod_cast (hN i).1) (hN i).2 (hcap _) (hmass _)
        (by omega)) (Nat.cast_nonneg _)
  have hdom' := (Rat.cast_lt (K := ℝ)).2 hdom
  simp only [Rat.cast_add, Rat.cast_pow, Rat.cast_sub, Rat.cast_sum, Rat.cast_mul,
    Rat.cast_natCast, Rat.cast_one] at hdom'
  simp only [sectorSum, Fin.sum_univ_succ, hm0, Nat.cast_one, one_mul] at h2
  linarith

/-- **The seed from a dominant sector.** Under `S.Valid`, the caps, the masses and
`Z_{2n₀} > 1`, the traces at `n₀ = q a` satisfy `(1 - p) Z_{n₀}² ≤ Z_{2n₀}` and
`Z_{2n₀} ≤ S.thermalBound`. -/
theorem seed (hS : S.Valid) (hcap : ∀ s v, S.cap s = some v → ∀ i ∈ sec s, |x i| < v)
    (hmass : ∀ s, ∑ i ∈ sec s, x i ^ S.q ≤ S.mass s)
    (h2 : 1 < S.sectorSum sec x (S.q * (2 * S.a))) :
    (1 - S.spatialDefect : ℝ) * S.sectorSum sec x (S.q * S.a) ^ 2 ≤
        S.sectorSum sec x (S.q * (2 * S.a)) ∧
      S.sectorSum sec x (S.q * (2 * S.a)) ≤ S.thermalBound := by
  classical
  obtain ⟨i₀, hi₀, hdom⟩ := exists_dominant hS hcap hmass h2
  have hV := hS
  obtain ⟨hq, hq0, ha, hell, hp1, hcap0, hm0, -, -, hN, -, hsp⟩ := hV
  obtain ⟨v₀, hv₀⟩ := Option.isSome_iff_exists.1 hcap0
  have hcapv : S.cap₀ = v₀ := by simp [cap₀, hv₀]
  set t₀ := x i₀ ^ S.q with ht₀
  have ht₀ell : (S.ell : ℝ) ^ S.q < t₀ := by
    rw [ht₀, ← hq.pow_abs (x i₀)]
    exact pow_lt_pow_left₀ hdom (by exact_mod_cast hell.le) hq0.ne'
  have ht₀cap : t₀ ≤ (S.cap₀ : ℝ) ^ S.q := by
    rw [hcapv]; exact pow_le_of_abs_le hq (hcap 0 v₀ hv₀ i₀ hi₀).le
  have hell0 : (0 : ℝ) ≤ (S.ell : ℝ) ^ S.q := pow_nonneg (by exact_mod_cast hell.le) _
  -- The rest of the spectrum: sector `0` without `i₀`, and the nontrivial sectors.
  set rest : ℕ → ℝ := fun n ↦ ∑ i ∈ (sec 0).erase i₀, x i ^ n +
    ∑ i : Fin m, (S.mult i.succ : ℝ) * ∑ j ∈ sec i.succ, x j ^ n
  have hsplit (n : ℕ) : S.sectorSum sec x n = x i₀ ^ n + rest n := by
    simp only [sectorSum, Fin.sum_univ_succ, hm0, Nat.cast_one, one_mul, rest,
      ← add_sum_erase (sec 0) (fun i ↦ x i ^ n) hi₀]
    ring
  have hrestJ {r : ℕ} (hr : 1 ≤ r) :
      ∑ i ∈ (sec 0).erase i₀, x i ^ (S.q * r) ≤ ((S.mass 0 : ℝ) - (S.ell : ℝ) ^ S.q) ^ r := by
    simp only [pow_mul]
    have hsp0 : x i₀ ^ S.q + ∑ i ∈ (sec 0).erase i₀, x i ^ S.q = ∑ i ∈ sec 0, x i ^ S.q :=
      add_sum_erase (sec 0) (fun i ↦ x i ^ S.q) hi₀
    have hrest0 : 0 ≤ ∑ i ∈ (sec 0).erase i₀, x i ^ S.q := sum_nonneg fun i _ ↦ hq.pow_nonneg _
    refine (sum_pow_le_pow_sum _ (fun i ↦ x i ^ S.q) (fun i _ ↦ hq.pow_nonneg _) hr).trans
      (pow_le_pow_left₀ hrest0 ?_ r)
    linarith [hmass 0]
  have hrest {r : ℕ} (hr : 1 ≤ r) : rest (S.q * r) ≤ S.restBound r := by
    have h1 := hrestJ hr
    have h2 : ∑ i : Fin m, (S.mult i.succ : ℝ) * ∑ j ∈ sec i.succ, x j ^ (S.q * r) ≤
        ∑ i : Fin m, (S.mult i.succ : ℝ) * S.sectorBound i.succ r :=
      sum_le_sum fun i _ ↦ mul_le_mul_of_nonneg_left
        (sum_le_sectorBound hq (by exact_mod_cast (hN i).1) (hN i).2 (hcap _) (hmass _) hr)
        (Nat.cast_nonneg _)
    simp only [rest, restBound, Rat.cast_add, Rat.cast_pow, Rat.cast_sub, Rat.cast_sum,
      Rat.cast_mul, Rat.cast_natCast]
    linarith
  have hrest_nonneg {r : ℕ} : 0 ≤ rest (S.q * r) := by
    have he : Even (S.q * r) := hq.mul_right r
    have h1 : 0 ≤ ∑ i ∈ (sec 0).erase i₀, x i ^ (S.q * r) := sum_nonneg fun i _ ↦ he.pow_nonneg _
    have h2 : 0 ≤ ∑ i : Fin m, (S.mult i.succ : ℝ) * ∑ j ∈ sec i.succ, x j ^ (S.q * r) :=
      sum_nonneg fun i _ ↦ mul_nonneg (Nat.cast_nonneg _) (sum_nonneg fun j _ ↦ he.pow_nonneg _)
    simp only [rest]
    linarith
  have hxr (r : ℕ) : x i₀ ^ (S.q * r) = t₀ ^ r := by rw [ht₀, pow_mul]
  have hRa := hrest (r := S.a) ha
  have hR2a := hrest (r := 2 * S.a) (by omega)
  refine ⟨?_, ?_⟩
  · -- `Z_{n₀} ≤ t₀ᵃ (1 + c)` with `c = R(a)/ℓ^{q a}`, and `(1 - p)(1 + c)² ≤ 1`.
    have hsp' : (1 - (S.spatialDefect : ℝ)) *
        (1 + (S.restBound S.a : ℝ) / ((S.ell : ℝ) ^ S.q) ^ S.a) ^ 2 ≤ 1 := by exact_mod_cast hsp
    set c : ℝ := (S.restBound S.a : ℝ) / ((S.ell : ℝ) ^ S.q) ^ S.a
    have hella : 0 < ((S.ell : ℝ) ^ S.q) ^ S.a := by
      have : (0 : ℝ) < S.ell := by exact_mod_cast hell
      positivity
    have hR : 0 ≤ (S.restBound S.a : ℝ) := hrest_nonneg.trans hRa
    have hc : (S.restBound S.a : ℝ) ≤ t₀ ^ S.a * c := by
      rw [mul_div_assoc', le_div_iff₀ hella, mul_comm (t₀ ^ S.a)]
      exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hell0 ht₀ell.le _) hR
    have ht₀0 : 0 ≤ t₀ := hell0.trans ht₀ell.le
    have hZ : S.sectorSum sec x (S.q * S.a) ≤ t₀ ^ S.a * (1 + c) := by
      rw [hsplit, hxr]; nlinarith
    have hZ' : 0 ≤ S.sectorSum sec x (S.q * S.a) := by
      rw [hsplit, hxr]; have := hrest_nonneg (r := S.a); positivity
    have hp : 0 ≤ 1 - (S.spatialDefect : ℝ) := by
      have : (S.spatialDefect : ℝ) ≤ 1 := by exact_mod_cast hp1
      linarith
    calc (1 - (S.spatialDefect : ℝ)) * S.sectorSum sec x (S.q * S.a) ^ 2
        ≤ (1 - S.spatialDefect) * (t₀ ^ S.a * (1 + c)) ^ 2 := by gcongr
      _ = t₀ ^ (2 * S.a) * ((1 - S.spatialDefect) * (1 + c) ^ 2) := by ring
      _ ≤ t₀ ^ (2 * S.a) * 1 := by gcongr
      _ ≤ S.sectorSum sec x (S.q * (2 * S.a)) := by
          rw [hsplit, hxr]; linarith [hrest_nonneg (r := 2 * S.a)]
  · rw [hsplit, hxr, thermalBound]
    have : t₀ ^ (2 * S.a) ≤ ((S.cap₀ : ℝ) ^ S.q) ^ (2 * S.a) :=
      pow_le_pow_left₀ (hell0.trans ht₀ell.le) ht₀cap _
    push_cast
    linarith

end SectorSeed

end CertifiedGaps
