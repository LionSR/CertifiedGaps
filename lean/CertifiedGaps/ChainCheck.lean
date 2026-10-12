/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Chain

/-!
# Rational chain certificates

A `ChainCertificate` is the finite rational data that `gapGT_of_chain` needs besides the seeds:
the chain `(s j, t j)`, the hand-over defect `ε` and the rate `r`. `ChainCertificate.Valid` is the
list of rational inequalities, which is decidable, so a concrete certificate is checked by
`decide`. `gapGT_of_chainCertificate` turns a valid certificate and the two seed inequalities into
the gap bound `log r / β₀`.

The data are produced by `lean/scripts/chain_certificate.jl` from the seed defects; that script
only searches, and Lean re-checks every inequality.
-/

open Real

namespace ThermalPurityGap

/-- The rational data of a chain certificate: upper bounds `s j, t j` on the defects at scales
`j = 0, …, J` (with `J + 1` the length of `s`), the hand-over defect `ε` and the rate `r`. -/
structure ChainCertificate where
  /-- Upper bounds on the spatial defects `1 - P_s(n₀ 2^j, β₀ 2^j)`. -/
  s : List ℚ
  /-- Upper bounds on the thermal defects `1 - P_τ(2 n₀ 2^j, β₀ 2^j)`. -/
  t : List ℚ
  /-- The rate: the gap bound is `log r / β₀`. -/
  r : ℚ
  /-- The common bound `≤ 1/100` on the defects at the last scale. -/
  ε : ℚ

namespace ChainCertificate

/-- The coupling defect `z(a, b) = 1 - (1 - b)(1 - a)²` over `ℚ`. -/
def zQ (a b : ℚ) : ℚ := 1 - (1 - b) * (1 - a) ^ 2

/-- The squaring bound `f(z) = z² / (2 (1 - z)²)` over `ℚ`. -/
def fQ (z : ℚ) : ℚ := z ^ 2 / (2 * (1 - z) ^ 2)

variable (c : ChainCertificate)

/-- The last scale `J`. -/
def J : ℕ := c.s.length - 1

/-- The spatial defect bound at scale `j`. -/
def sAt (j : ℕ) : ℚ := c.s.getD j 0

/-- The thermal defect bound at scale `j`. -/
def tAt (j : ℕ) : ℚ := c.t.getD j 0

/-- The rational conditions of `gapGT_of_chain`. -/
def Valid : Prop :=
  0 < c.r ∧ 0 < c.ε ∧ c.ε ≤ 1 / 100 ∧ c.sAt c.J ≤ c.ε ∧ c.tAt c.J ≤ c.ε ∧
    c.r ^ 2 ^ c.J * (5 * c.ε) ≤ 1 ∧
    (∀ j < c.J + 1, c.sAt j ≤ 1 ∧ c.tAt j ≤ 1) ∧
    (∀ j < c.J, zQ (c.sAt j) (c.tAt j) < 1 ∧ fQ (zQ (c.sAt j) (c.tAt j)) ≤ c.sAt (j + 1) ∧
      zQ (c.tAt j) (c.sAt (j + 1)) < 1 ∧ fQ (zQ (c.tAt j) (c.sAt (j + 1))) ≤ c.tAt (j + 1) ∧
      c.r ^ 2 ^ j * zQ (c.sAt j) (c.tAt j) < 1 - zQ (c.sAt j) (c.tAt j))

instance : Decidable c.Valid := by unfold Valid; infer_instance

theorem cast_zQ (a b : ℚ) : couplingDefect a b = (zQ a b : ℝ) := by
  simp [couplingDefect, zQ]

theorem cast_fQ (z : ℚ) : squaringDefect z = (fQ z : ℝ) := by
  simp [squaringDefect, fQ]

end ChainCertificate

variable {ι : ℕ → Type*} [∀ L, Fintype (ι L)] [∀ L, Nonempty (ι L)]

/-- **The two-purity criterion from a rational chain certificate.** Under the moment inequalities
above an even `n₀ > 0`, a valid certificate `c` and seeds `P_s(n₀, β₀) ≥ 1 - s 0`,
`P_τ(2n₀, β₀) ≥ 1 - t 0` give `Δ_L > log r / β₀` for every even `L ≥ n₀`. -/
theorem gapGT_of_chainCertificate (E : (L : ℕ) → ι L → ℝ) {n₀ : ℕ} (hn₀ : 0 < n₀)
    (hn₀e : Even n₀) {β₀ : ℝ} (hβ₀ : 0 < β₀) (hZ : HasSpatialMomentBounds E n₀)
    (c : ChainCertificate) (hc : c.Valid)
    (hs : 1 - (c.sAt 0 : ℝ) ≤ spatialPurity E n₀ β₀)
    (ht : 1 - (c.tAt 0 : ℝ) ≤ thermalPurity E (2 * n₀) β₀)
    {L : ℕ} (hL : n₀ ≤ L) (hLe : Even L) : GapGT (E L) (log (c.r : ℝ) / β₀) := by
  obtain ⟨hr, hε, hε₁, hsJ, htJ, htail, hrange, hstep⟩ := hc
  open ChainCertificate in
  refine gapGT_of_chain E hn₀ hn₀e hβ₀ hZ (fun j ↦ c.sAt j) (fun j ↦ c.tAt j) c.J
    (r := c.r) (ε := c.ε) (by exact_mod_cast hr) (fun j hj ↦ ?_) (fun j hj ↦ ?_) (fun j hj ↦ ?_)
    (by exact_mod_cast hε)
    (by rw [show (1 / 100 : ℝ) = ((1 / 100 : ℚ) : ℝ) by norm_num]; exact_mod_cast hε₁)
    (by exact_mod_cast hsJ)
    (by exact_mod_cast htJ) (by exact_mod_cast htail) hs ht hL hLe
  · obtain ⟨h1, h2⟩ := hrange j (by omega)
    exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩
  · obtain ⟨h1, h2, h3, h4, -⟩ := hstep j hj
    simp only [cast_zQ, cast_fQ]
    exact ⟨by exact_mod_cast h1, by exact_mod_cast h2, by exact_mod_cast h3,
      by exact_mod_cast h4⟩
  · obtain ⟨-, -, -, -, h5⟩ := hstep j hj
    simp only [cast_zQ]
    exact_mod_cast h5

end ThermalPurityGap
