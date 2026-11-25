/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Ising.Seed
import SeedKernel.Seeds

/-!
# The Ising ring at any inverse temperature

The low spectra `lowSpectrum6` and `lowSpectrum12` of the rings of `6` and `12` sites do not
depend on `β`. Through the general seeds of `SeedKernel.Seeds` they give the ring seeds at every
rational `β₀` for which a decidable rational condition holds (`RingCheck`), and so a gap for
every chain certificate that passes it (`ring_gap_of_check`). This file also records the
rational data (Ritz values, radius, thresholds) in the form the general seeds take.
-/

open CertifiedGaps Real

namespace SeedKernel.Ising

/-- The eigenvalue radius `10⁻¹⁰` of the residual checks, as a rational. -/
def radQ : ℚ := 1 / 10 ^ 10

/-- The Ritz value `lam / (20 S)` of a scaled integer candidate. -/
def ritz (lam : ℤ) : ℚ := lam / (20 * scale)

theorem cast_ritz (lam : ℤ) : (ritz lam : ℝ) = lam / (20 * scale) := by simp [ritz]

theorem cast_radQ : (radQ : ℝ) = rad := by norm_num [radQ, rad]

/-- The lower bound on the third eigenvalue of the 6-site ring. -/
def ringRho6 : ℚ := (-74 + 5 * 6 : ℤ) / 20 - 6 / 4

/-- The lower bound on the third eigenvalue of the 12-site ring. -/
def ringRho12 : ℚ := (-230 + 5 * 12 : ℤ) / 20 - 12 / 4

/-- `lowSpectrum6` with rational data. -/
theorem ring6 : LowSpectrum (isingModel.spectrum 6) ((ritz lam6_1 : ℝ) - radQ)
    ((ritz lam6_1 : ℝ) + radQ) ((ritz lam6_2 - radQ : ℚ) : ℝ) (ringRho6 : ℝ) := by
  convert lowSpectrum6 using 2 <;> push_cast [cast_ritz, cast_radQ, ringRho6] <;> norm_num

/-- `lowSpectrum12` with rational data. -/
theorem ring12 : LowSpectrum (isingModel.spectrum 12) ((ritz lam12_1 : ℝ) - radQ)
    ((ritz lam12_1 : ℝ) + radQ) ((ritz lam12_2 - radQ : ℚ) : ℝ) (ringRho12 : ℝ) := by
  convert lowSpectrum12 using 2 <;> push_cast [cast_ritz, cast_radQ, ringRho12] <;> norm_num

/-- The finite conditions of the two ring seeds at `β` with `k` Taylor terms: `P_s(6, β) ≥ 1 - s`
and `P_τ(12, β) ≥ 1 - t`. Decided by `decide +kernel`. -/
abbrev RingCheck (k : ℕ) (β s t : ℚ) : Prop :=
  SpatialCheck k (2 ^ 6) β (ritz lam6_1) (ritz lam12_1) radQ (ritz lam6_2 - radQ) ringRho6 s ∧
    ThermalCheck k (2 ^ 12) β (ritz lam12_1) radQ (ritz lam12_2 - radQ) ringRho12 t

/-- **The ring seeds at any `β`.** -/
theorem ring_seeds_of_check {k : ℕ} {β s t : ℚ} (h : RingCheck k β s t) :
    1 - (s : ℝ) ≤ isingModel.spatialPurity 6 β ∧
      1 - (t : ℝ) ≤ isingModel.thermalPurity 12 β := by
  refine ⟨spatialPurity_ge_of_check isingModel (L := 6) ring6 (near_of_lowSpectrum ring12)
    (by rw [Fintype.card_fin]; exact h.1), thermalPurity_ge_of_check isingModel ring12
    (by rw [Fintype.card_fin]; exact h.2)⟩

/-- **The Ising ring gap at any `β`.** For a valid chain certificate `c` whose seeds pass
`RingCheck` at `β > 0`, every even ring of `L ≥ 6` sites has gap `Δ_L > log r / β`. -/
theorem ring_gap_of_check (c : ChainCertificate) (hc : c.Valid) {k : ℕ} {β : ℚ} (hβ : 0 < β)
    (h : RingCheck k β (c.sAt 0) (c.tAt 0)) {L : ℕ} (hL : 6 ≤ L) (hLe : Even L) :
    GapGT (isingModel.spectrum L) (log (c.r : ℝ) / β) := by
  obtain ⟨hs, ht⟩ := ring_seeds_of_check h
  exact isingModel.gapGT_of_chainCertificate (n₀ := 6) (by norm_num) (by decide)
    (by exact_mod_cast hβ) hc ⟨rfl, rfl, rfl⟩ hs ht hL hLe

end SeedKernel.Ising
