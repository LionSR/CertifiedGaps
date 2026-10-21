/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.IsingOpen.Bridge
import SeedKernel.Comparison
import SeedKernel.IsingOpen.GeneratedChunks
import SeedKernel.Ising.Checks

/-!
# An unconditional gap for the free Ising chain

The longitudinal Ising chain `-∑ ZZ - (1/4) ∑ X - (1/20) ∑ Z` with open boundary conditions
(free ends). The kernel checks of `Generated.lean` and `GeneratedChunks.lean` (integer trial
vectors on the free chains of `2`, `8` and `12` sites, classical counts and comparison weights)
give low spectra of these chains (`lowSpectrum2`, `near8`, `lowSpectrum12`), which do not depend
on `β`. The general seeds of `SeedKernel.Seeds` turn them, with the ring spectra of
`SeedKernel.Ising.Checks`, into a gap theorem at any inverse temperature:

* `free_gap_of_check`: for a valid open-chain certificate with `N₀ = 3`, `K₀ = 4` whose seeds pass
  the decidable rational condition `FreeCheck` at `β > 0`, every free chain of `n ≥ 12 · 2^m`
  sites has gap `Δ_n > log(1/(10x)) / (β 2^m)`;
* `ising_free_gap`: at `β₀ = 97/4`, `Δ_n > log(10¹²/3801) / (97/2) ≈ 0.3998` for `n ≥ 24`;
* `ising_free_gap_long`: `Δ_n > log(10⁴³/1609) / 194 ≈ 0.4723` for `n ≥ 96`.

The bound on the eigenvalues above the lowest two comes from the weighted comparison
`SeedKernel.le_re_quadratic_of_weights` through `freeH_weights`, since the crude bound
`E(σ) - g n` of the ring fails at the ends.
-/

open Matrix CertifiedGaps Real SeedKernel.Ising

set_option maxRecDepth 100000

namespace SeedKernel.IsingOpen

/-- The free chain on `k + 1` sites is Hermitian. -/
theorem freeH_isHermitian (k : ℕ) : (freeH (1 / 4) (1 / 20) k).IsHermitian :=
  isingModel.openHamiltonian_isHermitian isingModel.freeBoundary (k + 1)

/-- The spectrum of the free chain on `k + 1` sites. -/
noncomputable abbrev spec (k : ℕ) : (Fin (k + 1) → Fin 2) → ℝ :=
  isingModel.openSpectrum isingModel.freeBoundary (k + 1)

/-! ### The low spectrum of a free chain from the kernel checks -/

/-- A residual check gives an eigenvalue within `rad` of `lam / (20 S)`. -/
theorem exists_near_of_check {k : ℕ} (v : ℕ) (lam : ℤ) (R : ℕ)
    (hdot : |((Ising.dot (k + 1) v v : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2)
    (hres : residualSq (k + 1) v lam ≤ R) (hR : (R : ℝ) * 10 ^ 18 ≤ (scale : ℝ) ^ 4) :
    ∃ i, |spec k i - lam / (20 * scale)| ≤ rad := by
  have hS : (0 : ℝ) < scale := by norm_num [scale]
  set K := (scale : ℝ) ^ 2 with hK
  have hKpos : 0 < K := by positivity
  have hres' : ((residualSq (k + 1) v lam : ℕ) : ℝ) ≤ R := by exact_mod_cast hres
  rw [IsingOpen.residualSq_eq, show (20 * (scale : ℝ)) ^ 2 = 400 * K by rw [hK]; ring] at hres'
  have h4 : (R : ℝ) ≤ K ^ 2 / 10 ^ 18 := by
    rw [le_div_iff₀ (by norm_num), hK, ← pow_mul]; exact hR
  refine exists_near_of_residual (freeH_isHermitian k) (Ising.vec v) (N := K * (1 - 1 / 10 ^ 12))
    (by norm_num [rad]) (by positivity) ?_ ?_
  · rw [normSqSum_eq_nsq, ← dot_self_eq]
    have := neg_abs_le (((Ising.dot (k + 1) v v : ℤ) : ℝ) - K)
    nlinarith
  · set X := normSqSum (freeH (1 / 4) (1 / 20) k *ᵥ Ising.vec v -
      ((lam / (20 * scale) : ℝ) : ℂ) • Ising.vec v)
    have hX : X ≤ K / (400 * 10 ^ 18) := by
      rw [le_div_iff₀ (by positivity)]
      have : 400 * K * X ≤ K ^ 2 / 10 ^ 18 := hres'.trans h4
      rw [le_div_iff₀ (by positivity)] at this
      nlinarith
    refine hX.trans ?_
    rw [rad, div_le_iff₀ (by positivity)]
    nlinarith

/-- The classical count and the weight checks bound the number of eigenvalues below
`ρ = rho / 20`. -/
theorem card_low_le {k : ℕ} (w : ℕ) (cut rho : ℤ) (hcount : classicalCount (k + 1) cut = 2)
    (hw : ∀ s < 2 ^ (k + 1), weightOK (k + 1) w cut rho s = true) :
    (Finset.univ.filter fun i ↦ spec k i < (rho : ℝ) / 20).card ≤ 2 := by
  classical
  refine le_trans ?_ ((IsingOpen.classicalCount_eq cut).symm.trans hcount).le
  refine card_eigenvalues_lt_le (freeH_isHermitian k) (lowSet k cut) _ fun x hx ↦ ?_
  exact SeedKernel.le_re_quadratic_of_weights (freeH_isHermitian k) (lowSet k cut)
    (fun σ ↦ (entry w (encode σ) : ℝ)) _ (fun σ hσ ↦ freeH_weights (1 / 4) (1 / 20) (by norm_num)
      _ _ _ (weights_spec w cut rho hw) σ hσ) x hx

/-- **The low spectrum of a free chain** from two residual checks, a classical count and the
weight checks. -/
theorem lowSpectrum_of_checks {k : ℕ} (u v w : ℕ) (a b cut rho : ℤ) (R : ℕ)
    (hdu : |((Ising.dot (k + 1) u u : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2)
    (hdv : |((Ising.dot (k + 1) v v : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2)
    (hres : residualSq (k + 1) u a + residualSq (k + 1) v b ≤ R)
    (hR : (R : ℝ) * 10 ^ 18 ≤ (scale : ℝ) ^ 4)
    (hab : (a : ℝ) / (20 * scale) + 2 * rad < b / (20 * scale))
    (hbρ : (b : ℝ) / (20 * scale) + rad < rho / 20)
    (hcount : classicalCount (k + 1) cut = 2)
    (hw : ∀ s < 2 ^ (k + 1), weightOK (k + 1) w cut rho s = true) :
    LowSpectrum (spec k) (a / (20 * scale) - rad) (a / (20 * scale) + rad)
      (b / (20 * scale) - rad) (rho / 20) := by
  classical
  exact lowSpectrum_of_near (exists_near_of_check (k := k) u a R hdu (by omega) hR)
    (exists_near_of_check (k := k) v b R hdv (by omega) hR) hab hbρ
    (card_low_le w cut rho hcount hw)

/-- A passing weight block covers its indices. -/
theorem weightBlock_spec {n w : ℕ} {cut rho : ℤ} {start count : ℕ}
    (h : weightBlock n w cut rho start count = true) :
    ∀ i < count, weightOK n w cut rho (start + i) = true := by
  simpa [weightBlock, List.all_eq_true, List.mem_range] using h

/-! ### The three free chains -/

/-- The free chain of `2` sites. -/
theorem lowSpectrum2 : LowSpectrum (spec 1)
    (lam2_1 / (20 * scale) - rad) (lam2_1 / (20 * scale) + rad)
    (lam2_2 / (20 * scale) - rad) ((rho2 : ℝ) / 20) :=
  lowSpectrum_of_checks (k := 1) u2_1 u2_2 w2 lam2_1 lam2_2 cut2 rho2
    (residualSq 2 u2_1 lam2_1 + residualSq 2 u2_2 lam2_2)
    (dot_real_bound (L := 2) (u := u2_1) dot2_1) (dot_real_bound (L := 2) (u := u2_2) dot2_2)
    le_rfl (res_real_bound res2) (by norm_num [lam2_1, lam2_2, scale, rad])
    (by norm_num [lam2_2, rho2, scale, rad]) count2
    (fun s hs ↦ by simpa using weightBlock_spec weights2 s hs)

/-- The free chain of `8` sites has an eigenvalue within `rad` of `lam8_1 / (20 S)`. -/
theorem near8 : ∃ i, |spec 7 i - lam8_1 / (20 * scale)| ≤ rad :=
  exists_near_of_check (k := 7) u8_1 lam8_1 _ (dot_real_bound (L := 8) (u := u8_1) dot8_1) le_rfl
    (res_real_bound res8)

/-- The block statistics of the 12-site chain are the whole-chain sums. -/
theorem stats12 :
    (Ising.dot 12 u12_1 u12_1, residualSq 12 u12_1 lam12_1) = sumBlocks statBlocks12_1 ∧
    (Ising.dot 12 u12_2 u12_2, residualSq 12 u12_2 lam12_2) = sumBlocks statBlocks12_2 := by
  rw [statBlocks12_1_eq, statBlocks12_2_eq]
  exact ⟨(sumBlocks_eq 12 _ _ 64 64 (by norm_num)).symm,
    (sumBlocks_eq 12 _ _ 64 64 (by norm_num)).symm⟩

theorem checks12 :
    (Ising.dot 12 u12_1 u12_1 - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat ∧
    (Ising.dot 12 u12_2 u12_2 - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat ∧
    (residualSq 12 u12_1 lam12_1 + residualSq 12 u12_2 lam12_2) * 10 ^ 18 ≤
      (scale ^ 4).toNat := by
  obtain ⟨h1, h2⟩ := stats12
  rw [statSum12_1, Prod.mk.injEq] at h1
  rw [statSum12_2, Prod.mk.injEq] at h2
  rw [h1.1, h2.1, h1.2, h2.2]
  decide +kernel

/-- The free chain of `12` sites. -/
theorem lowSpectrum12 : LowSpectrum (spec 11)
    (lam12_1 / (20 * scale) - rad) (lam12_1 / (20 * scale) + rad)
    (lam12_2 / (20 * scale) - rad) ((rho12 : ℝ) / 20) := by
  refine lowSpectrum_of_checks (k := 11) u12_1 u12_2 w12 lam12_1 lam12_2 cut12 rho12
    (residualSq 12 u12_1 lam12_1 + residualSq 12 u12_2 lam12_2) ?_ ?_ ?_ ?_ ?_ ?_ count12 ?_
    <;> try simp only [Nat.reduceAdd]
  · exact dot_real_bound (L := 12) (u := u12_1) checks12.1
  · exact dot_real_bound (L := 12) (u := u12_2) checks12.2.1
  · exact le_rfl
  · exact res_real_bound checks12.2.2
  · norm_num [lam12_1, lam12_2, scale, rad]
  · norm_num [lam12_2, rho12, scale, rad]
  · exact fun s hs ↦ weightBlocks_all 12 w12 cut12 rho12 64 64
      (by rw [← weightBlocks12_eq]; exact weightBlocks12_all) s (hs.trans_eq (by norm_num))

/-! ### The seeds at any inverse temperature -/

/-- `lowSpectrum2` with rational data. -/
theorem open2 : LowSpectrum (isingModel.openSpectrum isingModel.freeBoundary 2)
    ((ritz lam2_1 : ℝ) - radQ) ((ritz lam2_1 : ℝ) + radQ) ((ritz lam2_2 - radQ : ℚ) : ℝ)
    ((rho2 / 20 : ℚ) : ℝ) := by
  convert lowSpectrum2 using 2 <;> push_cast [cast_ritz, cast_radQ] <;> ring

/-- `near8` with rational data. -/
theorem open8 : ∃ i, |isingModel.openSpectrum isingModel.freeBoundary 8 i - (ritz lam8_1 : ℝ)| ≤
    (radQ : ℝ) := by
  simpa only [cast_ritz, cast_radQ] using near8

/-- `lowSpectrum12` with rational data. -/
theorem open12 : LowSpectrum (isingModel.openSpectrum isingModel.freeBoundary 12)
    ((ritz lam12_1 : ℝ) - radQ) ((ritz lam12_1 : ℝ) + radQ) ((ritz lam12_2 - radQ : ℚ) : ℝ)
    ((rho12 / 20 : ℚ) : ℝ) := by
  convert lowSpectrum12 using 2 <;> push_cast [cast_ritz, cast_radQ] <;> ring

/-- All finite conditions of the free chain at `β` with `k` Taylor terms, for a certificate with
`N₀ = 3`, `K₀ = 4` and one ring period in the boundary seed (`r = 1`): the two ring seeds
(`RingCheck`), the boundary trace `G₂ (d(u₀) Z₆) ≤ E₀ G₈` and the open thermal seed
`P_τ(12, β) ≥ 1 - η₀`. A decidable proposition about rationals. -/
abbrev FreeCheck (k : ℕ) (β : ℚ) (c : OpenCertificate) : Prop :=
  RingCheck k β c.u₀ c.u₀ ∧
    TraceCheck k (2 ^ 2) (2 ^ 6) 1 β (ritz lam2_1) (ritz lam2_2 - radQ) (rho2 / 20)
      (ritz Ising.lam6_1) (ritz Ising.lam6_2 - radQ) ringRho6 (ritz lam8_1) radQ
      (OpenCertificate.tailQ c.u₀) (c.EAt 0) ∧
    ThermalCheck k (2 ^ 12) β (ritz lam12_1) radQ (ritz lam12_2 - radQ) (rho12 / 20) (c.ηAt 0)

/-- **The free Ising chain at any `β`.** For a valid open-chain certificate `c` with
`N₀ = 3`, `K₀ = 4` whose seeds pass `FreeCheck` at `β > 0`, every free chain of
`n ≥ 12 · 2^m` sites has a unique ground state and gap `Δ_n > log(1/(10x)) / (β 2^m)`. -/
theorem free_gap_of_check (c : OpenCertificate) (hc : c.Valid) (hN : c.N₀ = 3) (hK : c.K₀ = 4)
    {k : ℕ} {β : ℚ} (hβ : 0 < β) (h : FreeCheck k β c) {n : ℕ} (hn : 12 * 2 ^ c.m ≤ n) :
    GapGT (isingModel.openSpectrum isingModel.freeBoundary n)
      (log (1 / (10 * (c.x : ℝ))) / (β * 2 ^ c.m)) := by
  obtain ⟨hring, htrace, hopen⟩ := h
  obtain ⟨hs, ht⟩ := ring_seeds_of_check hring
  refine isingModel.freeGapGT_of_traces rfl c hc (by exact_mod_cast hβ) ?_ ?_ (r := 1)
    (by rw [hN, hK]) ?_ ?_ (by rw [hN]; omega)
  · rw [hN]; exact hs
  · rw [hN]; exact ht
  · rw [hN, hK]
    exact openTrace_le_of_check isingModel isingModel.freeBoundary (n₁ := 2) (L := 6)
      (n₂ := 8) open2 ring6 open8 (by rw [Fintype.card_fin]; exact htrace)
  · rw [hN]
    exact openThermalPurity_ge_of_check isingModel isingModel.freeBoundary open12
      (by rw [Fintype.card_fin]; exact hopen)

/-! ### The gap at `β₀ = 97/4` -/

/-- The open-chain certificate of the free Ising chain at `β₀ = 97/4` with one base step, found
by `lean/scripts/open_certificate.jl 3 4 15e-7 11e-7 12e-9 97/4 1`. -/
def isingFreeCertificate : OpenCertificate where
  N₀ := 3
  K₀ := 4
  u₀ := 15e-7
  C := 4501e-3
  q := [1718e-7]
  E := [11e-7, 3801e-13]
  η := [12e-9, 1144e-15]
  x := 3801e-13

theorem isingFreeCertificate_valid : isingFreeCertificate.Valid := by decide +kernel

/-- The same seeds with three base steps (`open_certificate.jl ... 97/4 3`): a longer minimal
length and a sharper bound. -/
def isingFreeCertificate3 : OpenCertificate where
  N₀ := 3
  K₀ := 4
  u₀ := 15e-7
  C := 4501e-3
  q := [1718e-7, 152e-13, 6926e-25]
  E := [11e-7, 3801e-13, 1158e-23, 1609e-44]
  η := [12e-9, 1144e-15, 6353e-26, 1379e-46]
  x := 1609e-44

theorem isingFreeCertificate3_valid : isingFreeCertificate3.Valid := by decide +kernel

/-- The finite conditions at `β₀ = 97/4`, with `40` Taylor terms (both certificates share the
seeds). -/
theorem isingFree_check : FreeCheck 40 (97 / 4) isingFreeCertificate := by decide +kernel

theorem isingFree_check3 : FreeCheck 40 (97 / 4) isingFreeCertificate3 := by decide +kernel

/-- **An unconditional gap with open boundary conditions.** The free Ising chain
`-∑_{j<n-1} Zⱼ Zⱼ₊₁ - (1/4) ∑ Xⱼ - (1/20) ∑ Zⱼ`: every chain of `n ≥ 24` sites has a unique
ground state and gap `Δ_n > log(10¹²/3801) / (97/2)` (`≈ 0.3998`). No hypothesis and no
`sorry`. -/
theorem ising_free_gap {n : ℕ} (hn : 24 ≤ n) :
    GapGT (isingModel.openSpectrum isingModel.freeBoundary n)
      (log (10 ^ 12 / 3801) / (97 / 2)) := by
  have h := free_gap_of_check isingFreeCertificate isingFreeCertificate_valid rfl rfl
    (by norm_num) isingFree_check (n := n) (by rw [show isingFreeCertificate.m = 1 from rfl]; omega)
  rw [show isingFreeCertificate.m = 1 from rfl] at h
  convert h using 2 <;> norm_num [isingFreeCertificate]

/-- **The sharper gap for long free chains**: `Δ_n > log(10⁴³/1609) / 194` (`≈ 0.4723`) for
every `n ≥ 96`. -/
theorem ising_free_gap_long {n : ℕ} (hn : 96 ≤ n) :
    GapGT (isingModel.openSpectrum isingModel.freeBoundary n)
      (log (10 ^ 43 / 1609) / 194) := by
  have h := free_gap_of_check isingFreeCertificate3 isingFreeCertificate3_valid rfl rfl
    (by norm_num) isingFree_check3 (n := n)
    (by rw [show isingFreeCertificate3.m = 3 from rfl]; omega)
  rw [show isingFreeCertificate3.m = 3 from rfl] at h
  convert h using 2 <;> norm_num [isingFreeCertificate3]

end SeedKernel.IsingOpen

#print axioms SeedKernel.IsingOpen.free_gap_of_check
#print axioms SeedKernel.IsingOpen.isingFree_check
#print axioms SeedKernel.IsingOpen.ising_free_gap
#print axioms SeedKernel.IsingOpen.ising_free_gap_long
