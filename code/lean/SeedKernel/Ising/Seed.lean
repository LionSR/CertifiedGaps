/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Ising.Bridge
import SeedKernel.Ising.Generated
import SeedKernel.Ising.ChunkCoverage
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# The unconditional seed of the longitudinal Ising chain at `g = 1/4`, `h = 1/20`

From the kernel checks of `Generated.lean` and `ChunkCoverage.lean` (two integer trial vectors
on the rings of `6` and `12` sites with small residuals, and the count of classical
configurations below a cutoff) this file proves, with no hypotheses,

* `seed_spatial`: `1 - 9064·10⁻⁹ ≤ P_s(6, 97/4)`,
* `seed_thermal`: `1 - 3141·10⁻⁸ ≤ P_τ(12, 97/4)`,
* `first_point_seed`: the first of the thirty conjunctions of `IsingLongitudinal.seeds`.

The spectral input on each ring is `LowSpectrum` (one eigenvalue in `[θ₀ - r, θ₀ + r]`, another
at least `θ₁ - r`, all others at least `ρ`), from the residual bound
`CertifiedGaps.Spectral.exists_eigenvalue_near` and counting by comparison with the classical
energies (`CertifiedGaps.Spectral.eigCount_le_card_of_vanishing` with `re_quadratic_ge`). The
partition functions are then bounded by `Z ≥ e^{-β(θ₀+r)}` and
`Z ≤ e^{-β(θ₀-r)} (1 + e^{-β(θ₁-θ₀)} + 2^L e^{-β(ρ-θ₀+r)})`.
-/

open Matrix CertifiedGaps Real

namespace SeedKernel.Ising

/-! ### Partition functions from a low spectrum -/

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

/-! ### The low spectrum of a ring from the kernel checks -/

/-- The model of the first grid point. -/
noncomputable abbrev isingModel : RPRing (Fin 2) (Fin 1) := IsingLongitudinal.model (1 / 4) (1 / 20)

/-- The eigenvalue radius certified by the residual checks. -/
noncomputable def rad : ℝ := 1 / 10 ^ 10

theorem natAbs_le_real {x : ℤ} {k c : ℕ} (h : x.natAbs * k ≤ c) : |(x : ℝ)| * k ≤ c := by
  have : ((x.natAbs * k : ℕ) : ℝ) ≤ c := by exact_mod_cast h
  rwa [Nat.cast_mul, Nat.cast_natAbs, Int.cast_abs] at this

/-- A residual check gives an eigenvalue within `rad` of `lam / (20 S)`. -/
theorem exists_near_of_check {L : ℕ} (v : ℕ) (lam : ℤ) (R : ℕ)
    (hdot : |((dot L v v : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2)
    (hres : residualSq L v lam ≤ R) (hR : (R : ℝ) * 10 ^ 18 ≤ (scale : ℝ) ^ 4) :
    ∃ i, |(isingModel.spectrum L i) - lam / (20 * scale)| ≤ rad := by
  have hS : (0 : ℝ) < scale := by norm_num [scale]
  have : Nonempty (Fin L → Fin 2) := ⟨fun _ ↦ 0⟩
  obtain ⟨i, hi⟩ := Spectral.exists_eigenvalue_near (isingModel.hamiltonian_isHermitian L)
    (vec v) (lam / (20 * scale))
  -- The generic coercion `ℝ → 𝕜` is `Complex.ofReal` at `𝕜 = ℂ` only up to unfolding; restate
  -- `hi` in the form of `residualSq_eq` so that `linarith` sees the same atom.
  replace hi : (isingModel.spectrum L i - lam / (20 * scale)) ^ 2 * Spectral.nsq (vec (L := L) v)
      ≤ Spectral.nsq (isingModel.hamiltonian L *ᵥ vec v -
        ((lam / (20 * scale) : ℝ) : ℂ) • vec v) := hi
  refine ⟨i, ?_⟩
  rw [← dot_self_eq] at hi
  have hres' : ((residualSq L v lam : ℕ) : ℝ) ≤ R := by exact_mod_cast hres
  rw [residualSq_eq] at hres'
  set d := isingModel.spectrum L i - lam / (20 * scale)
  set N := ((dot L v v : ℤ) : ℝ)
  set K := (scale : ℝ) ^ 2 with hK
  have hKpos : 0 < K := by positivity
  have hN : K * (1 - 1 / 10 ^ 12) ≤ N := by
    have := neg_abs_le (N - K)
    nlinarith
  have hd : d ^ 2 * N * (400 * K) ≤ R := by
    have h1 := mul_le_mul_of_nonneg_right hi (sq_nonneg (20 * (scale : ℝ)))
    have h2 : (20 * (scale : ℝ)) ^ 2 = 400 * K := by rw [hK]; ring
    rw [h2] at h1 hres'
    exact h1.trans (by linarith)
  have hd' : d ^ 2 * (1 - 1 / 10 ^ 12) * 400 * K ^ 2 ≤ K ^ 2 / 10 ^ 18 := by
    have h3 : d ^ 2 * (K * (1 - 1 / 10 ^ 12)) * (400 * K) ≤ d ^ 2 * N * (400 * K) := by
      gcongr
    have h4 : (R : ℝ) ≤ K ^ 2 / 10 ^ 18 := by
      rw [le_div_iff₀ (by norm_num)]; rw [hK]; linarith [hR]
    nlinarith
  have hd2 : d ^ 2 ≤ rad ^ 2 := by
    have hK2 : 0 < K ^ 2 := by positivity
    have h5 : d ^ 2 * (1 - 1 / 10 ^ 12) * 400 ≤ 1 / 10 ^ 18 := by
      have := hd'
      rw [show K ^ 2 / 10 ^ 18 = 1 / 10 ^ 18 * K ^ 2 by ring] at this
      nlinarith
    rw [rad]
    nlinarith [sq_nonneg d]
  exact abs_le_of_sq_le_sq' hd2 (by norm_num [rad]) |> fun h ↦ abs_le.2 h

/-- The classical count bounds the number of eigenvalues below `ρ = cut / 20 - L / 4`: the form
of `H` is at least that of the classical energies shifted by `-L/4` (`re_quadratic_ge`), so
counting by vanishing coordinates applies to the configurations below the cut. -/
theorem card_low_le {L : ℕ} (cut : ℤ) (hcount : classicalCount L cut = 2) :
    (Finset.univ.filter fun i ↦ isingModel.spectrum L i < cut / 20 - L / 4).card ≤ 2 := by
  set B := Finset.univ.filter fun σ : Fin L → Fin 2 ↦ diagonal L (encode σ) < cut
  have hB : classicalCount L cut = B.card := by
    rw [classicalCount, show (fun a s ↦ if diagonal L s < cut then a + 1 else a) =
        (fun a s ↦ a + if diagonal L s < cut then 1 else 0) from
      funext₂ fun a s ↦ by split_ifs <;> rfl, foldl_add_eq_sum, zero_add, ← sum_encode,
      Finset.card_filter]
  refine le_trans ?_ (hB.symm.trans hcount).le
  refine Spectral.eigCount_le_card_of_vanishing (isingModel.hamiltonian_isHermitian L) B
    fun x hx ↦ ?_
  refine le_trans ?_ (re_quadratic_ge (1 / 4) (1 / 20) (by norm_num) x)
  rw [Spectral.nsq_eq_sum, Finset.mul_sum]
  refine Finset.sum_le_sum fun σ _ ↦ ?_
  by_cases hσ : σ ∈ B
  · simp [hx σ hσ]
  · refine mul_le_mul_of_nonneg_right ?_ (sq_nonneg _)
    have hd : ¬ (diagonal L (encode σ) < cut) := by simpa [B] using hσ
    have h20 := diagonal_encode σ
    have : (cut : ℝ) ≤ diagonal L (encode σ) := by exact_mod_cast not_lt.1 hd
    linarith

/-- **The low spectrum of a ring** from two residual checks and a classical count. -/
theorem lowSpectrum_of_checks {L : ℕ} (u v : ℕ) (a b cut : ℤ) (R : ℕ)
    (hdu : |((dot L u u : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2)
    (hdv : |((dot L v v : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2)
    (hres : residualSq L u a + residualSq L v b ≤ R) (hR : (R : ℝ) * 10 ^ 18 ≤ (scale : ℝ) ^ 4)
    (hab : (a : ℝ) / (20 * scale) + 2 * rad < b / (20 * scale) - 0)
    (hbρ : (b : ℝ) / (20 * scale) + rad < cut / 20 - L / 4)
    (hcount : classicalCount L cut = 2) :
    LowSpectrum (isingModel.spectrum L) (a / (20 * scale) - rad) (a / (20 * scale) + rad)
      (b / (20 * scale) - rad) (cut / 20 - L / 4) := by
  classical
  obtain ⟨i₀, h₀⟩ := exists_near_of_check (L := L) u a R hdu (by omega) hR
  obtain ⟨i₁, h₁⟩ := exists_near_of_check (L := L) v b R hdv (by omega) hR
  have h₀' := abs_le.1 h₀
  have h₁' := abs_le.1 h₁
  have hne : i₀ ≠ i₁ := by rintro rfl; linarith
  refine ⟨i₀, i₁, hne, by linarith, by linarith, by linarith, fun i hi₀ hi₁ ↦ ?_⟩
  by_contra hlt
  push Not at hlt
  have hsub : ({i₀, i₁, i} : Finset _) ⊆
      Finset.univ.filter fun j ↦ isingModel.spectrum L j < cut / 20 - L / 4 := by
    intro j hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    rcases hj with rfl | rfl | rfl <;> linarith
  have h3 : ({i₀, i₁, i} : Finset _).card = 3 := by
    rw [Finset.card_insert_of_notMem (by simp [hne, Ne.symm hi₀]),
      Finset.card_pair (Ne.symm hi₁)]
  have := (Finset.card_le_card hsub).trans (card_low_le cut hcount)
  omega

/-! ### The two rings -/

theorem chunkStats_eq_sum (n u v : ℕ) (a b : ℤ) (start count : ℕ) :
    chunkStats n u v a b start count = ∑ i ∈ Finset.range count, ((entry u (start + i) *
      entry u (start + i), entry v (start + i) * entry v (start + i),
      entry u (start + i) * entry v (start + i),
      (scale * applyH n u (start + i) - a * entry u (start + i)).natAbs ^ 2,
      (scale * applyH n v (start + i) - b * entry v (start + i)).natAbs ^ 2) : Stats) := by
  have hadd : ∀ acc x : Stats, addStats acc x = acc + x := fun _ _ ↦ rfl
  unfold chunkStats
  simp only [hadd]
  rw [show ((0, 0, 0, 0, 0) : Stats) = 0 from rfl, foldl_add_eq_sum, zero_add]

theorem sumStats_map_range (g : ℕ → Stats) (N : ℕ) :
    sumStats ((List.range N).map g) = ∑ b ∈ Finset.range N, g b := by
  have hadd : ∀ acc x : Stats, addStats acc x = acc + x := fun _ _ ↦ rfl
  rw [sumStats, List.foldl_map]
  simp only [hadd]
  rw [show ((0, 0, 0, 0, 0) : Stats) = 0 from rfl, foldl_add_eq_sum, zero_add]

theorem sum_blocks {M : Type*} [AddCommMonoid M] (F : ℕ → M) (c k : ℕ) :
    ∑ b ∈ Finset.range k, ∑ i ∈ Finset.range c, F (c * b + i) =
      ∑ s ∈ Finset.range (c * k), F s := by
  induction k with
  | zero => simp
  | succ k ih => rw [Finset.sum_range_succ, ih, Nat.mul_succ, Finset.sum_range_add]

/-- Chunked checks are the whole-ring sums: `k` blocks of `c` consecutive basis states, `c k = 2ⁿ`.
Proved for symbolic vectors so that the kernel never unfolds the concrete witnesses. -/
theorem sumStats_chunks (n u v : ℕ) (a b : ℤ) (c k : ℕ) (hck : c * k = 2 ^ n) :
    sumStats ((List.range k).map fun i ↦ chunkStats n u v a b (c * i) c) =
      (dot n u u, dot n v v, dot n u v, residualSq n u a, residualSq n v b) := by
  rw [sumStats_map_range]
  simp only [chunkStats_eq_sum]
  rw [sum_blocks (fun s ↦ ((entry u s * entry u s, entry v s * entry v s, entry u s * entry v s,
      (scale * applyH n u s - a * entry u s).natAbs ^ 2,
      (scale * applyH n v s - b * entry v s).natAbs ^ 2) : Stats)), hck]
  simp only [dot, residualSq, foldl_add_eq_sum, zero_add]
  ext <;> simp [Prod.fst_sum, Prod.snd_sum]

/-- The chunked checks of the 12-site ring are the whole-ring sums. -/
theorem stats12 : sumStats ((List.range 64).map fun b ↦
      chunkStats 12 compact12_1 compact12_2 lam12_1 lam12_2 (64 * b) 64) =
    (dot 12 compact12_1 compact12_1, dot 12 compact12_2 compact12_2,
      dot 12 compact12_1 compact12_2, residualSq 12 compact12_1 lam12_1,
      residualSq 12 compact12_2 lam12_2) :=
  sumStats_chunks 12 _ _ _ _ 64 64 (by norm_num)

theorem s2_toNat : ((scale ^ 2).toNat : ℝ) = (scale : ℝ) ^ 2 := by
  rw [show (scale ^ 2).toNat = 2 ^ 104 by decide]; norm_num [scale]

theorem s4_toNat : ((scale ^ 4).toNat : ℝ) = (scale : ℝ) ^ 4 := by
  rw [show (scale ^ 4).toNat = 2 ^ 208 by decide]; norm_num [scale]

theorem dot_real_bound {L u : ℕ}
    (h : (dot L u u - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat) :
    |((dot L u u : ℤ) : ℝ) - (scale : ℝ) ^ 2| * 10 ^ 12 ≤ (scale : ℝ) ^ 2 := by
  have := natAbs_le_real h
  rw [s2_toNat] at this
  push_cast at this ⊢
  linarith

theorem res_real_bound {R : ℕ} (h : R * 10 ^ 18 ≤ (scale ^ 4).toNat) :
    (R : ℝ) * 10 ^ 18 ≤ (scale : ℝ) ^ 4 := by
  have : ((R * 10 ^ 18 : ℕ) : ℝ) ≤ ((scale ^ 4).toNat : ℕ) := by exact_mod_cast h
  rw [s4_toNat] at this
  push_cast at this ⊢
  linarith

theorem dot6_1 : (dot 6 compact6_1 compact6_1 - scale ^ 2).natAbs * 10 ^ 12 ≤
    (scale ^ 2).toNat := by decide +kernel
theorem dot6_2 : (dot 6 compact6_2 compact6_2 - scale ^ 2).natAbs * 10 ^ 12 ≤
    (scale ^ 2).toNat := by decide +kernel
theorem res6 : (residualSq 6 compact6_1 lam6_1 + residualSq 6 compact6_2 lam6_2) * 10 ^ 18 ≤
    (scale ^ 4).toNat := by decide +kernel

/-- The six-site ring of the first grid point. -/
theorem lowSpectrum6 : LowSpectrum (isingModel.spectrum 6)
    (lam6_1 / (20 * scale) - rad) (lam6_1 / (20 * scale) + rad)
    (lam6_2 / (20 * scale) - rad) ((-74 + 5 * 6 : ℤ) / 20 - (6 : ℕ) / 4) :=
  lowSpectrum_of_checks (L := 6) compact6_1 compact6_2 lam6_1 lam6_2 (-74 + 5 * 6)
    (residualSq 6 compact6_1 lam6_1 + residualSq 6 compact6_2 lam6_2)
    (dot_real_bound (L := 6) (u := compact6_1) dot6_1)
    (dot_real_bound (L := 6) (u := compact6_2) dot6_2) le_rfl (res_real_bound res6)
    (by norm_num [lam6_1, lam6_2, scale, rad]) (by norm_num [lam6_2, scale, rad]) anchor6

/-- What `goodStats` asserts, for any statistics (no concrete sums are unfolded). -/
theorem goodStats_spec {a b c : ℤ} {d e : ℕ} (h : goodStats (a, b, c, d, e) = true) :
    (a - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat ∧
    (b - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat ∧
    (d + e) * 10 ^ 18 ≤ (scale ^ 4).toNat := by
  simp only [goodStats, Bool.and_eq_true, decide_eq_true_eq] at h
  exact ⟨h.1.1.1, h.1.1.2, h.2⟩

/-- The checked statistics of the 12-site ring, unpacked. -/
theorem checks12 :
    (dot 12 compact12_1 compact12_1 - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat ∧
    (dot 12 compact12_2 compact12_2 - scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat ∧
    (residualSq 12 compact12_1 lam12_1 + residualSq 12 compact12_2 lam12_2) * 10 ^ 18 ≤
      (scale ^ 4).toNat :=
  goodStats_spec (stats12 ▸ direct12_canonical)

/-- The twelve-site ring of the first grid point. -/
theorem lowSpectrum12 : LowSpectrum (isingModel.spectrum 12)
    (lam12_1 / (20 * scale) - rad) (lam12_1 / (20 * scale) + rad)
    (lam12_2 / (20 * scale) - rad) ((-230 + 5 * 12 : ℤ) / 20 - (12 : ℕ) / 4) :=
  lowSpectrum_of_checks (L := 12) compact12_1 compact12_2 lam12_1 lam12_2 (-230 + 5 * 12)
    (residualSq 12 compact12_1 lam12_1 + residualSq 12 compact12_2 lam12_2)
    (dot_real_bound (L := 12) (u := compact12_1) checks12.1)
    (dot_real_bound (L := 12) (u := compact12_2) checks12.2.1) le_rfl
    (res_real_bound checks12.2.2)
    (by norm_num [lam12_1, lam12_2, scale, rad]) (by norm_num [lam12_2, scale, rad]) anchor12

/-! ### The seed inequalities -/

theorem exp_neg_fourteen_le : exp (-14) ≤ 1 / 2.7182818283 ^ 14 := by
  rw [exp_neg, ← one_div, show (14 : ℝ) = (14 : ℕ) by norm_num, ← exp_one_pow]
  gcongr
  exact exp_one_gt_d9.le

/-- `Z ≤ A (1 + T)` with `A = e^{-β a₀}` and `T` bounded by powers of `e^{-14}`. -/
theorem partition_upper {L : ℕ} {a₀ b₀ a₁ ρ β : ℝ} (hβ : 0 ≤ β)
    (h : LowSpectrum (isingModel.spectrum L) a₀ b₀ a₁ ρ) :
    isingModel.partitionFunction L β ≤ exp (-(β * a₀)) *
      (1 + (exp (-(β * (a₁ - a₀))) + 2 ^ L * exp (-(β * (ρ - a₀))))) := by
  classical
  rw [RPRing.partitionFunction_eq_sum]
  refine (sum_exp_le hβ h).trans (le_of_eq ?_)
  simp only [Fintype.card_fun, Fintype.card_fin, Nat.cast_pow, Nat.cast_ofNat, mul_add,
    mul_one, ← exp_add, mul_left_comm (exp (-(β * a₀))) (2 ^ L)]
  ring_nf

theorem partition_lower {L : ℕ} {a₀ b₀ a₁ ρ β : ℝ} (hβ : 0 ≤ β)
    (h : LowSpectrum (isingModel.spectrum L) a₀ b₀ a₁ ρ) :
    exp (-(β * b₀)) ≤ isingModel.partitionFunction L β := by
  rw [RPRing.partitionFunction_eq_sum]
  exact sum_exp_ge hβ h

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

/-- **The spatial seed of the first grid point**, unconditionally. -/
theorem seed_spatial :
    1 - (9064e-9 : ℝ) ≤ (IsingLongitudinal.model (1/4) (1/20)).spatialPurity 6 (97/4) := by
  set β : ℝ := 97 / 4
  have hβ : (0 : ℝ) ≤ β := by norm_num [β]
  have hu := partition_upper hβ lowSpectrum6
  have hl := partition_lower hβ lowSpectrum12
  set θ := (lam6_1 : ℝ) / (20 * scale)
  set T := exp (-(β * ((lam6_2 : ℝ) / (20 * scale) - rad - (θ - rad)))) +
    2 ^ 6 * exp (-(β * (((-74 + 5 * 6 : ℤ) : ℝ) / 20 - ((6 : ℕ) : ℝ) / 4 - (θ - rad))))
  set q : ℝ := 1 / 2.7182818283 ^ 14
  have hq : exp (-14) ≤ q := exp_neg_fourteen_le
  have hT : T ≤ q + 2 ^ 6 * q ^ 4 := by
    have e1 : exp (-(β * ((lam6_2 : ℝ) / (20 * scale) - rad - (θ - rad)))) ≤ exp (-14) :=
      exp_le_exp.2 (by norm_num [β, θ, lam6_1, lam6_2, scale, rad])
    have e2 : exp (-(β * (((-74 + 5 * 6 : ℤ) : ℝ) / 20 - ((6 : ℕ) : ℝ) / 4 - (θ - rad)))) ≤
        exp (-14) ^ 4 := by
      rw [← exp_nat_mul]
      exact exp_le_exp.2 (by norm_num [β, θ, lam6_1, scale, rad])
    have hq0 : 0 ≤ exp (-14) := (exp_pos _).le
    have e2' : exp (-14) ^ 4 ≤ q ^ 4 := pow_le_pow_left₀ hq0 hq 4
    have e3 := mul_le_mul_of_nonneg_left (e2.trans e2') (by norm_num : (0 : ℝ) ≤ 2 ^ 6)
    simp only [T]
    exact add_le_add (e1.trans hq) e3
  have hT0 : 0 ≤ T := by positivity
  have hfinal := ratio_ge (Z₁ := isingModel.partitionFunction 6 β)
    (Z₂ := isingModel.partitionFunction 12 β) (exp_pos (-(β * (θ - rad))))
    (by rw [RPRing.partitionFunction_eq_sum]; positivity) hu ?_ hT0 hT
    (x := β * ((lam12_1 : ℝ) / (20 * scale) + rad - 2 * (θ - rad)))
  · refine le_trans ?_ hfinal
    norm_num [q, θ, β, lam6_1, lam12_1, scale, rad]
  · refine le_trans (le_of_eq ?_) hl
    rw [← exp_nat_mul, ← exp_add]
    · congr 1; push_cast; ring

/-- **The thermal seed of the first grid point**, unconditionally. -/
theorem seed_thermal :
    1 - (3141e-8 : ℝ) ≤ (IsingLongitudinal.model (1/4) (1/20)).thermalPurity 12 (97/4) := by
  set β : ℝ := 97 / 4
  have hβ : (0 : ℝ) ≤ β := by norm_num [β]
  have hu := partition_upper hβ lowSpectrum12
  have hl := partition_lower (β := 2 * β) (by positivity) lowSpectrum12
  set θ := (lam12_1 : ℝ) / (20 * scale)
  set T := exp (-(β * ((lam12_2 : ℝ) / (20 * scale) - rad - (θ - rad)))) +
    2 ^ 12 * exp (-(β * (((-230 + 5 * 12 : ℤ) : ℝ) / 20 - ((12 : ℕ) : ℝ) / 4 - (θ - rad))))
  set q : ℝ := 1 / 2.7182818283 ^ 14
  have hq : exp (-14) ≤ q := exp_neg_fourteen_le
  have hq0 : 0 ≤ exp (-14) := (exp_pos _).le
  have hq2 : exp (-14) ^ 2 ≤ q ^ 2 := pow_le_pow_left₀ hq0 hq 2
  have hT : T ≤ q ^ 2 + 2 ^ 12 * q ^ 2 := by
    have e1 : exp (-(β * ((lam12_2 : ℝ) / (20 * scale) - rad - (θ - rad)))) ≤ exp (-14) ^ 2 := by
      rw [← exp_nat_mul]
      exact exp_le_exp.2 (by norm_num [β, θ, lam12_1, lam12_2, scale, rad])
    have e2 : exp (-(β * (((-230 + 5 * 12 : ℤ) : ℝ) / 20 - ((12 : ℕ) : ℝ) / 4 - (θ - rad)))) ≤
        exp (-14) ^ 2 := by
      rw [← exp_nat_mul]
      exact exp_le_exp.2 (by norm_num [β, θ, lam12_1, scale, rad])
    have e3 := mul_le_mul_of_nonneg_left (e2.trans hq2) (by norm_num : (0 : ℝ) ≤ 2 ^ 12)
    simp only [T]
    exact add_le_add (e1.trans hq2) e3
  have hT0 : 0 ≤ T := by positivity
  have hfinal := ratio_ge (Z₁ := isingModel.partitionFunction 12 β)
    (Z₂ := isingModel.partitionFunction 12 (2 * β)) (exp_pos (-(β * (θ - rad))))
    (by rw [RPRing.partitionFunction_eq_sum]; positivity) hu ?_ hT0 hT
    (x := 4 * β * rad)
  · refine le_trans ?_ hfinal
    norm_num [q, β, rad]
  · refine le_trans (le_of_eq ?_) hl
    rw [← exp_nat_mul, ← exp_add]
    · congr 1; push_cast; ring

/-- The first of the thirty grid points of `IsingLongitudinal.points`. -/
abbrev firstPoint : IsingLongitudinal.Point :=
  IsingLongitudinal.points[0]'(by decide +kernel)

/-- **The first grid point of `IsingLongitudinal.seeds`**, with no hypotheses. -/
theorem first_point_seed :
    1 - (firstPoint.cert.sAt 0 : ℝ) ≤
      (IsingLongitudinal.model firstPoint.g firstPoint.h).spatialPurity 6 firstPoint.β₀ ∧
    1 - (firstPoint.cert.tAt 0 : ℝ) ≤
      (IsingLongitudinal.model firstPoint.g firstPoint.h).thermalPurity 12 firstPoint.β₀ := by
  constructor
  · convert seed_spatial using 1 <;>
      norm_num [firstPoint, IsingLongitudinal.points, ChainCertificate.sAt]
  · convert seed_thermal using 1 <;>
      norm_num [firstPoint, IsingLongitudinal.points, ChainCertificate.tAt]

/-- The first grid point lies in `IsingLongitudinal.points`. -/
theorem firstPoint_mem : firstPoint ∈ IsingLongitudinal.points := List.getElem_mem _

/-- **An unconditional spectral gap.** The Ising chain `-∑ ZZ - (1/4) ∑ X - (1/20) ∑ Z`: every
even ring of `L ≥ 6` sites has a unique ground state and gap `Δ_L > log 13480 / (97/4)`
(`≈ 0.3921`). No hypothesis and no `sorry`: the seed is `first_point_seed`. -/
theorem gap_first_point {L : ℕ} (hL : 6 ≤ L) (hLe : Even L) :
    GapGT ((IsingLongitudinal.model firstPoint.g firstPoint.h).spectrum L)
      (log (firstPoint.cert.r : ℝ) / firstPoint.β₀) :=
  (IsingLongitudinal.model firstPoint.g firstPoint.h).gapGT_of_chainCertificate (n₀ := 6)
    (by norm_num) (by decide) (by exact_mod_cast IsingLongitudinal.points_β₀_pos _ firstPoint_mem)
    (IsingLongitudinal.points_valid _ firstPoint_mem) ⟨rfl, rfl, rfl⟩ first_point_seed.1
    first_point_seed.2 hL hLe

/-- The same theorem with the numbers written out. -/
theorem ising_gap {L : ℕ} (hL : 6 ≤ L) (hLe : Even L) :
    GapGT ((IsingLongitudinal.model (1 / 4) (1 / 20)).spectrum L) (log 13480 / (97 / 4)) := by
  have h := gap_first_point hL hLe
  norm_num [firstPoint, IsingLongitudinal.points] at h ⊢
  exact h

end SeedKernel.Ising

#print axioms SeedKernel.Ising.seed_spatial
#print axioms SeedKernel.Ising.seed_thermal
#print axioms SeedKernel.Ising.ising_gap
