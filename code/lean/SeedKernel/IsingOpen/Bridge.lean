/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Ising.Bridge
import SeedKernel.IsingOpen.Checker
import SeedKernel.IsingOpen.Hamiltonian

/-!
# From the integer checker to the free Ising chain

The open-chain counterpart of `SeedKernel.Ising.Bridge`, with the encoding `encode` of
configurations as indices `s < 2^n` and the packed vectors `vec v` of that file:

* `20 E(σ)` is `diagonal n (encode σ)` (`diagonal_encode`), and `20 (H x)(σ)` is
  `applyH n v (encode σ)` for `x = vec v` (`applyH_encode`);
* `residualSq` is the squared residual (`residualSq_eq`), `classicalCount` the size of the set
  `B` of configurations below the cutoff (`classicalCount_eq`), and `weightOK` at every index is
  the hypothesis of the weighted comparison `le_re_quadratic_of_weights` (`weights_spec`);
* block sums are whole-chain sums (`sumBlocks_eq`, `weightBlocks_all`).
-/

open Matrix Complex CertifiedGaps SeedKernel.Ising

namespace SeedKernel.IsingOpen

variable {k : ℕ}

/-! ### The classical energies and the Hamiltonian -/

theorem diagonal_encode (σ : Fin (k + 1) → Fin 2) :
    (diagonal (k + 1) (encode σ) : ℝ) = 20 * openE (1 / 20) σ := by
  rw [diagonal, Nat.add_sub_cancel,
    show (fun a j ↦ a - 20 * z (encode σ) j * z (encode σ) (j + 1)) =
      (fun a j ↦ a + -(20 * z (encode σ) j * z (encode σ) (j + 1))) from
      funext₂ fun _ _ ↦ by ring,
    foldl_add_eq_sum, foldl_add_eq_sum, zero_add, zero_add,
    ← Fin.sum_univ_eq_sum_range (fun j ↦ -(20 * z (encode σ) j * z (encode σ) (j + 1))),
    ← Fin.sum_univ_eq_sum_range (fun j ↦ z (encode σ) j)]
  push_cast
  have hc (j : Fin k) : (z (encode σ) j : ℝ) = zR (σ j.castSucc) := by
    rw [← z_encode, Fin.val_castSucc]
  have hs (j : Fin k) : (z (encode σ) (j + 1) : ℝ) = zR (σ j.succ) := by
    rw [← z_encode, Fin.val_succ]
  simp only [z_encode, hc, hs, openE, Finset.sum_neg_distrib]
  simp only [mul_assoc, ← Finset.mul_sum]
  ring

/-- `20 (H x)(σ) = applyH n v (encode σ)` at `g = 1/4`, `h = 1/20`. -/
theorem applyH_encode (v : ℕ) (σ : Fin (k + 1) → Fin 2) :
    ((applyH (k + 1) v (encode σ) : ℤ) : ℂ) = 20 * (freeH (1 / 4) (1 / 20) k *ᵥ vec v) σ := by
  rw [openHamiltonian_mulVec_apply, applyH, foldl_add_eq_sum, zero_add,
    ← Fin.sum_univ_eq_sum_range (fun j ↦ entry v (encode σ ^^^ (1 <<< j)))]
  simp only [Ising.vec, ← encode_flipAt]
  push_cast
  rw [show ((diagonal (k + 1) (encode σ) : ℤ) : ℂ) =
      ((diagonal (k + 1) (encode σ) : ℝ) : ℂ) by norm_cast, IsingOpen.diagonal_encode]
  push_cast
  ring

/-! ### The checker's sums -/

/-- The residual of `vec v` against `θ = lam / (20 S)`: `‖H x - θ x‖² = residualSq / (20 S)²`. -/
theorem residualSq_eq (v : ℕ) (lam : ℤ) :
    ((residualSq (k + 1) v lam : ℕ) : ℝ) = (20 * (scale : ℝ)) ^ 2 * normSqSum
      (freeH (1 / 4) (1 / 20) k *ᵥ vec v - ((lam / (20 * scale) : ℝ) : ℂ) • vec v) := by
  rw [residualSq, foldl_add_eq_sum, zero_add, ← sum_encode, normSqSum, Finset.mul_sum]
  push_cast
  refine Finset.sum_congr rfl fun σ _ ↦ ?_
  set t : ℤ := scale * applyH (k + 1) v (encode σ) - lam * entry v (encode σ)
  have key : (t : ℂ) = ((20 * scale : ℝ) : ℂ) *
      ((freeH (1 / 4) (1 / 20) k *ᵥ vec v -
        ((lam / (20 * scale) : ℝ) : ℂ) • vec v : (Fin (k + 1) → Fin 2) → ℂ) σ) := by
    simp only [t, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
    push_cast
    rw [IsingOpen.applyH_encode]
    simp only [Ising.vec]
    have hS' : (scale : ℂ) ≠ 0 := by norm_num [scale]
    field_simp
  have h1 : ((t.natAbs : ℕ) : ℝ) ^ 2 = normSq (t : ℂ) := by
    rw [Nat.cast_natAbs, Int.cast_abs, sq_abs, ← Complex.ofReal_intCast, normSq_ofReal, sq]
  rw [h1, key, map_mul, normSq_ofReal, ← sq]
  push_cast
  rfl

/-- The configurations whose classical energy lies below the cutoff `cut / 20`. -/
def lowSet (k : ℕ) (cut : ℤ) : Finset (Fin (k + 1) → Fin 2) :=
  Finset.univ.filter fun σ ↦ diagonal (k + 1) (encode σ) < cut

theorem classicalCount_eq (cut : ℤ) : classicalCount (k + 1) cut = (lowSet k cut).card := by
  rw [classicalCount, show (fun a s ↦ if diagonal (k + 1) s < cut then a + 1 else a) =
      (fun a s ↦ a + if diagonal (k + 1) s < cut then 1 else 0) from
    funext₂ fun a s ↦ by split_ifs <;> rfl, foldl_add_eq_sum, zero_add, ← sum_encode,
    lowSet, Finset.card_filter]

theorem flipSum_encode (w : ℕ) (cut : ℤ) (σ : Fin (k + 1) → Fin 2) :
    ((flipSum (k + 1) w cut (encode σ) : ℤ) : ℝ) =
      ∑ j, (if flipAt j σ ∈ lowSet k cut then 0 else (entry w (encode (flipAt j σ)) : ℝ)) := by
  rw [flipSum, show (fun a j ↦ if diagonal (k + 1) (encode σ ^^^ (1 <<< j)) < cut then a else
      a + entry w (encode σ ^^^ (1 <<< j))) = (fun a j ↦ a +
      if diagonal (k + 1) (encode σ ^^^ (1 <<< j)) < cut then 0 else
        entry w (encode σ ^^^ (1 <<< j))) from funext₂ fun a j ↦ by split_ifs <;> simp,
    foldl_add_eq_sum, zero_add, ← Fin.sum_univ_eq_sum_range (fun j ↦
      if diagonal (k + 1) (encode σ ^^^ (1 <<< j)) < cut then 0 else
        entry w (encode σ ^^^ (1 <<< j)))]
  push_cast
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [← encode_flipAt]
  simp only [lowSet, Finset.mem_filter, Finset.mem_univ, true_and]

/-- **The weight checks give the comparison hypothesis.** If `weightOK` holds at every index, the
weights `w(σ) = entry w (encode σ)` satisfy the hypothesis of `le_re_quadratic_of_weights` with
`B = lowSet k cut` and `ρ = rho / 20`. -/
theorem weights_spec (w : ℕ) (cut rho : ℤ)
    (hw : ∀ s < 2 ^ (k + 1), weightOK (k + 1) w cut rho s = true) (σ : Fin (k + 1) → Fin 2)
    (hσ : σ ∉ lowSet k cut) :
    0 < (entry w (encode σ) : ℝ) ∧
      1 / 4 * ∑ j, (if flipAt j σ ∈ lowSet k cut then 0 else (entry w (encode (flipAt j σ)) : ℝ))
        ≤ (entry w (encode σ) : ℝ) * (openE (1 / 20) σ - rho / 20) := by
  have h := hw _ (encode_lt σ)
  have hlow : ¬ diagonal (k + 1) (encode σ) < cut := by simpa [lowSet] using hσ
  simp only [weightOK, hlow, decide_false, Bool.false_or, Bool.and_eq_true,
    decide_eq_true_eq] at h
  obtain ⟨hpos, hle⟩ := h
  have hle' : ((5 * flipSum (k + 1) w cut (encode σ) : ℤ) : ℝ) ≤
      ((entry w (encode σ) * (diagonal (k + 1) (encode σ) - rho) : ℤ) : ℝ) := by
    exact_mod_cast hle
  push_cast at hle'
  rw [flipSum_encode, IsingOpen.diagonal_encode] at hle'
  refine ⟨by exact_mod_cast hpos, ?_⟩
  linarith

/-! ### Blocks -/

theorem sum_range_blocks {M : Type*} [AddCommMonoid M] (F : ℕ → M) (c b : ℕ) :
    ∑ i ∈ Finset.range b, ∑ j ∈ Finset.range c, F (c * i + j) =
      ∑ s ∈ Finset.range (c * b), F s := by
  induction b with
  | zero => simp
  | succ b ih => rw [Finset.sum_range_succ, ih, Nat.mul_succ, Finset.sum_range_add]

theorem blockStats_eq_sum (n v : ℕ) (lam : ℤ) (start count : ℕ) :
    blockStats n v lam start count = ∑ i ∈ Finset.range count,
      ((entry v (start + i) * entry v (start + i),
        (scale * applyH n v (start + i) - lam * entry v (start + i)).natAbs ^ 2) : ℤ × ℕ) := by
  have hadd : ∀ (acc : ℤ × ℕ) (i : ℕ), (acc.1 + entry v (start + i) * entry v (start + i),
      acc.2 + (scale * applyH n v (start + i) - lam * entry v (start + i)).natAbs ^ 2) =
      acc + ((entry v (start + i) * entry v (start + i),
        (scale * applyH n v (start + i) - lam * entry v (start + i)).natAbs ^ 2) : ℤ × ℕ) :=
    fun _ _ ↦ rfl
  unfold blockStats
  simp only [hadd]
  rw [show ((0, 0) : ℤ × ℕ) = 0 from rfl, foldl_add_eq_sum, zero_add]

theorem sumBlocks_map_range (f : ℕ → ℤ × ℕ) (N : ℕ) :
    sumBlocks ((List.range N).map f) = ∑ b ∈ Finset.range N, f b := by
  have hadd : ∀ a x : ℤ × ℕ, (a.1 + x.1, a.2 + x.2) = a + x := fun _ _ ↦ rfl
  rw [sumBlocks, List.foldl_map]
  simp only [hadd]
  rw [show ((0, 0) : ℤ × ℕ) = 0 from rfl, foldl_add_eq_sum, zero_add]

/-- Block statistics are the whole-chain squared norm and squared residual, for `c k = 2ⁿ`.
Proved for a symbolic vector so that the kernel never unfolds the concrete witness. -/
theorem sumBlocks_eq (n v : ℕ) (lam : ℤ) (c b : ℕ) (hcb : c * b = 2 ^ n) :
    sumBlocks ((List.range b).map fun i ↦ blockStats n v lam (c * i) c) =
      (dot n v v, residualSq n v lam) := by
  rw [sumBlocks_map_range]
  simp only [blockStats_eq_sum]
  rw [sum_range_blocks (fun s ↦ ((entry v s * entry v s,
      (scale * applyH n v s - lam * entry v s).natAbs ^ 2) : ℤ × ℕ)), hcb]
  simp only [dot, residualSq, foldl_add_eq_sum, zero_add]
  ext <;> simp [Prod.fst_sum, Prod.snd_sum]

/-- Passing weight checks on `b` blocks of `c` indices cover every index below `c b`. -/
theorem weightBlocks_all (n w : ℕ) (cut rho : ℤ) (c b : ℕ)
    (h : ((List.range b).map fun i ↦ weightBlock n w cut rho (c * i) c).all id = true) :
    ∀ s < c * b, weightOK n w cut rho s = true := by
  intro s hs
  have hc : 0 < c := Nat.pos_of_ne_zero fun h0 ↦ by simp [h0] at hs
  have hb : s / c < b := (Nat.div_lt_iff_lt_mul hc).2 (by rwa [mul_comm])
  simp only [List.all_map, List.all_eq_true, List.mem_range, Function.comp_apply, id] at h
  have hblock := h _ hb
  simp only [weightBlock, List.all_eq_true, List.mem_range] at hblock
  have := hblock (s % c) (Nat.mod_lt _ hc)
  rwa [Nat.div_add_mod] at this

end SeedKernel.IsingOpen
