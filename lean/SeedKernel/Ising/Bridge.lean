/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Ising.Checker
import SeedKernel.Ising.Chunked
import SeedKernel.Ising.Hamiltonian

/-!
# From the integer checker to the Ising Hamiltonian

The kernel checks of `Checker.lean` and `Chunked.lean` work on `Nat` basis indices `s < 2^L`
(bit `j` of `s` is the spin at site `j`) and on integer vectors packed into one `Nat`. This
file identifies them with the matrix `IsingLongitudinal.model (1/4) (1/20)`:

* `encode σ` is the index of a configuration (`finFunctionFinEquiv`); flipping spin `j` is
  `xor 2^j` (`encode_flipAt`);
* `20 E(σ)` is `diagonal L (encode σ)` (`diagonal_encode`), and `20 (H x)(σ)` is
  `applyH L v (encode σ)` for `x = vec v` (`applyH_encode`);
* the `foldl` sums of the checker are sums over all configurations (`dot_eq`, `residualSq_eq`,
  `classicalCount_eq`, `chunkStats_sum`).
-/

open Matrix Complex CertifiedGaps

namespace SeedKernel.Ising

variable {L : ℕ}

/-! ### Folds as sums -/

theorem foldl_add_eq_sum {β : Type*} [AddCommMonoid β] (f : ℕ → β) (b : β) (N : ℕ) :
    (List.range N).foldl (fun a s ↦ a + f s) b = b + ∑ s ∈ Finset.range N, f s := by
  induction N generalizing b with
  | zero => simp
  | succ N ih => rw [List.range_succ, List.foldl_append, ih, Finset.sum_range_succ]; simp [add_assoc]

/-! ### Encoding configurations as integers -/

/-- The index `∑ⱼ σⱼ 2ʲ` of a configuration. -/
def encode (σ : Fin L → Fin 2) : ℕ := (finFunctionFinEquiv σ : ℕ)

theorem encode_lt (σ : Fin L → Fin 2) : encode σ < 2 ^ L := (finFunctionFinEquiv σ).2

theorem encode_div_mod (σ : Fin L → Fin 2) (j : Fin L) : encode σ / 2 ^ (j : ℕ) % 2 = σ j := by
  have h := finFunctionFinEquiv_symm_apply_val (finFunctionFinEquiv σ) j
  rw [Equiv.symm_apply_apply] at h
  exact h.symm

theorem testBit_encode (σ : Fin L → Fin 2) (i : ℕ) :
    (encode σ).testBit i = decide (∃ h : i < L, σ ⟨i, h⟩ = 1) := by
  by_cases hi : i < L
  · rw [Nat.testBit_eq_decide_div_mod_eq, encode_div_mod σ ⟨i, hi⟩]
    simp only [hi, exists_true_left]
    congr 1
    exact propext ⟨fun h ↦ Fin.ext h, fun h ↦ by rw [h]; rfl⟩
  · rw [Nat.testBit_lt_two_pow ((encode_lt σ).trans_le (Nat.pow_le_pow_right (by norm_num)
      (by omega)))]
    simp [hi]

theorem encode_injective : Function.Injective (encode (L := L)) :=
  fun _ _ h ↦ finFunctionFinEquiv.injective (Fin.ext h)

/-- Flipping spin `j` is `xor 2ʲ` on indices. -/
theorem encode_flipAt (j : Fin L) (σ : Fin L → Fin 2) :
    encode (flipAt j σ) = encode σ ^^^ (1 <<< (j : ℕ)) := by
  apply Nat.eq_of_testBit_eq
  intro i
  rw [Nat.testBit_xor, Nat.one_shiftLeft, Nat.testBit_two_pow, testBit_encode, testBit_encode]
  by_cases hi : i < L
  · by_cases hij : (j : ℕ) = i
    · have hj : (⟨i, hi⟩ : Fin L) = j := Fin.ext hij.symm
      simp only [hi, exists_true_left, hij, decide_true, hj, flipAt, Function.update_self]
      generalize σ j = a
      fin_cases a <;> rfl
    · have hj : (⟨i, hi⟩ : Fin L) ≠ j := fun h ↦ hij (by rw [← h])
      simp [hi, hij, flipAt, Function.update_of_ne hj]
  · have : (j : ℕ) ≠ i := by omega
    simp [hi, this]

/-- Sums over configurations are sums over indices `s < 2^L`. -/
theorem sum_encode {β : Type*} [AddCommMonoid β] (F : ℕ → β) :
    ∑ σ : Fin L → Fin 2, F (encode σ) = ∑ s ∈ Finset.range (2 ^ L), F s := by
  rw [← Fin.sum_univ_eq_sum_range]
  exact Equiv.sum_comp finFunctionFinEquiv (fun s : Fin (2 ^ L) ↦ F s)

/-! ### The classical energies and the Hamiltonian -/

theorem z_encode (σ : Fin L → Fin 2) (j : Fin L) : (z (encode σ) j : ℝ) = zR (σ j) := by
  rw [z, Nat.shiftRight_eq_div_pow, encode_div_mod]
  generalize σ j = a
  fin_cases a <;> norm_num [zR]

theorem finRotate_val (j : Fin L) : (finRotate L j : ℕ) = ((j : ℕ) + 1) % L := by
  rcases L with _ | L
  · exact j.elim0
  rw [finRotate_succ_apply, Fin.val_add]
  rcases L with _ | L
  · simp [Fin.val_zero]
  · simp

theorem diagonal_encode (σ : Fin L → Fin 2) :
    (diagonal L (encode σ) : ℝ) = 20 * diagE (1 / 20) σ := by
  rw [diagonal, show (fun a j ↦ a - 20 * z (encode σ) j * z (encode σ) ((j + 1) % L) -
      z (encode σ) j) = (fun a j ↦ a + (-(20 * z (encode σ) j * z (encode σ) ((j + 1) % L)) -
      z (encode σ) j)) from funext₂ fun _ _ ↦ by ring, foldl_add_eq_sum, zero_add,
    ← Fin.sum_univ_eq_sum_range (fun j ↦ (-(20 * z (encode σ) j * z (encode σ) ((j + 1) % L)) -
      z (encode σ) j))]
  push_cast
  simp only [diagE]
  have hrot (j : Fin L) : (z (encode σ) ((j + 1) % L) : ℝ) = zR (σ (finRotate L j)) := by
    rw [← z_encode, finRotate_val]
  simp only [z_encode, hrot, Finset.sum_sub_distrib, Finset.sum_neg_distrib]
  simp only [mul_assoc, ← Finset.mul_sum]
  ring

/-- The integer vector stored in the packed `Nat` `v`, as a vector on configurations. -/
def vec (v : ℕ) (σ : Fin L → Fin 2) : ℂ := ((entry v (encode σ) : ℤ) : ℂ)

/-- `20 (H x)(σ) = applyH L v (encode σ)` at `g = 1/4`, `h = 1/20`. -/
theorem applyH_encode (v : ℕ) (σ : Fin L → Fin 2) :
    ((applyH L v (encode σ) : ℤ) : ℂ) =
      20 * ((IsingLongitudinal.model (1 / 4) (1 / 20)).hamiltonian L *ᵥ vec v) σ := by
  rw [hamiltonian_mulVec_apply, applyH, foldl_add_eq_sum, zero_add,
    ← Fin.sum_univ_eq_sum_range (fun j ↦ entry v (encode σ ^^^ (1 <<< j)))]
  simp only [vec, ← encode_flipAt]
  push_cast
  rw [show ((diagonal L (encode σ) : ℤ) : ℂ) = ((diagonal L (encode σ) : ℝ) : ℂ) by norm_cast,
    diagonal_encode]
  push_cast
  ring

/-! ### The checker's sums -/

theorem dot_self_eq (v : ℕ) :
    ((dot L v v : ℤ) : ℝ) = normSqSum (vec (L := L) v) := by
  rw [dot, foldl_add_eq_sum, zero_add, ← sum_encode]
  push_cast
  simp [normSqSum, vec, normSq_apply]

/-- The residual of `vec v` against `θ = lam / (20 S)`: `‖H x - θ x‖² = residualSq / (20 S)²`. -/
theorem residualSq_eq (v : ℕ) (lam : ℤ) :
    ((residualSq L v lam : ℕ) : ℝ) = (20 * (scale : ℝ)) ^ 2 * normSqSum
      ((IsingLongitudinal.model (1 / 4) (1 / 20)).hamiltonian L *ᵥ vec v -
        ((lam / (20 * scale) : ℝ) : ℂ) • vec v) := by
  rw [residualSq, foldl_add_eq_sum, zero_add, ← sum_encode, normSqSum, Finset.mul_sum]
  push_cast
  refine Finset.sum_congr rfl fun σ _ ↦ ?_
  have hS : (scale : ℝ) ≠ 0 := by norm_num [scale]
  set t : ℤ := scale * applyH L v (encode σ) - lam * entry v (encode σ)
  have key : (t : ℂ) = ((20 * scale : ℝ) : ℂ) *
      (((IsingLongitudinal.model (1 / 4) (1 / 20)).hamiltonian L *ᵥ vec v -
        ((lam / (20 * scale) : ℝ) : ℂ) • vec v : (Fin L → Fin 2) → ℂ) σ) := by
    simp only [t, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
    push_cast
    rw [applyH_encode]
    simp only [vec]
    have hS' : (scale : ℂ) ≠ 0 := by norm_num [scale]
    field_simp
  have h1 : ((t.natAbs : ℕ) : ℝ) ^ 2 = normSq (t : ℂ) := by
    rw [Nat.cast_natAbs, Int.cast_abs, sq_abs, ← Complex.ofReal_intCast, normSq_ofReal, sq]
  rw [h1, key, map_mul, normSq_ofReal, ← sq]
  push_cast
  rfl

end SeedKernel.Ising
