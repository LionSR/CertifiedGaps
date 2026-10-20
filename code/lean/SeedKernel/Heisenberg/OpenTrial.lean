/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.OpenTrialBound
import SeedKernel.Heisenberg.TrialBlocks
import SeedKernel.Heisenberg.TrialIntertwiner
import SeedKernel.Heisenberg.OpenTrialSpec
import SeedKernel.Heisenberg.OpenTrialBlocks1
import SeedKernel.Heisenberg.OpenTrialBlocks2

/-!
# The opposite-field trial bounds (kernel-checked)

OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (appendix "Exact trial vectors
for opposite fields", Proposition `apptr:bound`) bounds the ground-state energy of the
opposite-field chain

  `H_m^{opp} = ∑ⱼ Sⱼ · Sⱼ₊₁ - (3/5) S₁ᶻ + (3/5) Sₘᶻ`

by `E₀ + a m < α`, `a = 700741/500000`, `α = 2021/5000`, at `m = 34, 78, 142, 240`, with the
open matrix-product trial vector `ψ(σ) = bᵀ A_{σ₁} ⋯ A_{σₘ} c`: the bond-dimension-26 tensors
`trialA` of the periodic trial and the endpoint columns `b`, `c` of the paper's endpoint table.
This file proves the four bounds (`oppositeField_trial_34`, …) from the kernel checks of
`OpenTrialBlocks1` and `OpenTrialBlocks2`.

The proof has three parts.

1. **Reduction.** The intertwiner `W` of the periodic trial (`T W = W B`, `O W = W C`,
   `V W = c 1`, kernel-checked in `TrialChecks*`) turns every contraction `x P y` of products
   `P` of `T` and `O` into `c⁻¹ (x W) Q (V y)` with `Q` block diagonal (`reduce`). The row
   vectors `L₀ W` and `L₀ G W` vanish outside the copy `k = Jᵦ` of each block, so each contraction
   is a sum of one contraction per block of size `mᵦ ≤ 32`.
2. **Block iterations.** `Params.run_spec` identifies the kernel's outputs with those block
   contractions.
3. **Arithmetic.** The stored outputs satisfy the four integer comparisons
   (`totals_lt`, `decide +kernel`), which give the hypothesis of
   `exists_oppositeFieldSpectrum_lt_of_openMPS`.

## Main results

* `oppositeField_trial_34`, `oppositeField_trial_78`, `oppositeField_trial_142`,
  `oppositeField_trial_240`: `∃ i, oppositeFieldSpectrum (3/5) m i < -a m + α`.
-/

namespace SeedKernel.OpenTrial

open Matrix SeedKernel.Trial CertifiedGaps CertifiedGaps.SpinOneHeisenberg
open scoped Kronecker

set_option maxRecDepth 100000

/-- Every block passes its kernel check. -/
theorem openBlockCheck_all : ∀ b : Fin 11, OpenBlockCheck b
  | ⟨0, _⟩ => openBlockCheck_0
  | ⟨1, _⟩ => openBlockCheck_1
  | ⟨2, _⟩ => openBlockCheck_2
  | ⟨3, _⟩ => openBlockCheck_3
  | ⟨4, _⟩ => openBlockCheck_4
  | ⟨5, _⟩ => openBlockCheck_5
  | ⟨6, _⟩ => openBlockCheck_6
  | ⟨7, _⟩ => openBlockCheck_7
  | ⟨8, _⟩ => openBlockCheck_8
  | ⟨9, _⟩ => openBlockCheck_9
  | ⟨10, _⟩ => openBlockCheck_10

/-! ### The trial vector -/

/-- The left endpoint `b` of the trial vector. -/
def openLeft (a : Fin 26) : ℤ := leftN a

/-- The right endpoint `c` of the trial vector. -/
def openRight (a : Fin 26) : ℤ := rightN a

/-- The transfer matrix `T`. -/
abbrev Tm : Matrix Pair Pair ℤ := MPS.transfer trialA
/-- The bond insertion `O`. -/
abbrev Om : Matrix Pair Pair ℤ := MPS.insertion exchangeZ trialA
/-- The field insertion `G = ∑ₛ mₛ Aₛ ⊗ Aₛ`. -/
abbrev Gm : Matrix Pair Pair ℤ := OpenMPS.siteInsertion szDiag trialA

/-! ### The intertwiner -/

/-- `T W = W B` (kernel-checked in `TrialChecks`). -/
theorem T_mul_W : Tm * Wm = Wm * Bm :=
  T_mul_Wm wIdxOk wAbsOk aBoundOk (fun b ↦ (blockCheck_all b).1) tB_offset checkT

/-- `O W = W C` (kernel-checked in `TrialChecksInsertion1`–`TrialChecksInsertion4`). -/
theorem O_mul_W : Om * Wm = Wm * Cm :=
  O_mul_Wm wIdxOk wAbsOk aBoundOk (fun b ↦ (blockCheck_all b).2.1) oB_offset checkO1 checkO2
    checkO3 checkO4

/-- `W V = c 1`: `V W = c 1` for square matrices, through `ℚ`. -/
theorem W_mul_V : Wm * Vm = (trialScale : ℤ) • (1 : Matrix Pair Pair ℤ) := by
  have hV := Vm_mul_Wm wIdxOk wAbsOk vIdxOk vAbsOk checkV
  set f := Int.castRingHom ℚ
  have hc : (trialScale : ℚ) ≠ 0 := by norm_num [trialScale]
  have h1 : ((trialScale : ℚ)⁻¹ • Vm.map f) * Wm.map f = 1 := by
    rw [Matrix.smul_mul, ← Matrix.map_mul, hV]
    ext i j
    simp only [Matrix.smul_apply, map_apply, Matrix.one_apply, smul_eq_mul]
    split_ifs <;> simp [f, hc]
  have h2 := (mul_eq_one_comm_of_card_eq Iota Pair ℚ card_Iota).mp h1
  ext i j
  have h3 := congrFun (congrFun h2 i) j
  rw [Matrix.mul_smul, ← Matrix.map_mul, Matrix.smul_apply, map_apply, smul_eq_mul] at h3
  have h4 : ((Wm * Vm) i j : ℚ) = (trialScale : ℚ) * (1 : Matrix Pair Pair ℚ) i j := by
    rw [← h3]; field_simp; rfl
  simp only [Matrix.smul_apply, Matrix.one_apply, smul_eq_mul] at h4 ⊢
  split_ifs at h4 ⊢ <;> exact_mod_cast h4

/-- **Reduction by the intertwiner**: `c (x P y) = (x W) Q (V y)` when `P W = W Q`. -/
theorem intertwine {P : Matrix Pair Pair ℤ} {Q : Matrix Iota Iota ℤ} (h : P * Wm = Wm * Q)
    (x y : Pair → ℤ) :
    (trialScale : ℤ) * (x ⬝ᵥ (P *ᵥ y)) = (x ᵥ* Wm) ⬝ᵥ (Q *ᵥ (Vm *ᵥ y)) := by
  have hy : (trialScale : ℤ) • y = Wm *ᵥ (Vm *ᵥ y) := by
    rw [mulVec_mulVec, W_mul_V, smul_mulVec, one_mulVec]
  rw [← smul_eq_mul, ← dotProduct_smul, ← mulVec_smul, hy, mulVec_mulVec, h, ← mulVec_mulVec,
    dotProduct_mulVec]

/-- `Tᵏ W = W Bᵏ`. -/
theorem pow_mul_W (k : ℕ) : Tm ^ k * Wm = Wm * Bm ^ k := by
  induction k with
  | zero => simp
  | succ k ih => rw [pow_succ, Matrix.mul_assoc, T_mul_W, ← Matrix.mul_assoc, ih,
      Matrix.mul_assoc, ← pow_succ]

/-- `Tʲ O Tᵏ W = W Bʲ C Bᵏ`. -/
theorem bond_mul_W (j k : ℕ) : Tm ^ j * Om * Tm ^ k * Wm = Wm * (Bm ^ j * Cm * Bm ^ k) := by
  rw [Matrix.mul_assoc, pow_mul_W, ← Matrix.mul_assoc, Matrix.mul_assoc _ Om, O_mul_W,
    ← Matrix.mul_assoc, pow_mul_W, Matrix.mul_assoc, Matrix.mul_assoc, Matrix.mul_assoc]

/-- `Bᵏ` is block diagonal with blocks `tᵦᵏ`. -/
theorem Bm_pow (k : ℕ) : Bm ^ k = copies bd bm fun b ↦ tMat b ^ k := copies_pow _ k

/-- `Bʲ C Bᵏ` is block diagonal with blocks `tᵦʲ oᵦ tᵦᵏ`. -/
theorem Bm_bond (j k : ℕ) :
    Bm ^ j * Cm * Bm ^ k = copies bd bm fun b ↦ tMat b ^ j * oMat b * tMat b ^ k := by
  rw [Bm_pow, Bm_pow, Cm, copies_mul, copies_mul]

/-! ### One copy per block -/

/-- The weight-zero copy `k = Jᵦ` of block `b`. -/
def jIdx (b : Fin 11) : Fin (bd b) := ⟨bJ b, by unfold bd; omega⟩

/-- The copy `Jᵦ` of a vector on the columns of `W`. -/
def copyJ (x : Iota → ℤ) (b : Fin 11) : Fin (bm b) → ℤ := fun c ↦ x ⟨b, (jIdx b, c)⟩

/-- A vector vanishing outside the copies `Jᵦ` pairs with a block matrix block by block. -/
theorem dotProduct_copies_mulVec (M : ∀ b, Matrix (Fin (bm b)) (Fin (bm b)) ℤ) (x y : Iota → ℤ)
    (hx : ∀ b (k : Fin (bd b)) c, k ≠ jIdx b → x ⟨b, (k, c)⟩ = 0) :
    x ⬝ᵥ (copies bd bm M *ᵥ y) = ∑ b, copyJ x b ⬝ᵥ (M b *ᵥ copyJ y b) := by
  rw [dotProduct, Fintype.sum_sigma]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [Fintype.sum_prod_type, Finset.sum_eq_single (jIdx b)]
  · simp only [dotProduct, copyJ]
    refine Finset.sum_congr rfl fun c _ ↦ ?_
    rw [mulVec_copies]
    rfl
  · intro k _ hk
    exact Finset.sum_eq_zero fun c _ ↦ by rw [hx b k c hk, zero_mul]
  · simp

/-- **Reduction to the blocks**: `c (x P y) = ∑ᵦ (x W)_J Mᵦ (V y)_J` when `P W = W copies M` and
`x W` vanishes outside the copies `Jᵦ`. -/
theorem reduce {P : Matrix Pair Pair ℤ} {M : ∀ b, Matrix (Fin (bm b)) (Fin (bm b)) ℤ}
    (h : P * Wm = Wm * copies bd bm M) (x y : Pair → ℤ)
    (hx : ∀ b (k : Fin (bd b)) c, k ≠ jIdx b → (x ᵥ* Wm) ⟨b, (k, c)⟩ = 0) :
    (trialScale : ℤ) * (x ⬝ᵥ (P *ᵥ y)) =
      ∑ b, copyJ (x ᵥ* Wm) b ⬝ᵥ (M b *ᵥ copyJ (Vm *ᵥ y) b) := by
  rw [intertwine h, dotProduct_copies_mulVec M _ _ hx]

/-! ### The reduced vectors from their lists -/

/-- The entries of a concatenation of sparse rows add. -/
theorem sparseEntry_append (g g' : List (ℕ × ℤ)) (j : ℕ) :
    sparseEntry (g ++ g') j = sparseEntry g j + sparseEntry g' j := by
  simp [sparseEntry]

/-- Dropping zero values does not change the entries of a sparse row. -/
theorem sparseEntry_filter (g : List (ℕ × ℤ)) (j : ℕ) :
    sparseEntry (g.filter fun t ↦ t.2 ≠ 0) j = sparseEntry g j := by
  induction g with
  | nil => rfl
  | cons t g ih =>
    by_cases h : t.2 = 0
    · rw [List.filter_cons_of_neg (by simp [h]), ih]
      simp [sparseEntry, h]
    · rw [List.filter_cons_of_pos (by simp [h])]
      simp only [sparseEntry, List.map_cons, List.sum_cons] at ih ⊢
      rw [ih]

/-- The entries of the sparse row `(j, v j)`, `j < n`. -/
theorem sparseEntry_range_map (v : ℕ → ℤ) (a : ℕ) :
    ∀ n, sparseEntry ((List.range n).map fun j ↦ (j, v j)) a = if a < n then v a else 0
  | 0 => by simp [sparseEntry]
  | n + 1 => by
    rw [List.range_succ, List.map_append, sparseEntry_append, sparseEntry_range_map v a n]
    simp only [sparseEntry, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
    by_cases h : a < n
    · simp [h, show a < n + 1 by omega, show n ≠ a by omega]
    · by_cases h' : a = n
      · subst h'; simp
      · simp [h, show ¬ a < n + 1 by omega, Ne.symm h']

/-- `b Aₛ` is represented by the sparse list `leftA s`. -/
theorem vecMul_trialA_eq (s : Fin 3) (a : Fin 26) :
    (openLeft ᵥ* trialA s) a = sparseEntry (leftA s) a := by
  rw [leftA, sparseEntry_filter, sparseEntry_range_map]
  simp only [a.isLt, ↓reduceIte, vecMul, dotProduct, openLeft]
  rw [show ((a : ℕ)) = ((a : Fin 26) : ℕ) from rfl, sparseDot_trialSparseT]
  exact Finset.sum_congr rfl fun e _ ↦ mul_comm _ _

/-- The labels of `leftA s` are labels of the virtual space. -/
theorem leftA_lt (s : Fin 3) : ∀ t ∈ leftA s, t.1 < 26 := by
  intro t ht
  simp only [leftA, List.mem_filter, List.mem_map, List.mem_range] at ht
  obtain ⟨⟨j, hj, rfl⟩, _⟩ := ht
  exact hj

/-- The labels of the left endpoint are labels of the virtual space. -/
theorem openLeftSparse_lt : ∀ t ∈ openLeftSparse, t.1 < 26 := by decide

/-- `Aₛ c` is the list `rightA s`. -/
theorem mulVec_trialA_eq (s : Fin 3) (a : Fin 26) :
    (trialA s *ᵥ openRight) a = (rightA s).getD a 0 := by
  rw [rightA, getD_ofFn26, sparseDot_trialSparse]
  rfl

/-- The weighted groups of `f ⊗ f` as a double sparse sum. -/
theorem groupsEntry_map (f : List (ℕ × ℤ)) (k : ℤ) (G : ℕ → List (ℕ × ℤ)) (j : ℕ) :
    groupsEntry (f.map fun t' ↦ (k * t'.2, G t'.1)) j =
      k * sparseDot f fun a' ↦ sparseEntry (G a') j := by
  induction f with
  | nil => simp [groupsEntry, sparseDot]
  | cons t f ih =>
    simp only [groupsEntry, sparseDot, List.map_cons, List.sum_cons, List.map_map] at ih ⊢
    rw [ih]
    ring

/-- The double sparse sum over `g ⊗ f` is the sum over the weighted groups. -/
theorem sparseDot_pairGroups (g f : List (ℕ × ℤ)) (wg : List (List (ℕ × ℤ))) (j : ℕ) :
    sparseDot g (fun a ↦ sparseDot f fun a' ↦ sparseEntry (wg.getD (a' + 26 * a) []) j) =
      groupsEntry (g.flatMap fun t ↦ f.map fun t' ↦ (t.2 * t'.2, wg.getD (t'.1 + 26 * t.1) []))
        j := by
  induction g with
  | nil => simp [groupsEntry, sparseDot]
  | cons t g ih =>
    rw [List.flatMap_cons]
    simp only [groupsEntry, List.map_append, List.sum_append] at ih ⊢
    rw [← ih]
    have := groupsEntry_map f t.2 (fun a' ↦ wg.getD (a' + 26 * t.1) []) j
    simp only [groupsEntry] at this
    rw [this]
    simp [sparseDot]

/-- **`(f ⊗ f) W` from its weighted groups.** -/
theorem vecMul_Wm_pair {f : Fin 26 → ℤ} {fS : List (ℕ × ℤ)}
    (hf : ∀ a : Fin 26, f a = sparseEntry fS a) (hfS : ∀ t ∈ fS, t.1 < 26) (b : Fin 11)
    (y : Fin (bd b) × Fin (bm b)) :
    (OpenMPS.pair f ᵥ* Wm) ⟨b, y⟩ = groupsEntry (pairGroups fS (wG b)) (localIdx y) := by
  set F : ℕ → ℤ := fun q ↦ sparseEntry ((wG b).getD q []) (localIdx y)
  have hin (a : ℕ) :
      ∑ a' : Fin 26, f a' * F (a' + 26 * a) = sparseDot fS fun a' ↦ F (a' + 26 * a) := by
    simp_rw [hf]
    exact sum_sparseEntry_mul (fun a' ↦ F (a' + 26 * a)) fS hfS
  have hout : ∑ a : Fin 26, f a * sparseDot fS (fun a' ↦ F (a' + 26 * a)) =
      sparseDot fS fun a ↦ sparseDot fS fun a' ↦ F (a' + 26 * a) := by
    simp_rw [hf]
    exact sum_sparseEntry_mul (fun a ↦ sparseDot fS fun a' ↦ F (a' + 26 * a)) fS hfS
  rw [pairGroups, ← sparseDot_pairGroups, ← hout]
  simp only [vecMul, dotProduct, Fintype.sum_prod_type, OpenMPS.pair, Wm, pairIdx]
  refine Finset.sum_congr rfl fun a _ ↦ ?_
  rw [← hin, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a' _ ↦ ?_
  simp only [F]
  ring

/-- `(b ⊗ b) G = (bA₂ ⊗ bA₂) - (bA₀ ⊗ bA₀)`. -/
theorem pair_vecMul_G (f : Fin 26 → ℤ) :
    OpenMPS.pair f ᵥ* Gm = OpenMPS.pair (f ᵥ* trialA 2) - OpenMPS.pair (f ᵥ* trialA 0) := by
  simp only [Gm, OpenMPS.siteInsertion, Fin.sum_univ_three, vecMul_add, vecMul_smul,
    OpenMPS.vecMul_pair_kronecker]
  ext p
  simp [szDiag, mz, OpenMPS.pair]
  ring

/-- `G (c ⊗ c) = (A₂c ⊗ A₂c) - (A₀c ⊗ A₀c)`. -/
theorem G_mulVec_pair (f : Fin 26 → ℤ) :
    Gm *ᵥ OpenMPS.pair f = OpenMPS.pair (trialA 2 *ᵥ f) - OpenMPS.pair (trialA 0 *ᵥ f) := by
  simp only [Gm, OpenMPS.siteInsertion, Fin.sum_univ_three, add_mulVec, smul_mulVec,
    OpenMPS.kronecker_mulVec_pair]
  ext p
  simp [szDiag, mz, OpenMPS.pair]
  ring

/-- The sum of negated values is the negated sum. -/
theorem sum_map_neg {X : Type*} (g : X → ℤ) :
    ∀ l : List X, (l.map fun x ↦ -g x).sum = -(l.map g).sum
  | [] => by simp
  | x :: l => by simp [sum_map_neg g l]; ring

/-- `L₀ W` from the groups of `groupsL`. -/
theorem L0W_eq (b : Fin 11) (y : Fin (bd b) × Fin (bm b)) :
    (OpenMPS.pair openLeft ᵥ* Wm) ⟨b, y⟩ = groupsEntry (groupsL b) (localIdx y) :=
  vecMul_Wm_pair (fun _ ↦ rfl) openLeftSparse_lt b y

/-- `L₀ G W` from the groups of `groupsG`. -/
theorem L0GW_eq (b : Fin 11) (y : Fin (bd b) × Fin (bm b)) :
    ((OpenMPS.pair openLeft ᵥ* Gm) ᵥ* Wm) ⟨b, y⟩ = groupsEntry (groupsG b) (localIdx y) := by
  rw [pair_vecMul_G, sub_vecMul, Pi.sub_apply, vecMul_Wm_pair (vecMul_trialA_eq 2) (leftA_lt 2),
    vecMul_Wm_pair (vecMul_trialA_eq 0) (leftA_lt 0), groupsG]
  simp only [groupsEntry, List.map_append, List.sum_append, List.map_map]
  rw [sub_eq_add_neg, ← sum_map_neg]
  congr 2
  ext q
  simp [neg_mul]

/-- Groups in the copy `J` give zero at an index of another copy. -/
theorem groupsEntry_eq_zero {m J : ℕ} {g : List (ℤ × List (ℕ × ℤ))} (hg : InCopy m J g)
    {j : ℕ} (hj : j / m ≠ J) : groupsEntry g j = 0 := by
  unfold groupsEntry
  refine List.sum_eq_zero fun z hz ↦ ?_
  simp only [List.mem_map] at hz
  obtain ⟨q, hq, rfl⟩ := hz
  have : sparseEntry q.2 j = 0 := by
    unfold sparseEntry
    refine List.sum_eq_zero fun z hz ↦ ?_
    simp only [List.mem_map] at hz
    obtain ⟨t, ht, rfl⟩ := hz
    have hne : t.1 ≠ j := fun e ↦ hj (by rw [← e]; exact hg q hq t ht)
    simp [hne]
  rw [this, mul_zero]

/-- The copy of a local column `c + m k` is `k`. -/
theorem localIdx_div {b : Fin 11} (y : Fin (bd b) × Fin (bm b)) : localIdx y / bm b = y.1 := by
  have := y.2.isLt
  unfold localIdx
  rw [Nat.add_mul_div_left _ _ (by omega), Nat.div_eq_of_lt this, zero_add]

/-- A vector given by groups in the copy `Jᵦ` vanishes on the other copies. -/
theorem copy_zero {x : Iota → ℤ} {g : Fin 11 → List (ℤ × List (ℕ × ℤ))}
    (hx : ∀ b y, x ⟨b, y⟩ = groupsEntry (g b) (localIdx y)) (hg : ∀ b, InCopy (bm b) (bJ b) (g b))
    (b : Fin 11) (k : Fin (bd b)) (c : Fin (bm b)) (hk : k ≠ jIdx b) : x ⟨b, (k, c)⟩ = 0 := by
  rw [hx, groupsEntry_eq_zero (hg b)]
  rw [localIdx_div]
  exact fun e ↦ hk (Fin.ext e)

/-- Entry `i` of `List.ofFn f` is `f i`. -/
theorem getD_ofFn {n : ℕ} (f : Fin n → ℤ) (i : Fin n) : (List.ofFn f).getD i 0 = f i := by
  rw [List.getD_eq_getElem _ _ (by simp), List.getElem_ofFn]

/-- The copy `Jᵦ` of a row vector given by groups is the list `rowCopy`. -/
theorem copyJ_row {x : Iota → ℤ} {g : List (ℤ × List (ℕ × ℤ))} (b : Fin 11)
    (hx : ∀ y, x ⟨b, y⟩ = groupsEntry g (localIdx y)) :
    copyJ x b = toVec (bm b) (rowCopy b g) := by
  ext c
  rw [copyJ, hx, toVec, rowCopy, getD_ofFn]
  rfl

/-- The copy `Jᵦ` of `V y` is the list `colCopy`. -/
theorem copyJ_col (b : Fin 11) (Y : ℕ → ℤ) :
    copyJ (Vm *ᵥ fun i ↦ Y (pairIdx i)) b = toVec (bm b) (colCopy b Y) := by
  ext c
  rw [copyJ, mulVec_Vm vIdxOk, toVec, colCopy, getD_ofFn]
  rfl

/-- `u = c ⊗ c` on the doubled index is `uN`. -/
theorem pair_openRight : OpenMPS.pair openRight = fun i ↦ uN (pairIdx i) := by
  ext i
  have h1 : ((i.2 : ℕ) + 26 * i.1) / 26 = i.1 := by omega
  have h2 : ((i.2 : ℕ) + 26 * i.1) % 26 = i.2 := by omega
  simp only [OpenMPS.pair, uN, pairIdx, h1, h2, openRight]

/-- `G u` on the doubled index is `guN`. -/
theorem G_pair_openRight : Gm *ᵥ OpenMPS.pair openRight = fun i ↦ guN (pairIdx i) := by
  ext i
  have h1 : ((i.2 : ℕ) + 26 * i.1) / 26 = i.1 := by omega
  have h2 : ((i.2 : ℕ) + 26 * i.1) % 26 = i.2 := by omega
  rw [G_mulVec_pair]
  simp only [Pi.sub_apply, OpenMPS.pair, guN, pairIdx, h1, h2, mulVec_trialA_eq]

/-! ### The block parameters -/

/-- The parameters of every block are well formed. -/
theorem params_good (b : Fin 11) : (params b).Good where
  t := (blockCheck_all b).1
  o := (blockCheck_all b).2.1
  l := List.length_ofFn
  l' := List.length_ofFn
  v := List.length_ofFn
  v' := List.length_ofFn

/-- The copy `Jᵦ` of `L₀ W` is the left vector `ℓ` of block `b`. -/
theorem L0W_copy (b : Fin 11) : copyJ (OpenMPS.pair openLeft ᵥ* Wm) b = (params b).lv :=
  copyJ_row b (L0W_eq b)

/-- The copy `Jᵦ` of `L₀ G W` is the left field vector `ℓ′` of block `b`. -/
theorem L0GW_copy (b : Fin 11) :
    copyJ ((OpenMPS.pair openLeft ᵥ* Gm) ᵥ* Wm) b = (params b).lv' :=
  copyJ_row b (L0GW_eq b)

/-- The copy `Jᵦ` of `V u` is the right vector `ν` of block `b`. -/
theorem Vu_copy (b : Fin 11) : copyJ (Vm *ᵥ OpenMPS.pair openRight) b = (params b).vv := by
  rw [pair_openRight]; exact copyJ_col b uN

/-- The copy `Jᵦ` of `V G u` is the right field vector `ν′` of block `b`. -/
theorem VGu_copy (b : Fin 11) :
    copyJ (Vm *ᵥ (Gm *ᵥ OpenMPS.pair openRight)) b = (params b).vv' := by
  rw [G_pair_openRight]; exact copyJ_col b guN

/-- `L₀ W` vanishes outside the copies `Jᵦ` (`InCopy`, kernel-checked). -/
theorem L0W_zero : ∀ b (k : Fin (bd b)) c, k ≠ jIdx b →
    (OpenMPS.pair openLeft ᵥ* Wm) ⟨b, (k, c)⟩ = 0 :=
  copy_zero (fun b y ↦ L0W_eq b y) fun b ↦ (openBlockCheck_all b).1

/-- `L₀ G W` vanishes outside the copies `Jᵦ` (`InCopy`, kernel-checked). -/
theorem L0GW_zero : ∀ b (k : Fin (bd b)) c, k ≠ jIdx b →
    ((OpenMPS.pair openLeft ᵥ* Gm) ᵥ* Wm) ⟨b, (k, c)⟩ = 0 :=
  copy_zero (fun b y ↦ L0GW_eq b y) fun b ↦ (openBlockCheck_all b).2.1

/-! ### The norm and the energy per block -/

/-- **The norm per block**: `c L₀ T^(k+1) u = ∑ᵦ ℓᵦ t^(k+1) νᵦ`. -/
theorem norm_reduce (k : ℕ) :
    (trialScale : ℤ) * OpenMPS.norm openLeft openRight trialA (k + 1) =
      ∑ b, ((params b).expected k).1 := by
  rw [OpenMPS.norm, reduce (M := fun b ↦ tMat b ^ (k + 1)) (by rw [pow_mul_W, Bm_pow]) _ _
    L0W_zero]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [L0W_copy, Vu_copy]
  exact dotProduct_mulVec _ _ _

/-- The bond energy per block. -/
theorem bond_reduce (k : ℕ) :
    (trialScale : ℤ) * OpenMPS.bondEnergy exchangeZ openLeft openRight trialA k =
      ∑ b, ∑ j ∈ Finset.range k, (params b).lv ⬝ᵥ
        (((params b).T ^ j * (params b).O * (params b).T ^ (k - 1 - j)) *ᵥ (params b).vv) := by
  rw [OpenMPS.bondEnergy, Finset.mul_sum, Finset.sum_comm,
    ← Fin.sum_univ_eq_sum_range (fun j ↦ ∑ b, (params b).lv ⬝ᵥ
      (((params b).T ^ j * (params b).O * (params b).T ^ (k - 1 - j)) *ᵥ (params b).vv)) k]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [reduce (M := fun b ↦ tMat b ^ (j : ℕ) * oMat b * tMat b ^ (k - 1 - j))
    (by rw [bond_mul_W, Bm_bond]) _ _ L0W_zero]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [L0W_copy, Vu_copy]
  rfl

/-- The first-site field per block. -/
theorem left_reduce (k : ℕ) :
    (trialScale : ℤ) * OpenMPS.leftField szDiag openLeft openRight trialA k =
      ∑ b, (params b).lv' ⬝ᵥ ((params b).T ^ k *ᵥ (params b).vv) := by
  rw [OpenMPS.leftField, ← mulVec_mulVec, dotProduct_mulVec,
    reduce (M := fun b ↦ tMat b ^ k) (by rw [pow_mul_W, Bm_pow]) _ _ L0GW_zero]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [L0GW_copy, Vu_copy]
  rfl

/-- The last-site field per block. -/
theorem right_reduce (k : ℕ) :
    (trialScale : ℤ) * OpenMPS.rightField szDiag openLeft openRight trialA k =
      ∑ b, (params b).lv ⬝ᵥ ((params b).T ^ k *ᵥ (params b).vv') := by
  rw [OpenMPS.rightField, ← mulVec_mulVec,
    reduce (M := fun b ↦ tMat b ^ k) (by rw [pow_mul_W, Bm_pow]) _ _ L0W_zero]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [L0W_copy, VGu_copy]
  rfl

/-- **The energy per block**: the shifted energy without the site shift,
`c (500000 ∑ⱼ L₀ Tʲ O T^(k-1-j) u - 300000 L₀ G Tᵏ u + 300000 L₀ Tᵏ G u)`, is the sum of the
block outputs `zₖ ν + 300000 ℓ tᵏ ν′`. -/
theorem energy_reduce (k : ℕ) :
    (trialScale : ℤ) * (500000 * OpenMPS.bondEnergy exchangeZ openLeft openRight trialA k -
      300000 * OpenMPS.leftField szDiag openLeft openRight trialA k +
      300000 * OpenMPS.rightField szDiag openLeft openRight trialA k) =
      ∑ b, (((params b).expected k).2.1 + 300000 * ((params b).expected k).2.2) := by
  rw [mul_add, mul_sub, mul_left_comm, mul_left_comm _ (300000 : ℤ), mul_left_comm _ (300000 : ℤ),
    bond_reduce, left_reduce, right_reduce, Finset.mul_sum, Finset.mul_sum, Finset.mul_sum,
    ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  simp only [Params.expected, Params.zVec, sub_dotProduct, smul_dotProduct, sum_dotProduct,
    smul_eq_mul, dotProduct_mulVec]

/-! ### The stored outputs -/

/-- The stored output of block `b` at the length number `i` (`m = 34, 78, 142, 240`). -/
def outAt (b : Fin 11) (i : ℕ) : ℤ × ℤ × ℤ := (openOutputs.getD b []).getD i (0, 0, 0)

/-- The stored outputs are the block contractions (`Params.run_spec`). -/
theorem outputs_eq (b : Fin 11) :
    outAt b 0 = (params b).expected 33 ∧ outAt b 1 = (params b).expected 77 ∧
      outAt b 2 = (params b).expected 141 ∧ outAt b 3 = (params b).expected 239 := by
  have h := (openBlockCheck_all b).2.2
  have hs := Params.run_spec (params_good b) (by rw [h])
  have hs' : openOutputs.getD b [] =
      [(params b).expected 33, (params b).expected 77, (params b).expected 141,
        (params b).expected 239] := by
    rw [h] at hs
    exact hs
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp only [outAt, hs'] <;> rfl

/-- `c ⟨ψ, ψ⟩` from the stored outputs. -/
def totalV (i : ℕ) : ℤ := ∑ b, (outAt b i).1

/-- `c 500000 ⟨ψ, (H^{opp} + a m) ψ⟩` from the stored outputs. -/
def totalH (i : ℕ) (m : ℤ) : ℤ :=
  ∑ b, ((outAt b i).2.1 + 300000 * (outAt b i).2.2 + 700741 * m * (outAt b i).1)

/-- **The four integer comparisons**, `totalH < 202100 totalV`. -/
theorem totals_lt : totalH 0 34 < 202100 * totalV 0 ∧ totalH 1 78 < 202100 * totalV 1 ∧
    totalH 2 142 < 202100 * totalV 2 ∧ totalH 3 240 < 202100 * totalV 3 := by
  decide +kernel

/-- The stored norms are positive, so each trial vector is nonzero. -/
theorem totals_pos : 0 < totalV 0 ∧ 0 < totalV 1 ∧ 0 < totalV 2 ∧ 0 < totalV 3 := by
  decide +kernel

/-! ### The four trial bounds -/

/-- The trial bound at `m = k + 1` from the stored outputs at the length number `i`. -/
theorem trial_of_outputs (k i : ℕ) (hout : ∀ b, outAt b i = (params b).expected k)
    (hlt : totalH i (k + 1) < 202100 * totalV i) :
    ∃ j, oppositeFieldSpectrum (3 / 5) (k + 1) j <
      -(700741 / 500000 : ℝ) * (k + 1 : ℕ) + 2021 / 5000 := by
  refine exists_oppositeFieldSpectrum_lt_of_openMPS openLeft openRight trialA k ?_
  have hc : (0 : ℤ) < trialScale := by decide
  have hV : (trialScale : ℤ) * OpenMPS.norm openLeft openRight trialA (k + 1) = totalV i := by
    rw [norm_reduce, totalV]
    exact Finset.sum_congr rfl fun b _ ↦ by rw [hout]
  have hE := energy_reduce k
  have hH : (trialScale : ℤ) * (500000 * OpenMPS.bondEnergy exchangeZ openLeft openRight trialA k -
      300000 * OpenMPS.leftField szDiag openLeft openRight trialA k +
      300000 * OpenMPS.rightField szDiag openLeft openRight trialA k +
      700741 * (k + 1) * OpenMPS.norm openLeft openRight trialA (k + 1)) = totalH i (k + 1) := by
    rw [mul_add, hE, mul_left_comm, hV, totalH, totalV, Finset.mul_sum, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun b _ ↦ by rw [hout]
  refine lt_of_mul_lt_mul_left ?_ hc.le
  rw [hH, mul_left_comm, hV]
  exact hlt

/-- **The opposite-field trial bound at `m = 34`**: `H_34^{opp}` has an eigenvalue below
`-34 a + α`. -/
theorem oppositeField_trial_34 :
    ∃ i, oppositeFieldSpectrum (3 / 5) 34 i < -(700741 / 500000 : ℝ) * 34 + 2021 / 5000 := by
  simpa using trial_of_outputs 33 0 (fun b ↦ (outputs_eq b).1) totals_lt.1

/-- **The opposite-field trial bound at `m = 78`**. -/
theorem oppositeField_trial_78 :
    ∃ i, oppositeFieldSpectrum (3 / 5) 78 i < -(700741 / 500000 : ℝ) * 78 + 2021 / 5000 := by
  simpa using trial_of_outputs 77 1 (fun b ↦ (outputs_eq b).2.1) totals_lt.2.1

/-- **The opposite-field trial bound at `m = 142`**. -/
theorem oppositeField_trial_142 :
    ∃ i, oppositeFieldSpectrum (3 / 5) 142 i < -(700741 / 500000 : ℝ) * 142 + 2021 / 5000 := by
  simpa using trial_of_outputs 141 2 (fun b ↦ (outputs_eq b).2.2.1) totals_lt.2.2.1

/-- **The opposite-field trial bound at `m = 240`**. -/
theorem oppositeField_trial_240 :
    ∃ i, oppositeFieldSpectrum (3 / 5) 240 i < -(700741 / 500000 : ℝ) * 240 + 2021 / 5000 := by
  simpa using trial_of_outputs 239 3 (fun b ↦ (outputs_eq b).2.2.2) totals_lt.2.2.2

end SeedKernel.OpenTrial
