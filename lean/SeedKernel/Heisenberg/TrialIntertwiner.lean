/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.TrialBound
import SeedKernel.Heisenberg.TrialCompute
import SeedKernel.Heisenberg.TrialReduction

/-!
# The intertwiner of the trial-state transfer matrix

The doubled virtual space `Fin 26 × Fin 26` (index `p = 26 a + b`) carries the transfer matrix
`T = mpsTransfer trialA` and the insertion `O = mpsInsertion trialA`. Both commute with the
virtual `SU(2)` action and with the swap of the two copies, so an integer change of basis `W`
(columns: lowered highest-weight vectors, one group per spin `J`, swap sign and weight) brings
them to the block matrices `B = copies d m t`, `C = copies d m o`: `dᵦ = 2Jᵦ + 1` copies of the
block `tᵦ` (resp. `oᵦ`). `V` is `c W⁻¹`.

This file defines `W`, `V`, `B`, `C` from the generated data and proves the identities
`T W = W B`, `O W = W C`, `V W = c 1` from the list comparisons of `TrialCompute`, which
`TrialChecks` runs in the kernel. Each identity is evaluated on one radix vector over the
columns (`radG`, digits of width `w` larger than the a priori entry bounds of `EntryBound`), so
one comparison of integers per packed row decides every entry (`radix_eq_zero`).

* `iotaEquiv`: the column index `ι = Σ b, Fin dᵦ × Fin mᵦ` numbered `0, …, 675`.
* `mulVec_Wm`, `mulVec_Vm`, `mulVec_copies`: sparse evaluation of `W`, `V`, `B` on vectors.
* `mulVec_T`, `mulVec_O`: `T` and `O` on a vectorized `26 × 26` matrix `M` are
  `M ↦ ∑ Aₛ M Aₛ'ᵀ` (Kronecker products).
* `rowPack_T`, `rowPack_O`: their packed rows are the lists `TPack`, `OPack`.
-/

namespace SeedKernel.Trial

open Matrix CertifiedGaps.SpinOneHeisenberg
open scoped Kronecker

/-! ### Index sets -/

/-- The doubled virtual index `(a, b)`. -/
abbrev Pair := Fin 26 × Fin 26

/-- `p = 26 a + b`. -/
def pairIdx (i : Pair) : ℕ := (i.2 : ℕ) + 26 * i.1

/-- Columns of `W`: block, weight (lowering power `k`), copy `c`. -/
abbrev Iota := Σ b : Fin 11, Fin (bd b) × Fin (bm b)

/-- Local column `c + mᵦ k` inside block `b`. -/
def localIdx {b : Fin 11} (y : Fin (bd b) × Fin (bm b)) : ℕ := (y.2 : ℕ) + bm b * y.1

theorem sum_bd_bm : ∑ b : Fin 11, bd b * bm b = 676 := by decide

/-- The columns numbered `0, …, 675`. -/
def iotaEquiv : Iota ≃ Fin 676 :=
  (Equiv.sigmaCongrRight fun _ ↦ finProdFinEquiv).trans
    (finSigmaFinEquiv.trans (finCongr sum_bd_bm))

theorem bStart_eq : ∀ b : Fin 11,
    ∑ i : Fin b, bd (Fin.castLE b.2.le i) * bm (Fin.castLE b.2.le i) = bStart b := by
  decide

theorem iotaEquiv_val (x : Iota) : (iotaEquiv x : ℕ) = bStart x.1 + localIdx x.2 := by
  rcases x with ⟨b, y⟩
  simp only [iotaEquiv, Equiv.trans_apply, Equiv.sigmaCongrRight_apply, finCongr_apply,
    Fin.val_cast, finSigmaFinEquiv_apply]
  rw [bStart_eq]
  rfl

theorem card_Pair : Fintype.card Pair = 676 := by simp

theorem card_Iota : Fintype.card Iota = Fintype.card Pair := by
  rw [Fintype.card_congr iotaEquiv]; simp

theorem card_Iota' : Fintype.card Iota = 676 := by rw [card_Iota, card_Pair]

/-! ### The matrices -/

/-- Block `tᵦ`. -/
def tMat (b : Fin 11) : Matrix (Fin (bm b)) (Fin (bm b)) ℤ := (tB b).mat (bm b)
/-- Block `oᵦ`. -/
def oMat (b : Fin 11) : Matrix (Fin (bm b)) (Fin (bm b)) ℤ := (oB b).mat (bm b)

/-- `B`: `dᵦ` copies of each `tᵦ`. -/
def Bm : Matrix Iota Iota ℤ := copies bd bm tMat
/-- `C`: `dᵦ` copies of each `oᵦ`. -/
def Cm : Matrix Iota Iota ℤ := copies bd bm oMat

/-- The intertwiner `W`. -/
def Wm : Matrix Pair Iota ℤ := fun i x ↦ sparseEntry ((wG x.1).getD (pairIdx i) []) (localIdx x.2)
/-- `V = c W⁻¹`. -/
def Vm : Matrix Iota Pair ℤ := fun x i ↦ sparseEntry ((vG x.1).getD (localIdx x.2) []) (pairIdx i)

/-! ### Sparse evaluation -/

theorem mem_getD_or_nil {α : Type*} {l : List (List α)} {n : ℕ} :
    l.getD n [] ∈ l ∨ l.getD n [] = [] := by
  rcases Nat.lt_or_ge n l.length with h | h
  · left; rw [List.getD_eq_getElem _ _ h]; exact List.getElem_mem h
  · right; simp [h]

theorem sum_prod_localIdx {b : Fin 11} (f : ℕ → ℤ) :
    ∑ y : Fin (bd b) × Fin (bm b), f (localIdx y) = ∑ j : Fin (bd b * bm b), f j :=
  (Equiv.sum_comp finProdFinEquiv (fun j : Fin (bd b * bm b) ↦ f j))

theorem sum_pair_pairIdx (f : ℕ → ℤ) : ∑ i : Pair, f (pairIdx i) = ∑ j : Fin 676, f j :=
  (Equiv.sum_comp finProdFinEquiv (fun j : Fin (26 * 26) ↦ f j))

theorem mulVec_Wm (hW : WIdxOk) (Z : Fin 11 → ℕ → ℤ) (i : Pair) :
    (Wm *ᵥ fun x ↦ Z x.1 (localIdx x.2)) i =
      ∑ b : Fin 11, sparseDot ((wG b).getD (pairIdx i) []) (Z b) := by
  simp only [mulVec, dotProduct, Wm, Fintype.sum_sigma]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [sum_prod_localIdx (fun j ↦ sparseEntry ((wG b).getD (pairIdx i) []) j * Z b j)]
  apply sum_sparseEntry_mul
  rcases mem_getD_or_nil (l := wG b) (n := pairIdx i) with h | h
  · exact hW b _ h
  · rw [h]; simp

theorem mulVec_Vm (hV : VIdxOk) (Y : ℕ → ℤ) (x : Iota) :
    (Vm *ᵥ fun i ↦ Y (pairIdx i)) x = sparseDot ((vG x.1).getD (localIdx x.2) []) Y := by
  simp only [mulVec, dotProduct, Vm]
  rw [sum_pair_pairIdx (fun j ↦ sparseEntry ((vG x.1).getD (localIdx x.2) []) j * Y j)]
  apply sum_sparseEntry_mul
  rcases mem_getD_or_nil (l := vG x.1) (n := localIdx x.2) with h | h
  · exact hV x.1 _ h
  · rw [h]; simp

theorem mulVec_copies {d m : Fin 11 → ℕ} (t : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ)
    (v : (Σ b, Fin (d b) × Fin (m b)) → ℤ) (b : Fin 11) (k : Fin (d b)) (c : Fin (m b)) :
    (copies d m t *ᵥ v) ⟨b, (k, c)⟩ = ∑ c' : Fin (m b), t b c c' * v ⟨b, (k, c')⟩ := by
  simp only [mulVec, dotProduct, copies, Fintype.sum_sigma, blockDiagonal'_apply]
  rw [Finset.sum_eq_single b]
  · simp only [↓reduceDIte, cast_eq, kroneckerMap_apply, Fintype.sum_prod_type, one_apply]
    rw [Finset.sum_eq_single k]
    · simp
    · intro k' _ hk; simp [Ne.symm hk]
    · simp
  · intro b' _ hb; simp [Ne.symm hb]
  · simp

/-! ### `T` and `O` on vectorized matrices -/

/-- `O = Y_z² + Y_u Y_d + Y_d Y_u` with `Y_z = A₂⊗A₂ - A₀⊗A₀`, `Y_u = A₁⊗A₀ + A₂⊗A₁`,
`Y_d = A₀⊗A₁ + A₁⊗A₂`: the exchange `SᶻSᶻ + (S⁺S⁻ + S⁻S⁺)/2` written as products of transfer
matrices with one spin operator inserted. -/
theorem mpsInsertion_eq_sum {D : Type*} [Fintype D] (A : Fin 3 → Matrix D D ℤ) :
    mpsInsertion A =
      (A 2 ⊗ₖ A 2 - A 0 ⊗ₖ A 0) * (A 2 ⊗ₖ A 2 - A 0 ⊗ₖ A 0) +
        (A 1 ⊗ₖ A 0 + A 2 ⊗ₖ A 1) * (A 0 ⊗ₖ A 1 + A 1 ⊗ₖ A 2) +
        (A 0 ⊗ₖ A 1 + A 1 ⊗ₖ A 2) * (A 1 ⊗ₖ A 0 + A 2 ⊗ₖ A 1) := by
  simp only [mpsInsertion, Fintype.sum_prod_type, Fin.sum_univ_three]
  simp only [exchangeZ, mz, Prod.mk.injEq, Fin.isValue]
  norm_num
  simp only [add_mul, mul_add, sub_mul, mul_sub, ← mul_kronecker_mul]
  abel

theorem mpsTransfer_eq_sum {D : Type*} (A : Fin 3 → Matrix D D ℤ) :
    mpsTransfer A = A 0 ⊗ₖ A 0 + A 1 ⊗ₖ A 1 + A 2 ⊗ₖ A 2 := by
  simp [mpsTransfer, Fin.sum_univ_three]

/-- `M ↦ A s M (A s')ᵀ`. -/
def sand (s s' : Fin 3) (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ :=
  trialA s * M * (trialA s')ᵀ

/-- `T` on vectorized matrices. -/
def PhiT (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ :=
  sand 0 0 M + sand 1 1 M + sand 2 2 M
/-- `Y_z` on vectorized matrices. -/
def Phiz (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ := sand 2 2 M - sand 0 0 M
/-- `Y_u` on vectorized matrices. -/
def Phiu (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ := sand 1 0 M + sand 2 1 M
/-- `Y_d` on vectorized matrices. -/
def Phid (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ := sand 0 1 M + sand 1 2 M
/-- `O` on vectorized matrices. -/
def PhiO (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ :=
  Phiz (Phiz M) + Phiu (Phid M) + Phid (Phiu M)

theorem kron_mulVec_sand (s s' : Fin 3) (M : Matrix (Fin 26) (Fin 26) ℤ) :
    (trialA s ⊗ₖ trialA s') *ᵥ vecM M = vecM (sand s s' M) := kronecker_mulVec _ _ _

theorem mulVec_T (M : Matrix (Fin 26) (Fin 26) ℤ) :
    mpsTransfer trialA *ᵥ vecM M = vecM (PhiT M) := by
  rw [mpsTransfer_eq_sum, add_mulVec, add_mulVec, kron_mulVec_sand, kron_mulVec_sand,
    kron_mulVec_sand, PhiT, vecM_add, vecM_add]

theorem mulVec_O (M : Matrix (Fin 26) (Fin 26) ℤ) :
    mpsInsertion trialA *ᵥ vecM M = vecM (PhiO M) := by
  have hz (N : Matrix (Fin 26) (Fin 26) ℤ) :
      (trialA 2 ⊗ₖ trialA 2 - trialA 0 ⊗ₖ trialA 0) *ᵥ vecM N = vecM (Phiz N) := by
    rw [sub_mulVec, kron_mulVec_sand, kron_mulVec_sand, Phiz, vecM_sub]
  have hu (N : Matrix (Fin 26) (Fin 26) ℤ) :
      (trialA 1 ⊗ₖ trialA 0 + trialA 2 ⊗ₖ trialA 1) *ᵥ vecM N = vecM (Phiu N) := by
    rw [add_mulVec, kron_mulVec_sand, kron_mulVec_sand, Phiu, vecM_add]
  have hd (N : Matrix (Fin 26) (Fin 26) ℤ) :
      (trialA 0 ⊗ₖ trialA 1 + trialA 1 ⊗ₖ trialA 2) *ᵥ vecM N = vecM (Phid N) := by
    rw [add_mulVec, kron_mulVec_sand, kron_mulVec_sand, Phid, vecM_add]
  rw [mpsInsertion_eq_sum, add_mulVec, add_mulVec, ← mulVec_mulVec, ← mulVec_mulVec,
    ← mulVec_mulVec, hz, hz, hd, hu, hu, hd, PhiO, vecM_add, vecM_add]

/-- The term `(A s A s') M (A t A t')ᵀ` of `O`. -/
def oterm (s s' t t' : Fin 3) (M : Matrix (Fin 26) (Fin 26) ℤ) : Matrix (Fin 26) (Fin 26) ℤ :=
  trialA s * trialA s' * M * (trialA t * trialA t')ᵀ

/-- `O` as twelve terms. -/
theorem PhiO_expand (M : Matrix (Fin 26) (Fin 26) ℤ) :
    PhiO M = oterm 2 2 2 2 M - oterm 2 0 2 0 M - oterm 0 2 0 2 M + oterm 0 0 0 0 M +
      oterm 1 0 0 1 M + oterm 1 1 0 2 M + oterm 2 0 1 1 M + oterm 2 1 1 2 M +
      oterm 0 1 1 0 M + oterm 0 2 1 1 M + oterm 1 1 2 0 M + oterm 1 2 2 1 M := by
  simp only [PhiO, Phiz, Phiu, Phid, sand, oterm, Matrix.mul_add, Matrix.mul_sub, Matrix.add_mul,
    Matrix.sub_mul, transpose_mul, Matrix.mul_assoc]
  abel

/-! ### Radix vectors restricted to a range of blocks -/

/-- `R^(position of x - start of block lo)` on the columns of the blocks `lo ≤ b < hi`, zero on
the other columns (`R = 2^w`). -/
def radG (w lo hi : ℕ) (x : Iota) : ℤ := maskW w lo hi x.1 (localIdx x.2)

theorem length_wG (b : Fin 11) : (wG b).length = 676 := length_decodeSparse _ _

theorem getD_flatW (F : Fin 11 → ℕ → ℤ) (p : ℕ) (hp : p < 676) :
    (flatW F).getD p 0 = ∑ b : Fin 11, sparseDot ((wG b).getD p []) (F b) := by
  rw [flatW, getD_vsumN 676 _ _ p hp, List.map_map, Fin.sum_univ_def]
  · congr 1
    apply List.map_congr_left
    intro b _
    simp only [Function.comp, blockVec]
    rw [show (0 : ℤ) = sparseDot [] (F b) from rfl, List.getD_map]
  · intro l hl
    simp only [List.mem_map] at hl
    obtain ⟨b, _, rfl⟩ := hl
    simp [blockVec, length_wG]

theorem length_flatW (F : Fin 11 → ℕ → ℤ) : (flatW F).length = 676 := by
  apply length_vsumN
  intro l hl
  simp only [List.mem_map] at hl
  obtain ⟨b, _, rfl⟩ := hl
  simp [blockVec, length_wG]

/-- Abbreviation for `26 × 26` list matrices as Mathlib matrices. -/
abbrev mat26 (Y : List (List ℤ)) : Matrix (Fin 26) (Fin 26) ℤ := toMat 26 26 Y

theorem sq_flat (F : Fin 11 → ℕ → ℤ) : Sq 26 (chunksZ 26 26 (flatW F)) :=
  sq_chunksZ 26 26 _ (by rw [length_flatW])

theorem sq_YG (w lo hi : ℕ) : Sq 26 (YG w lo hi) := sq_flat _

theorem sq_ZG (w lo hi : ℕ) (X : Fin 11 → OM) : Sq 26 (ZG w lo hi X) := sq_flat _

theorem vecM_flat (F : Fin 11 → ℕ → ℤ) :
    vecM (mat26 (chunksZ 26 26 (flatW F))) =
      fun i : Pair ↦ ∑ b : Fin 11, sparseDot ((wG b).getD (pairIdx i) []) (F b) := by
  ext i
  simp only [vecM, mat26, toMat_apply]
  rw [lget_chunksZ 26 26 _ _ _ i.1.isLt i.2.isLt, getD_flatW _ _ (by omega), pairIdx,
    Nat.add_comm]

theorem mulVec_Wm_radG (hW : WIdxOk) (w lo hi : ℕ) :
    Wm *ᵥ radG w lo hi = vecM (mat26 (YG w lo hi)) := by
  rw [YG, vecM_flat]
  ext i
  exact mulVec_Wm hW (maskW w lo hi) i

theorem packN_eq_sum (w : ℕ) : ∀ (m : ℕ) (l : List ℕ), l.length = m →
    (packN w l : ℤ) = ∑ j : Fin m, (l.getD j 0 : ℤ) * ((2 ^ (w * j) : ℕ) : ℤ)
  | 0, [], _ => by simp [packN]
  | m + 1, x :: l, h => by
    rw [packN_cons, Fin.sum_univ_succ]
    push_cast
    rw [packN_eq_sum w m l (by simpa using h)]
    simp only [Fin.val_zero, Nat.mul_zero, pow_zero, List.getD_cons_zero, Fin.val_succ,
      List.getD_cons_succ, Nat.mul_succ, pow_add]
    push_cast
    rw [Finset.mul_sum]
    congr 1
    · ring
    · refine Finset.sum_congr rfl fun j _ ↦ ?_; ring

theorem ubT_eq {m : ℕ} {X : OM} (hX : X.Good m) (w : ℕ) (c : Fin m) :
    ubT w X m c = ∑ c' : Fin m, X.mat m c c' * ((2 ^ (w * c') : ℕ) : ℤ) := by
  have hc : (c : ℕ) < X.2.length := by rw [hX.1]; exact c.isLt
  have hrow : (X.2.getD c []).length = m := by
    rw [List.getD_eq_getElem _ _ hc]; exact (hX.2.2 _ (List.getElem_mem hc)).1
  rw [ubT, packN_eq_sum w m _ hrow, packN_eq_sum w m _ (by simp), Finset.mul_sum,
    ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun c' _ ↦ ?_
  have hc' : (c' : ℕ) < (X.2.getD c []).length := by rw [hrow]; exact c'.isLt
  have hrep : (List.replicate m 1).getD c' 0 = 1 := by
    rw [List.getD_eq_getElem _ _ (by simp)]; simp
  have hc'' : (c' : ℕ) < X.2[(c : ℕ)].length := by rw [← List.getD_eq_getElem _ _ hc]; exact hc'
  simp only [OM.mat, toMatB, toMat_apply]
  rw [lget_shiftL _ _ _ _ hc hc'', sh, hrep, List.getD_eq_getElem _ _ hc,
    List.getD_eq_getElem _ _ hc'']
  push_cast; ring

theorem mulVec_copies_radG {X : Fin 11 → OM} (hX : ∀ b, (X b).Good (bm b)) (w lo hi : ℕ) :
    copies bd bm (fun b ↦ (X b).mat (bm b)) *ᵥ radG w lo hi =
      fun x ↦ zBG w lo hi X x.1 (localIdx x.2) := by
  ext ⟨b, k, c⟩
  rw [mulVec_copies, zBG, localIdx]
  simp only [radG, maskW, localIdx]
  split_ifs with hb
  · have hm : 0 < bm b := Nat.lt_of_le_of_lt (Nat.zero_le _) c.isLt
    have h1 : ((c : ℕ) + bm b * k) / bm b = k := by
      rw [Nat.add_mul_div_left _ _ hm, Nat.div_eq_of_lt c.isLt, zero_add]
    have h2 : ((c : ℕ) + bm b * k) % bm b = c := by
      rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt c.isLt]
    simp only [h1, h2]
    rw [ubT_eq (hX b), Finset.mul_sum]
    refine Finset.sum_congr rfl fun c' _ ↦ ?_
    rw [show w * (bStart b - gStart lo + ((c' : ℕ) + bm b * k)) =
      w * (bStart b - gStart lo + bm b * k) + w * c' by ring, pow_add]
    push_cast; ring
  · simp

theorem mulVec_Wm_copies_radG (hW : WIdxOk) {X : Fin 11 → OM} (hX : ∀ b, (X b).Good (bm b))
    (w lo hi : ℕ) : Wm *ᵥ (copies bd bm (fun b ↦ (X b).mat (bm b)) *ᵥ radG w lo hi) =
      vecM (mat26 (ZG w lo hi X)) := by
  rw [mulVec_copies_radG hX, ZG, vecM_flat]
  ext i
  exact mulVec_Wm hW (zBG w lo hi X) i

/-! ### Packed rows -/

/-- The packed row `∑_b f(a, b) S^b`. -/
def rowPack (w : ℕ) (f : Pair → ℤ) (a : Fin 26) : ℤ := ∑ b : Fin 26, f (a, b) * powS w b

theorem rowPack_add (w : ℕ) (f g : Pair → ℤ) (a : Fin 26) :
    rowPack w (f + g) a = rowPack w f a + rowPack w g a := by
  simp [rowPack, add_mul, Finset.sum_add_distrib]

theorem rowPack_sub (w : ℕ) (f g : Pair → ℤ) (a : Fin 26) :
    rowPack w (f - g) a = rowPack w f a - rowPack w g a := by
  simp [rowPack, sub_mul, Finset.sum_sub_distrib]

/-- Packed row of a sandwich `X M Yᵀ`: `∑_c X(a, c) ∑_d M(c, d) ∑_b Y(b, d) S^b`. -/
theorem rowPack_sandwich (w : ℕ) (X M Y : Matrix (Fin 26) (Fin 26) ℤ) (a : Fin 26) :
    rowPack w (vecM (X * M * Yᵀ)) a =
      ∑ c : Fin 26, X a c * ∑ d : Fin 26, M c d * ∑ b : Fin 26, Y b d * powS w b := by
  simp only [rowPack, vecM, mul_apply, transpose_apply, Finset.sum_mul, Finset.mul_sum]
  calc ∑ b : Fin 26, ∑ d : Fin 26, ∑ c : Fin 26, X a c * M c d * Y b d * powS w b
      = ∑ b : Fin 26, ∑ c : Fin 26, ∑ d : Fin 26, X a c * M c d * Y b d * powS w b :=
        Finset.sum_congr rfl fun _ _ ↦ Finset.sum_comm
    _ = ∑ c : Fin 26, ∑ b : Fin 26, ∑ d : Fin 26, X a c * M c d * Y b d * powS w b :=
        Finset.sum_comm
    _ = ∑ c : Fin 26, ∑ d : Fin 26, ∑ b : Fin 26, X a c * M c d * Y b d * powS w b :=
        Finset.sum_congr rfl fun _ _ ↦ Finset.sum_comm
    _ = _ := by
        refine Finset.sum_congr rfl fun c _ ↦ Finset.sum_congr rfl fun d _ ↦
          Finset.sum_congr rfl fun b _ ↦ ?_
        ring

theorem getD_ofFn26 (f : Fin 26 → ℤ) (i : Fin 26) : (List.ofFn f).getD i 0 = f i := by
  rw [List.getD_eq_getElem _ _ (by simp), List.getElem_ofFn]

theorem sparseDot_trialSparse (s : Fin 3) (a : Fin 26) (F : ℕ → ℤ) :
    sparseDot ((trialSparse s).getD a []) F = ∑ c : Fin 26, trialA s a c * F c := by
  rw [trialSparse, List.getD_eq_getElem _ _ (by simp), List.getElem_ofFn, sparseDot_filter,
    sparseDot_ofFn]
  rfl

theorem sparseDot_trialSparseT (s : Fin 3) (d : Fin 26) (F : ℕ → ℤ) :
    sparseDot ((trialSparseT s).getD d []) F = ∑ e : Fin 26, trialA s e d * F e := by
  rw [trialSparseT, List.getD_eq_getElem _ _ (by simp), List.getElem_ofFn, sparseDot_filter,
    sparseDot_ofFn' (fun b : Fin 26 ↦ trialEntry s b d)]
  rfl

theorem length_gamL (s : Fin 3) (g : List ℤ) : (gamL s g).length = 26 := by simp [gamL]
theorem length_lcombL (s : Fin 3) (v : List ℤ) : (lcombL s v).length = 26 := by simp [lcombL]
theorem length_sigS (w : ℕ) : (sigS w).length = 26 := by simp [sigS]

theorem gamL_getD (s : Fin 3) (g : List ℤ) (d : Fin 26) :
    (gamL s g).getD d 0 = ∑ e : Fin 26, trialA s e d * g.getD e 0 := by
  rw [gamL, getD_ofFn26, sparseDot_trialSparseT]

theorem lcombL_getD (s : Fin 3) (v : List ℤ) (a : Fin 26) :
    (lcombL s v).getD a 0 = ∑ c : Fin 26, trialA s a c * v.getD c 0 := by
  rw [lcombL, getD_ofFn26, sparseDot_trialSparse]

theorem sigS_getD (w : ℕ) (b : Fin 26) : (sigS w).getD b 0 = powS w b := by
  rw [sigS, getD_ofFn26]

theorem prowL_getD {Y : List (List ℤ)} (hY : Sq 26 Y) {g : List ℤ} (hg : g.length = 26)
    (c : Fin 26) : (prowL Y g).getD c 0 = ∑ d : Fin 26, mat26 Y c d * g.getD d 0 := by
  have hc : (c : ℕ) < Y.length := by rw [hY.1]; exact c.isLt
  have hrow := hY.2 _ (List.getElem_mem hc)
  rw [prowL, List.getD_eq_getElem _ _ (by simpa using hc), List.getElem_map,
    ldot_eq_sum 26 _ _ hrow hg]
  refine Finset.sum_congr rfl fun d _ ↦ ?_
  simp only [mat26, toMat_apply, lget, List.getD_eq_getElem _ _ hc]

theorem gam_sigS (s : Fin 3) (w : ℕ) (d : Fin 26) :
    (gamL s (sigS w)).getD d 0 = ∑ b : Fin 26, trialA s b d * powS w b := by
  rw [gamL_getD]
  simp only [sigS_getD]

theorem gam_gam_sigS (t t' : Fin 3) (w : ℕ) (d : Fin 26) :
    (gamL t' (gamL t (sigS w))).getD d 0 =
      ∑ b : Fin 26, (trialA t * trialA t') b d * powS w b := by
  rw [gamL_getD]
  simp only [gam_sigS, mul_apply, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun b _ ↦ Finset.sum_congr rfl fun e _ ↦ ?_
  ring

theorem lcomb_lcomb (s s' : Fin 3) (v : List ℤ) (a : Fin 26) :
    (lcombL s (lcombL s' v)).getD a 0 =
      ∑ c : Fin 26, (trialA s * trialA s') a c * v.getD c 0 := by
  rw [lcombL_getD]
  simp only [lcombL_getD, mul_apply, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun c _ ↦ Finset.sum_congr rfl fun e _ ↦ ?_
  ring

theorem rowPack_sand {Y : List (List ℤ)} (hY : Sq 26 Y) (w : ℕ) (s : Fin 3) (a : Fin 26) :
    rowPack w (vecM (sand s s (mat26 Y))) a =
      (lcombL s (prowL Y (gamL s (sigS w)))).getD a 0 := by
  rw [sand, rowPack_sandwich, lcombL_getD]
  refine Finset.sum_congr rfl fun c _ ↦ ?_
  rw [prowL_getD hY (length_gamL _ _)]
  simp only [gam_sigS]

/-- **Packed rows of `T`.** -/
theorem rowPack_T {Y : List (List ℤ)} (hY : Sq 26 Y) (w : ℕ) (a : Fin 26) :
    rowPack w (mpsTransfer trialA *ᵥ vecM (mat26 Y)) a = (TPack w Y).getD a 0 := by
  rw [mulVec_T, PhiT, vecM_add, vecM_add, rowPack_add, rowPack_add, rowPack_sand hY,
    rowPack_sand hY, rowPack_sand hY, TPack,
    getD_vadd _ _ (by rw [length_vadd, length_lcombL, length_lcombL]; rfl) (length_lcombL _ _)
      a a.isLt,
    getD_vadd _ _ (length_lcombL _ _) (length_lcombL _ _) a a.isLt]

theorem rowPack_oterm {Y : List (List ℤ)} (hY : Sq 26 Y) (w : ℕ) (s s' t t' : Fin 3)
    (a : Fin 26) :
    rowPack w (vecM (oterm s s' t t' (mat26 Y))) a = (OTerm w Y s s' t t').getD a 0 := by
  rw [oterm, rowPack_sandwich, OTerm, lcomb_lcomb]
  refine Finset.sum_congr rfl fun c _ ↦ ?_
  rw [prowL_getD hY (length_gamL _ _)]
  simp only [gam_gam_sigS]

/-- **Packed rows of `O`.** -/
theorem rowPack_O {Y : List (List ℤ)} (hY : Sq 26 Y) (w : ℕ) (a : Fin 26) :
    rowPack w (mpsInsertion trialA *ᵥ vecM (mat26 Y)) a = (OPack w Y).getD a 0 := by
  have hl (s s' t t' : Fin 3) : (OTerm w Y s s' t t').length = 26 := length_lcombL _ _
  have hn (s s' t t' : Fin 3) : (vneg (OTerm w Y s s' t t')).length = 26 := by
    rw [length_vneg, hl]
  have hv {u v : List ℤ} (hu : u.length = 26) (hv : v.length = 26) :
      (vadd u v).length = 26 := by
    rw [length_vadd, hu, hv]; rfl
  have ga {u v : List ℤ} (hu : u.length = 26) (hv : v.length = 26) :
      (vadd u v).getD a 0 = u.getD a 0 + v.getD a 0 := getD_vadd u v hu hv a a.isLt
  rw [mulVec_O, PhiO_expand]
  simp only [vecM_add, vecM_sub, rowPack_add, rowPack_sub, rowPack_oterm hY]
  simp only [OPack]
  rw [ga (hv (hv (hv (hv (hv (hv (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hv (hv (hv (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hv (hv (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _))
      (hl _ _ _ _),
    ga (hv (hv (hv (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _))
      (hl _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _)) (hl _ _ _ _)) (hl _ _ _ _),
    ga (hv (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _)) (hl _ _ _ _),
    ga (hv (hl _ _ _ _) (hn _ _ _ _)) (hn _ _ _ _),
    ga (hl _ _ _ _) (hn _ _ _ _), getD_vneg, getD_vneg]
  ring

/-- Packed rows of a `26 × 26` list matrix. -/
theorem rowPack_vecM {Z : List (List ℤ)} (hZ : Sq 26 Z) (w : ℕ) (a : Fin 26) :
    rowPack w (vecM (mat26 Z)) a = (ZPack w Z).getD a 0 := by
  have ha : (a : ℕ) < Z.length := by rw [hZ.1]; exact a.isLt
  have hrow := hZ.2 _ (List.getElem_mem ha)
  rw [ZPack, List.getD_eq_getElem _ _ (by simpa using ha), List.getElem_map,
    ldot_eq_sum 26 _ _ hrow (length_sigS w), rowPack]
  refine Finset.sum_congr rfl fun b _ ↦ ?_
  rw [sigS_getD]
  simp only [vecM, mat26, toMat_apply, lget, List.getD_eq_getElem _ _ ha]

/-- Radix uniqueness with an injective choice of exponents. -/
theorem radix_eq_zero_inj {ι : Type*} [Fintype ι] (R : ℤ) (hR : 0 < R) (e : ι → ℕ)
    (he : Function.Injective e) (d : ι → ℤ) (hd : ∀ i, |d i| < R)
    (h : ∑ i, d i * R ^ e i = 0) (i : ι) : d i = 0 := by
  classical
  set N := Finset.univ.sup e + 1
  have hlt (j : ι) : e j < N := Nat.lt_succ_of_le (Finset.le_sup (Finset.mem_univ j))
  set g : ι → Fin N := fun j ↦ ⟨e j, hlt j⟩ with hg
  set d' : Fin N → ℤ := fun k ↦ ∑ j ∈ Finset.univ.filter (fun j ↦ g j = k), d j with hd'
  have hfib (j : ι) : d' (g j) = d j := by
    have : Finset.univ.filter (fun j' ↦ g j' = g j) = {j} := by
      ext j'
      simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton, hg,
        Fin.mk.injEq]
      exact ⟨fun hj ↦ he hj, fun hj ↦ hj ▸ rfl⟩
    simp only [hd', this, Finset.sum_singleton]
  have hbound (k : Fin N) : |d' k| < R := by
    by_cases hk : ∃ j, g j = k
    · obtain ⟨j, rfl⟩ := hk
      rw [hfib]
      exact hd j
    · push Not at hk
      have : Finset.univ.filter (fun j ↦ g j = k) = ∅ := by
        ext j
        simp [hk j]
      simp only [hd', this, Finset.sum_empty, abs_zero]
      exact hR
  have hsum : ∑ k : Fin N, d' k * R ^ (k : ℕ) = 0 := by
    rw [← h]
    simp only [hd', Finset.sum_mul]
    rw [← Finset.sum_fiberwise_of_maps_to (s := Finset.univ) (t := Finset.univ) (g := g)
      (fun j _ ↦ Finset.mem_univ (g j)) (fun j ↦ d j * R ^ e j)]
    refine Finset.sum_congr rfl fun k _ ↦ Finset.sum_congr rfl fun j hj ↦ ?_
    rw [(Finset.mem_filter.1 hj).2.symm]
  have := radix_eq_zero R hR N d' hbound hsum (g i)
  rwa [hfib] at this

theorem bStart_add_le : ∀ hi : ℕ, hi ≤ 11 → ∀ b : Fin 11, (b : ℕ) < hi →
    bStart b + bd b * bm b ≤ gStart hi := by
  decide

theorem gStart_le : ∀ lo : ℕ, lo ≤ 11 → ∀ b : Fin 11, lo ≤ (b : ℕ) → gStart lo ≤ bStart b := by
  decide

theorem localIdx_lt {b : Fin 11} (y : Fin (bd b) × Fin (bm b)) : localIdx y < bd b * bm b := by
  have := (finProdFinEquiv y).isLt
  simpa [localIdx] using this

/-- **Radix uniqueness on packed rows.** If every packed row of `D` on the radix vector of the
blocks `lo ≤ b < hi` vanishes (slots of `L = start hi - start lo` columns), the columns of those
blocks of `D` vanish. -/
theorem eq_zero_of_rowPack (D : Matrix Pair Iota ℤ) (w lo hi : ℕ) (K : ℤ)
    (hD : EntryBound D K) (hK : K < 2 ^ w) (hlo : lo ≤ 11) (hhi : hi ≤ 11)
    (h : ∀ a, rowPack ((gStart hi - gStart lo) * w) (D *ᵥ radG w lo hi) a = 0) (i : Pair)
    (x : Iota) (hx : lo ≤ (x.1 : ℕ) ∧ (x.1 : ℕ) < hi) : D i x = 0 := by
  classical
  obtain ⟨a, b⟩ := i
  set L := gStart hi - gStart lo with hL
  set R : ℤ := 2 ^ w
  have hR : 0 < R := by positivity
  let P : Fin 26 × Iota → Prop := fun y ↦ lo ≤ (y.2.1 : ℕ) ∧ (y.2.1 : ℕ) < hi
  -- Position of a column of the group, counted from the first column of block `lo`.
  have hpos (y : Fin 26 × Iota) (hy : P y) :
      gStart lo ≤ bStart y.2.1 ∧ bStart y.2.1 + localIdx y.2.2 < gStart hi := by
    refine ⟨gStart_le lo hlo y.2.1 hy.1, ?_⟩
    have := bStart_add_le hi hhi y.2.1 hy.2
    have := localIdx_lt y.2.2
    omega
  let e : {y // P y} → ℕ := fun y ↦ (bStart y.1.2.1 - gStart lo + localIdx y.1.2.2) + L * y.1.1
  have he : Function.Injective e := by
    rintro ⟨⟨b1, x1⟩, h1⟩ ⟨⟨b2, x2⟩, h2⟩ heq
    obtain ⟨p1a, p1b⟩ : gStart lo ≤ bStart x1.1 ∧ bStart x1.1 + localIdx x1.2 < gStart hi :=
      hpos _ h1
    obtain ⟨p2a, p2b⟩ : gStart lo ≤ bStart x2.1 ∧ bStart x2.1 + localIdx x2.2 < gStart hi :=
      hpos _ h2
    simp only [e] at heq
    have hm1 : bStart x1.1 - gStart lo + localIdx x1.2 < L := by omega
    have hm2 : bStart x2.1 - gStart lo + localIdx x2.2 < L := by omega
    have hb : (b1 : ℕ) = b2 := by
      have := congrArg (· / L) heq
      rwa [Nat.add_mul_div_left _ _ (by omega), Nat.add_mul_div_left _ _ (by omega),
        Nat.div_eq_of_lt hm1, Nat.div_eq_of_lt hm2, Nat.zero_add, Nat.zero_add] at this
    have hp : bStart x1.1 + localIdx x1.2 = bStart x2.1 + localIdx x2.2 := by
      rw [hb] at heq; omega
    have hx : x1 = x2 := by
      apply iotaEquiv.injective
      apply Fin.ext
      rw [iotaEquiv_val, iotaEquiv_val, hp]
    subst hx
    rw [Subtype.mk.injEq, Prod.mk.injEq]
    exact ⟨Fin.ext hb, rfl⟩
  let F : Fin 26 × Iota → ℤ := fun y ↦
    D (a, y.1) y.2 * R ^ (bStart y.2.1 - gStart lo + localIdx y.2.2 + L * y.1)
  have hsum : ∑ y : {y // P y}, D (a, y.1.1) y.1.2 * R ^ e y = 0 := by
    have h1 : ∑ y : {y // P y}, D (a, y.1.1) y.1.2 * R ^ e y =
        ∑ y ∈ Finset.univ.filter P, F y :=
      (Finset.sum_subtype (Finset.univ.filter P) (fun y ↦ by simp) F).symm
    have h2 : ∑ y ∈ Finset.univ.filter P, F y = rowPack (L * w) (D *ᵥ radG w lo hi) a := by
      rw [Finset.sum_filter, Fintype.sum_prod_type, rowPack]
      refine Finset.sum_congr rfl fun b' _ ↦ ?_
      simp only [mulVec, dotProduct, Finset.sum_mul]
      refine Finset.sum_congr rfl fun y _ ↦ ?_
      simp only [F, radG, maskW, powS, R]
      by_cases hy : lo ≤ (y.1 : ℕ) ∧ (y.1 : ℕ) < hi
      · simp only [hy, and_self, ↓reduceIte]
        push_cast
        rw [← pow_mul, mul_assoc, ← pow_add]
        congr 2
        ring
      · simp only [hy, ↓reduceIte, mul_zero, zero_mul]
    rw [h1, h2, h a]
  have := radix_eq_zero_inj R hR e he (fun y ↦ D (a, y.1.1) y.1.2)
    (fun y ↦ lt_of_le_of_lt (hD _ _) hK) hsum ⟨(b, x), hx⟩
  simpa using this

/-! ### Entry bounds -/

theorem OM.entryBound {m : ℕ} {X : OM} (hX : X.Good m) : EntryBound (X.mat m) X.1 := by
  intro i j
  have hi : (i : ℕ) < X.2.length := by rw [hX.1]; exact i.isLt
  have hrow := (hX.2.2 _ (List.getElem_mem hi))
  have hj : (j : ℕ) < X.2[(i : ℕ)].length := by rw [hrow.1]; exact j.isLt
  have hx := hrow.2 _ (List.getElem_mem hj)
  simp only [OM.mat, toMatB, toMat_apply]
  rw [lget_shiftL _ _ _ _ hi hj, sh, abs_le]
  constructor <;> omega

theorem EntryBound.copies {d m : Fin 11 → ℕ} {t : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ}
    {K : ℤ} (hK : 0 ≤ K) (h : ∀ b, EntryBound (t b) K) : EntryBound (copies d m t) K := by
  rintro ⟨b, k, c⟩ ⟨b', k', c'⟩
  simp only [SeedKernel.Trial.copies, blockDiagonal'_apply]
  split_ifs with hb
  · subst hb
    simp only [cast_eq, kroneckerMap_apply, one_apply]
    split_ifs
    · simpa using h b c c'
    · simpa using hK
  · simpa using hK

theorem entryBound_copies_OM {X : Fin 11 → OM} {K : ℕ} (hX : ∀ b, (X b).Good (bm b))
    (hK : ∀ b, (X b).1 ≤ K) : EntryBound (copies bd bm fun b ↦ (X b).mat (bm b)) K :=
  EntryBound.copies (by positivity) fun b ↦ (OM.entryBound (hX b)).mono (by exact_mod_cast hK b)


theorem abs_sparseEntry_getD_le {l : List (List (ℕ × ℤ))} {K : ℤ} (hK : 0 ≤ K)
    (h : ∀ g ∈ l, (g.map fun t ↦ |t.2|).sum ≤ K) (n j : ℕ) : |sparseEntry (l.getD n []) j| ≤ K := by
  refine (abs_sparseEntry_le _ _).trans ?_
  rcases mem_getD_or_nil (l := l) (n := n) with hm | hm
  · exact h _ hm
  · rw [hm]; simpa using hK

theorem entryBound_Wm (h : WAbsOk) : EntryBound Wm (2 ^ 11) :=
  fun i x ↦ abs_sparseEntry_getD_le (by norm_num) (h x.1) _ _

theorem entryBound_Vm (h : VAbsOk) : EntryBound Vm (2 ^ 18) :=
  fun x i ↦ abs_sparseEntry_getD_le (by norm_num) (h x.1) _ _


theorem entryBound_T (h : ABoundOk) : EntryBound (mpsTransfer trialA) (3 * (640 * 640)) := by
  have hA (s : Fin 3) : EntryBound (trialA s) 640 := fun i j ↦ h s i j
  rw [mpsTransfer_eq_sum]
  have := (((hA 0).kronecker (hA 0)).add ((hA 1).kronecker (hA 1))).add
    ((hA 2).kronecker (hA 2))
  exact this.mono (by norm_num)

theorem entryBound_O (h : ABoundOk) :
    EntryBound (mpsInsertion trialA) (3 * (676 * ((2 * (640 * 640)) * (2 * (640 * 640))))) := by
  have hA (s : Fin 3) : EntryBound (trialA s) 640 := fun i j ↦ h s i j
  have hk (s s' : Fin 3) : EntryBound (trialA s ⊗ₖ trialA s') (640 * 640) :=
    (hA s).kronecker (hA s')
  have hz := (hk 2 2).sub (hk 0 0)
  have hu := (hk 1 0).add (hk 2 1)
  have hd := (hk 0 1).add (hk 1 2)
  rw [mpsInsertion_eq_sum]
  have := ((hz.mul hz).add (hu.mul hd)).add (hd.mul hu)
  simp only [Fintype.card_prod, Fintype.card_fin] at this
  exact this.mono (by norm_num)


/-! ### The intertwining identities from their list checks -/

/-- The radix vector of all blocks is the plain radix vector `x ↦ R^(position of x)`. -/
def radI (w : ℕ) (x : Iota) : ℤ := ((2 ^ (w * (iotaEquiv x : ℕ)) : ℕ) : ℤ)

theorem radG_all (w : ℕ) : radG w 0 11 = radI w := by
  ext x
  have hx : 0 ≤ (x.1 : ℕ) ∧ (x.1 : ℕ) < 11 := ⟨Nat.zero_le _, x.1.isLt⟩
  simp only [radG, maskW, radI, iotaEquiv_val, hx, and_self, ↓reduceIte]
  rfl

theorem eq_zero_of_mulVec_radI {α : Type*} (D : Matrix α Iota ℤ) (w : ℕ) (K : ℤ)
    (hD : EntryBound D K) (hK : K < 2 ^ w) (h : D *ᵥ radI w = 0) : D = 0 := by
  ext i x
  have hrow : ∑ y, D i y * radI w y = 0 := by
    have := congrFun h i; simpa [mulVec, dotProduct] using this
  set R : ℤ := 2 ^ w
  have hR : 0 < R := by positivity
  have hsum : ∑ j : Fin 676, D i (iotaEquiv.symm j) * R ^ (j : ℕ) = 0 := by
    rw [← hrow, ← Equiv.sum_comp iotaEquiv.symm]
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [radI, Equiv.apply_symm_apply]
    push_cast [R]; rw [← pow_mul]
  have := radix_eq_zero R hR 676 (fun j ↦ D i (iotaEquiv.symm j))
    (fun j ↦ lt_of_le_of_lt (hD _ _) hK) hsum (iotaEquiv x)
  simpa using this

/-- From one packed comparison on the blocks `lo ≤ b < hi` to the vanishing of those columns of
`E W - W X`. -/
theorem columns_eq_of_pack (E : Matrix Pair Pair ℤ) {X : Fin 11 → OM}
    (hX : ∀ b, (X b).Good (bm b)) (hW : WIdxOk) (w lo hi : ℕ) (hlo : lo ≤ 11) (K : ℤ)
    (hD : EntryBound (E * Wm - Wm * copies bd bm (fun b ↦ (X b).mat (bm b))) K)
    (hK : K < 2 ^ w)
    (hhi : hi ≤ 11)
    (hpack : ∀ a, rowPack ((gStart hi - gStart lo) * w) (E *ᵥ vecM (mat26 (YG w lo hi))) a =
      (ZPack ((gStart hi - gStart lo) * w) (ZG w lo hi X)).getD a 0)
    (i : Pair) (x : Iota) (hx : lo ≤ (x.1 : ℕ) ∧ (x.1 : ℕ) < hi) :
    (E * Wm - Wm * copies bd bm (fun b ↦ (X b).mat (bm b))) i x = 0 := by
  refine eq_zero_of_rowPack _ w lo hi K hD hK hlo hhi (fun a ↦ ?_) i x hx
  rw [sub_mulVec, ← mulVec_mulVec, ← mulVec_mulVec, mulVec_Wm_radG hW,
    mulVec_Wm_copies_radG hW hX, rowPack_sub, hpack, rowPack_vecM (sq_ZG _ _ _ _), sub_self]

/-- **`T W = W B`**, given the kernel checks of `TrialChecks` (one group of `676` columns of
width `42`, slots of `676 · 42 = 28392` bits). -/
theorem T_mul_Wm (hW : WIdxOk) (hWa : WAbsOk) (hA : ABoundOk)
    (hB : ∀ b, (tB b).Good (bm b)) (hBo : ∀ b, (tB b).1 ≤ 2 ^ 19)
    (hc : TPack 28392 (YG 42 0 11) = ZPack 28392 (ZG 42 0 11 tB)) :
    mpsTransfer trialA * Wm = Wm * Bm := by
  have hD := ((entryBound_T hA).mul (entryBound_Wm hWa)).sub
    ((entryBound_Wm hWa).mul (entryBound_copies_OM hB hBo))
  rw [card_Pair, card_Iota'] at hD
  have hu : (gStart 11 - gStart 0) * 42 = 28392 := by decide
  rw [← sub_eq_zero]
  ext i x
  exact columns_eq_of_pack _ hB hW 42 0 11 (by norm_num) _ hD (by norm_num) (le_refl _)
    (fun a ↦ by rw [hu, rowPack_T (sq_YG _ _ _), hc]) i x ⟨Nat.zero_le _, x.1.isLt⟩

/-- **`O W = W C`**, given the kernel checks of `TrialChecksInsertion1`–`TrialChecksInsertion4`
(columns of width `71` in four groups: the blocks `0–3`, `4`, `5–6` and `7–10`, with `194`, `130`,
`180` and `172` columns). -/
theorem O_mul_Wm (hW : WIdxOk) (hWa : WAbsOk) (hA : ABoundOk)
    (hC : ∀ b, (oB b).Good (bm b)) (hCo : ∀ b, (oB b).1 ≤ 2 ^ 37)
    (hc1 : OPack 13774 (YG 71 0 4) = ZPack 13774 (ZG 71 0 4 oB))
    (hc2 : OPack 9230 (YG 71 4 5) = ZPack 9230 (ZG 71 4 5 oB))
    (hc3 : OPack 12780 (YG 71 5 7) = ZPack 12780 (ZG 71 5 7 oB))
    (hc4 : OPack 12212 (YG 71 7 11) = ZPack 12212 (ZG 71 7 11 oB)) :
    mpsInsertion trialA * Wm = Wm * Cm := by
  have hD := ((entryBound_O hA).mul (entryBound_Wm hWa)).sub
    ((entryBound_Wm hWa).mul (entryBound_copies_OM hC hCo))
  rw [card_Pair, card_Iota'] at hD
  rw [← sub_eq_zero]
  ext i x
  have hgroup (lo hi u : ℕ) (hlo : lo ≤ 11) (hhi : hi ≤ 11) (hu : (gStart hi - gStart lo) * 71 = u)
      (hc : OPack u (YG 71 lo hi) = ZPack u (ZG 71 lo hi oB))
      (hx : lo ≤ (x.1 : ℕ) ∧ (x.1 : ℕ) < hi) :
      (mpsInsertion trialA * Wm - Wm * Cm) i x = 0 :=
    columns_eq_of_pack _ hC hW 71 lo hi hlo _ hD (by norm_num) hhi
      (fun a ↦ by rw [hu, rowPack_O (sq_YG _ _ _), hc]) i x hx
  have hb := x.1.isLt
  rcases Nat.lt_or_ge (x.1 : ℕ) 4 with h4 | h4
  · exact hgroup 0 4 13774 (by norm_num) (by norm_num) (by decide) hc1 ⟨Nat.zero_le _, h4⟩
  rcases Nat.lt_or_ge (x.1 : ℕ) 5 with h5 | h5
  · exact hgroup 4 5 9230 (by norm_num) (by norm_num) (by decide) hc2 ⟨h4, h5⟩
  rcases Nat.lt_or_ge (x.1 : ℕ) 7 with h7 | h7
  · exact hgroup 5 7 12780 (by norm_num) (by norm_num) (by decide) hc3 ⟨h5, h7⟩
  · exact hgroup 7 11 12212 (by norm_num) (by norm_num) (by decide) hc4 ⟨h7, hb⟩

/-- **`V W = c 1`**, given the kernel checks of `TrialChecks`. -/
theorem Vm_mul_Wm (hW : WIdxOk) (hWa : WAbsOk) (hV : VIdxOk) (hVa : VAbsOk) (hc : VCheck 39) :
    Vm * Wm = (trialScale : ℤ) • 1 := by
  have h1 : EntryBound ((trialScale : ℤ) • (1 : Matrix Iota Iota ℤ)) 161280 := by
    intro i j
    simp only [Matrix.smul_apply, one_apply, smul_eq_mul]
    split_ifs <;> simp [trialScale]
  have hD := ((entryBound_Vm hVa).mul (entryBound_Wm hWa)).sub h1
  rw [card_Pair] at hD
  rw [← sub_eq_zero]
  refine eq_zero_of_mulVec_radI _ 39 _ hD (by norm_num) ?_
  rw [sub_mulVec, ← mulVec_mulVec, ← radG_all, mulVec_Wm_radG hW]
  ext x
  have hY : vecM (mat26 (YG 39 0 11)) =
      fun i : Pair ↦ (fun p ↦ lget (YG 39 0 11) (p / 26) (p % 26)) (pairIdx i) := by
    ext i
    have h1 : ((i.2 : ℕ) + 26 * i.1) / 26 = i.1 := by omega
    have h2 : ((i.2 : ℕ) + 26 * i.1) % 26 = i.2 := by omega
    simp only [vecM, mat26, toMat_apply, pairIdx, h1, h2]
  have hv : (Vm *ᵥ vecM (mat26 (YG 39 0 11))) x = sparseDot ((vG x.1).getD (localIdx x.2) [])
      (fun p ↦ lget (YG 39 0 11) (p / 26) (p % 26)) := by
    rw [hY]; exact mulVec_Vm hV (fun p ↦ lget (YG 39 0 11) (p / 26) (p % 26)) x
  rw [Pi.sub_apply, hv, smul_mulVec, one_mulVec, Pi.smul_apply, smul_eq_mul, radG_all, radI,
    iotaEquiv_val]
  have := hc x.1 (localIdx x.2) (by
    have := (finProdFinEquiv x.2).isLt; simpa [localIdx] using this)
  rw [this, Pi.zero_apply, sub_self]

end SeedKernel.Trial
