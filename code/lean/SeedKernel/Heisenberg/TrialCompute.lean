/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialSparse
import SeedKernel.Heisenberg.TrialPowers
import SeedKernel.Heisenberg.TrialRecipe

/-!
# Kernel computations of the intertwiner checks

The computable side of `TrialIntertwiner`: the decoded data and the lists of integers that
`TrialChecks` compares with `decide +kernel`. This file imports no physics, so that the kernel
checks run with a small environment.

A vector on the doubled index `p = 26 a + c` evaluated on a radix vector is a `26 × 26` list
matrix of radix numbers (`YG`, `ZG`). The comparisons are made on packed rows: entry `(a, b)`
enters row `a` with weight `S^b`, `S = 2^u` (`TPack`, `OPack`, `ZPack`), where the slot width
`u = L w` holds the `L` columns of width `w`. The radix vector may be restricted to a range of
blocks `lo ≤ b < hi` (with exponents counted from the first column of block `lo`, and
`L = start hi - start lo`), which keeps the numbers small when a check is split.
-/

namespace SeedKernel.Trial

open CertifiedGaps.SpinOneHeisenberg

/-! ### Blocks and data -/

/-- Spin of block `b`. -/
def bJ (b : Fin 11) : ℕ := blockJ.getD b 0
/-- Multiplicity (size of `tᵦ`) of block `b`. -/
def bm (b : Fin 11) : ℕ := blockDim.getD b 0
/-- Number of copies `2J + 1` of block `b`. -/
def bd (b : Fin 11) : ℕ := 2 * bJ b + 1
/-- First column of block `b`. -/
def bStart (b : Fin 11) : ℕ := blockStart.getD b 0
/-- First column of block `lo` (`676` past the last block). -/
def gStart (lo : ℕ) : ℕ := blockStart.getD lo 676

/-- Block `tᵦ` as an offset matrix. -/
def tB (b : Fin 11) : OM := tOM b
/-- Block `oᵦ` as an offset matrix. -/
def oB (b : Fin 11) : OM := oOM b

instance : Inhabited SparseCode := ⟨⟨0, (0, 0), (0, 0), (0, 0, 0)⟩⟩

/-- The groups of `W` restricted to block `b`, one per doubled index. -/
def wG (b : Fin 11) : List (List (ℕ × ℤ)) := decodeSparse 676 (wCode.getD b default)
/-- The groups of `V` restricted to block `b`, one per local column. -/
def vG (b : Fin 11) : List (List (ℕ × ℤ)) := decodeSparse (bd b * bm b) (vCode.getD b default)

/-! ### Well-formedness and bounds of the data -/

/-- Indices of `W`'s groups lie inside their block. -/
def WIdxOk : Prop := ∀ b : Fin 11, ∀ g ∈ wG b, ∀ t ∈ g, t.1 < bd b * bm b
/-- Indices of `V`'s groups are doubled indices. -/
def VIdxOk : Prop := ∀ b : Fin 11, ∀ g ∈ vG b, ∀ t ∈ g, t.1 < 676
/-- Each group of `W` has absolute sum at most `2^11`. -/
def WAbsOk : Prop := ∀ b : Fin 11, ∀ g ∈ wG b, (g.map fun t ↦ |t.2|).sum ≤ 2 ^ 11
/-- Each group of `V` has absolute sum at most `2^18`. -/
def VAbsOk : Prop := ∀ b : Fin 11, ∀ g ∈ vG b, (g.map fun t ↦ |t.2|).sum ≤ 2 ^ 18
/-- Every entry of the recipe is at most `640` in absolute value. -/
def ABoundOk : Prop := ∀ s : Fin 3, ∀ i j : Fin 26, |trialEntry s i j| ≤ 640

instance : Decidable WIdxOk := by unfold WIdxOk; infer_instance
instance : Decidable VIdxOk := by unfold VIdxOk; infer_instance
instance : Decidable WAbsOk := by unfold WAbsOk; infer_instance
instance : Decidable VAbsOk := by unfold VAbsOk; infer_instance
instance : Decidable ABoundOk := by unfold ABoundOk; infer_instance

/-! ### `W` and `W B` on a radix vector -/

/-- The values of block `b`'s groups of `W` against `F`, one per doubled index. -/
def blockVec (b : Fin 11) (F : ℕ → ℤ) : List ℤ := (wG b).map fun g ↦ sparseDot g F

/-- `∑ᵦ` of the block values, as a flat list over the doubled index. -/
def flatW (F : Fin 11 → ℕ → ℤ) : List ℤ :=
  vsumN 676 ((List.finRange 11).map fun b ↦ blockVec b (F b))

/-- The radix weights of block `b` restricted to the blocks `lo ≤ b < hi`: local column `j`
gets `R^(start b - start lo + j)`, `R = 2^w`. -/
def maskW (w lo hi : ℕ) (b : Fin 11) (j : ℕ) : ℤ :=
  if lo ≤ (b : ℕ) ∧ (b : ℕ) < hi then ((2 ^ (w * (bStart b - gStart lo + j)) : ℕ) : ℤ) else 0

/-- `W` on the restricted radix vector, as a `26 × 26` list matrix. -/
def YG (w lo hi : ℕ) : List (List ℤ) := chunksZ 26 26 (flatW (maskW w lo hi))

/-- Row `c` of an offset matrix on the radix vector `(R^c')`: the packed row minus the offset
times the packed all-ones row. -/
def ubT (w : ℕ) (X : OM) (m c : ℕ) : ℤ :=
  (packN w (X.2.getD c []) : ℤ) - X.1 * (packN w (List.replicate m 1) : ℤ)

/-- The block matrix applied to the restricted radix vector, at local column `j` of block `b`. -/
def zBG (w lo hi : ℕ) (X : Fin 11 → OM) (b : Fin 11) (j : ℕ) : ℤ :=
  if lo ≤ (b : ℕ) ∧ (b : ℕ) < hi then
    ((2 ^ (w * (bStart b - gStart lo + bm b * (j / bm b))) : ℕ) : ℤ) * ubT w (X b) (bm b) (j % bm b)
  else 0

/-- `W B` (or `W C`) on the restricted radix vector, as a `26 × 26` list matrix. -/
def ZG (w lo hi : ℕ) (X : Fin 11 → OM) : List (List ℤ) := chunksZ 26 26 (flatW (zBG w lo hi X))

/-! ### Packed rows -/

/-- `S^b`, `S = 2^u`. -/
def powS (u b : ℕ) : ℤ := ((2 ^ (u * b) : ℕ) : ℤ)

/-- `(S^b)_b`. -/
def sigS (u : ℕ) : List ℤ := List.ofFn fun b : Fin 26 ↦ powS u b

/-- `(A s)ᵀ g` for a vector `g` of length `26`. -/
def gamL (s : Fin 3) (g : List ℤ) : List ℤ :=
  List.ofFn fun d : Fin 26 ↦ sparseDot ((trialSparseT s).getD d []) fun e ↦ g.getD e 0

/-- `(A s) v` for a vector `v` of length `26`. -/
def lcombL (s : Fin 3) (v : List ℤ) : List ℤ :=
  List.ofFn fun a : Fin 26 ↦ sparseDot ((trialSparse s).getD a []) fun c ↦ v.getD c 0

/-- `Y g`: one dot product per row. -/
def prowL (Y : List (List ℤ)) (g : List ℤ) : List ℤ := Y.map fun r ↦ ldot r g

/-- Packed rows of `T` applied to the vectorized list matrix `Y`. -/
def TPack (u : ℕ) (Y : List (List ℤ)) : List ℤ :=
  vadd (vadd (lcombL 0 (prowL Y (gamL 0 (sigS u)))) (lcombL 1 (prowL Y (gamL 1 (sigS u)))))
    (lcombL 2 (prowL Y (gamL 2 (sigS u))))

/-- The term `(A s A s') Y (A t A t')ᵀ` of `O`, on packed rows. -/
def OTerm (u : ℕ) (Y : List (List ℤ)) (s s' t t' : Fin 3) : List ℤ :=
  lcombL s (lcombL s' (prowL Y (gamL t' (gamL t (sigS u)))))

/-- Packed rows of `O` applied to the vectorized list matrix `Y`: the twelve terms of
`Y_z² + Y_u Y_d + Y_d Y_u`. -/
def OPack (u : ℕ) (Y : List (List ℤ)) : List ℤ :=
  vadd (vadd (vadd (vadd (vadd (vadd (vadd (vadd (vadd (vadd (vadd
    (OTerm u Y 2 2 2 2) (vneg (OTerm u Y 2 0 2 0))) (vneg (OTerm u Y 0 2 0 2)))
    (OTerm u Y 0 0 0 0)) (OTerm u Y 1 0 0 1)) (OTerm u Y 1 1 0 2)) (OTerm u Y 2 0 1 1))
    (OTerm u Y 2 1 1 2)) (OTerm u Y 0 1 1 0)) (OTerm u Y 0 2 1 1)) (OTerm u Y 1 1 2 0))
    (OTerm u Y 1 2 2 1)

/-- Packed rows of a `26 × 26` list matrix. -/
def ZPack (u : ℕ) (Z : List (List ℤ)) : List ℤ := Z.map fun r ↦ ldot r (sigS u)

/-! ### `V W` on a radix vector -/

/-- The list check of `V W = c 1` on the radix vector of width `w`. -/
def VCheck (w : ℕ) : Prop :=
  ∀ b : Fin 11, ∀ r < bd b * bm b,
    sparseDot ((vG b).getD r []) (fun p ↦ lget (YG w 0 11) (p / 26) (p % 26)) =
      trialScale * ((2 ^ (w * (bStart b + r)) : ℕ) : ℤ)

instance (w : ℕ) : Decidable (VCheck w) := by unfold VCheck; infer_instance

end SeedKernel.Trial
