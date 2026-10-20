/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialData

/-!
# The integer matrix-product trial state

The trial state of the spin-one seed is the periodic matrix-product state with bond dimension
`26` built from 23 integers (`models/spin1/src/trials.jl`, `trial_matrices()`, which matches the
upstream evidence `matrix-data.json`). The virtual space is a sum of eight multiplets `(l, a)` of
spin `l / 2` (four of spin `1/2`, three of spin `3/2`, one of spin `5/2`); the label `(l, a, d)`
with `d = 2m ∈ {-l, -l + 2, …, l}` runs over their states in the Julia order. The tensor maps
`(l, a, d)` to `(k, b, d + 2(s - 1))` with a reduced matrix element `p` (eight diagonal and
fifteen off-diagonal integers) times an integer Clebsch–Gordan-type factor.

* `trialA s`: the recipe as a Mathlib matrix (`s : Fin 3` labels `m = s - 1`).
* `trialA_eq_code`: it equals the matrices written by the Julia generator (`TrialData`).
* `trialSparse s`: its nonzero entries row by row, used by the kernel computations.
-/

namespace CertifiedGaps.SpinOneHeisenberg

open Matrix

/-- The multiplets `(l, a)`: `l = 2j` and the copy index `a`. -/
def trialMultiplets : List (ℕ × ℕ) :=
  [(1, 0), (1, 1), (1, 2), (1, 3), (3, 0), (3, 1), (3, 2), (5, 0)]

/-- The 26 virtual labels `(l, a, d)`, `d = -l, -l + 2, …, l`, in the Julia order. -/
def trialLabels : List (ℕ × ℕ × ℤ) :=
  trialMultiplets.flatMap fun x ↦ (List.range (x.1 + 1)).map fun i ↦ (x.1, x.2, 2 * (i : ℤ) - x.1)

/-- The diagonal reduced elements, one per multiplet. -/
def trialDiagonal : List ℤ := [-289, -84, -9, 33, -22, 27, 71, -11]

/-- The pairs of multiplets `(x, y)` with `l_y = l_x + 2`, in the Julia order. -/
def trialPairs : List ((ℕ × ℕ) × (ℕ × ℕ)) :=
  trialMultiplets.flatMap fun x ↦ (trialMultiplets.filter fun y ↦ y.1 = x.1 + 2).map fun y ↦ (x, y)

/-- The off-diagonal reduced elements, one per pair of `trialPairs`. -/
def trialOffDiagonal : List ℤ :=
  [-11, -35, -28, -23, 16, -49, -6, -25, -20, -34, 9, -32, 1, 2, -16]

/-- The reduced matrix element between multiplets `x = (l, a)` and `y = (k, b)`. -/
def trialReduced (x y : ℕ × ℕ) : ℤ :=
  if y.1 = x.1 then (if y = x then trialDiagonal.getD (trialMultiplets.idxOf x) 0 else 0)
  else if y.1 = x.1 + 2 then trialOffDiagonal.getD (trialPairs.idxOf (x, y)) 0
  else if x.1 = y.1 + 2 then -trialOffDiagonal.getD (trialPairs.idxOf (y, x)) 0
  else 0

/-- The Clebsch–Gordan-type factor for `l → k` at `d`, physical label `s` (`m = s - 1`). -/
def trialFactor (l k : ℕ) (d : ℤ) (s : Fin 3) : ℤ :=
  let L : ℤ := l
  let r : ℤ := (d + L) / 2
  if k = l + 2 then ![2 * (L + 1 - r) * (L + 2 - r), 2 * (r + 1) * (L - r + 1), (r + 1) * (r + 2)] s
  else if k = l then ![2 * (L - r + 1), d, -(r + 1)] s
  else ![2, -2, 1] s

/-- Entry `(i, j)` of the tensor with physical label `s`. -/
def trialEntry (s : Fin 3) (i j : ℕ) : ℤ :=
  let x := trialLabels.getD i (0, 0, 0)
  let y := trialLabels.getD j (0, 0, 0)
  if y.2.2 = x.2.2 + 2 * ((s : ℤ) - 1) then
    trialReduced (x.1, x.2.1) (y.1, y.2.1) * trialFactor x.1 y.1 x.2.2 s
  else 0

/-- **The trial tensor** `A s : Matrix (Fin 26) (Fin 26) ℤ` of the spin-one seed,
`s ↔ m = s - 1`. -/
def trialA (s : Fin 3) : Matrix (Fin 26) (Fin 26) ℤ := Matrix.of fun i j ↦ trialEntry s i j

/-- The recipe as rows of a list matrix. -/
def trialRows (s : Fin 3) : List (List ℤ) :=
  List.ofFn fun i : Fin 26 ↦ List.ofFn fun j : Fin 26 ↦ trialEntry s i j

theorem toMat_trialRows (s : Fin 3) : SeedKernel.Trial.toMat 26 26 (trialRows s) = trialA s := by
  ext i j
  have h1 : (trialRows s).getD i [] = List.ofFn fun j : Fin 26 ↦ trialEntry s i j := by
    rw [trialRows, List.getD_eq_getElem _ _ (by rw [List.length_ofFn]; exact i.isLt),
      List.getElem_ofFn]
  rw [SeedKernel.Trial.toMat_apply, SeedKernel.Trial.lget, h1,
    List.getD_eq_getElem _ _ (by rw [List.length_ofFn]; exact j.isLt), List.getElem_ofFn]
  rfl

/-- The nonzero entries of row `i` of `trialA s`, as `(column, value)`. -/
def trialSparse (s : Fin 3) : List (List (ℕ × ℤ)) :=
  List.ofFn fun i : Fin 26 ↦
    (List.ofFn fun j : Fin 26 ↦ ((j : ℕ), trialEntry s i j)).filter fun t ↦ t.2 ≠ 0

/-- The nonzero entries of column `d` of `trialA s`, as `(row, value)`. -/
def trialSparseT (s : Fin 3) : List (List (ℕ × ℤ)) :=
  List.ofFn fun d : Fin 26 ↦
    (List.ofFn fun b : Fin 26 ↦ ((b : ℕ), trialEntry s b d)).filter fun t ↦ t.2 ≠ 0

end CertifiedGaps.SpinOneHeisenberg
