/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialCompute
import SeedKernel.Heisenberg.OpenTrialData

/-!
# Kernel computations of the open-chain trial bound

The computable side of `OpenTrial`: the lists that `OpenTrialBlocks` evaluates with
`decide +kernel`. This file imports no physics, so that the kernel checks run with a small
environment.

The norm and the energy of the opposite-field trial vector `ψ(σ) = bᵀ A_{σ₁} ⋯ A_{σₘ} c` are
contractions `L₀ ⋯ u` (`L₀ = b ⊗ b`, `u = c ⊗ c`) of the transfer matrix `T`, the bond insertion
`O` and the field insertion `G`. The intertwiner `W` of the periodic trial (`T W = W B`,
`O W = W C`, `V W = c 1`) reduces them to the blocks `tᵦ`, `oᵦ`: `L₀ W` and `L₀ G W` live on the
copy `k = Jᵦ` of each block (their supporting groups of `W` lie in that copy, `InCopy`), so each
block contributes the row iteration of `blockRun` on vectors of length `mᵦ ≤ 32`:

  `x_r = ℓ t^r`,  `z₁ = -300000 ℓ′`,  `z_{r+1} = z_r t + 500000 x_{r-1} o`,

with `ℓ`, `ℓ′` the copy `Jᵦ` of `L₀ W`, `L₀ G W` and `ν`, `ν′` that of `V u`, `V G u`. At each
length `m` it records `(x_m ν, z_m ν, x_{m-1} ν′)`.

The products by `t` and `o` use one big multiplication per entry (`vecmat`): the row is stored
with the offset `β = max |xᵢ|` and packed in radix `2^w`, and each column of the matrix is packed
reversed once per block (`colsOf`). The width `w` is fixed per block, so the run also checks the
lane condition `n² (2β) (2 off) < 2^w` of every product (`vecmatOk`).
-/

namespace SeedKernel.OpenTrial

open SeedKernel.Trial CertifiedGaps.SpinOneHeisenberg

/-! ### Packed vector-matrix products -/

/-- The largest absolute value of a list of integers. -/
def maxAbs : List ℤ → ℕ
  | [] => 0
  | x :: xs => max x.natAbs (maxAbs xs)

/-- The `n` columns of the rows `Y`, each packed reversed in radix `2^w` with its sum. -/
def colsOf (w : ℕ) (Y : List (List ℕ)) (n : ℕ) : List (ℕ × ℕ) :=
  (List.range n).map fun j ↦ (packN w (colN Y j).reverse, sumN (colN Y j))

/-- The row `x` times the matrix stored with offset `off` whose packed columns are `cols`: the
row is stored with the offset `β = maxAbs x`, and every entry is one packed product
(`entryB`). -/
def vecmat (w off n : ℕ) (cols : List (ℕ × ℕ)) (x : List ℤ) : List ℤ :=
  let β := maxAbs x
  let xs := x.map fun v ↦ (v + (β : ℤ)).toNat
  let px := packN w xs
  let sx := sumN xs
  cols.map fun c ↦ (entryB w n β off px sx c.1 c.2 : ℤ) - ((n * β * off : ℕ) : ℤ)

/-- The lane condition of `vecmat`: every convolution coefficient fits in a digit. -/
def vecmatOk (w off n : ℕ) (x : List ℤ) : Bool := n * n * (2 * maxAbs x) * (2 * off) < 2 ^ w

/-! ### The block iteration -/

/-- The lengths `m` of the four trial bounds. -/
def lengths : List ℕ := [34, 78, 142, 240]

/-- The parameters of one block iteration. -/
structure Params where
  /-- The block size. -/
  n : ℕ
  /-- The block `t`, stored with an offset. -/
  t : OM
  /-- The block `o`, stored with an offset. -/
  o : OM
  /-- The lane width of the products by `t`. -/
  wt : ℕ
  /-- The lane width of the products by `o`. -/
  wo : ℕ
  /-- The left vector `ℓ` (copy `J` of `L₀ W`). -/
  l : List ℤ
  /-- The left field vector `ℓ′` (copy `J` of `L₀ G W`). -/
  l' : List ℤ
  /-- The right vector `ν` (copy `J` of `V u`). -/
  v : List ℤ
  /-- The right field vector `ν′` (copy `J` of `V G u`). -/
  v' : List ℤ

namespace Params

variable (P : Params)

/-- The packed columns of `t`. -/
def ct : List (ℕ × ℕ) := colsOf P.wt P.t.2 P.n

/-- The packed columns of `o`. -/
def co : List (ℕ × ℕ) := colsOf P.wo P.o.2 P.n

/-- `x ↦ x t`. -/
def mt (x : List ℤ) : List ℤ := vecmat P.wt P.t.1 P.n P.ct x

/-- `x ↦ x o`. -/
def mo (x : List ℤ) : List ℤ := vecmat P.wo P.o.1 P.n P.co x

/-- One step `(x_{r-1}, x_r, z_r) ↦ (x_r, x_{r+1}, z_{r+1})`. -/
def step (s : List ℤ × List ℤ × List ℤ) : List ℤ × List ℤ × List ℤ :=
  (s.2.1, P.mt s.2.1, vadd (P.mt s.2.2) ((P.mo s.1).map (500000 * ·)))

/-- The lane conditions of the three products of `step`. -/
def stepOk (s : List ℤ × List ℤ × List ℤ) : Bool :=
  vecmatOk P.wt P.t.1 P.n s.2.1 && vecmatOk P.wt P.t.1 P.n s.2.2 &&
    vecmatOk P.wo P.o.1 P.n s.1

/-- The outputs `(x_r ν, z_r ν, x_{r-1} ν′)` of a state. -/
def out (s : List ℤ × List ℤ × List ℤ) : ℤ × ℤ × ℤ :=
  (ldot s.2.1 P.v, ldot s.2.2 P.v, ldot s.1 P.v')

/-- `k` steps from the state `s` at `r`, recording the outputs at the `lengths` and the lane
conditions. -/
def loop : ℕ → ℕ → List ℤ × List ℤ × List ℤ → List (ℤ × ℤ × ℤ) → Bool →
    List (ℤ × ℤ × ℤ) × Bool
  | 0, _, _, acc, ok => (acc, ok)
  | k + 1, r, s, acc, ok =>
    loop k (r + 1) (P.step s) (if lengths.contains r then acc ++ [P.out s] else acc)
      (ok && P.stepOk s)

/-- The initial state `(ℓ, ℓ t, -300000 ℓ′)` at `r = 1`. -/
def init : List ℤ × List ℤ × List ℤ := (P.l, P.mt P.l, P.l'.map (-300000 * ·))

/-- The block iteration: `240` states, the outputs at `m = 34, 78, 142, 240`, and whether every
lane condition holds. -/
def run : List (ℤ × ℤ × ℤ) × Bool :=
  P.loop 240 1 P.init [] (vecmatOk P.wt P.t.1 P.n P.l)

end Params

/-! ### The reduced endpoint vectors -/

/-- The left endpoint as a function of the label. -/
def leftN (a : ℕ) : ℤ := sparseEntry openLeftSparse a

/-- The right endpoint as a function of the label. -/
def rightN (a : ℕ) : ℤ := sparseEntry openRightSparse a

/-- `b Aₛ` as a sparse list (column `j` of `Aₛ` against `b`). -/
def leftA (s : Fin 3) : List (ℕ × ℤ) :=
  ((List.range 26).map fun j ↦ (j, sparseDot ((trialSparseT s).getD j []) leftN)).filter
    fun t ↦ t.2 ≠ 0

/-- `Aₛ c` as a list of length `26`. -/
def rightA (s : Fin 3) : List ℤ :=
  List.ofFn fun i : Fin 26 ↦ sparseDot ((trialSparse s).getD i []) rightN

/-- The groups of `W` met by the pair vector `f ⊗ f` of a sparse `f`, with their weights
`fₐ f_{a'}`: the group of the doubled index `26 a + a'`. -/
def pairGroups (f : List (ℕ × ℤ)) (wg : List (List (ℕ × ℤ))) : List (ℤ × List (ℕ × ℤ)) :=
  f.flatMap fun t ↦ f.map fun t' ↦ (t.2 * t'.2, wg.getD (t'.1 + 26 * t.1) [])

/-- The weighted groups of `L₀ W` in block `b`. -/
def groupsL (b : Fin 11) : List (ℤ × List (ℕ × ℤ)) := pairGroups openLeftSparse (wG b)

/-- The weighted groups of `L₀ G W = (bA₂ ⊗ bA₂) W - (bA₀ ⊗ bA₀) W` in block `b`. -/
def groupsG (b : Fin 11) : List (ℤ × List (ℕ × ℤ)) :=
  pairGroups (leftA 2) (wG b) ++ (pairGroups (leftA 0) (wG b)).map fun q ↦ (-q.1, q.2)

/-- Entry `j` of the row vector with weighted groups `g`. -/
def groupsEntry (g : List (ℤ × List (ℕ × ℤ))) (j : ℕ) : ℤ :=
  (g.map fun q ↦ q.1 * sparseEntry q.2 j).sum

/-- Every group lies in the copy `J` of a block of size `m`. -/
def InCopy (m J : ℕ) (g : List (ℤ × List (ℕ × ℤ))) : Prop := ∀ q ∈ g, ∀ t ∈ q.2, t.1 / m = J

instance (m J : ℕ) (g : List (ℤ × List (ℕ × ℤ))) : Decidable (InCopy m J g) := by
  unfold InCopy; infer_instance

/-- The copy `J` of a row vector given by weighted groups. -/
def rowCopy (b : Fin 11) (g : List (ℤ × List (ℕ × ℤ))) : List ℤ :=
  List.ofFn fun c : Fin (bm b) ↦ groupsEntry g (c + bm b * bJ b)

/-- `u = c ⊗ c` on the doubled index `p = 26 a + a'`. -/
def uN (p : ℕ) : ℤ := rightN (p / 26) * rightN (p % 26)

/-- `G u = A₂c ⊗ A₂c - A₀c ⊗ A₀c` on the doubled index. -/
def guN (p : ℕ) : ℤ :=
  (rightA 2).getD (p / 26) 0 * (rightA 2).getD (p % 26) 0 -
    (rightA 0).getD (p / 26) 0 * (rightA 0).getD (p % 26) 0

/-- The copy `J` of `V Y` in block `b`. -/
def colCopy (b : Fin 11) (Y : ℕ → ℤ) : List ℤ :=
  List.ofFn fun c : Fin (bm b) ↦ sparseDot ((vG b).getD (c + bm b * bJ b) []) Y

/-- The parameters of the iteration of block `b`. -/
def params (b : Fin 11) : Params where
  n := bm b
  t := tB b
  o := oB b
  wt := (openWidths.getD b (0, 0)).1
  wo := (openWidths.getD b (0, 0)).2
  l := rowCopy b (groupsL b)
  l' := rowCopy b (groupsG b)
  v := colCopy b uN
  v' := colCopy b guN

/-- The kernel obligation of block `b`: the groups of `L₀ W` and `L₀ G W` lie in the copy `Jᵦ`,
and the iteration passes its lane conditions and reproduces the stored outputs. -/
def OpenBlockCheck (b : Fin 11) : Prop :=
  InCopy (bm b) (bJ b) (groupsL b) ∧ InCopy (bm b) (bJ b) (groupsG b) ∧
    (params b).run = (openOutputs.getD b [], true)

instance (b : Fin 11) : Decidable (OpenBlockCheck b) := by
  unfold OpenBlockCheck; infer_instance

end SeedKernel.OpenTrial
