import SeedKernel.Ising.Checker

/-!
Exact arithmetic checks for trial vectors and comparison weights of the longitudinal Ising chain
with free ends, `-∑_{j<n-1} Z_j Z_{j+1} - (1/4) ∑ X_j - (1/20) ∑ Z_j` on `n` sites. This module
imports no spectral axioms; `Bridge.lean` identifies the checks with the Lean Hamiltonian.

Basis states are `s < 2^n` (bit `j` is the spin at site `j`) and vectors are packed as in
`SeedKernel.Ising.Checker` (`entry`). Everything is scaled by `20`, so that `20 H` has integer
entries: the diagonal is `diagonal n s` and each spin flip carries `-5`.
-/
namespace SeedKernel.IsingOpen

open SeedKernel.Ising

/-- The diagonal `20 E(s) = -20 ∑_{j<n-1} z_j z_{j+1} - ∑_j z_j` of `20 H` on the open chain. -/
def diagonal (n s : Nat) : Int :=
  (List.range (n - 1)).foldl (fun a j => a - 20 * z s j * z s (j + 1)) 0 -
    (List.range n).foldl (fun a j => a + z s j) 0

/-- `(20 H v)(s)` for the packed vector `v`. -/
def applyH (n v s : Nat) : Int :=
  diagonal n s * entry v s -
    5 * (List.range n).foldl (fun a j => a + entry v (s ^^^ (1 <<< j))) 0

/-- `‖S · 20 H v - lam v‖²` for the scale `S = scale`. -/
def residualSq (n v : Nat) (lam : Int) : Nat :=
  (List.range (2^n)).foldl (fun a s =>
    a + (scale * applyH n v s - lam * entry v s).natAbs ^ 2) 0

/-- Number of basis states whose diagonal lies strictly below the cutoff. -/
def classicalCount (n : Nat) (cut : Int) : Nat :=
  (List.range (2^n)).foldl (fun a s => if diagonal n s < cut then a + 1 else a) 0

/-- The weighted neighbour sum `∑_j w(s ⊕ 2^j)` over the flips that stay above the cutoff. -/
def flipSum (n w : Nat) (cut : Int) (s : Nat) : Int :=
  (List.range n).foldl (fun a j =>
    let t := s ^^^ (1 <<< j)
    if diagonal n t < cut then a else a + entry w t) 0

/-- The comparison condition at one basis state `s` above the cutoff: the weight is positive
and `5 ∑_j w(s ⊕ 2^j) ≤ w(s) (20 E(s) - rho)`, that is `g ∑ w(flip)/w(s) ≤ E(s) - rho/20`. -/
def weightOK (n w : Nat) (cut rho : Int) (s : Nat) : Bool :=
  decide (diagonal n s < cut) ||
    (decide (0 < entry w s) && decide (5 * flipSum n w cut s ≤ entry w s * (diagonal n s - rho)))

/-- The comparison condition on the block `start, …, start + count - 1`. -/
def weightBlock (n w : Nat) (cut rho : Int) (start count : Nat) : Bool :=
  (List.range count).all fun i => weightOK n w cut rho (start + i)

/-- Squared norm and squared residual of one vector on a block of basis states. -/
def blockStats (n v : Nat) (lam : Int) (start count : Nat) : Int × Nat :=
  (List.range count).foldl (fun acc i =>
    let s := start + i
    let x := entry v s
    (acc.1 + x * x, acc.2 + (scale * applyH n v s - lam * x).natAbs ^ 2)) (0, 0)

/-- The sum of block statistics. -/
def sumBlocks (l : List (Int × Nat)) : Int × Nat :=
  l.foldl (fun a x => (a.1 + x.1, a.2 + x.2)) (0, 0)

end SeedKernel.IsingOpen
