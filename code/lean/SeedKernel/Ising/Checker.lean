import Lean

/-!
Exact arithmetic checks for two fixed-point trial vectors of the longitudinal Ising chain.
This module imports no spectral axioms. A successful check is an arithmetic theorem only;
the spectral enclosure theorem and RPRing Hamiltonian identification are separate obligations.
-/
namespace SeedKernel.Ising

def scale : Int := 2^52
def bias : Int := 2^54

def entry (v s : Nat) : Int := Int.ofNat ((v >>> (56*s)) % 2^56) - bias

def z (s j : Nat) : Int := 1 - 2 * Int.ofNat ((s >>> j) % 2)

def diagonal (n s : Nat) : Int :=
  (List.range n).foldl (fun a j => a - 20 * z s j * z s ((j+1)%n) - z s j) 0

def applyH (n v s : Nat) : Int :=
  diagonal n s * entry v s -
    5 * (List.range n).foldl (fun a j => a + entry v (s ^^^ (1 <<< j))) 0

def dot (n u v : Nat) : Int :=
  (List.range (2^n)).foldl (fun a s => a + entry u s * entry v s) 0

def residualSq (n v : Nat) (lam : Int) : Nat :=
  (List.range (2^n)).foldl (fun a s =>
    a + (scale * applyH n v s - lam * entry v s).natAbs ^ 2) 0

def checkVectors (n u v : Nat) (a b : Int) : Bool :=
  let s2 := (scale^2).toNat
  (dot n u u - scale^2).natAbs * 10^12 ≤ s2 &&
  (dot n v v - scale^2).natAbs * 10^12 ≤ s2 &&
  (dot n u v).natAbs * 10^12 ≤ s2 &&
  (residualSq n u a + residualSq n v b) * 10^18 ≤ s2^2 &&
  a + scale < b && b + scale < (-(21 : Int)*n+82-5*n)*scale

/-- Number of diagonal entries strictly below a proposed classical cutoff. -/
def classicalCount (n : Nat) (cut : Int) : Nat :=
  (List.range (2^n)).foldl (fun a s => if diagonal n s < cut then a+1 else a) 0

end SeedKernel.Ising
