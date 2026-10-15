import SeedKernel.Ising.Checker
import SeedKernel.Packed

/-! Arithmetic implementation of packed Gram/residual checks on `2^n` lanes of `width` bits,
built on the shared lane primitives of `SeedKernel.Packed`: `ones` is `repunit`, `middle` the
top `lane`, and `rawDot` the shared reversal-and-multiply `dot`. The generic lane-decoding proof
is NOT supplied here; this is explicitly an additional gap before spectral use. -/
namespace SeedKernel.Ising.Packed

def width : Nat := 288
def whole (n : Nat) : Nat := (1 <<< (width * 2^n)) - 1
def ones (n : Nat) : Nat := repunit width (2^n)
def lowMask (n j : Nat) : Nat :=
  ((1 <<< (width * 2^j)) - 1) * (whole n / ((1 <<< (width * 2^(j+1))) - 1))
def highMask (n j : Nat) : Nat := whole n - lowMask n j

/-- Expand contiguous 56-bit slots to 288-bit slots in `n` masked shifts. -/
def widen (n u : Nat) : Nat :=
  ((List.range n).reverse).foldl (fun a j =>
    let low := ((1 <<< (56*2^j))-1) *
      (whole n / ((1 <<< (width*2^(j+1)))-1))
    (a &&& low) + ((a &&& (low <<< (56*2^j))) <<< ((width-56)*2^j))) u
/-- Flip spin `j`: swap the lane blocks whose index bit `j` is 0 and 1. -/
def flip (n j u : Nat) : Nat :=
  ((u &&& lowMask n j) <<< (width * 2^j)) +
  ((u &&& highMask n j) >>> (width * 2^j))
def middle (n u : Nat) : Nat := lane width (2^n-1) u
def sumDigits (n u : Nat) : Nat := middle n (u * ones n)
def rawDot (n u v : Nat) : Nat := SeedKernel.dot width (2^n) u v
def centeredDot (n u v : Nat) (c : Int) : Int :=
  Int.ofNat (rawDot n u v) - c * (sumDigits n u + sumDigits n v) + c*c*(2^n : Nat)

def hPacked (n u : Nat) : Int :=
  let pos := (List.range n).foldl (fun a j =>
    a + 40 * (u &&& (highMask n j ^^^ highMask n ((j+1)%n))) +
      2 * (u &&& highMask n j)) 0
  let neg := (List.range n).foldl (fun a j => a + flip n j u) 0
  Int.ofNat pos - (21*n : Nat)*Int.ofNat u - 5*Int.ofNat neg

def residual (n u : Nat) (lam : Int) : Nat :=
  (scale * (hPacked n u - bias * hPacked n (ones n)) -
    lam * (Int.ofNat u - bias * ones n) + (2^128 : Int)*ones n).toNat
def resNorm (n u : Nat) (lam : Int) : Int :=
  let r := residual n u lam
  centeredDot n r r (2^128)

def check (n u v : Nat) (a b : Int) : Bool :=
  let s2 := (scale^2).toNat
  (centeredDot n u u bias - scale^2).natAbs * 10^12 ≤ s2 &&
  (centeredDot n v v bias - scale^2).natAbs * 10^12 ≤ s2 &&
  (centeredDot n u v bias).natAbs * 10^12 ≤ s2 &&
  0 ≤ resNorm n u a && 0 ≤ resNorm n v b &&
  (resNorm n u a + resNorm n v b).toNat * 10^18 ≤ s2^2 &&
  a + scale < b && b + scale < (-(21 : Int)*n+82-5*n)*scale

end SeedKernel.Ising.Packed
