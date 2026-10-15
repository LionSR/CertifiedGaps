import SeedKernel.Packed

/-! # Packed spin-1 Heisenberg ring

The spin-1 ring `H = Σⱼ Sⱼ·Sⱼ₊₁` acting on full `3ⁿ`-lane packed vectors (lane index
`Σⱼ sⱼ 3ʲ`, `sⱼ = mⱼ + 1`), built from the shared lane primitives of `SeedKernel.Packed`.
General correctness obligations are in `Spin1/Obligations.lean`. No native execution is used in
proofs. -/
namespace SeedKernel.Spin1

/-- All bits of all lanes whose j-th ternary digit is t. O(1) bignum recipe. -/
def digitMask (w n j t : Nat) : Nat :=
  let p := 3^j
  let period := 3*p
  let repeats := ((1 <<< (w*3^n)) - 1) / ((1 <<< (w*period)) - 1)
  (repeats * ((1 <<< (w*p)) - 1)) <<< (w*t*p)

/-- H in the standard Sz basis, closing pi twist optional; n >= 2 is required. -/
def hApply (w n : Nat) (twist : Bool) (x : Signed) : Signed := Id.run do
  let mut y : Signed := ⟨0,0⟩
  for j in List.range n do
    let k := (j+1)%n
    let j0 := digitMask w n j 0
    let j1 := digitMask w n j 1
    let j2 := digitMask w n j 2
    let k0 := digitMask w n k 0
    let k1 := digitMask w n k 1
    let k2 := digitMask w n k 2
    let plus := (j0 &&& k0) ||| (j2 &&& k2)
    let minus := (j0 &&& k2) ||| (j2 &&& k0)
    y := y.add (x.select plus)
    y := y.add ((x.select minus).negate)
    let up := ((x.select ((j0 ||| j1) &&& (k1 ||| k2))).move w (3^j) (3^k))
    let down := ((x.select ((j1 ||| j2) &&& (k0 ||| k1))).move w (3^k) (3^j))
    let off := up.add down
    y := y.add (if twist && j+1 == n then off.negate else off)
  return y

/-- nI minus the Marshall conjugate of the untwisted H, for even n. Nonnegative entries. -/
def filterApply (w n x : Nat) : Nat := Id.run do
  let mut y := 0
  for j in List.range n do
    let k := (j+1)%n
    let j0 := digitMask w n j 0
    let j1 := digitMask w n j 1
    let j2 := digitMask w n j 2
    let k0 := digitMask w n k 0
    let k1 := digitMask w n k 1
    let k2 := digitMask w n k 2
    let opposite := (j0 &&& k2) ||| (j2 &&& k0)
    y := y + (x &&& (j1 ||| k1)) + 2*(x &&& opposite)
    y := y + move w (3^j) (3^k) (x &&& ((j0 ||| j1) &&& (k1 ||| k2)))
    y := y + move w (3^k) (3^j) (x &&& ((j1 ||| j2) &&& (k0 ||| k1)))
  return y

def neel (n : Nat) : Nat :=
  ((List.range n).filter fun j => j%2 == 1).foldl (fun s j => s+2*3^j) 0

/-- Regenerated integer trial: an exact power filter; no float input or stored vector. -/
def trial (w n degree : Nat) : Nat :=
  (List.range degree).foldl (fun x _ => filterApply w n x) (1 <<< (w*neel n))

/-- Integer witnesses G=V'V, A=V'(nI-H_M)V, C=((nI-H_M)V)'((nI-H_M)V). -/
def moments (w n degree : Nat) : Nat × Nat × Nat :=
  let v := trial w n degree
  let av := filterApply w n v
  (dot w (3^n) v v, dot w (3^n) v av, dot w (3^n) av av)

/-- H Rayleigh quotient <= -num/den, and squared residual <= (epsNum/epsDen)^2.
    These are arithmetic witnesses only until the packing and spectral bridges are proved. -/
def ritzCheck (w n degree num den epsNum epsDen : Nat) : Bool :=
  let (g,a,c) := moments w n degree
  g > 0 && den*a ≥ (den*n+num)*g &&
    epsDen^2*(c*g-a*a) ≤ epsNum^2*g*g && a*a ≤ c*g

end SeedKernel.Spin1
