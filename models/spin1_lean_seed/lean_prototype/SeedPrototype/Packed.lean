import Std

/-! Executable packed-Nat arithmetic. General correctness obligations are in Obligations.lean.
    Checked examples do not replace those obligations. No native execution is used in proofs. -/
namespace SeedPrototype

/-- Radix is 2^w; lane zero is least significant. -/
def lane (w i x : Nat) : Nat := (x >>> (w * i)) % (1 <<< w)

def repunit (w d : Nat) : Nat := ((1 <<< (w*d)) - 1) / ((1 <<< w) - 1)

def pack (w : Nat) (xs : List Nat) : Nat :=
  xs.foldr (fun x acc => x + (acc <<< w)) 0

/-- All bits of all lanes whose j-th ternary digit is t. O(1) bignum recipe. -/
def digitMask (w n j t : Nat) : Nat :=
  let p := 3^j
  let period := 3*p
  let repeats := ((1 <<< (w*3^n)) - 1) / ((1 <<< (w*period)) - 1)
  (repeats * ((1 <<< (w*p)) - 1)) <<< (w*t*p)

/-- Select and translate a set of source lanes by a fixed base-3 index displacement. -/
def move (w a b x : Nat) : Nat :=
  if b ≤ a then x <<< (w*(a-b)) else x >>> (w*(b-a))

structure Signed where
  pos : Nat
  neg : Nat
  deriving DecidableEq, Repr

def Signed.value (w i : Nat) (x : Signed) : Int :=
  (lane w i x.pos : Int) - (lane w i x.neg : Int)

def Signed.add (x y : Signed) : Signed := ⟨x.pos+y.pos, x.neg+y.neg⟩
def Signed.negate (x : Signed) : Signed := ⟨x.neg,x.pos⟩
def Signed.select (m : Nat) (x : Signed) : Signed := ⟨x.pos &&& m,x.neg &&& m⟩
def Signed.move (w a b : Nat) (x : Signed) : Signed :=
  ⟨SeedPrototype.move w a b x.pos, SeedPrototype.move w a b x.neg⟩

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

/-- Per-lane floor(x/2^s), for s <= w; masks inter-lane contamination. -/
def floorUnsigned (w d s x : Nat) : Nat :=
  (x >>> s) &&& (((1 <<< (w-s))-1)*repunit w d)

/-- Signed floor with a 2^biasBits offset. Preconditions: s<=biasBits<w,
    each negative limb <=2^biasBits and each positive limb+2^biasBits < 2^w. -/
def floorSigned (w d s biasBits : Nat) (x : Signed) : Signed :=
  let bias := (1 <<< biasBits)*repunit w d
  ⟨floorUnsigned w d s (x.pos+bias-x.neg),
    (1 <<< (biasBits-s))*repunit w d⟩

/-- Reverse d lanes via ceil(log2 d) butterfly swaps, padded to a power of two. -/
def reverseLanes (w d x : Nat) : Nat := Id.run do
  let levels := Nat.log2 (d-1)+1
  let size := 1 <<< levels
  let mut y := x
  for level in List.range levels do
    let g := 1 <<< level
    let low := (((1 <<< (w*size))-1) / ((1 <<< (w*2*g))-1)) * ((1 <<< (w*g))-1)
    y := ((y &&& low) <<< (w*g)) ||| ((y >>> (w*g)) &&& low)
  return y >>> (w*(size-d))

/-- One multiplication, provided every convolution coefficient is below the radix.
    For lanes <=M, d*M^2 < 2^w suffices. Reversal is O(log d) bitwise operations. -/
def dot (w d x y : Nat) : Nat :=
  lane w (d-1) (x*reverseLanes w d y)

def signedDot (w d : Nat) (x y : Signed) : Int :=
  (dot w d x.pos y.pos : Int) + (dot w d x.neg y.neg : Int) -
  (dot w d x.pos y.neg : Int) - (dot w d x.neg y.pos : Int)

def gram (w d : Nat) (vs : List Nat) : List (List Nat) :=
  vs.map fun x => vs.map fun y => dot w d x y

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

end SeedPrototype
