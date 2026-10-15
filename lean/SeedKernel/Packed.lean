import Std

/-! # Packed `Nat` lanes

A vector of `d` small naturals is stored as one big `Nat` in radix `2^w` (lane zero least
significant), so that the Lean kernel, which evaluates `Nat` arithmetic through GMP, performs a
whole vector operation in a few big-integer steps. This module holds the lane primitives shared
by the spin-1 (`SeedKernel.Spin1`) and Ising (`SeedKernel.Ising`) checkers.

Everything here is executable arithmetic, checked by `decide +kernel` only (never
`native_decide`). Correctness of the primitives against a lane-by-lane specification is an open
obligation (`SeedKernel/Obligations.lean`); checked examples do not replace it. -/
namespace SeedKernel

/-- Radix is 2^w; lane zero is least significant. -/
def lane (w i x : Nat) : Nat := (x >>> (w * i)) % (1 <<< w)

/-- Every one of the `d` lanes equal to one. -/
def repunit (w d : Nat) : Nat := ((1 <<< (w*d)) - 1) / ((1 <<< w) - 1)

def pack (w : Nat) (xs : List Nat) : Nat :=
  xs.foldr (fun x acc => x + (acc <<< w)) 0

/-- Select and translate a set of source lanes by a fixed index displacement. -/
def move (w a b x : Nat) : Nat :=
  if b ≤ a then x <<< (w*(a-b)) else x >>> (w*(b-a))

/-- A signed packed vector as a difference of two unsigned ones. -/
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
  ⟨SeedKernel.move w a b x.pos, SeedKernel.move w a b x.neg⟩

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

/-- Hoist reversal once per vector; k^2 dot products then use k^2 big multiplications. -/
def gramPrepared (w d : Nat) (vs : List Nat) : List (List Nat) :=
  let backwards := vs.map (reverseLanes w d)
  vs.map fun x => backwards.map fun y => lane w (d-1) (x*y)

/-- Repeat a mask by doubling. This avoids huge integer divisions in reversal masks. -/
def reversalMask (w levels level : Nat) : Nat :=
  let g := 1 <<< level
  (List.range (levels-level-1)).foldl
    (fun m j => m ||| (m <<< (w*2*g*(1 <<< j)))) ((1 <<< (w*g))-1)

def reverseFast (w d x : Nat) : Nat := Id.run do
  let levels := Nat.log2 (d-1)+1
  let size := 1 <<< levels
  let mut y := x
  for level in List.range levels do
    let g := 1 <<< level
    let low := reversalMask w levels level
    y := ((y &&& low) <<< (w*g)) ||| ((y >>> (w*g)) &&& low)
  return y >>> (w*(size-d))

def dotFast (w d x y : Nat) : Nat := lane w (d-1) (x*reverseFast w d y)

end SeedKernel
