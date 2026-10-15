import SeedPrototype.Packed
namespace SeedPrototype

structure BondMasks where
  zero : Nat
  opposite : Nat
  up : Nat
  down : Nat
  a : Nat
  b : Nat
  deriving Inhabited

/-- Hoist repeated denominators and digit masks out of every matvec/filter step. -/
def prepare (w n : Nat) : List BondMasks :=
  let allBits := (1 <<< (w*3^n))-1
  let masks := ((List.range n).map fun j =>
    let p := 3^j
    let x := allBits / ((1 <<< (w*3*p))-1) * ((1 <<< (w*p))-1)
    (x,x <<< (w*p),x <<< (w*2*p))).toArray
  (List.range n).map fun j =>
    let k := (j+1)%n
    let (j0,j1,j2) := masks[j]!
    let (k0,k1,k2) := masks[k]!
    {zero := j1 ||| k1, opposite := (j0 &&& k2) ||| (j2 &&& k0),
     up := (j0 ||| j1) &&& (k1 ||| k2), down := (j1 ||| j2) &&& (k0 ||| k1),
     a := 3^j, b := 3^k}

def filterPrepared (w : Nat) (bs : List BondMasks) (x : Nat) : Nat :=
  bs.foldl (fun y m => y + (x &&& m.zero) + 2*(x &&& m.opposite) +
    move w m.a m.b (x &&& m.up) + move w m.b m.a (x &&& m.down)) 0

def momentsPrepared (w n degree : Nat) : Nat × Nat × Nat :=
  let bs := prepare w n
  let v := (List.range degree).foldl (fun x _ => filterPrepared w bs x)
    (1 <<< (w*neel n))
  let av := filterPrepared w bs v
  (dot w (3^n) v v, dot w (3^n) v av, dot w (3^n) av av)

end SeedPrototype

namespace SeedPrototype
/-- Hoist reversal once per vector; k^2 dot products then use k^2 big multiplications. -/
def gramPrepared (w d : Nat) (vs : List Nat) : List (List Nat) :=
  let backwards := vs.map (reverseLanes w d)
  vs.map fun x => backwards.map fun y => lane w (d-1) (x*y)
end SeedPrototype

namespace SeedPrototype
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
end SeedPrototype
