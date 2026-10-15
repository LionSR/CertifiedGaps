import SeedKernel.Spin1.Packed

/-! Prepared-mask variants of the spin-1 filter: hoist repeated denominators and digit masks
out of every matvec/filter step. The generic fast reversal and prepared Gram matrix live in
`SeedKernel.Packed`. -/
namespace SeedKernel.Spin1

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

end SeedKernel.Spin1
