import LeanSeedSmall.Checker

namespace LeanSeedSmall

abbrev Stats := Int × Int × Int × Nat × Nat

def addStats (a b : Stats) : Stats :=
  (a.1+b.1, a.2.1+b.2.1, a.2.2.1+b.2.2.1,
    a.2.2.2.1+b.2.2.2.1, a.2.2.2.2+b.2.2.2.2)

/-- Direct scalar evaluation of a consecutive block of computational-basis states. -/
def chunkStats (n u v : Nat) (a b : Int) (start count : Nat) : Stats :=
  (List.range count).foldl (fun acc i =>
    let s := start+i
    let x := entry u s
    let y := entry v s
    addStats acc (x*x, y*y, x*y,
      (scale*applyH n u s-a*x).natAbs^2,
      (scale*applyH n v s-b*y).natAbs^2)) (0,0,0,0,0)

def sumStats (l : List Stats) : Stats := l.foldl addStats (0,0,0,0,0)

def goodStats (s : Stats) : Bool :=
  (s.1-scale^2).natAbs * 10^12 ≤ (scale^2).toNat &&
  (s.2.1-scale^2).natAbs * 10^12 ≤ (scale^2).toNat &&
  s.2.2.1.natAbs * 10^12 ≤ (scale^2).toNat &&
  (s.2.2.2.1+s.2.2.2.2)*10^18 ≤ (scale^4).toNat

end LeanSeedSmall
