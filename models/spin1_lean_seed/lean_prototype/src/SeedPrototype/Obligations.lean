import SeedPrototype.PackingProofs
import SeedPrototype.Optimized

/-! Precise open proof obligations, expressed as definitions of propositions, NOT axioms.
    Nothing in the project assumes these propositions. See REPORT.md for their status. -/
namespace SeedPrototype

def DotCorrect : Prop := ∀ (w d m x y : Nat),
  0 < w → 0 < d → x < 2^(w*d) → y < 2^(w*d) →
  (∀ i < d, lane w i x ≤ m ∧ lane w i y ≤ m) → d*m*m < 2^w →
  dot w d x y = ((List.range d).map (fun i => lane w i x * lane w i y)).sum

def ReverseCorrect : Prop := ∀ (w d x : Nat),
  0 < w → 0 < d → x < 2^(w*d) →
  ∀ i < d, lane w i (reverseLanes w d x) = lane w (d-1-i) x

def FastReverseCorrect : Prop := ∀ (w d x : Nat),
  0 < w → 0 < d → x < 2^(w*d) →
  ∀ i < d, lane w i (reverseFast w d x) = lane w (d-1-i) x

def MaskCorrect : Prop := ∀ (w n j t i : Nat),
  0 < w → j < n → t < 3 → i < 3^n →
  lane w i (digitMask w n j t) = if i / 3^j % 3 = t then 2^w-1 else 0

def FloorCorrect : Prop := ∀ (w d s biasBits : Nat) (x : Signed),
  s ≤ biasBits → biasBits < w →
  x.pos < 2^(w*d) → x.neg < 2^(w*d) →
  (∀ i < d, lane w i x.neg ≤ 2^biasBits ∧ lane w i x.pos + 2^biasBits < 2^w) →
  ∀ i < d, Signed.value w i (floorSigned w d s biasBits x) =
    Signed.value w i x / (2^s : Nat)

def magnetization (n s : Nat) : Int :=
  ((List.range n).map fun j => (s / 3^j % 3 : Int)-1).sum

/-- Independent row specification of the Sz Hamiltonian, including a closing pi twist. -/
def rowApply (n : Nat) (twist : Bool) (v : Nat → Int) (s : Nat) : Int :=
  ((List.range n).map fun j =>
    let k := (j+1)%n
    let dj := s/3^j%3
    let dk := s/3^k%3
    let sign : Int := if twist && j+1 == n then -1 else 1
    ((dj : Int)-1)*((dk : Int)-1)*v s +
    (if dj < 2 && 0 < dk then sign*v (s+3^j-3^k) else 0) +
    (if 0 < dj && dk < 2 then sign*v (s+3^k-3^j) else 0)).sum

/-- Bounds cover each limb and every partial row sum (at most 3n contributions). -/
def HApplyCorrect : Prop := ∀ (w n m : Nat) (twist : Bool) (x : Signed),
  2 ≤ n → 0 < w → x.pos < 2^(w*3^n) → x.neg < 2^(w*3^n) →
  (∀ i < 3^n, lane w i x.pos ≤ m ∧ lane w i x.neg ≤ m) → 3*n*m < 2^w →
  ∀ i < 3^n, Signed.value w i (hApply w n twist x) =
    rowApply n twist (fun j => Signed.value w j x) i

/-- The even-ring Marshall phase is (-1)^(sum of ternary digits on odd sites). -/
def marshall (n s : Nat) : Int :=
  if (((List.range n).filter fun j => j%2 == 1).map fun j => s/3^j%3).sum % 2 == 0
  then 1 else -1

def FilterCorrect : Prop := ∀ (w n m x : Nat),
  2 ≤ n → n%2 = 0 → 0 < w → x < 2^(w*3^n) →
  (∀ i < 3^n, lane w i x ≤ m) → 4*n*m < 2^w →
  ∀ i < 3^n, (lane w i (filterApply w n x) : Int) =
    (n : Int)*(lane w i x : Int) -
    marshall n i * rowApply n false (fun j => marshall n j * (lane w j x : Int)) i

end SeedPrototype
