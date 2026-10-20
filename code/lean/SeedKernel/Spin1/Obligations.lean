import SeedKernel.Obligations
import SeedKernel.Spin1.Packed

/-! Precise open proof obligations for the packed spin-1 operators, expressed as definitions of
propositions, NOT axioms. Nothing in the library assumes these propositions. -/
namespace SeedKernel.Spin1

def MaskCorrect : Prop := ∀ (w n j t i : Nat),
  0 < w → j < n → t < 3 → i < 3^n →
  lane w i (digitMask w n j t) = if i / 3^j % 3 = t then 2^w-1 else 0

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

end SeedKernel.Spin1
