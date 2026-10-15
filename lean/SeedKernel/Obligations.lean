import SeedKernel.Packed

/-! Precise open proof obligations for the shared lane primitives, expressed as definitions of
propositions, NOT axioms. Nothing in the library assumes these propositions. The spin-1
specific obligations are in `Spin1/Obligations.lean`. -/
namespace SeedKernel

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

def FloorCorrect : Prop := ∀ (w d s biasBits : Nat) (x : Signed),
  s ≤ biasBits → biasBits < w →
  x.pos < 2^(w*d) → x.neg < 2^(w*d) →
  (∀ i < d, lane w i x.neg ≤ 2^biasBits ∧ lane w i x.pos + 2^biasBits < 2^w) →
  ∀ i < d, Signed.value w i (floorSigned w d s biasBits x) =
    Signed.value w i x / (2^s : Nat)

end SeedKernel
