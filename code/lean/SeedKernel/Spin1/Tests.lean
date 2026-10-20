import SeedKernel.Spin1.Obligations
set_option maxRecDepth 100000
set_option maxHeartbeats 0

/-! Kernel-checked comparisons of the packed spin-1 operators with the independent row
specification `rowApply`, on small rings. -/
namespace SeedKernel.Spin1

def referenceTest (n : Nat) (twist : Bool) : Bool :=
  let p := pack 32 ((List.range (3^n)).map fun i => i%7)
  let q := pack 32 ((List.range (3^n)).map fun i => (i+3)%5)
  let x : Signed := ⟨p,q⟩
  let y := hApply 32 n twist x
  (List.range (3^n)).all fun i => Signed.value 32 i y ==
    rowApply n twist (fun j => Signed.value 32 j x) i

def filterTest (n : Nat) : Bool :=
  let x := pack 32 ((List.range (3^n)).map fun i => i%7)
  let y := filterApply 32 n x
  (List.range (3^n)).all fun i => (lane 32 i y : Int) ==
    (n : Int)*(lane 32 i x : Int) - marshall n i *
      rowApply n false (fun j => marshall n j * (lane 32 j x : Int)) i

theorem signed_h_cases :
    [2,3,4].all (fun n => referenceTest n false && referenceTest n true) = true := by
  decide +kernel

theorem filter_cases : [2,4].all filterTest = true := by decide +kernel

end SeedKernel.Spin1
