import SeedPrototype.Obligations
import SeedPrototype.Optimized
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace SeedPrototype

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

theorem signed_floor_case :
    let x : Signed := ⟨pack 8 [3,0,17,0],pack 8 [0,3,0,17]⟩
    let y := floorSigned 8 4 1 5 x
    (List.range 4).map (fun i => Signed.value 8 i y) = [1,-2,8,-9] := by
  decide +kernel

theorem reverse_case : reverseLanes 8 5 (pack 8 [1,2,3,4,5]) =
    pack 8 [5,4,3,2,1] := by decide +kernel

theorem gram_case : gram 16 3 [pack 16 [1,2,3],pack 16 [4,0,1],pack 16 [2,1,0]] =
    [[14,7,4],[7,17,8],[4,8,5]] := by decide +kernel

-- A concrete counterexample if the convolution carry hypothesis is dropped.
theorem dot_carry_counterexample : dot 2 2 (pack 2 [3,3]) (pack 2 [3,3]) = 0 := by
  decide +kernel

theorem fast_reverse_case : reverseFast 8 5 (pack 8 [1,2,3,4,5]) =
    pack 8 [5,4,3,2,1] := by decide +kernel

theorem gram_prepared_case : gramPrepared 16 3 [pack 16 [1,2,3],pack 16 [4,0,1]] =
    gram 16 3 [pack 16 [1,2,3],pack 16 [4,0,1]] := by decide +kernel

#print axioms signed_h_cases
#print axioms signed_floor_case
#print axioms gram_case
end SeedPrototype
