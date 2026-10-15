import SeedKernel.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0

/-! Kernel-checked examples of the shared lane primitives (`decide +kernel`). They are
regression checks, not proofs of the obligations in `SeedKernel/Obligations.lean`. -/
namespace SeedKernel

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

end SeedKernel
