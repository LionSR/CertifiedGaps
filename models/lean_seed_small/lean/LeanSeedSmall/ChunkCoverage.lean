import LeanSeedSmall.GeneratedChunks

namespace LeanSeedSmall

/-- Tie the generator-emitted list to a fixed enumeration. A repeated or omitted chunk
cannot pass this equality, even if all its individual arithmetic equalities are correct. -/
theorem allChunks12_eq_range : allChunks12 =
    (List.range 64).map (fun b =>
      chunkStats 12 compact12_1 compact12_2 lam12_1 lam12_2 (64*b) 64) := by rfl

/-- Direct checked bounds with a canonical, generator-independent enumeration of blocks. -/
theorem direct12_canonical : goodStats (sumStats
    ((List.range 64).map (fun b =>
      chunkStats 12 compact12_1 compact12_2 lam12_1 lam12_2 (64*b) 64))) = true := by
  rw [← allChunks12_eq_range]
  exact direct12_chunked

#print axioms allChunks12_eq_range
#print axioms direct12_canonical
end LeanSeedSmall
