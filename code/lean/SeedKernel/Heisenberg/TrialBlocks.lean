/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialPowers
import SeedKernel.Heisenberg.TrialChecksInsertion4

/-!
# Kernel evaluation of the block traces

One `decide +kernel` per block evaluates the packed power chain of `TrialPowers` and compares
the four traces with the stored values. The blocks are checked one after another so that the
peak memory is that of the largest block (`m = 32`). The import of `TrialChecksInsertion4` only
orders the kernel runs of the trial files, so that a parallel build holds one of them at a time.
-/

namespace SeedKernel.Trial

set_option maxRecDepth 100000

theorem blockCheck_0 : BlockCheck 0 := by decide +kernel
theorem blockCheck_1 : BlockCheck 1 := by decide +kernel
theorem blockCheck_2 : BlockCheck 2 := by decide +kernel
theorem blockCheck_3 : BlockCheck 3 := by decide +kernel
theorem blockCheck_4 : BlockCheck 4 := by decide +kernel
theorem blockCheck_5 : BlockCheck 5 := by decide +kernel
theorem blockCheck_6 : BlockCheck 6 := by decide +kernel
theorem blockCheck_7 : BlockCheck 7 := by decide +kernel
theorem blockCheck_8 : BlockCheck 8 := by decide +kernel
theorem blockCheck_9 : BlockCheck 9 := by decide +kernel
theorem blockCheck_10 : BlockCheck 10 := by decide +kernel

/-- Every block passes its kernel check. -/
theorem blockCheck_all : ∀ b : Fin 11, BlockCheck b
  | ⟨0, _⟩ => blockCheck_0
  | ⟨1, _⟩ => blockCheck_1
  | ⟨2, _⟩ => blockCheck_2
  | ⟨3, _⟩ => blockCheck_3
  | ⟨4, _⟩ => blockCheck_4
  | ⟨5, _⟩ => blockCheck_5
  | ⟨6, _⟩ => blockCheck_6
  | ⟨7, _⟩ => blockCheck_7
  | ⟨8, _⟩ => blockCheck_8
  | ⟨9, _⟩ => blockCheck_9
  | ⟨10, _⟩ => blockCheck_10

end SeedKernel.Trial
