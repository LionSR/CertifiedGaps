/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.Data.Rat.Defs

/-!
# Endpoint arithmetic of the Ising seed

The four rational log-partition endpoints of `SeedBridge.LogEnclosures` and the two exact
inequalities between them that the seed bridge needs.
-/

namespace SeedKernel.Ising

def logSmallUpper : ℚ := 2726798329802173 / 17592186044416
def logLargeLower : ℚ := 5453596500152007 / 17592186044416
def logLargeUpper : ℚ := 2726798388197837 / 8796093022208
def logDoubleLower : ℚ := 2726798250076001 / 4398046511104

-- These check the endpoint arithmetic, not the claims that these endpoints enclose log Z.
theorem spatial_arithmetic :
    -(9064 / 10^9 : ℚ) ≤ logLargeLower - 2 * logSmallUpper := by decide +kernel
theorem thermal_arithmetic :
    -(3141 / 10^8 : ℚ) ≤ logDoubleLower - 2 * logLargeUpper := by decide +kernel

#print axioms spatial_arithmetic
#print axioms thermal_arithmetic

end SeedKernel.Ising
