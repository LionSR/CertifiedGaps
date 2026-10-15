import Mathlib.Data.Rat.Defs

namespace LeanSeedSmall

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

end LeanSeedSmall
