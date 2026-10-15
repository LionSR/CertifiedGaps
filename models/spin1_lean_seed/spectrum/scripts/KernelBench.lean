import Lean

set_option maxRecDepth 1000000
set_option maxHeartbeats 0

-- These are kernel checks, not native_decide. Run each size in a separate process
-- via kernel_bench.jl to measure startup separately and avoid timing elaboration
-- of decimal literals. Dense patterned operands prevent a power-of-two shortcut.
def packedA (bits : Nat) : Nat := ((1 <<< bits) - 1) / 3
def packedB (bits : Nat) : Nat := ((1 <<< bits) - 1) / 5

example : let a := packedA 4000000; let b := packedB 4000000
    (a * b) / a = b := by decide +kernel

#print axioms packedA
