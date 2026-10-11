/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Hamiltonian

/-!
# Trotter slices of the Heisenberg ring

The finite approximants used to prove the spatial moment inequalities. One Trotter slice of
width `τ` is the ordered product of the first-order bond factors `1 - τ hⱼ`, with
`hⱼ = Sⱼ · Sⱼ₊₁`: first over the even bonds `j`, then over the odd bonds.

* `bond L j` is `hⱼ`, so `hamiltonian L = ∑ⱼ bond L j`.
* `bondOrder L` lists the bonds, even ones first.
* `trotterSlice L τ` is `∏_{j ∈ bondOrder L} (1 - τ hⱼ)`.
* `trotterTrace M β L` is `Re Tr (trotterSlice L (β/M))^M`, the approximant of
  `partitionFunction L β`.
-/

open Matrix

namespace SpinOneHeisenberg

/-- The bond term `hⱼ = Sⱼ · Sⱼ₊₁` of the ring. -/
def bond (L : ℕ) (j : Fin L) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  bondOp j (finRotate L j) exchange

theorem hamiltonian_eq_sum_bond (L : ℕ) : hamiltonian L = ∑ j, bond L j := rfl

/-- The bonds of the ring in Trotter order: even bonds, then odd bonds. -/
def bondOrder (L : ℕ) : List (Fin L) :=
  (List.finRange L).filter (fun j ↦ decide (Even j.val)) ++
    (List.finRange L).filter (fun j ↦ !decide (Even j.val))

/-- One first-order Trotter slice `∏_{j ∈ bondOrder L} (1 - τ hⱼ)`. -/
def trotterSlice (L : ℕ) (τ : ℝ) : Matrix (Fin L → Fin 3) (Fin L → Fin 3) ℂ :=
  ((bondOrder L).map fun j ↦ 1 - (τ : ℂ) • bond L j).prod

/-- The Trotter approximant `Re Tr (trotterSlice L (β/M))^M` of `Z_L(β)`. -/
noncomputable def trotterTrace (M : ℕ) (β : ℝ) (L : ℕ) : ℝ :=
  (trace (trotterSlice L (β / M) ^ M)).re

end SpinOneHeisenberg
