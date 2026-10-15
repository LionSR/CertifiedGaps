/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Hamiltonian

/-!
# Real spin-one factors

In the Cartesian basis `(S^α)_{s s'} = -i ε_{α s s'}` of `SpinOneHeisenberg.spinOne`, the
matrices `G_α = i S^α = ε_α` are real and antisymmetric, and `-S · S = ∑_α G_α ⊗ G_α`. They are
the bond factors of every spin-one model (`SpinOneLargeD`, `SpinOneDimerized`).
-/

open Matrix

namespace CertifiedGaps.SpinOne

/-- The real bond factors `G_α = ε_α = i S^α`, with entries `ε_{α s s'}`. -/
def G (α : Fin 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  of fun s s' ↦ (SpinOneHeisenberg.leviCivita α s s' : ℝ)

/-- The factors `G_α` are antisymmetric. -/
theorem G_transpose (α : Fin 3) : (G α)ᵀ = -G α := by
  ext s s'; simp [G, SpinOneHeisenberg.leviCivita_swap α s s']

end CertifiedGaps.SpinOne
