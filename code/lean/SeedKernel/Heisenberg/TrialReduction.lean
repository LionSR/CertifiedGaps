/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.Data.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.SemiringInverse
import Mathlib.Data.Rat.Cast.Defs

/-!
# Trace reduction by an integer intertwiner

If an integer matrix `W` intertwines `T` and `O` with `B` and `C` (`T W = W B`, `O W = W C`) and
has an integer left inverse up to a nonzero scalar (`V W = c 1`, square), then
`tr (T^k) = tr (B^k)` and `tr (T^k O) = tr (B^k C)`. With `B` and `C` block diagonal with blocks
`1 ⊗ tᵦ` and `1 ⊗ oᵦ`, the traces become `∑ᵦ dᵦ tr (tᵦ^k)` and `∑ᵦ dᵦ tr (tᵦ^k oᵦ)`.
-/

namespace SeedKernel.Trial

open Matrix
open scoped Kronecker

section Intertwiner

variable {I ι : Type*} [Fintype I] [Fintype ι] [DecidableEq I] [DecidableEq ι]

theorem pow_mul_intertwiner {T : Matrix I I ℚ} {W : Matrix I ι ℚ} {B : Matrix ι ι ℚ}
    (h : T * W = W * B) : ∀ k : ℕ, T ^ k * W = W * B ^ k
  | 0 => by simp
  | k + 1 => by
    rw [pow_succ, Matrix.mul_assoc, h, ← Matrix.mul_assoc, pow_mul_intertwiner h k,
      Matrix.mul_assoc, ← pow_succ]

/-- **Trace reduction.** An integer intertwiner with an integer left inverse up to a nonzero
scalar transports the traces of `T^k` and `T^k O` to those of `B^k` and `B^k C`. -/
theorem trace_eq_of_intertwiner (T O : Matrix I I ℤ) (W : Matrix I ι ℤ) (V : Matrix ι I ℤ)
    (B C : Matrix ι ι ℤ) (c : ℤ) (hc : c ≠ 0) (hcard : Fintype.card ι = Fintype.card I)
    (hT : T * W = W * B) (hO : O * W = W * C) (hV : V * W = c • 1) (k : ℕ) :
    trace (T ^ k) = trace (B ^ k) ∧ trace (T ^ k * O) = trace (B ^ k * C) := by
  set f := Int.castRingHom ℚ
  have hTq : T.map f * W.map f = W.map f * B.map f := by rw [← Matrix.map_mul, ← Matrix.map_mul, hT]
  have hOq : O.map f * W.map f = W.map f * C.map f := by rw [← Matrix.map_mul, ← Matrix.map_mul, hO]
  set V' : Matrix ι I ℚ := (c : ℚ)⁻¹ • V.map f
  have hcq : (c : ℚ) ≠ 0 := by exact_mod_cast hc
  have hVW : V' * W.map f = 1 := by
    have : (V * W).map f = (c : ℚ) • (1 : Matrix ι ι ℚ) := by
      ext i j
      rw [hV]
      simp only [Matrix.map_apply, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
      split_ifs <;> simp [f]
    rw [Matrix.smul_mul, ← Matrix.map_mul, this, smul_smul, inv_mul_cancel₀ hcq, one_smul]
  have hWV : W.map f * V' = 1 := (mul_eq_one_comm_of_card_eq ι I ℚ hcard).mp hVW
  have htr {X : Matrix I I ℤ} : ((trace X : ℤ) : ℚ) = trace (X.map f) :=
    (AddMonoidHom.map_trace (Int.castRingHom ℚ).toAddMonoidHom X)
  have htr' {X : Matrix ι ι ℤ} : ((trace X : ℤ) : ℚ) = trace (X.map f) :=
    (AddMonoidHom.map_trace (Int.castRingHom ℚ).toAddMonoidHom X)
  have hpowT : (T ^ k).map f = T.map f ^ k := Matrix.map_pow T f k
  have hpowB : (B ^ k).map f = B.map f ^ k := Matrix.map_pow B f k
  have key := pow_mul_intertwiner hTq k
  constructor
  · have : trace (T.map f ^ k) = trace (B.map f ^ k) := by
      calc trace (T.map f ^ k) = trace (T.map f ^ k * (W.map f * V')) := by rw [hWV, mul_one]
        _ = trace ((W.map f * B.map f ^ k) * V') := by rw [← Matrix.mul_assoc, key]
        _ = trace (V' * (W.map f * B.map f ^ k)) := trace_mul_comm _ _
        _ = trace (B.map f ^ k) := by rw [← Matrix.mul_assoc, hVW, Matrix.one_mul]
    exact_mod_cast (by rw [htr, htr', hpowT, hpowB, this] :
      ((trace (T ^ k) : ℤ) : ℚ) = ((trace (B ^ k) : ℤ) : ℚ))
  · have : trace (T.map f ^ k * O.map f) = trace (B.map f ^ k * C.map f) := by
      calc trace (T.map f ^ k * O.map f)
          = trace (T.map f ^ k * O.map f * (W.map f * V')) := by rw [hWV, mul_one]
        _ = trace ((W.map f * (B.map f ^ k * C.map f)) * V') := by
            rw [← Matrix.mul_assoc, Matrix.mul_assoc _ (O.map f), hOq, ← Matrix.mul_assoc, key,
              Matrix.mul_assoc, Matrix.mul_assoc, Matrix.mul_assoc, Matrix.mul_assoc]
        _ = trace (V' * (W.map f * (B.map f ^ k * C.map f))) := trace_mul_comm _ _
        _ = trace (B.map f ^ k * C.map f) := by rw [← Matrix.mul_assoc, hVW, Matrix.one_mul]
    have e1 : (T ^ k * O).map f = T.map f ^ k * O.map f := by rw [Matrix.map_mul, hpowT]
    have e2 : (B ^ k * C).map f = B.map f ^ k * C.map f := by rw [Matrix.map_mul, hpowB]
    exact_mod_cast (by rw [htr, htr', e1, e2, this] :
      ((trace (T ^ k * O) : ℤ) : ℚ) = ((trace (B ^ k * C) : ℤ) : ℚ))

end Intertwiner

section Blocks

variable {β : Type*} [Fintype β] [DecidableEq β] {d m : β → ℕ}

/-- The block-diagonal matrix with blocks `1 ⊗ tᵦ` (`dᵦ` copies of `tᵦ`). -/
def copies (d m : β → ℕ) (t : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ) :
    Matrix (Σ b, Fin (d b) × Fin (m b)) (Σ b, Fin (d b) × Fin (m b)) ℤ :=
  blockDiagonal' fun b ↦ (1 : Matrix (Fin (d b)) (Fin (d b)) ℤ) ⊗ₖ t b

theorem copies_mul (t o : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ) :
    copies d m t * copies d m o = copies d m fun b ↦ t b * o b := by
  rw [copies, copies, ← blockDiagonal'_mul, copies]
  congr 1; ext1 b
  rw [← mul_kronecker_mul, one_mul]

theorem copies_pow (t : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ) :
    ∀ k : ℕ, copies d m t ^ k = copies d m fun b ↦ t b ^ k
  | 0 => by
    rw [pow_zero, copies]
    simp only [pow_zero, one_kronecker_one]
    rw [← blockDiagonal'_one]
    rfl
  | k + 1 => by rw [pow_succ, copies_pow t k, copies_mul, ← funext fun b ↦ pow_succ (t b) k]

theorem trace_copies (t : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ) :
    trace (copies d m t) = ∑ b, (d b : ℤ) * trace (t b) := by
  rw [copies, trace_blockDiagonal']
  congr 1; ext b
  rw [trace_kronecker, trace_one, Fintype.card_fin]

/-- Traces of powers of the block matrices: `dᵦ` copies of each block. -/
theorem trace_copies_pow (t o : ∀ b, Matrix (Fin (m b)) (Fin (m b)) ℤ) (k : ℕ) :
    trace (copies d m t ^ k) = ∑ b, (d b : ℤ) * trace (t b ^ k) ∧
      trace (copies d m t ^ k * copies d m o) = ∑ b, (d b : ℤ) * trace (t b ^ k * o b) := by
  rw [copies_pow, copies_mul, trace_copies, trace_copies]
  exact ⟨rfl, rfl⟩

end Blocks

end SeedKernel.Trial
