/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialData

/-!
# Traces of powers of one reduced block

For one block `t` (size `n`) with insertion `o`, the four numbers

    tr (t^72),  tr (t^70 o),  tr (t^120),  tr (t^118 o)

are computed by repeated squaring with the packed products of `SeedKernel.Trial.lmulB`:
`t², t⁴, …, t³²`, then `t³⁶ = t³² t⁴`, `t²o`, `t³⁴o = t³² (t²o)`, `t⁴⁸, t⁵⁶, t⁶⁰` and
`t⁵⁸o = t⁵⁶ (t²o)`, and each trace is read from the diagonal of one more product
(`tr (t^72) = tr (t³⁶ t³⁶)`, `tr (t^70 o) = tr (t³⁶ (t³⁴o))`, and so on).

`blockTraces_spec` proves, for every block satisfying the shape and offset bounds `Good`, that
these computed numbers are the Mathlib traces of the Mathlib matrix powers. The kernel only
evaluates the chain, once per block (`TrialBlocks*.lean`).
-/

namespace SeedKernel.Trial

open Matrix

/-- An integer matrix stored as `(β, rows)`: entries `x - β` for the naturals `x` of `rows`. -/
abbrev OM := ℕ × List (List ℕ)

/-- The `n × n` integer matrix of an offset matrix. -/
def OM.mat (n : ℕ) (X : OM) : Matrix (Fin n) (Fin n) ℤ := toMatB n n X.1 X.2

/-- Product of `n × n` offset matrices; the offset of the result is `n β_X β_Y`. -/
def omul (n : ℕ) (X Y : OM) : OM := (n * X.1 * Y.1, lmulB X.1 Y.1 X.2 Y.2 n)

/-- `tr (X Y)` from the diagonal of the packed product. -/
def otrace (n : ℕ) (X Y : OM) : ℤ := (ltraceB X.1 Y.1 X.2 Y.2 n : ℤ) - n * (n * X.1 * Y.1)

/-- Shape and offset bound of an `n × n` offset matrix: `n ≥ 1` rows of length `n`, entries at
most twice the offset (so that every stored entry lies in `[-β, β]`). -/
def OM.Good (n : ℕ) (X : OM) : Prop :=
  X.2.length = n ∧ X.2 ≠ [] ∧ ∀ r ∈ X.2, r.length = n ∧ ∀ x ∈ r, x ≤ 2 * X.1

instance (n : ℕ) (X : OM) : Decidable (X.Good n) := by unfold OM.Good; infer_instance

theorem omul_good {n : ℕ} {X Y : OM} (hX : X.Good n) (hY : Y.Good n) : (omul n X Y).Good n := by
  obtain ⟨hXn, hXne, hXr⟩ := hX
  obtain ⟨hYn, hYne, hYr⟩ := hY
  refine ⟨by simp [omul, length_lmulB, hXn], ?_, ?_⟩
  · simp only [omul]
    intro h
    have := length_lmulB X.1 Y.1 n X.2 Y.2
    rw [h, hXn] at this
    exact hXne (List.length_eq_zero_iff.mp (by rw [hXn]; exact this.symm))
  · have := lmulB_shape X.1 Y.1 n X.2 Y.2 hYne (by rw [hYn]; exact hXr) hYr
    simpa [omul, hYn] using this

theorem omul_mat {n : ℕ} {X Y : OM} (hX : X.Good n) (hY : Y.Good n) :
    (omul n X Y).mat n = X.mat n * Y.mat n :=
  toMatB_lmulB X.1 Y.1 X.2 Y.2 hY.1 hY.2.1 hX.2.2 hY.2.2

theorem otrace_eq {n : ℕ} {X Y : OM} (hX : X.Good n) (hY : Y.Good n) :
    otrace n X Y = trace (X.mat n * Y.mat n) :=
  (trace_toMatB_mul X.1 Y.1 X.2 Y.2 hX.1 hY.1 hY.2.1 hX.2.2 hY.2.2).symm

section Chain

variable (n : ℕ) (t o : OM)

/-- `t²`. -/
def P2 : OM := omul n t t
/-- `t⁴`. -/
def P4 : OM := omul n (P2 n t) (P2 n t)
/-- `t⁸`. -/
def P8 : OM := omul n (P4 n t) (P4 n t)
/-- `t¹⁶`. -/
def P16 : OM := omul n (P8 n t) (P8 n t)
/-- `t³²`. -/
def P32 : OM := omul n (P16 n t) (P16 n t)
/-- `t³⁶`. -/
def P36 : OM := omul n (P32 n t) (P4 n t)
/-- `t⁴⁸`. -/
def P48 : OM := omul n (P32 n t) (P16 n t)
/-- `t⁵⁶`. -/
def P56 : OM := omul n (P48 n t) (P8 n t)
/-- `t⁶⁰`. -/
def P60 : OM := omul n (P56 n t) (P4 n t)
/-- `t² o`. -/
def Q2 : OM := omul n (P2 n t) o
/-- `t³⁴ o`. -/
def Q34 : OM := omul n (P32 n t) (Q2 n t o)
/-- `t⁵⁸ o`. -/
def Q58 : OM := omul n (P56 n t) (Q2 n t o)

/-- `(tr (t^72), tr (t^70 o), tr (t^120), tr (t^118 o))`, computed by packed products. -/
def blockTraces : ℤ × ℤ × ℤ × ℤ :=
  (otrace n (P36 n t) (P36 n t), otrace n (P36 n t) (Q34 n t o),
    otrace n (P60 n t) (P60 n t), otrace n (P60 n t) (Q58 n t o))

end Chain

/-- **The computed traces are the traces of the matrix powers.** -/
theorem blockTraces_spec {n : ℕ} {t o : OM} (ht : t.Good n) (ho : o.Good n) :
    blockTraces n t o = (trace (t.mat n ^ 72), trace (t.mat n ^ 70 * o.mat n),
      trace (t.mat n ^ 120), trace (t.mat n ^ 118 * o.mat n)) := by
  -- Goodness and the matrix of every stage of the chain.
  have g2 : (P2 n t).Good n := omul_good ht ht
  have g4 : (P4 n t).Good n := omul_good g2 g2
  have g8 : (P8 n t).Good n := omul_good g4 g4
  have g16 : (P16 n t).Good n := omul_good g8 g8
  have g32 : (P32 n t).Good n := omul_good g16 g16
  have g36 : (P36 n t).Good n := omul_good g32 g4
  have g48 : (P48 n t).Good n := omul_good g32 g16
  have g56 : (P56 n t).Good n := omul_good g48 g8
  have g60 : (P60 n t).Good n := omul_good g56 g4
  have q2 : (Q2 n t o).Good n := omul_good g2 ho
  have q34 : (Q34 n t o).Good n := omul_good g32 q2
  have q58 : (Q58 n t o).Good n := omul_good g56 q2
  set T := t.mat n
  have m2 : (P2 n t).mat n = T ^ 2 := by rw [P2, omul_mat ht ht, sq]
  have m4 : (P4 n t).mat n = T ^ 4 := by rw [P4, omul_mat g2 g2, m2, ← pow_add]
  have m8 : (P8 n t).mat n = T ^ 8 := by rw [P8, omul_mat g4 g4, m4, ← pow_add]
  have m16 : (P16 n t).mat n = T ^ 16 := by rw [P16, omul_mat g8 g8, m8, ← pow_add]
  have m32 : (P32 n t).mat n = T ^ 32 := by rw [P32, omul_mat g16 g16, m16, ← pow_add]
  have m36 : (P36 n t).mat n = T ^ 36 := by rw [P36, omul_mat g32 g4, m32, m4, ← pow_add]
  have m48 : (P48 n t).mat n = T ^ 48 := by rw [P48, omul_mat g32 g16, m32, m16, ← pow_add]
  have m56 : (P56 n t).mat n = T ^ 56 := by rw [P56, omul_mat g48 g8, m48, m8, ← pow_add]
  have m60 : (P60 n t).mat n = T ^ 60 := by rw [P60, omul_mat g56 g4, m56, m4, ← pow_add]
  have n2 : (Q2 n t o).mat n = T ^ 2 * o.mat n := by rw [Q2, omul_mat g2 ho, m2]
  have n34 : (Q34 n t o).mat n = T ^ 34 * o.mat n := by
    rw [Q34, omul_mat g32 q2, m32, n2, ← mul_assoc, ← pow_add]
  have n58 : (Q58 n t o).mat n = T ^ 58 * o.mat n := by
    rw [Q58, omul_mat g56 q2, m56, n2, ← mul_assoc, ← pow_add]
  simp only [blockTraces, otrace_eq g36 g36, otrace_eq g36 q34, otrace_eq g60 g60,
    otrace_eq g60 q58, m36, n34, m60, n58, ← mul_assoc, ← pow_add]

/-- Block `b` of the reduced transfer matrix, `tᵦ`, decoded from `TrialData`. -/
def tOM (b : ℕ) : OM :=
  ((tCode.getD b (0, 0, 0)).1, decodeSquare (blockDim.getD b 0) (tCode.getD b (0, 0, 0)))

/-- Block `b` of the reduced energy insertion, `oᵦ`, decoded from `TrialData`. -/
def oOM (b : ℕ) : OM :=
  ((oCode.getD b (0, 0, 0)).1, decodeSquare (blockDim.getD b 0) (oCode.getD b (0, 0, 0)))

/-- The stored traces of block `b`, in the order of `blockTraces`. -/
def storedTraces (b : ℕ) : ℤ × ℤ × ℤ × ℤ :=
  ((blockTrace72.getD b (0, 0)).1, (blockTrace72.getD b (0, 0)).2,
    (blockTrace120.getD b (0, 0)).1, (blockTrace120.getD b (0, 0)).2)

/-- The kernel obligation of block `b`: the data are well formed and the packed chain
reproduces the stored traces. -/
def BlockCheck (b : ℕ) : Prop :=
  (tOM b).Good (blockDim.getD b 0) ∧ (oOM b).Good (blockDim.getD b 0) ∧
    blockTraces (blockDim.getD b 0) (tOM b) (oOM b) = storedTraces b

instance (b : ℕ) : Decidable (BlockCheck b) := by unfold BlockCheck; infer_instance

end SeedKernel.Trial
