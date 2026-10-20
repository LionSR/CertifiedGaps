/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import SeedKernel.Heisenberg.TrialBlocks
import SeedKernel.Heisenberg.TrialIntertwiner

/-!
# The trial-state energy bound of the spin-one seed (kernel-checked)

For the integer matrix-product trial state `trialA` (bond dimension `26`), with transfer matrix
`T = mpsTransfer trialA` and energy insertion `O = mpsInsertion trialA` (`TrialBound.lean`),

    tr (T^72) > 0,   500000 · tr (T^70 O) < -700741 · tr (T^72),

and the same with `120`, `118`. Since `⟨ψ|ψ⟩ = tr T^N` and `⟨ψ|Sⱼ·Sⱼ₊₁|ψ⟩ = tr (T^(N-2) O)`, the
energy per bond of the trial state on rings of `72` and `120` sites is below `-1.401482`.

The proof has three parts, all checked by the kernel (`decide +kernel`, no native code):

1. **Reduction** (`trace_trial_blocks`). An integer intertwiner `W` with `T W = W B`,
   `O W = W C` and `V W = c 1` (`TrialIntertwiner`, checked in `TrialChecks` and
   `TrialChecksInsertion1`–`TrialChecksInsertion4`) reduces the traces over the 676-dimensional
   doubled virtual space to eleven blocks `tᵦ` (spin `J ≤ 5` and swap sign), each counted
   `2J + 1` times (`TrialReduction`).
2. **Block powers** (`trace_block`). Packed big-integer products compute `tr (tᵦ^N)` and
   `tr (tᵦ^(N-2) oᵦ)` exactly (`TrialPowers`, `TrialBlocks`).
3. **Arithmetic.** The weighted sums equal the stored integers `trialV72`, `trialU72`, … (the
   upstream certificate values), which satisfy the inequalities.
-/

namespace CertifiedGaps.SpinOneHeisenberg

open Matrix SeedKernel.Trial

/-- The recipe matrices are the matrices written by the Julia generator, which asserts them equal
to the upstream evidence `matrix-data.json`. -/
theorem trialA_eq_code :
    trialA 0 = toMatB 26 26 trialMatrixCode0.1 (decodeSquare 26 trialMatrixCode0) ∧
      trialA 1 = toMatB 26 26 trialMatrixCode1.1 (decodeSquare 26 trialMatrixCode1) ∧
      trialA 2 = toMatB 26 26 trialMatrixCode2.1 (decodeSquare 26 trialMatrixCode2) := by
  obtain ⟨h0, h1, h2⟩ := trialRows_eq_code
  refine ⟨?_, ?_, ?_⟩
  · rw [← toMat_trialRows, h0]; rfl
  · rw [← toMat_trialRows, h1]; rfl
  · rw [← toMat_trialRows, h2]; rfl

theorem trialScale_ne_zero : (trialScale : ℤ) ≠ 0 := by decide

/-- **Reduction to the blocks.** -/
theorem trace_trial_blocks (k : ℕ) :
    trace (mpsTransfer trialA ^ k) = ∑ b : Fin 11, (bd b : ℤ) * trace (tMat b ^ k) ∧
      trace (mpsTransfer trialA ^ k * mpsInsertion trialA) =
        ∑ b : Fin 11, (bd b : ℤ) * trace (tMat b ^ k * oMat b) := by
  have hB : ∀ b : Fin 11, (tB b).Good (bm b) := fun b ↦ (blockCheck_all b).1
  have hC : ∀ b : Fin 11, (oB b).Good (bm b) := fun b ↦ (blockCheck_all b).2.1
  have hT := T_mul_Wm wIdxOk wAbsOk aBoundOk hB tB_offset checkT
  have hO := O_mul_Wm wIdxOk wAbsOk aBoundOk hC oB_offset checkO1 checkO2 checkO3 checkO4
  have hV := Vm_mul_Wm wIdxOk wAbsOk vIdxOk vAbsOk checkV
  obtain ⟨h1, h2⟩ := trace_eq_of_intertwiner _ _ Wm Vm Bm Cm _ trialScale_ne_zero card_Iota
    hT hO hV k
  obtain ⟨h3, h4⟩ := trace_copies_pow (d := bd) (m := bm) tMat oMat k
  exact ⟨h1.trans h3, h2.trans h4⟩

/-- **The block traces** equal the stored integers. -/
theorem trace_block (b : Fin 11) :
    trace (tMat b ^ 72) = (storedTraces b).1 ∧ trace (tMat b ^ 70 * oMat b) = (storedTraces b).2.1 ∧
      trace (tMat b ^ 120) = (storedTraces b).2.2.1 ∧
        trace (tMat b ^ 118 * oMat b) = (storedTraces b).2.2.2 := by
  obtain ⟨ht, ho, hc⟩ := blockCheck_all b
  have h := blockTraces_spec ht ho
  rw [hc] at h
  refine ⟨?_, ?_, ?_, ?_⟩ <;> rw [h] <;> rfl

theorem trace_T_72 : trace (mpsTransfer trialA ^ 72) = trialV72 := by
  have h (b : Fin 11) : (bd b : ℤ) * trace (tMat b ^ 72) = (bd b : ℤ) * (storedTraces b).1 := by
    rw [(trace_block b).1]
  rw [(trace_trial_blocks 72).1, Finset.sum_congr rfl fun b _ ↦ h b]
  decide +kernel

theorem trace_TO_70 : trace (mpsTransfer trialA ^ 70 * mpsInsertion trialA) = trialU72 := by
  have h (b : Fin 11) :
      (bd b : ℤ) * trace (tMat b ^ 70 * oMat b) = (bd b : ℤ) * (storedTraces b).2.1 := by
    rw [(trace_block b).2.1]
  rw [(trace_trial_blocks 70).2, Finset.sum_congr rfl fun b _ ↦ h b]
  decide +kernel

theorem trace_T_120 : trace (mpsTransfer trialA ^ 120) = trialV120 := by
  have h (b : Fin 11) :
      (bd b : ℤ) * trace (tMat b ^ 120) = (bd b : ℤ) * (storedTraces b).2.2.1 := by
    rw [(trace_block b).2.2.1]
  rw [(trace_trial_blocks 120).1, Finset.sum_congr rfl fun b _ ↦ h b]
  decide +kernel

theorem trace_TO_118 : trace (mpsTransfer trialA ^ 118 * mpsInsertion trialA) = trialU120 := by
  have h (b : Fin 11) :
      (bd b : ℤ) * trace (tMat b ^ 118 * oMat b) = (bd b : ℤ) * (storedTraces b).2.2.2 := by
    rw [(trace_block b).2.2.2]
  rw [(trace_trial_blocks 118).2, Finset.sum_congr rfl fun b _ ↦ h b]
  decide +kernel

/-- **Trial-state bound at `N = 72`.** -/
theorem trial_check_72 :
    0 < trace (mpsTransfer trialA ^ 72) ∧
    500000 * trace (mpsTransfer trialA ^ 70 * mpsInsertion trialA) <
      -700741 * trace (mpsTransfer trialA ^ 72) := by
  rw [trace_T_72, trace_TO_70]
  decide +kernel

/-- **Trial-state bound at `N = 120`.** -/
theorem trial_check_120 :
    0 < trace (mpsTransfer trialA ^ 120) ∧
    500000 * trace (mpsTransfer trialA ^ 118 * mpsInsertion trialA) <
      -700741 * trace (mpsTransfer trialA ^ 120) := by
  rw [trace_T_120, trace_TO_118]
  decide +kernel

end CertifiedGaps.SpinOneHeisenberg
