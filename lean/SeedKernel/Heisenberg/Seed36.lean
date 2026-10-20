/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.Seed36
import SeedKernel.Heisenberg.Trial
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# The spin-one gap from the trace table alone

The kernel-checked trial state (`trial_check_72`: OpenAI's integer matrix-product state of bond
dimension 26 on 72 sites has `500000 Tr (T⁷⁰ O) < -700741 Tr T⁷²`) discharges the trial-state input
of the seed at `(36, 49/4)`. What remains is the trace table at `n = 4, …, 12`
(`TraceTable`, Appendix A of the certificate note):

* `seed36_of_table`: the starting defects `1 - P_s(36, 49/4) ≤ 0.047916`,
  `1 - P_τ(72, 49/4) ≤ 0.181521` of `Main.lean`;
* `gap_of_table`: every even `L ≥ 36` has a unique ground state and
  `Δ_L > log(1.05882)/(49/4) ≈ 0.00467`;
* `gap_log_twenty_of_table`: OpenAI's eq. (1.4), `Δ_L > (log 20)/784` for every even `L ≥ 2304`.
-/

open Real Matrix

namespace CertifiedGaps.SpinOneHeisenberg

/-- The trial-state input holds: `e^{-βκ·72} Z₇₂(β) > 1` for every `β > 0`. -/
theorem one_lt_shiftedZ_seventytwo {β : ℝ} (hβ : 0 < β) : 1 < shiftedZ 72 β :=
  one_lt_shiftedZ_of_mps trialA trial_check_72.2 hβ

/-- **The seed at `(36, 49/4)` from the trace table.** -/
theorem seed36_of_table (htab : TraceTable) :
    1 - 47916 / 10 ^ 6 ≤ spatialPurity 36 (49 / 4) ∧
      1 - 181521 / 10 ^ 6 ≤ thermalPurity 72 (49 / 4) :=
  seed36_of_table_mps htab trialA trial_check_72.2

/-- **The gap for every even `L ≥ 36` from the trace table.** -/
theorem gap_of_table (htab : TraceTable) {L : ℕ} (hL : 36 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log (105882 / 100000) / (49 / 4)) :=
  gap_thirtysix_of_table_mps htab trialA trial_check_72.2 hL hLe

/-- **OpenAI's eq. (1.4) from the trace table**: `Δ_L > (log 20)/784` for every even `L ≥ 2304`,
which the sharper `gap_of_table` implies. -/
theorem gap_log_twenty_of_table (htab : TraceTable) {L : ℕ} (hL : 2304 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log 20 / 784) := by
  refine (gap_of_table htab (by omega) hLe).mono ?_
  rw [div_le_div_iff₀ (by norm_num) (by norm_num)]
  -- `log 20 ≤ 3` from `e > 2.718`, and `log 1.05882 ≥ 1/20` from `e^{1/20} < 1/(1 - 1/20)`.
  have h20 : log 20 ≤ 3 := by
    rw [log_le_iff_le_exp (by norm_num)]
    have e1 := Real.exp_one_gt_d9
    have : exp 3 = exp 1 ^ 3 := by rw [← exp_nat_mul]; norm_num
    rw [this]
    nlinarith [pow_le_pow_left₀ (by norm_num : (0 : ℝ) ≤ 2.7182818283) e1.le 3]
  have hr : (1 : ℝ) / 20 ≤ log (105882 / 100000) := by
    rw [le_log_iff_exp_le (by norm_num)]
    have := Real.exp_bound_div_one_sub_of_interval' (x := (1 : ℝ) / 20) (by norm_num)
      (by norm_num)
    linarith [show (1 : ℝ) / (1 - 1 / 20) ≤ 105882 / 100000 by norm_num]
  nlinarith

end CertifiedGaps.SpinOneHeisenberg
