/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryData
import CertifiedGaps.Heisenberg.BoundaryInit
import CertifiedGaps.Trace.Frobenius

/-!
# The fifteen short-chain thermal traces of the boundary-field paper

The appendix "Certified short-chain thermal moments" of OpenAI's *A boundary-field gap for the
spin-one Heisenberg chain* (2026) initializes the boundary argument from fifteen thermal traces
of open spin-one chains of `n = 2, …, 6` sites at `b₀ = 49/4`:

  `H_{n,g} = ∑_{j < n-1} Sⱼ · Sⱼ₊₁ - (3/5) S₀ᶻ + (3/5) R_g`,  `R_{±1} = ± S_{n-1}ᶻ`,
  `R_0 = S_{n-1}ˣ`,  `Z_{n,g} = Tr exp[-b₀ (H_{n,g} + a n)]`,  `a = 700741/500000`.

This file defines these objects, relates them to the opposite-field chain of
`Heisenberg/OpenChain.lean`, states the rational enclosures of the fifteen traces as the
proposition `BoundaryTraceTable`, and derives from it the two conclusions of the paper's
Proposition appth:bounds that do not involve the transfer operator:

* `sixth_of_table`: `Y₆(b₀) < .027`, in the form of `BoundaryInit.ThermalData.sixth`;
* `moment_polynomial_lt_of_table`: `I(6) - .28 I(5) - .3404 I(4) + .0504 I(3) + .0324 I(2)
  < .003745`, with `I(n) = (Z_{n,1} + Z_{n,-1} + 4 Z_{n,0})/6` (eq. (appth:polynomial)).

`Heisenberg/BoundaryThermal.lean` proves `BoundaryTraceTable` from the paper's integer table
(`boundaryTraceTable_of_totals`). The reduction of the leading-weight bound
`c₀ P_*(λ₀)² ≤ I(6) - .28 I(5) - …` (the two-boundary transfer identity and the averaging over
the rotation group, eqs. (appth:average), (appth:polynomial)) is not part of these files.

## Conventions

The spin matrices are those of `Local/SpinOne.lean` (Cartesian basis, `S 0 = Sˣ`,
`S 2 = Sᶻ`), as in `endFieldHamiltonian`. The right field is indexed by an integer `g`; only
`g ∈ {1, -1, 0}` is used. With these conventions `H_{n,1}` is the opposite-field chain
`endFieldHamiltonian (3/5) (-(3/5)) n` and `H_{n,-1}` the equal-field chain, so that
`Z_{n,1} = Y_n(b₀)` is the shifted opposite-field trace (`trace_one`), as the paper states.

## Differences from the paper

* The paper prints only the integer table (appth:integer-table), the enclosure
  `.0266319440 < Z_{6,1} < .0266319607` and the five centers of `I(n)`. The fifteen intervals of
  `traceLo`, `traceHi` are the paper's rational enclosures (appth:partition-interval), with the
  exponential enclosed by a Taylor polynomial of degree `40` instead of `80`, rounded outward to
  ten significant digits. They imply the paper's enclosure of `Z_{6,1}` (`trace_six_one_mem`)
  and the polynomial bound `≤ .0037440319` (paper: `< .0037442021`).
-/

open Matrix

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryTraces

/-! ### The fifteen Hamiltonians and traces -/

/-- The right-boundary operator `R_g`: `g Sᶻ` for `g ≠ 0` and `Sˣ` for `g = 0`. -/
noncomputable def rightField (g : ℤ) : Matrix (Fin 3) (Fin 3) ℂ :=
  if g = 0 then SpinOne.S 0 else (g : ℂ) • SpinOne.S 2

/-- The boundary Hamiltonian `H_{n,g} = ∑_{j < n-1} Sⱼ · Sⱼ₊₁ - (3/5) S₀ᶻ + (3/5) R_g` on `n`
sites (`0` for `n = 0`). -/
noncomputable def hamiltonian (g : ℤ) : (n : ℕ) → Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℂ
  | 0 => 0
  | n + 1 => ∑ j : Fin n, bondOp j.castSucc j.succ SpinOne.exchange +
      siteOp 0 (-((3 / 5 : ℝ) : ℂ) • SpinOne.S 2) +
        siteOp (Fin.last n) (((3 / 5 : ℝ) : ℂ) • rightField g)

/-- The shift `a = 700741/500000` per site. -/
noncomputable def shift : ℝ := 700741 / 500000

/-- The shifted thermal trace `Z_{n,g} = Tr exp[-b₀ (H_{n,g} + a n)]` at `b₀ = 49/4`. -/
noncomputable def trace (n : ℕ) (g : ℤ) : ℝ :=
  gibbsTrace (hamiltonian g n + ((shift * n : ℝ) : ℂ) • 1) (49 / 4)

/-- The right-field average `I(n) = (Z_{n,1} + Z_{n,-1} + 4 Z_{n,0})/6` (eq. (appth:average)). -/
noncomputable def moment (n : ℕ) : ℝ := (trace n 1 + trace n (-1) + 4 * trace n 0) / 6

/-- `H_{n,1}` is the opposite-field chain. -/
theorem hamiltonian_one (n : ℕ) :
    hamiltonian 1 n = endFieldHamiltonian (3 / 5) (-(3 / 5)) n := by
  cases n with
  | zero => rfl
  | succ n =>
    simp only [hamiltonian, endFieldHamiltonian, rightField, one_ne_zero, ite_false,
      Int.cast_one, one_smul]
    congr 2
    push_cast; rw [neg_neg]

/-- `H_{n,-1}` is the equal-field chain. -/
theorem hamiltonian_neg_one (n : ℕ) :
    hamiltonian (-1) n = endFieldHamiltonian (3 / 5) (3 / 5) n := by
  cases n with
  | zero => rfl
  | succ n =>
    simp only [hamiltonian, endFieldHamiltonian, rightField, neg_eq_zero, one_ne_zero,
      ite_false, Int.cast_neg, Int.cast_one, neg_smul, one_smul, smul_neg]

/-- A scalar shift multiplies the trace by `e^{-βa n}`. -/
theorem trace_eq (n : ℕ) (g : ℤ) :
    trace n g = BoundaryInit.shiftFactor (49 / 4) ^ n * gibbsTrace (hamiltonian g n) (49 / 4) := by
  rw [trace, Trace.gibbsTrace_add_smul_one, BoundaryInit.shiftFactor_pow, shift]
  congr 2; ring

/-- **`Z_{n,1} = Y_n(b₀)`**: the trace with `g = 1` is the shifted opposite-field trace. -/
theorem trace_one (n : ℕ) :
    trace n 1 = BoundaryInit.shiftFactor (49 / 4) ^ n *
      oppositeFieldPartitionFunction (3 / 5) n (49 / 4) := by
  rw [trace_eq, hamiltonian_one]
  rfl

/-- `Z_{n,-1}` is the shifted equal-field trace. -/
theorem trace_neg_one (n : ℕ) :
    trace n (-1) = BoundaryInit.shiftFactor (49 / 4) ^ n *
      gibbsTrace (endFieldHamiltonian (3 / 5) (3 / 5) n) (49 / 4) := by
  rw [trace_eq, hamiltonian_neg_one]

/-! ### The table -/

open BoundaryThermal (lengths fields)

/-- The lower ends of the enclosures of `Z_{n,g}`, `g = 1, -1, 0`. -/
def traceLo : ℕ → ℤ → ℚ
  | 2, g => if g = 1 then 1319009172 / 10 ^ 10 else if g = -1 then 5387274464 / 10 ^ 14
      else 7557416221 / 10 ^ 12
  | 3, g => if g = 1 then 2455126214 / 10 ^ 14 else if g = -1 then 5408387857 / 10 ^ 11
      else 4526316255 / 10 ^ 12
  | 4, g => if g = 1 then 4606108509 / 10 ^ 11 else if g = -1 then 9570975498 / 10 ^ 14
      else 6363464067 / 10 ^ 12
  | 5, g => if g = 1 then 2761131291 / 10 ^ 13 else if g = -1 then 3166071767 / 10 ^ 11
      else 5989192150 / 10 ^ 12
  | 6, g => if g = 1 then 2663194405 / 10 ^ 11 else if g = -1 then 5308152809 / 10 ^ 13
      else 6608376591 / 10 ^ 12
  | _, _ => 0

/-- The upper ends of the enclosures of `Z_{n,g}`, `g = 1, -1, 0`. -/
def traceHi : ℕ → ℤ → ℚ
  | 2, g => if g = 1 then 1319009206 / 10 ^ 10 else if g = -1 then 5387281149 / 10 ^ 14
      else 7557417013 / 10 ^ 12
  | 3, g => if g = 1 then 2455126600 / 10 ^ 14 else if g = -1 then 5408387876 / 10 ^ 11
      else 4526316309 / 10 ^ 12
  | 4, g => if g = 1 then 4606109163 / 10 ^ 11 else if g = -1 then 9571005288 / 10 ^ 14
      else 6363466496 / 10 ^ 12
  | 5, g => if g = 1 then 2761131724 / 10 ^ 13 else if g = -1 then 3166071814 / 10 ^ 11
      else 5989192353 / 10 ^ 12
  | 6, g => if g = 1 then 2663196067 / 10 ^ 11 else if g = -1 then 5308176266 / 10 ^ 13
      else 6608384864 / 10 ^ 12
  | _, _ => 0

/-- **The table of the fifteen traces**: `traceLo n g ≤ Z_{n,g} ≤ traceHi n g` for
`n = 2, …, 6` and `g = 1, -1, 0`. -/
def BoundaryTraceTable : Prop :=
  ∀ n ∈ lengths, ∀ g ∈ fields, (traceLo n g : ℝ) ≤ trace n g ∧ trace n g ≤ traceHi n g

/-! ### The two consequences -/

section Consequences

variable (htab : BoundaryTraceTable)
include htab

/-- The paper's enclosure `.0266319440 < Z_{6,1} < .0266319607` (eq. (appth:Y6-interval)). -/
theorem trace_six_one_mem :
    (266319440 / 10 ^ 10 : ℝ) < trace 6 1 ∧ trace 6 1 < 266319607 / 10 ^ 10 := by
  obtain ⟨h1, h2⟩ := htab 6 (by decide) 1 (by decide)
  norm_num [traceLo, traceHi] at h1 h2
  constructor <;> linarith

/-- **`Y₆(b₀) < .027`**, in the form of `BoundaryInit.ThermalData.sixth`. -/
theorem sixth_of_table :
    BoundaryInit.shiftFactor (49 / 4) ^ 6 * oppositeFieldPartitionFunction (3 / 5) 6 (49 / 4) <
      27 / 1000 := by
  rw [← trace_one]
  have := (trace_six_one_mem htab).2
  linarith

/-- **The leading-weight polynomial** (eq. (appth:polynomial)):
`I(6) - .28 I(5) - .3404 I(4) + .0504 I(3) + .0324 I(2) < .003745`. -/
theorem moment_polynomial_lt_of_table :
    moment 6 - 28 / 100 * moment 5 - 3404 / 10000 * moment 4 + 504 / 10000 * moment 3 +
      324 / 10000 * moment 2 < 3745 / 10 ^ 6 := by
  have e : ∀ n ∈ lengths, ∀ g ∈ fields,
      (traceLo n g : ℝ) ≤ trace n g ∧ trace n g ≤ traceHi n g := htab
  obtain ⟨a1, a2⟩ := e 2 (by decide) 1 (by decide)
  obtain ⟨a3, a4⟩ := e 2 (by decide) (-1) (by decide)
  obtain ⟨a5, a6⟩ := e 2 (by decide) 0 (by decide)
  obtain ⟨b1, b2⟩ := e 3 (by decide) 1 (by decide)
  obtain ⟨b3, b4⟩ := e 3 (by decide) (-1) (by decide)
  obtain ⟨b5, b6⟩ := e 3 (by decide) 0 (by decide)
  obtain ⟨c1, c2⟩ := e 4 (by decide) 1 (by decide)
  obtain ⟨c3, c4⟩ := e 4 (by decide) (-1) (by decide)
  obtain ⟨c5, c6⟩ := e 4 (by decide) 0 (by decide)
  obtain ⟨d1, d2⟩ := e 5 (by decide) 1 (by decide)
  obtain ⟨d3, d4⟩ := e 5 (by decide) (-1) (by decide)
  obtain ⟨d5, d6⟩ := e 5 (by decide) 0 (by decide)
  obtain ⟨f1, f2⟩ := e 6 (by decide) 1 (by decide)
  obtain ⟨f3, f4⟩ := e 6 (by decide) (-1) (by decide)
  obtain ⟨f5, f6⟩ := e 6 (by decide) 0 (by decide)
  norm_num [traceLo, traceHi] at a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 b6
  norm_num [traceLo, traceHi] at c1 c2 c3 c4 c5 c6 d1 d2 d3 d4 d5 d6
  norm_num [traceLo, traceHi] at f1 f2 f3 f4 f5 f6
  simp only [moment]
  linarith

end Consequences

end BoundaryTraces

end CertifiedGaps.SpinOneHeisenberg
