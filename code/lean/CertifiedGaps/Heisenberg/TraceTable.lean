/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Scheme

/-!
# The short-ring trace table of the spin-one chain

The spin-one seeds start from a table of ring partition functions at `b = 49/4`: for
`n = 4, …, 12`, the traces `Z_n` and the traces `Z_n(P)` with the closing bond twisted by the
`π` rotation `P` about `z`, both shifted by `κ` per bond (the source-1 table of
`models/spin1/results/source1.json`; eq. (c-twist) of the certificate note). This file records
the table centers and, from the generic sector bounds (`TraceBounds.momentBounds`), the eigenvalue
caps of the lean_seeds design C (`models/lean_seeds/results/seed_design/certificate.json`,
design `real_even`), which uses only the even lengths with trace error `10⁻³`. The seed at
`(36, 49/4)` (`Seed36.lean`) uses all nine lengths.
-/

namespace CertifiedGaps.SpinOneHeisenberg

open PolynomialFilter

/-- The centers of the table of the shifted traces `Z_n(49/4)` for `n = 4, …, 12`. -/
def traceCenter : ℕ → ℚ
  | 4 => 124885379 / 1000000
  | 5 => 117 / 40000
  | 6 => 6435307 / 500000
  | 7 => 54069 / 1000000
  | 8 => 4639283 / 1000000
  | 9 => 39091 / 200000
  | 10 => 2654823 / 1000000
  | 11 => 377629 / 1000000
  | 12 => 1900547 / 1000000
  | _ => 0

/-- The centers of the table of the shifted twisted traces `Z_n(49/4, P)` for
`n = 4, …, 12`. -/
def twistedCenter : ℕ → ℚ
  | 4 => 65757 / 1000000
  | 5 => 1480371 / 500000
  | 6 => 286599 / 1000000
  | 7 => 960207 / 500000
  | 8 => 25507 / 50000
  | 9 => 757423 / 500000
  | 10 => 675419 / 1000000
  | 11 => 82159 / 62500
  | 12 => 157481 / 200000
  | _ => 0

/-- The lengths of the table. -/
def tableLengths : List ℕ := [4, 5, 6, 7, 8, 9, 10, 11, 12]

/-! ### Design C: caps from the even lengths -/

/-- The even lengths of design C. -/
def tableLengthsDesignC : List ℕ := [4, 6, 8, 10, 12]

/-- The table of design C: the centers `± 10⁻³` at the even lengths, untwisted and twisted. -/
def tableDesignC : TraceBounds 1 :=
  TraceBounds.centered ![traceCenter, twistedCenter] tableLengthsDesignC (1 / 1000)

/-- The factor `q` of the invariant-sector filter of design C. -/
def qJDesignC : List ℚ := [4359 / 500000, 0, -32529 / 100000, 0, 1]

/-- The factor `q` of the nontrivial-sector filter of design C. -/
def qNDesignC : List ℚ := [11397 / 1000000, 0, -48363 / 125000, 0, 1]

/-- The invariant-sector filter of design C, `X⁴ (X⁴ - 0.32529 X² + 0.008718)²`. -/
def filterJDesignC : List ℚ := filterPoly 2 qJDesignC

/-- The nontrivial-sector filter of design C, `X⁴ (X⁴ - 0.386904 X² + 0.011397)²`. -/
def filterNDesignC : List ℚ := filterPoly 2 qNDesignC

/-- The filter bounds agree with the stored design-C certificate. -/
theorem pairingBound_filterJDesignC :
    pairingBound filterJDesignC (tableDesignC.loS SectorScheme.d2 0)
        (tableDesignC.hiS SectorScheme.d2 0) =
      some (95825389346247003 / 200000000000000000) := by
  decide +kernel

theorem pairingBound_filterNDesignC :
    pairingBound filterNDesignC (tableDesignC.loS SectorScheme.d2 1)
        (tableDesignC.hiS SectorScheme.d2 1) =
      some (101418917085866523 / 2000000000000000000) := by
  decide +kernel

/-- **Caps of design C.** `D₂` sector data whose traces obey the even part of the table within
`10⁻³` have `|x i| < 1.0018483311` in the invariant sector and `|x i| < 0.8722877194` in a
nontrivial sector. -/
theorem caps_designC {ι : Type*} {sec : Fin (SectorScheme.d2.m + 1) → Finset ι} {x : ι → ℝ}
    (h : tableDesignC.Holds (SectorScheme.d2.trace sec x)) :
    (∀ i ∈ sec 0, |x i| < (10018483311 / 10000000000 : ℚ)) ∧
      ∀ i ∈ sec 1, |x i| < (4361438597 / 5000000000 : ℚ) := by
  have hM := TraceBounds.momentBounds SectorScheme.d2_valid h
  exact ⟨fun i hi ↦ abs_lt_of_filter_bounds _ x (eval_filterPoly_nonneg _ _) (hM 0)
      pairingBound_filterJDesignC (by decide +kernel) hi,
    fun i hi ↦ abs_lt_of_filter_bounds _ x (eval_filterPoly_nonneg _ _) (hM 1)
      pairingBound_filterNDesignC (by decide +kernel) hi⟩

end CertifiedGaps.SpinOneHeisenberg
