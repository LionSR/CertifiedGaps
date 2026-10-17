/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Lean.Meta.Basic

/-!
# The `@[gap_record]` attribute

`@[gap_record]` marks a theorem as a leaderboard entry. `lake exe leaderboard_audit` lists every
theorem carrying the tag, checks that its statement is a `Leaderboard.Record.Unconditional` or a
`Leaderboard.Record.Conditional`, and checks its axioms.

The attribute lives in its own module for two reasons: an attribute cannot be used in the module
that declares it, and the audit executable links this module natively, so it must depend on Lean
core only (linking Mathlib into an executable would compile all of Mathlib to C).
-/

namespace Leaderboard

open Lean Meta

/-- The two record heads. Spelled as names, not imported, to keep this module free of Mathlib. -/
def recordHeads : List Name := [`Leaderboard.Record.Unconditional, `Leaderboard.Record.Conditional]

/-- The record head of a statement: the head constant after unfolding reducible definitions only,
so that an `abbrev` of a record is accepted but no model definition is unfolded. -/
def recordHead? (type : Expr) : MetaM (Option Name) := do
  let type ← whnfR (← instantiateMVars type)
  match type.getAppFn.constName? with
  | some n => return if recordHeads.contains n then some n else none
  | none => return none

/-- The tag of leaderboard records. Only theorems whose statement is a record can carry it. -/
initialize gapRecordAttr : TagAttribute ←
  registerTagAttribute `gap_record "a spectral-gap leaderboard record (see leaderboard/README.md)"
    (validate := fun declName => do
      let info ← getConstInfo declName
      unless info.isTheorem do
        throwError "@[gap_record]: '{declName}' must be a theorem"
      let head ← (recordHead? info.type).run'
      if head.isNone then
        throwError "@[gap_record]: the statement of '{declName}' must be \
          `Leaderboard.Record.Unconditional …` or `Leaderboard.Record.Conditional …`")

end Leaderboard
