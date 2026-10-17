/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Lean

/-!
# `lake exe leaderboard_audit [out.json]`

Lists every record of the `Leaderboard.Records.*` modules as JSON (see `leaderboard/README.md`).
The executable links Lean core only and loads the library from its `.olean` files at run time,
so it never compiles Mathlib to C. Numbers are read by kernel reduction, never by compiled code,
which an `implemented_by` could make disagree with the statement.
-/

open Lean Meta

namespace Leaderboard.Audit

/-- The axioms a record may use. -/
def allowedAxioms : List Name := [``propext, ``Classical.choice, ``Quot.sound]

/-- A record head and the positions of its arguments. -/
structure Head where
  /-- The definition, `Leaderboard.Record.*`. -/
  name : Name
  /-- The `kind` field of the JSON. -/
  kind : String
  /-- The number of arguments, implicit ones included. -/
  arity : Nat
  /-- The positions of `spectrum`, `L₀`, `r`, `β`. -/
  spectrum : Nat
  L₀ : Nat
  r : Nat
  β : Nat
  /-- The positions of `n₀`, `s`, `t`, for a seeded record. -/
  seeds : Option (Nat × Nat × Nat)
  /-- The positions of the interval ends `a`, `b`, for a record on a parameter interval. -/
  interval : Option (Nat × Nat) := none

/-- The record heads of `Leaderboard/Record.lean`, with their argument positions. -/
def heads : List Head := [
  -- `Unconditional {d} spectrum L₀ r β`
  { name := `Leaderboard.Record.Unconditional, kind := "unconditional", arity := 5,
    spectrum := 1, L₀ := 2, r := 3, β := 4, seeds := none },
  -- `Seeded {d} [Fintype d] [DecidableEq d] spectrum L₀ r β n₀ s t`
  { name := `Leaderboard.Record.Seeded, kind := "conditional", arity := 10,
    spectrum := 3, L₀ := 4, r := 5, β := 6, seeds := some (7, 8, 9) },
  -- `SeededOn {d} [Fintype d] [DecidableEq d] family a b L₀ r β n₀ s t`
  { name := `Leaderboard.Record.SeededOn, kind := "conditional", arity := 12,
    spectrum := 3, L₀ := 6, r := 7, β := 8, seeds := some (9, 10, 11), interval := some (4, 5) }]

/-- The record head of a statement and its arguments, after unfolding reducible definitions
only (so that no model definition is unfolded). -/
def recordHead? (type : Expr) : MetaM (Option (Head × Array Expr)) := do
  let type ← whnfR (← instantiateMVars type)
  let some n := type.getAppFn.constName? | return none
  let some h := heads.find? (·.name == n) | return none
  let args := type.getAppArgs
  return if args.size == h.arity then some (h, args) else none

/-- Weak head normal form by the kernel, which ignores `implemented_by` and transparency. -/
def kernelWhnf (e : Expr) : MetaM Expr := do
  match Kernel.whnf (← getEnv) {} e with
  | .ok e => return e
  | .error ex => throwError "kernel reduction failed: {ex.toMessageData {}}"

/-- A closed natural-number term, reduced by the kernel to a literal. -/
partial def evalNat (e : Expr) : MetaM Nat := do
  match ← kernelWhnf e with
  | .lit (.natVal n) => return n
  | .const ``Nat.zero _ => return 0
  | .app (.const ``Nat.succ _) a => return (← evalNat a) + 1
  | e => throwError "not a natural-number literal after kernel reduction: {e}"

/-- A closed integer term, reduced by the kernel to `Int.ofNat n` or `Int.negSucc n`. -/
def evalInt (e : Expr) : MetaM Int := do
  let e ← kernelWhnf e
  match e.getAppFn.constName?, e.getAppArgs with
  | some ``Int.ofNat, #[n] => return (← evalNat n)
  | some ``Int.negSucc, #[n] => return -((← evalNat n : Int) + 1)
  | _, _ => throwError "not an integer literal after kernel reduction: {e}"

/-- An exact rational as `"num"` or `"num/den"`. -/
def ratString (q : Rat) : String :=
  if q.den == 1 then s!"{q.num}" else s!"{q.num}/{q.den}"

/-- A closed rational term, reduced by the kernel to the structure `Rat.mk' num den _ _`. The
structure invariant makes `num/den` the reduced form. -/
def evalRat (e : Expr) : MetaM String := do
  let e ← kernelWhnf e
  match e.getAppFn.constName?, e.getAppArgs with
  | some ``Rat.mk', #[num, den, _, _] => return ratString (mkRat (← evalInt num) (← evalNat den))
  | _, _ => throwError "not a rational literal after kernel reduction: {e}"

/-- A numeral of `ℕ`, `ℤ`, `ℚ` or `ℝ` built from literals, `+`, `-`, `*`, `/` and casts. Real
numbers cannot be reduced by the kernel, so model parameters are read from their syntax. -/
partial def evalNumeral (e : Expr) : MetaM Rat := do
  let e := e.consumeMData
  match e.getAppFn.constName?, e.getAppArgs with
  | some ``OfNat.ofNat, #[_, .lit (.natVal n), _] => return n
  | some ``OfScientific.ofScientific, #[_, _, .lit (.natVal m), .const b _, .lit (.natVal x)] =>
    return if b == ``Bool.true then mkRat m (10 ^ x) else (m * 10 ^ x : Nat)
  | some ``HAdd.hAdd, #[_, _, _, _, a, b] => return (← evalNumeral a) + (← evalNumeral b)
  | some ``HSub.hSub, #[_, _, _, _, a, b] => return (← evalNumeral a) - (← evalNumeral b)
  | some ``HMul.hMul, #[_, _, _, _, a, b] => return (← evalNumeral a) * (← evalNumeral b)
  | some ``HDiv.hDiv, #[_, _, _, _, a, b] => return (← evalNumeral a) / (← evalNumeral b)
  | some ``Neg.neg, #[_, _, a] => return -(← evalNumeral a)
  | some ``Nat.cast, #[_, _, a] | some ``Int.cast, #[_, _, a] | some `Rat.cast, #[_, _, a] =>
    evalNumeral a
  | _, _ => throwError "a model parameter must be a numeral, not {e}"

/-- The explicit arguments of an application, which are the model parameters. -/
def explicitArgs (e : Expr) : MetaM (Array Expr) := do
  let info ← getFunInfo e.getAppFn
  let args := e.getAppArgs
  return (args.zip (info.paramInfo.map (·.isExplicit))).filterMap fun (a, x) ↦
    if x then some a else none

/-- The spectrum of a record as `{model, frame, params, printed}`. For `R.spectrum` of a
reflection-positive ring `R = M p₁ … pₖ` the model is `M` and the frame `spectrum`; otherwise
(`physicalSpectrum p₁ … pₖ`, `SpinOneHeisenberg.spectrum`) the model is the head constant and
the frame its last name component. For the family `fun x ↦ …` of a record on an interval, the
parameter `x` is listed as `"x"`. `printed` is the full term without notation. -/
def spectrumJson (e : Expr) : MetaM Json := do
  let e ← instantiateMVars e
  let opts := ({} : Options) |>.set `pp.notation false |>.set `pp.fieldNotation false
    |>.set `pp.fullNames true
  let printed := (← withOptions (fun _ ↦ opts) (ppExpr e)).pretty 100000
  lambdaTelescope e fun xs e ↦ do
    let param (a : Expr) : MetaM String := do
      if xs.contains a then return "x"
      return ratString (← evalNumeral a)
    let (frame, m) := match e.getAppFn.constName? with
      | some `CertifiedGaps.RPRing.spectrum => ("spectrum", e.getAppArgs.back!)
      | some n => (n.componentsRev.head!.toString, e)
      | none => ("", e)
    let some model := m.getAppFn.constName? |
      throwError "the spectrum must be a constant applied to numerals, not {e}"
    let params ← (← explicitArgs m).mapM param
    return Json.mkObj [("model", toJson model.toString), ("frame", toJson frame),
      ("params", toJson params), ("printed", toJson printed)]

/-- The JSON entry of the record `n` of module `mod`, whose statement has head `h` applied to
`args`. Throws on a disallowed axiom or on arguments that are not literals. -/
def recordJson (n mod : Name) (h : Head) (args : Array Expr) : MetaM Json := do
  let axioms ← collectAxioms n
  let bad := axioms.filter (!allowedAxioms.contains ·)
  unless bad.isEmpty do
    throwError (s!"disallowed axioms {bad.toList} (`sorry` adds sorryAx, `native_decide` adds \
      Lean.ofReduceBool or `*._native.native_decide.ax_*`; only {allowedAxioms} are allowed)")
  let seeds ← match h.seeds with
    | some (n₀, s, t) =>
      pure [("n0", toJson (← evalNat args[n₀]!)), ("s", toJson (← evalRat args[s]!)),
        ("t", toJson (← evalRat args[t]!))]
    | none => pure [("n0", Json.null), ("s", Json.null), ("t", Json.null)]
  let interval ← match h.interval with
    | some (a, b) => pure (toJson [← evalRat args[a]!, ← evalRat args[b]!])
    | none => pure Json.null
  return Json.mkObj <| [
    ("name", toJson n.toString), ("module", toJson mod.toString), ("kind", toJson h.kind),
    ("spectrum", ← spectrumJson args[h.spectrum]!), ("L0", toJson (← evalNat args[h.L₀]!)),
    ("r", toJson (← evalRat args[h.r]!)), ("beta", toJson (← evalRat args[h.β]!))] ++ seeds ++
    [("interval", interval), ("axioms", toJson (axioms.map (·.toString) |>.qsort (· < ·)))]

/-- Every record of the `Leaderboard.Records.*` modules, sorted by name, and every problem. -/
def run : MetaM (Array (Name × Json) × Array String) := do
  let env ← getEnv
  let mut entries := #[]
  let mut errs := #[]
  for h : idx in [0:env.header.moduleNames.size] do
    let mod := env.header.moduleNames[idx]
    unless (`Leaderboard.Records).isPrefixOf mod do continue
    for n in env.header.moduleData[idx]!.constNames do
      let some info := env.find? n | continue
      unless info.isTheorem && !n.isInternalDetail do continue
      let some (head, args) ← recordHead? info.type | continue
      try entries := entries.push (n, ← recordJson n mod head args)
      catch e => errs := errs.push s!"{n}: {← e.toMessageData.toString}"
  if entries.isEmpty then errs := errs.push "no records found in Leaderboard.Records.*"
  return (entries.qsort (·.1.toString < ·.1.toString), errs)

end Leaderboard.Audit

open Leaderboard.Audit in
/-- `leaderboard_audit [out.json]`: writes the records to `out.json` (stdout if `-` or absent),
or writes nothing and exits 1. -/
unsafe def main (args : List String) : IO UInt32 := do
  let out := args.headD "-"
  initSearchPath (← findSysroot)
  enableInitializersExecution
  let env ← importModules #[{ module := `Leaderboard }] {} (trustLevel := 1024)
    (loadExts := true)
  let ctx : Core.Context := { fileName := "<leaderboard_audit>", fileMap := default }
  let (entries, errs) ← (run.run' {} {}).toIO' ctx { env }
  unless errs.isEmpty do
    IO.eprintln "leaderboard_audit: FAILED"
    for e in errs do IO.eprintln s!"  {e}"
    return 1
  let text := (Json.arr (entries.map (·.2))).pretty ++ "\n"
  if out == "-" then
    IO.print text
    return 0
  if let some dir := (System.FilePath.mk out).parent then
    if dir.toString != "" then IO.FS.createDirAll dir
  IO.FS.writeFile out text
  IO.eprintln s!"leaderboard_audit: {entries.size} records written to {out}"
  return 0
