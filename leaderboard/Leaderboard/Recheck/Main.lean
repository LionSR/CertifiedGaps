/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Lean

/-!
# `lake exe leaderboard_recheck [--leanchecker]`

Re-checks in the Lean kernel every declaration of the `Leaderboard` library **and of every module
of this repository that it imports** (`CertifiedGaps.*`, `SeedKernel.*`), from the `.olean`
files. `lake build` already sends each declaration to the kernel, but a file can switch that off
(`set_option debug.skipKernelTC true`) or write declarations into the environment from
metaprograms; replaying the compiled declarations detects both.

The check is the one of `leanchecker` (the toolchain's upstreamed `lean4checker`): replay the
declarations through `Lean.Kernel.Environment.replay` on top of an environment that contains
only their imports. `leanchecker` re-imports Mathlib once per module, which for the 50-odd
modules here takes about ten minutes; this executable imports the external dependencies
(Mathlib, Lean, …) once and replays all repository modules together, which is the same kernel
work. External packages are not replayed: they are pinned upstream releases.

With `--leanchecker` it instead runs the toolchain's `leanchecker` on each repository module in
turn, as an independent cross-check (slow).
-/

open Lean

namespace Leaderboard.Recheck

/-- The root module whose import closure is checked. -/
def rootModule : Name := `Leaderboard

/-- The modules of this repository in the import closure of `rootModule`, in import order: those
whose `.olean` lies in this package's build directory rather than in `.lake/packages`. -/
def repoModules (env : Environment) (buildDir : System.FilePath) : IO (Array Name) := do
  let mut mods := #[]
  for m in env.header.moduleNames do
    let olean ← IO.FS.realPath (← findOLean m)
    if olean.toString.startsWith (buildDir.toString ++ "/") then
      mods := mods.push m
  return mods

/-- The declarations of module `m`, read from its `.olean` parts as `leanchecker` does: the last
(most private) part contains all of them. -/
def moduleConstants (m : Name) : IO (Array ConstantInfo) := do
  let mFile ← findOLean m
  let mut fnames := #[mFile]
  let sFile := OLeanLevel.server.adjustFileName mFile
  if ← sFile.pathExists then
    fnames := fnames.push sFile
    let pFile := OLeanLevel.private.adjustFileName mFile
    if ← pFile.pathExists then fnames := fnames.push pFile
  let parts ← readModuleDataParts fnames
  let some (data, _) := parts.back? | throw <| IO.userError s!"no module data for {m}"
  return data.constants

/-- Replays every declaration of `mods` on top of their external imports; returns the number of
declarations replayed. -/
def replayAll (env : Environment) (mods : Array Name) : IO Nat := do
  let modSet := mods.foldl (·.insert ·) (∅ : NameSet)
  -- The external modules imported directly by a repository module; their closure is the base.
  let mut base : Array Import := #[]
  let mut seen : NameSet := ∅
  for m in mods do
    let some idx := env.getModuleIdx? m | throw <| IO.userError s!"module {m} not loaded"
    for imp in env.header.moduleData[idx]!.imports do
      unless modSet.contains imp.module || seen.contains imp.module do
        seen := seen.insert imp.module
        base := base.push imp
  IO.println s!"leaderboard_recheck: importing {base.size} external modules"
  let baseEnv ← importModules base {}
  let mut newConstants : Std.HashMap Name ConstantInfo := {}
  for m in mods do
    for ci in ← moduleConstants m do
      newConstants := newConstants.insert ci.name ci
  IO.println s!"leaderboard_recheck: replaying {newConstants.size} declarations of \
    {mods.size} modules"
  discard <| baseEnv.toKernelEnv.replay newConstants
  return newConstants.size

/-- Runs the toolchain's `leanchecker` on each module in turn; `true` if it accepts all. -/
def leancheckerAll (mods : Array Name) : IO Bool := do
  let exe := (← findSysroot) / "bin" / ("leanchecker" ++ System.FilePath.exeExtension)
  let mut ok := true
  for m in mods do
    IO.println s!"leanchecker {m}"
    let out ← IO.Process.output { cmd := exe.toString, args := #[m.toString] }
    unless out.exitCode == 0 do
      IO.eprintln s!"{out.stdout}{out.stderr}"
      ok := false
  return ok

end Leaderboard.Recheck

open Leaderboard.Recheck in
/-- `leaderboard_recheck [--leanchecker]`. -/
def main (args : List String) : IO UInt32 := do
  initSearchPath (← findSysroot)
  -- This executable lives in `<package>/.lake/build/bin`.
  let buildDir ← IO.FS.realPath <| ((← IO.appPath).parent.bind (·.parent)).getD "."
  let env ← importModules #[{ module := rootModule }] {}
  let mods ← repoModules env buildDir
  if mods.isEmpty then
    IO.eprintln s!"leaderboard_recheck: no repository modules found under {buildDir}"
    return 1
  try
    if args.contains "--leanchecker" then
      unless ← leancheckerAll mods do
        IO.eprintln "leaderboard_recheck: FAILED"
        return 1
      IO.println s!"leaderboard_recheck: leanchecker accepts all {mods.size} modules"
    else
      let n ← replayAll env mods
      IO.println s!"leaderboard_recheck: the kernel accepts all {n} declarations of \
        {mods.size} modules"
    return 0
  catch e =>
    IO.eprintln s!"leaderboard_recheck: FAILED\n{e}"
    return 1
