import LeanChecker

/-
Run with `KERNEL_CHECK_MODULE=Module.Name lake env lean lean/scripts/KernelCheck.lean`.

Lean 4.36.0-rc1's leanchecker crashes locally when replaying a second module in one process.
Its CLI also expands namespace prefixes, so invoking it once per module is not exact selection.
Call the unchanged bundled replay routine once instead: it checks all declaration data parts
with the same kernel and imported environment. A separate process bounds each replay's lifetime.
Using #eval avoids redefining the main function imported from LeanChecker; errors still fail Lean.
-/
#eval show IO Unit from do
  let some module ← IO.getEnv "KERNEL_CHECK_MODULE"
    | throw <| IO.userError "KERNEL_CHECK_MODULE is required"
  if module.toName.isAnonymous then
    throw <| IO.userError "KERNEL_CHECK_MODULE must name a module"
  Lean.initSearchPath (← Lean.findSysroot)
  IO.println s!"replaying {module}"
  replayFromImports module.toName
