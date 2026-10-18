## Summary

<!-- What does this pull request add or change? For a record: the model, parameter point, bound,
lengths covered, and the claimed tier. -->

## Checklist for a new or improved record

- [ ] Row added or updated in `home_page/_data/records.yml`; an improved bound moves the previous
      entry into `history`
- [ ] `julia --project=tools/leaderboard tools/leaderboard/check_records.jl` passes
- [ ] Tier claimed (delete the others): `lean-unconditional` (Unconditional Lean 4 theorem) /
      `lean-conditional` (Lean 4 theorem from certified seed) / `certified`
      (Computer-certified bound). Open targets and challenges carry no tier
- [ ] Computer certificate (if any): `models/<slug>/` follows `AGENTS.md`, `results/certificate*.json`
      committed, and `check.jl` passes (run time: ...)
- [ ] Lean record (if any): named `Leaderboard.Records.<Model>.<point>`, its `L0`, `r`, and `beta`
      match the YAML row, and it follows `leaderboard/README.md`
      (CI enforces axiom and kernel checks)
- [ ] Root `README.md` summary table updated (per `AGENTS.md`) for a new or updated model

## For reviewers

- [ ] The Hamiltonian, parameter point, and unit cell in the row match the model directory and,
      for a Lean record, the spectrum defined in the theorem
- [ ] The trusted arithmetic base in the model's `README.md` satisfies `AGENTS.md` and the assigned
      verification tier matches the theorem proving the displayed bound
