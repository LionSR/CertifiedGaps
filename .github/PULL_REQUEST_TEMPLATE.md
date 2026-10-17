## Summary

<!-- What does this pull request add or change? For a record: the model, parameter point, bound,
lengths covered, and the claimed tier. -->

## Checklist for a new or improved record

- [ ] Row added or updated in `home_page/_data/records.yml`; an improved bound moves the old
      one into `history`
- [ ] `julia --project=tools/leaderboard tools/leaderboard/check_records.jl` passes
- [ ] Tier claimed: `lean-unconditional` / `lean-conditional` / `certified` / `unreviewed` /
      `assumed` / `target` (delete the others)
- [ ] Computer certificate (if any): `models/<slug>/` follows `AGENTS.md`, `results/certificate*.json`
      committed, and `check.jl` passes (run time: ...)
- [ ] Lean record (if any): name `Leaderboard.Records.<Model>.<point>`, its `L0`, `r` and `beta`
      match the YAML row, and it follows `leaderboard/README.md`
      (CI enforces its axiom and kernel rules)
- [ ] Root `README.md` summary table updated (per `AGENTS.md`) for a new or updated model

## For reviewers

- [ ] The Hamiltonian, parameter point and unit cell in the row match the model directory and,
      for a Lean record, the spectrum named in the theorem
- [ ] The trusted base in the model's `README.md` is acceptable and the tier is right
