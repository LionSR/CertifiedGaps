# Contributing

Contributions are welcome: a gap certificate for a new model, a sharper bound for an existing
one, a Lean record that upgrades a computer-certified bound, or a proposed target. The full guide
is the [contribute page](https://LionSR.github.io/haldane-gap-diagrams/contribute.html) of the
website (source: [`home_page/contribute.md`](home_page/contribute.md)); the results are on the
[leaderboard](https://LionSR.github.io/haldane-gap-diagrams/leaderboard.html).

## In short

1. **A computer certificate** (tier *certified*) lives in `models/<slug>/` and follows
   [`AGENTS.md`](AGENTS.md): an exact reflection-positivity identity, rigorous seed enclosures,
   the certified recurrence, and a standalone `check.jl`. Add a row to the summary table of the
   root [`README.md`](README.md).
2. **A Lean record** (tiers *Lean, unconditional* and *Lean, conditional*) is a theorem tagged
   `@[gap_record]` of type `Leaderboard.Record.Unconditional` or `Leaderboard.Record.Conditional`
   under `leaderboard/Records/`, using only the standard axioms. See
   [`leaderboard/README.md`](leaderboard/README.md).
3. **Every record** has a row in [`home_page/_data/records.yml`](home_page/_data/records.yml)
   (schema in its header). Validate it with

   ```sh
   julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
   julia --project=tools/leaderboard tools/leaderboard/check_records.jl
   ```

4. **A target** is proposed with the "Propose a target model" issue form; a problem with a
   published record is reported with the "Report a problem with a record" form.

## What continuous integration checks

* `blueprint.yml` (every pull request): the Lean build including `lake build Leaderboard`; the
  record audit `lake exe leaderboard_audit` (axioms) and kernel re-check `lake exe
  leaderboard_recheck`; `tools/leaderboard/check_records.jl` against the generated Lean records;
  the blueprint and its declaration check; and the site build.
* `certificates.yml` (pull requests that touch `models/`): the fast `check.jl` of every changed
  model directory.

Reviewers then check that the Hamiltonian in each row is the one the certificate and the Lean
theorem are about, that the trusted base is acceptable, and that the claimed tier is right.
