# Contributing

Contributions to the [leaderboard of certified spectral gaps](https://sirui-lu.com/CertifiedGaps/leaderboard.html)
are submitted via pull request:

* a computer-certified gap bound for a new model or parameter point,
* a tighter certified bound or smaller threshold length $L_0$ for an existing record, or
* a Lean 4 formalization (seeded or unconditional) of a certified bound.

See the [contribution guide](https://sirui-lu.com/CertifiedGaps/contribute.html)
(source: [`home_page/contribute.md`](home_page/contribute.md)) for full details and the
[open challenges](https://sirui-lu.com/CertifiedGaps/leaderboard.html#challenges) for concrete
targets.

## Summary of requirements

1. **Computer certificate (`code/models/<slug>/`):** Follow [`AGENTS.md`](AGENTS.md). Provide an exact
   reflection-positivity factorization, rigorous enclosures of the seed purity defects at
   $(n_0, \beta_0)$, the gap bound from the doubling recurrence, and a standalone checker
   `check.jl`. Update the summary table in the root [`README.md`](README.md).
2. **Lean 4 theorem (`code/leaderboard/`):** Follow [`code/leaderboard/README.md`](code/leaderboard/README.md).
   Proofs may depend only on Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`).
3. **Leaderboard entry (`home_page/_data/records.yml`):** Add or update the row in
   [`home_page/_data/records.yml`](home_page/_data/records.yml) and run the validator locally:

   ```sh
   julia --project=code/tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
   julia --project=code/tools/leaderboard code/tools/leaderboard/check_records.jl
   ```

4. **Issues:** Open an issue to coordinate work on an open problem or use the
   [target model form](https://github.com/LionSR/CertifiedGaps/issues/new?template=target.yml) to
   propose a new model.

## CI and review

Every pull request runs automated CI checks that build the Lean library, audit all theorem axioms
and parameters (`leaderboard_audit`), replay proofs in the Lean kernel (`leanchecker` when Lean
files change and weekly), cross-check every `records.yml` entry against its certificate and Lean 4
theorem (`check_records.jl`), and run `check.jl` for each modified model. A maintainer then
verifies that the Hamiltonian, unit-cell convention, trusted arithmetic base, and verification tier
match across the model directory, the Lean 4 theorem, and the leaderboard entry.
