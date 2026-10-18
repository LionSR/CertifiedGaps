---
title: Contribute a result
description: How to submit a new or improved spectral gap certificate or Lean 4 theorem via pull request
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign repo_blob = repo_url | append: '/blob/' | append: site.branch | append: '/' %}

Contributions to the [leaderboard]({{ '/leaderboard.html' | relative_url }}) are submitted as pull
requests to [{{ site.repository }}]({{ repo_url }}). A submission can be:

* a computer-certified spectral gap bound for a new model or parameter point,
* a tighter lower bound or smaller threshold length $$L_0$$ for an existing record, or
* a Lean 4 formalization (seeded or unconditional) of a certified bound.

See the [open challenges]({{ '/leaderboard.html' | relative_url }}#challenges) for concrete targets.
Every record certifies that for all ring lengths $$L \ge L_0$$ of the stated parity, the
Hamiltonian $$H_L$$ has a unique ground state and a spectral gap $$\Delta_L = E_1(L) - E_0(L)$$ at
least an explicit rational or logarithmic bound.

## 1. Computer-certified gap bound

Computer certificates live in `models/<slug>/` and follow the repository conventions in
[`AGENTS.md`]({{ repo_blob }}AGENTS.md). Every certificate requires four components:

1. **Reflection-positivity proof.** Explicit matrices $$C_\alpha$$ with exact rational or
   closed-form algebraic entries satisfying
   $$-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$$ (or
   $$-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$$ with an internal cell reflection
   $$F$$), verified by an exact algebraic identity, together with a real symmetric on-site term
   $$o$$.
2. **Seed enclosures.** Rigorous rational upper bounds
   $$\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$$ and
   $$\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$$ at seed $$(n_0, \beta_0)$$, derived from
   exact rational/integer arithmetic, Arb ball arithmetic, interval arithmetic, or verified
   symmetry-sector eigenvalue enclosures (`CertifiedGaps.VerifiedSpectra`) with explicit a priori
   IEEE binary64 rounding-error bounds. Unverified floating-point numbers, Lanczos/ARPACK runs,
   and DMRG are not accepted as certificates.
3. **Gap bound.** `CertifiedGaps.certified_gap_bound` (interval arithmetic) or
   `CertifiedGaps.Recurrence.rational_gap_chain` (exact rational certificate) applied to
   $$(\varepsilon_{\mathrm{s}}, \varepsilon_\tau, \beta_0)$$.
4. **Standalone checker (`check.jl`).** Recomputes steps 1–3 from the committed certificate in a
   few minutes (and supports `--full` for a full spectral replay) and exits non-zero on any
   mismatch.

## 2. Lean 4 theorem

The Lean library `Leaderboard` in `leaderboard/` packages one theorem per Lean-verified record; see
[`leaderboard/README.md`]({{ repo_blob }}leaderboard/README.md) for full instructions.

* **Unconditional records** (`Record.Unconditional`) discharge both the finite seed bounds and the
  doubling recurrence inside the Lean kernel with no hypotheses.
* **Seeded records** (`Record.Seeded`, `Record.SeededOn`, `Record.TwoSeeded`) take as hypotheses
  only the seed purity inequalities at $$(n_0, \beta_0)$$, which are established by the companion
  computer certificate in `models/<slug>/`. Running `tools/leaderboard/lean_chains.jl` generates
  the 14-step rational doubling chain checked in the Lean kernel via `decide +kernel`.
* Proofs may depend only on Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`);
  `sorry`, custom `axiom` declarations, and `native_decide` are rejected by CI.

## 3. Proposing a target model

To propose a model that is not yet certified, open an issue using the
[target model form]({{ repo_url }}/issues/new?template=target.yml) with the Hamiltonian, a
reflection-positive basis (if known), any floating-point estimate of an admissible seed
$$(n_0, \beta_0)$$, and the benchmark gap.

## 4. Updating the leaderboard data

Add or update the record in
[`home_page/_data/records.yml`]({{ repo_blob }}home_page/_data/records.yml) (schema documented in
the file header) and run the validator locally before opening a pull request:

```sh
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl
```

## Automated CI verification

CI workflows run the following checks according to which paths a pull request modifies, as well as
on weekly scheduled runs:

| Check | Trigger | What it verifies |
|---|---|---|
| `lake build` and `leaderboard_audit` | Every PR and push to `main` (`blueprint.yml`) | Compiles all Lean 4 theorems, verifies that only standard axioms are used, reduces record parameters to exact kernel literals, and exports `lean_records.json` |
| `leanchecker` | PRs/pushes modifying Lean code (`lean/`, `leaderboard/`, toolchain files) and weekly (`kernel-check.yml`) | Replays all declarations in `Leaderboard` and the repository modules it imports through the Lean kernel in a fresh process |
| `check_records.jl` | Every PR and push to `main` (`blueprint.yml`), and on certificate runs (`certificates.yml`) | Verifies that every row in `records.yml` matches its committed certificate (seed, defects, and exact bound), that all displayed decimals are rounded down, and that every Lean entry matches `lean_records.json` |
| `models/<slug>/check.jl` | PRs/pushes modifying non-documentation files in `models/<slug>/` or `julia/`, and weekly across all models (`certificates.yml`) | Reruns the standalone certificate checker for each modified model (or all models when `julia/` changes or on the weekly schedule) |

## Maintainer review

Because automated checks cannot verify physical modeling intent, a human reviewer checks that:

* the Hamiltonian, coupling parameters, boundary conditions, and unit-cell convention in
  `records.yml` match both the model implementation in `models/<slug>/` and the Lean 4 spectrum
  definition,
* the reflection-positivity factorization is exact and the arithmetic trusted base documented in
  `models/<slug>/README.md` satisfies [`AGENTS.md`]({{ repo_blob }}AGENTS.md),
* the assigned verification tier matches the theorem proving the displayed bound, and
* the numerical benchmark gap cites its source.
