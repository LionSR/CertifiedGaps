---
title: Contribute a result
description: How to add a proven gap to the leaderboard, or improve one, by pull request
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign repo_blob = repo_url | append: '/blob/' | append: site.branch | append: '/' %}

Anyone can add a result to the [leaderboard]({{ '/leaderboard.html' | relative_url }}). A result
can be a proven gap for a new model, a larger proven gap for a known point, or a stronger proof
of a known bound, such as a Lean theorem for a bound that so far has a computer proof. You submit
it as a pull request to [{{ site.repository }}]({{ repo_url }}). Automatic checks run on every
pull request, and then a person reviews it. The
[open challenges]({{ '/leaderboard.html' | relative_url }}#challenges) list concrete next steps.

Each result states that for every ring length $$L \ge L_0$$ of the stated parity, the
Hamiltonian $$H_L$$ has a unique ground state and a gap $$\Delta_L$$ at least an exact number.

## A computer-checked proof

The repository's [`AGENTS.md`]({{ repo_blob }}AGENTS.md) defines what a computer proof must
contain. It goes in a new directory `models/<name>/` with four parts.

1. **Reflection positivity.** Exact matrices $$C_\alpha$$ with
   $$-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$$, checked by an exact algebraic
   identity, and a real symmetric on-site term.
2. **The starting scale.** Rigorous rational upper bounds on $$1 - P_{\mathrm{s}}(n_0, \beta_0)$$
   and $$1 - P_\tau(2n_0, \beta_0)$$, from exact, interval, or floating-point arithmetic with an
   explicit bound on the rounding error. Unverified floating-point numbers, Lanczos runs and DMRG
   do not count.
3. **The gap bound.** The repository's recurrence applied to these two numbers.
4. **A checker program** `check.jl` that recomputes everything from the stored certificate in a
   few minutes and fails on any mismatch.

## A Lean theorem

The Lean library `Leaderboard` in the directory `leaderboard/` holds one theorem per result. Its
[README]({{ repo_blob }}leaderboard/README.md) explains how to add one. The theorem either has no
assumptions, or assumes only the two inequalities at its starting scale (at each of two starting
scales, for a proof with two), which a computer proof in `models/` then establishes. For a row
with such a computer proof, `tools/leaderboard/lean_chains.jl` writes the chain of rational
numbers that the Lean proof checks. The theorem may use only Lean's standard axioms, so no
`sorry` and no `native_decide`. The automatic checks build the library, list the axioms of every
theorem and write the file `lean_records.json` that the leaderboard reads. A reviewer checks that
the Hamiltonian in the theorem is the physical model of the leaderboard row. No program can check
that.

## A proposed target

If the method should apply to a model that is not on the leaderboard, open an issue with the
[target form]({{ repo_url }}/issues/new?template=target.yml). Give the Hamiltonian, a
reflection-positive basis if you know one, any floating-point estimate of the starting scale
$$(n_0, \beta_0)$$, and the known gap.

## The leaderboard row

Each result has a row in [`home_page/_data/records.yml`]({{ repo_blob }}home_page/_data/records.yml),
whose header documents every field. The row gives the model, the parameters, the ring lengths,
the exact bound with a decimal rounded down, the starting scale, the paths of the certificate and
checker, the Lean theorem if there is one, the GitHub handles of the contributors, and earlier
bounds for the same point. Check the file before you open the pull request:

```sh
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl
```

## What the automatic checks do

| check | what it establishes |
|---|---|
| Lean build and audit | every Lean theorem compiles and uses only the standard axioms |
| kernel re-check (when Lean code changes, and weekly) | every proof passes a second, independent check in the Lean kernel |
| `check_records.jl` | every row matches its certificate (starting scale, defects and bound, compared exactly), every decimal is rounded down, and every Lean row matches its theorem |
| `check.jl` of each changed model | the computer proof recomputes to the stored certificate |

Every number on the leaderboard is rechecked this way on each update and every week.

## What a reviewer checks

* The Hamiltonian, the parameters and the unit cell in the row match the model directory and the
  Lean theorem.
* The reflection-positivity identity is exact, and the arithmetic the proof relies on, stated in
  the model's README, is acceptable.
* The badge is right: a Lean badge only when the Lean theorem proves the bound shown.
* The numerical comparison gap names its source.
