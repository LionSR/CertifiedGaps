---
title: Contributing
description: How to submit a new or improved spectral gap certificate or Lean 4 theorem via pull request
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign repo_blob = repo_url | append: '/blob/' | append: site.branch | append: '/' %}

New results come in as pull requests to [{{ site.repository }}]({{ repo_url }}). A pull request
may certify a gap for a new model or coupling, improve an existing bound or lower its $$L_0$$, or
prove an existing bound in Lean 4. The
[open problems]({{ '/leaderboard.html#challenges' | relative_url }}) on the leaderboard are good
places to start.

## A computer-certified bound

Each model lives in its own directory `models/<slug>/` and follows the conventions in
[`AGENTS.md`]({{ repo_blob }}AGENTS.md). A certificate has four parts.

1. A proof of reflection positivity: explicit matrices $$C_\alpha$$, with exact rational or
   algebraic entries, such that $$-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$$ (or
   $$-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$$ for a cell with an internal
   reflection $$F$$), checked as an exact identity, and a real symmetric on-site term $$o$$.
2. Rigorous upper bounds $$\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$$ and
   $$\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$$ at the seed. They may come from exact
   arithmetic, from ball or interval arithmetic, or from verified enclosures of the ring spectra
   (`CertifiedGaps.VerifiedSpectra`) with explicit bounds on floating-point rounding. Unverified
   floating-point numbers, Lanczos runs and DMRG do not count.
3. The gap bound itself, from `CertifiedGaps.certified_gap_bound` (interval arithmetic) or
   `CertifiedGaps.Recurrence.rational_gap_chain` (exact rationals) applied to
   $$(\varepsilon_{\mathrm{s}}, \varepsilon_\tau, \beta_0)$$.
4. A checker, `check.jl`, that recomputes the first three parts from the stored certificate in a
   few minutes, replays the full spectra with `--full`, and fails on any mismatch.

## A Lean 4 theorem

Lean theorems for the leaderboard live in the library `Leaderboard` in `leaderboard/`; its
[README]({{ repo_blob }}leaderboard/README.md) has the details. An unconditional record
(`Record.Unconditional`) proves both the seed bounds and the doubling steps inside the Lean kernel.
A seeded record (`Record.Seeded`, `Record.SeededOn`, `Record.TwoSeeded`) assumes only the two seed
inequalities at $$(n_0, \beta_0)$$, which the computer certificate supplies;
`tools/leaderboard/lean_chains.jl` generates its 14 rational doubling steps, which the kernel
checks with `decide +kernel`. Proofs may use only Lean's standard axioms (`propext`,
`Classical.choice`, `Quot.sound`); CI rejects `sorry`, new axioms and `native_decide`.

## A new model

To suggest a model that has no certificate yet, open an issue with the
[target model form]({{ repo_url }}/issues/new?template=target.yml). Give the Hamiltonian, a
reflection-positive basis if you know one, any floating-point estimate of a working seed
$$(n_0, \beta_0)$$, and the expected gap.

## The leaderboard entry

Add or update the entry in
[`home_page/_data/records.yml`]({{ repo_blob }}home_page/_data/records.yml), whose header
documents the fields, and run the validator before opening the pull request:

```sh
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl
```

## What CI checks

| Check | When it runs | What it does |
|---|---|---|
| `lake build`, `leaderboard_audit` | every pull request and push to `main` (`blueprint.yml`) | builds the Lean theorems, confirms they use only the standard axioms, and writes `lean_records.json` |
| `leanchecker` | changes to Lean code, and weekly (`kernel-check.yml`) | replays `Leaderboard` and everything it imports through the Lean kernel in a fresh process |
| `check_records.jl` | every pull request and push to `main`, and with certificate runs | matches each entry of `records.yml` against its certificate and its Lean theorem, and checks that every displayed decimal is rounded down |
| `models/<slug>/check.jl` | changes to a model or to `julia/`, and weekly for all models (`certificates.yml`) | re-runs the model's checker |

## What a reviewer checks

Some things no program can check. A maintainer reads each pull request to confirm that the
Hamiltonian, couplings, boundary conditions and unit cell in `records.yml` agree with the code in
`models/<slug>/` and with the Lean definition of the spectrum; that the reflection-positivity
identity is exact and the trusted arithmetic described in the model's README meets
[`AGENTS.md`]({{ repo_blob }}AGENTS.md); and that the numerical estimate of the gap cites its
source.
