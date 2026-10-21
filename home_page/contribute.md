---
title: Contributing
description: How to submit a new or improved spectral gap certificate or Lean 4 theorem via pull request
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign repo_blob = repo_url | append: '/blob/' | append: site.branch | append: '/' %}

Send new results as pull requests to [{{ site.repository }}]({{ repo_url }}). A pull request can
do one of these tasks:

* certify a gap bound for a new model or a new coupling,
* sharpen an existing bound or lower its threshold $$L_0$$,
* prove an existing bound in Lean 4.

The [open problems]({{ '/leaderboard.html#challenges' | relative_url }}) on the leaderboard are
good tasks to start with.

## Certify a bound by computer

Put each model in its own directory `code/models/<slug>/` following
[`AGENTS.md`]({{ repo_blob }}AGENTS.md). A certificate has four parts.

1. Reflection positivity. Give matrices $$C_\alpha$$ with exact rational or algebraic entries such
   that $$-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$$, or
   $$-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$$ for a cell with an internal
   reflection $$F$$. Check this identity exactly. The on-site term $$o$$ must be real and
   symmetric.
2. Seed bounds. Calculate rigorous upper bounds on the defects
   $$p = 1 - P_{\mathrm S}(n_0, \beta_0)$$ and $$q = 1 - P_{\mathrm T}(2n_0, \beta_0)$$. Use exact
   arithmetic, ball or interval arithmetic, or verified enclosures of ring spectra
   (`CertifiedGaps.VerifiedSpectra`) with explicit bounds on rounding errors. Numerical
   values without proven error bounds, Lanczos runs and DMRG are not sufficient.
3. Gap bound. Apply `CertifiedGaps.certified_gap_bound` (interval arithmetic) or
   `CertifiedGaps.Recurrence.rational_gap_chain` (exact rationals) to the bounds on $$p$$ and
   $$q$$ and to $$\beta_0$$.
4. Checker. Write `check.jl`. It recomputes parts 1 to 3 from the stored certificate in
   a few minutes. With `--full`, it also replays the full spectra. It exits with an error
   on any mismatch.

## Prove a bound in Lean 4

Put the theorem in the library `Leaderboard` in `code/leaderboard/`. Its
[README]({{ repo_blob }}code/leaderboard/README.md) gives the details.

* An unconditional record (`Record.Unconditional`) proves the seed bounds and the doubling steps
  in the Lean kernel.
* A seeded record (`Record.Seeded`, `Record.SeededOn`, `Record.TwoSeeded`) uses the two seed
  inequalities at $$(n_0, \beta_0)$$ as hypotheses. The computer certificate supplies these
  inequalities. Run `code/tools/leaderboard/lean_chains.jl` to generate the 14 exact doubling steps. The
  kernel checks them with `decide +kernel`.

Use only the standard axioms of Lean (`propext`, `Classical.choice`, `Quot.sound`). CI rejects
`sorry`, new axioms and `native_decide`.

## Propose a new model

Open an issue with the [target model form]({{ repo_url }}/issues/new?template=target.yml). Give
the Hamiltonian, a reflection-positive basis if you know one, a numerical estimate of a seed
$$(n_0, \beta_0)$$ if you have one, and the expected gap.

## Update the leaderboard

Add or change the entry in
[`home_page/_data/records.yml`]({{ repo_blob }}home_page/_data/records.yml). The header of the
file describes the fields. Run the validator before you open the pull request:

```sh
julia --project=code/tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=code/tools/leaderboard code/tools/leaderboard/check_records.jl
```

## Automatic checks

| Check | When it runs | What it does |
|---|---|---|
| `lake build`, `leaderboard_audit` | each pull request and each push to `main` (`blueprint.yml`) | builds the Lean theorems, checks that they use only the standard axioms, and writes `lean_records.json` |
| `leanchecker` | changes to Lean code, and each week (`kernel-check.yml`) | replays `Leaderboard` and all its imports in the Lean kernel, in a fresh process |
| `check_records.jl` | each pull request and each push to `main`, and with certificate runs | compares each entry of `records.yml` with its certificate and its Lean theorem, and checks that each displayed decimal is rounded down |
| `code/models/<slug>/check.jl` | changes to a model or to `code/julia/`, and each week for all models (`certificates.yml`) | reruns the checker of each changed model |

## Human review

Programs cannot check all parts of a submission. A maintainer reads each pull request and checks
these items:

* The Hamiltonian, couplings, boundary conditions and unit cell in `records.yml` agree with the
  code in `code/models/<slug>/` and with the Lean definition of the spectrum.
* The reflection-positivity identity is exact.
* The trusted arithmetic in the README of the model follows [`AGENTS.md`]({{ repo_blob }}AGENTS.md).
* The numerical estimate of the gap has a source.
