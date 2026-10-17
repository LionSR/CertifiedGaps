---
title: Contribute a record
description: How to submit a new certified spectral gap, or improve one, by pull request
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign repo_blob = repo_url | append: '/blob/' | append: site.branch | append: '/' %}

The [leaderboard]({{ '/leaderboard.html' | relative_url }}) is open: anyone can submit a new
model, a sharper bound for an existing one, or a stronger verification of an existing bound
(for example, a Lean record for a result that so far has only a computer certificate). Every
submission is a pull request against
[`{{ site.repository }}`]({{ repo_url }}), checked by continuous integration and then by a
reviewer.

A record states: *for every ring length $$L \ge L_0$$ of the stated parity, the Hamiltonian
$$H_L$$ has a unique ground state and a spectral gap $$\Delta_L \ge b$$*, with $$b$$ an exact
number. There are four ways to contribute.

## (a) A rigorous computer certificate (tier *certified*)

Follow [`AGENTS.md`]({{ repo_blob }}AGENTS.md), which defines a gap certificate precisely. In
short, a new directory `models/<slug>/` must contain:

1. **Reflection positivity:** exact matrices $$C_\alpha$$ with
   $$-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$$ (or the internally reflected form
   for multi-site cells), verified by an exact algebraic identity, and a real symmetric on-site
   term.
2. **Seed enclosures:** rigorous rational upper bounds
   $$\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$$ and
   $$\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$$ from exact, interval, or verified
   floating-point arithmetic with an explicit rounding analysis. A floating-point number without
   an error bound, a Lanczos run, or DMRG is not a certificate.
3. **The gap bound:** `CertifiedGaps.certified_gap_bound` or
   `CertifiedGaps.Recurrence.rational_gap_chain` applied to $$(\varepsilon_{\mathrm{s}}, \varepsilon_\tau, \beta_0)$$.
4. **A standalone checker** `check.jl` that recomputes all of this from the stored certificate
   in a few minutes, supports `--full`, and fails loudly on any mismatch.

The directory follows the standard layout (`README.md` with the required headings, `IDEAS.md`,
`certify.jl`, `check.jl`, `test/runtests.jl`, `results/certificate*.json`), and the summary
table in the root `README.md` gets a row.

## (b) A Lean record (tiers *Lean, unconditional* and *Lean, conditional*)

The Lean library `Leaderboard` (top-level directory `leaderboard/`; see its
[README]({{ repo_blob }}leaderboard/README.md)) holds one typed theorem per record:

* add a theorem tagged `@[gap_record]` of type `Leaderboard.Record.Unconditional` (no
  hypotheses) or `Leaderboard.Record.Conditional` (the two seed inequalities as hypotheses)
  under `leaderboard/Records/`;
* use only the standard axioms (`propext`, `Classical.choice`, `Quot.sound`): no `sorry`, no
  `native_decide`, no new `axiom`, no `implemented_by` or `extern`, no
  `set_option debug.skipKernelTC`;
* add the record's name to its row in `records.yml` (field `lean`, with the same `L0`, `r` and
  `beta` as the theorem).

Continuous integration builds the library, audits the axioms of every record, re-checks it in
the kernel, and writes `lean_records.json`, from which the leaderboard takes the Lean tier.
**Reviewers check the Hamiltonian:** the audit proves that the theorem is about the spectrum it
names, but only a human can confirm that this spectrum is the physical model in the row.

## (c) Conditional certification

The most common route combines (a) and (b): a Lean theorem proves the gap from the two seed
inequalities (tier *Lean, conditional*), and a computer certificate in `models/<slug>/`
establishes those inequalities. The seeds named in the Lean hypothesis must be exactly the
defects stored in the certificate, and the reviewer checks that they match.

## (d) Proposing a target

If you know a model where the method should apply, open an issue with the
[target template]({{ repo_url }}/issues/new?template=target.yml): the Hamiltonian, a
reflection-positive basis, any floating-point estimate of the seed scale $$(n_0, \beta_0)$$, and
the known gap. Targets appear in the “Open targets” table.

## Adding or updating a row

Edit [`home_page/_data/records.yml`]({{ repo_blob }}home_page/_data/records.yml); its header
documents the schema. Each row has an `id`, a `family`, the `parameter`, the covered `lengths`,
the exact `bound` (either `{r, beta}`, meaning $$\Delta_L > \log(r)/\beta$$, or `{value}`) with
a `decimal` rounded **down**, the `seed`, the `tier`, the `certificate` and `check` paths, and
optionally a `lean` record, a `numerical_gap` for comparison, `notes`, and a `history` list of
earlier bounds. When a bound improves, move the old one into `history`. Then run

```sh
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl
```

## What continuous integration checks

| check | workflow |
|---|---|
| the Lean build, including `lake build Leaderboard` | `blueprint.yml` |
| `lake exe leaderboard_audit`: axioms of every record, written to `lean_records.json` | `blueprint.yml` |
| `lake exe leaderboard_recheck`: kernel re-check of every record | `blueprint.yml` |
| `check_records.jl`: schema and tiers, paths, rounded-down decimals (exact arithmetic), agreement with the Lean records, and a scan of `leaderboard/` for `sorry`, `native_decide`, `axiom`, `implemented_by`, `extern`, `debug.skipKernelTC` | `blueprint.yml` |
| the blueprint, the declaration check and the site build | `blueprint.yml` |
| `check.jl` (fast mode) of every model directory a pull request changes | `certificates.yml` |

## What reviewers check

* The Hamiltonian, the parameter point and the unit cell in the row match the model directory
  and, for Lean records, the spectrum named in the theorem.
* The reflection-positivity identity is exact, and the trusted base stated in the model's
  `README.md` is acceptable.
* The claimed tier is the right one (a certificate without independent review is *unreviewed*;
  imported estimates are *assumed*).
* The numerical comparison gap has a stated source and is labeled as a numerical estimate.
