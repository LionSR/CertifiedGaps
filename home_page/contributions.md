---
title: Community contributions
description: Contributor credit and verification status for additions to the CertifiedGaps repository
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}

The [leaderboard]({{ '/leaderboard.html' | relative_url }}) includes repository contributions
beyond the results presented in the paper. This page credits those additions separately; it does
not extend the manuscript's claims. The linked model documentation and certificates specify the
Hamiltonians, length conventions, arithmetic assumptions, and reproduction commands.

## Maximilian Lutz

[Maximilian Lutz (lm1909)](https://github.com/lm1909) contributed the following work, merged on
November 11–12, 2026. The free-fermion PRs disclose work done with Claude Code; see each PR for its
implementation and test report. A merged PR is not by itself a mathematical verification claim.

| Contribution | Results and scope | Source |
|---|---|---|
| Coupled seeds and rigorous Rayleigh–Ritz bounds | Reusable low-eigenvalue upper bounds give partition-function lower bounds on the doubled ring. Dimerized spin-½ Heisenberg bounds improve at δ = 1/2, 3/10, 1/4, with a new certificate at δ = 1/5. | [PR #107]({{ repo_url }}/pull/107), [model]({{ repo_url }}/tree/main/code/models/dimerized_heisenberg) |
| Two-leg Ising ladder | Improved g = 3 bound; new certificates at g = 11/4 and 4. Covers every even length L ≥ 4 rungs. | [PR #108]({{ repo_url }}/pull/108), [model]({{ repo_url }}/tree/main/code/models/two_leg_tfim) |
| Odd seeds in Lean | Proves transfer-operator positivity and the odd-seed criterion for the dimerized spin-1 chain, including its δ = 7/10 record. | [PR #110]({{ repo_url }}/pull/110), [Lean source]({{ repo_url }}/tree/main/code/leaderboard/Leaderboard/OddSeed) |
| Free-fermion examples | Computer-certified TFIM, XY-field, and cluster-chain bounds from interval-enclosed closed-form partition functions. Each model verifies reflection positivity, including an explicit cell-frame change for the cluster chain. | [PR #111]({{ repo_url }}/pull/111), [#112]({{ repo_url }}/pull/112), [#113]({{ repo_url }}/pull/113) |
| Spin-1 small seed | $$\Delta_L > 0.039036761$$ for every even $$L \ge 60$$, from twisted traces on rings of at most 12 sites plus a stored length-120 integer MPS trial. | [PR #116]({{ repo_url }}/pull/116), [model]({{ repo_url }}/tree/main/code/models/spin1_small_seed) |

## What has been verified

* **Computer-certified bounds are not automatically Lean-verified bounds.** The new coupled-seed,
  free-fermion, and small-seed results are computer-certified. Their checkers use the exact or
  interval arithmetic and rounding assumptions documented by each model.
* **The odd-seed Lean result is conditional on its finite seed inequalities.** Lean checks the
  transfer-operator argument and doubling chain; the external computer certificate supplies the
  finite seed. This is not an end-to-end kernel check of the seed.
* **The improved coupled-seed bounds do not replace the older Lean seed records.** The dimerized
  spin-½ and g = 3 Ising-ladder entries retain their earlier n₀ = 4 Lean certificates. Their newer
  n₀ = 6 computer-certified bounds have a separate verification status.
* **The small-seed trace-table kernel check remains open.** Its rational replay and audited finite
  inputs yield a computer certificate, not an unconditional Lean proof. The small-seed bound is
  not an improvement over the project's 0.0729 spin-1 bound: its contribution is the smaller
  finite trace input, and it also uses the length-120 MPS trial.
* **Fast checking and full replay are different.** Coupled-seed fast checks use stored, hash-pinned
  Ritz bounds; `--full` recomputes them. PR #107 reports a pre-existing cross-platform replay issue
  for older n₀ = 4 certificates, while its new n₀ = 6 full replays passed. The linked PRs preserve
  the detailed test reports and limitations.

For current numerical bounds, length thresholds, certificate links, and contributor handles, use
[the per-model leaderboard]({{ '/leaderboard.html' | relative_url }}). No new numerical calculation
or proof is claimed by this website summary.
