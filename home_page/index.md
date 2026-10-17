---
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign records = site.data.records.records %}

A quantum spin chain is *gapped* when its lowest energy level stays separated from the next one
however long the chain grows. A computer can diagonalize only short chains, so a gap for every
length needs an argument that turns one finite computation into a statement about all lengths.

The **spacetime bootstrap** (the *two-purity doubling criterion*) is that argument. For a
reflection-positive ring it compares the partition functions $$Z_L(\beta) = \operatorname{Tr} e^{-\beta H_L}$$
through the **spatial purity** $$P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$$ and the
**thermal purity** $$P_\tau(n, \beta) = Z_n(2\beta)/Z_n(\beta)^2$$. If both are close enough to 1
at one *seed* scale $$(n_0, \beta_0)$$, doubling length and inverse temperature together keeps
them close to 1 at every scale, and every even ring $$L \ge n_0$$ has a unique ground state and a
gap above an explicit constant. We certify the seeds by computer in exact or interval arithmetic
and prove the step from seed to gap in Lean 4.

## Headline results

* **Unconditional, in Lean:** the Ising ring
  $$-\sum_j Z_j Z_{j+1} - \tfrac14 \sum_j X_j - \tfrac1{20} \sum_j Z_j$$ has a unique ground
  state and $$\Delta_L > \log(13480)/(97/4) \approx 0.392$$ for every even $$L \ge 6$$; the seed
  is checked in the Lean kernel.
* **The spin-1 Heisenberg chain** (Haldane's conjecture): $$\Delta_L > \log(20)/784 \approx 0.00382$$
  for every even $$L \ge 2304$$, proved in Lean from the seeds $$P_{\mathrm{s}}(2304, 784) \ge 0.99$$
  and $$P_\tau(4608, 784) \ge 0.99$$ of the
  [OpenAI preprint](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026),
  which this repository reproduces in exact arithmetic; an exact computation not yet reviewed
  gives $$\Delta_L > 0.0729$$ for every even $$L \ge 36$$.
* **{{ records | size }} leaderboard entries** for eight certified model families and open
  targets: see the [leaderboard]({{ '/leaderboard.html' | relative_url }}).

## Read

* [Explainer]({{ '/explainer/' | relative_url }}): an illustrated, undergraduate-level walk
  through the proof (its films are rendered locally).
* Notes (PDF): [main note]({{ '/notes/thermal-mpo-gap-note.pdf' | relative_url }}) ·
  [scope]({{ '/notes/thermal-mpo-gap-scope.pdf' | relative_url }}) ·
  [certificate]({{ '/notes/thermal-mpo-gap-certificate.pdf' | relative_url }}) ·
  [2+1 dimensions]({{ '/notes/thermal-mpo-gap-2plus1d.pdf' | relative_url }})

## Lean

* [Blueprint]({{ '/blueprint/' | relative_url }}) ([PDF]({{ '/blueprint.pdf' | relative_url }})) ·
  [dependency graph]({{ '/blueprint/dep_graph_document.html' | relative_url }}) ·
  [API docs]({{ '/docs/' | relative_url }})

## Take part

* [Leaderboard]({{ '/leaderboard.html' | relative_url }}) ·
  [how to contribute]({{ '/contribute.html' | relative_url }}) ·
  [GitHub repository]({{ repo_url }})

Every number carries a verification tier: a Lean theorem (unconditional, or conditional on two
computer-certified seed inequalities), a rigorous computer certificate with a standalone checker,
or a floating-point estimate for comparison only. The rules are in
[`AGENTS.md`]({{ repo_url }}/blob/{{ site.branch }}/AGENTS.md).
