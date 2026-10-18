---
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign data = site.data.records %}

In 1983, Haldane conjectured that the spin-1 antiferromagnetic Heisenberg chain
$$H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1}$$ has a uniform spectral gap
$$\Delta_L = E_1(L) - E_0(L) > 0$$ in the thermodynamic limit [[1](#ref-1), [2](#ref-2)]. In
October 2026, OpenAI released two preprints proving uniform gap bounds for the periodic spin-1
ring [[3](#ref-3)] and for open spin-1 chains with endpoint fields [[4](#ref-4)] by combining
finite-scale numerical certificates with a coupled spatial-and-thermal doubling argument. This
project reformulates that mechanism as a general **spacetime bootstrap** (*two-purity doubling
criterion*) using continuous-time tensor networks, reproduces the spin-1 finite certificates in
exact integer and rational arithmetic, formalizes the analytic proofs in Lean 4 [[7](#ref-7)], and
applies the method to certify uniform spectral gaps across eight families of one-dimensional
quantum spin chains and ladders.

For a periodic ring Hamiltonian $$H_L$$ of $$L$$ unit cells with partition function
$$Z_L(\beta) = \operatorname{Tr} e^{-\beta H_L}$$, define the **spatial purity**
$$P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$$ and the **thermal purity**
$$P_\tau(n, \beta) = Z_n(2\beta)/Z_n(\beta)^2$$. When the nearest-neighbor bond interaction is
reflection positive [[6](#ref-6)], the spatial transfer operator $$X_\beta$$
($$Z_L(\beta) = \operatorname{Tr} X_\beta^L$$) is self-adjoint, so $$P_{\mathrm{s}}$$ and
$$P_\tau$$ place symmetric constraints on the transfer spectrum in space and the energy spectrum
in imaginary time. If the seed purity defects $$\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$$
and $$\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$$ lie in the admissible region
($$\max(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) < 0.12731$$) at a single even seed length
$$n_0$$ and inverse temperature $$\beta_0$$, a doubling recurrence simultaneously doubles $$n$$
and $$\beta$$ across all dyadic scales, proving that every even ring $$L \ge n_0$$ has a unique
ground state and an explicit lower bound on $$\Delta_L$$. Each record pairs a computer-checked seed
certificate at $$(n_0, \beta_0)$$ with a Lean 4 theorem covering the step from the seed to all
lengths $$L \ge n_0$$.

## Headline results

* **Spin-1 Heisenberg chain** (Haldane gap [[1](#ref-1), [2](#ref-2)]): $$\Delta_L > 0.0729$$ for
  every even ring length $$L \ge 36$$ (numerical bulk gap $$\approx 0.4105$$ [[5](#ref-5)]). The
  seed inequalities at $$(n_0, \beta_0) = (60, 833/32)$$ and $$(36, 49/4)$$ are certified in
  exact rational arithmetic from the finite trace and variational inputs of the OpenAI periodic
  preprint [[3](#ref-3)], reproduced here and independently audited; Lean 4 proves the bound from
  those two seeds. We also formalize in Lean 4 both periodic bounds from the preprint
  ($$\Delta_L > (\log 20)/784 \approx 0.00382$$ for even $$L \ge 2304$$ and
  $$\Delta_L > \tfrac{4}{105}\log\tfrac{80}{79} \approx 0.000479$$ for even $$L \ge 60$$) and the
  open spin-1 chain with endpoint fields $$h = 3/5$$ [[4](#ref-4)]
  ($$\Delta_n > (\log 10)/392 \approx 0.00587$$ for odd $$n \ge 1921$$).
* **Unconditional Lean 4 theorem:** For the Ising chain in transverse and longitudinal fields,
  $$H_L = -\sum_j Z_j Z_{j+1} - \tfrac14 \sum_j X_j - \tfrac1{20} \sum_j Z_j$$, Lean 4 proves
  $$\Delta_L > 0.3922$$ for every even $$L \ge 6$$ with no external hypotheses; both the finite
  seed bounds at $$(n_0, \beta_0) = (6, 97/4)$$ and the doubling recurrence are checked inside
  the Lean kernel.
* **{{ data.records | where_exp: "r", "r.tier != 'target'" | size }} certified gap records**
  across eight model families—including spin-1 and spin-1/2 dimerized chains, two-leg Heisenberg
  and Ising ladders, single-ion anisotropy, and staggered fields—all except one accompanied by a
  Lean 4 theorem from the certified seed. Browse the full table on the
  [leaderboard]({{ '/leaderboard.html' | relative_url }}).

## News

{%- comment -%}
  Hand-written items from _data/news.yml, merged with the leaderboard rows that no item already
  announces (by pull request or row id); rows set in one pull request form one item. Items are
  sorted as "date|html" strings joined by ␞, a character that does not occur in the data.
{%- endcomment -%}
{%- assign news = site.data.news -%}
{%- assign seen_prs = "," -%}{%- assign seen_ids = "," -%}
{%- for n in news -%}
  {%- for p in n.prs %}{% assign seen_prs = seen_prs | append: p | append: "," %}{% endfor -%}
  {%- for i in n.records %}{% assign seen_ids = seen_ids | append: i | append: "," %}{% endfor -%}
{%- endfor -%}
{%- assign items = "" -%}
{%- for n in news -%}
  {%- if n.link -%}
    {%- assign first = n.link | slice: 0 -%}
    {%- if first == "/" %}{% assign href = n.link | relative_url %}{% else %}{% assign href = n.link %}{% endif -%}
  {%- else -%}
    {%- assign href = repo_url | append: "/pull/" | append: n.prs[0] -%}
  {%- endif -%}
  {%- capture item %}{{ n.date }}|<a href="{{ href }}">{{ n.text }}</a>{% endcapture -%}
  {%- assign items = items | append: item | append: "␞" -%}
{%- endfor -%}
{%- assign groups = data.records | where_exp: "r", "r.tier != 'target'" | group_by: "pr" -%}
{%- for g in groups -%}
  {%- for r in g.items -%}
    {%- capture id_key %},{{ r.id }},{% endcapture %}{% capture pr_key %},{{ r.pr }},{% endcapture -%}
    {%- if seen_ids contains id_key or r.pr and seen_prs contains pr_key %}{% continue %}{% endif -%}
    {%- if g.name != "" and forloop.index > 1 %}{% continue %}{% endif -%}
    {%- assign best = r -%}
    {%- if g.name != "" -%}
      {%- for x in g.items %}{% assign a = x.bound.decimal | times: 1.0 %}{% assign b = best.bound.decimal | times: 1.0 %}{% if a > b %}{% assign best = x %}{% endif %}{% endfor -%}
    {%- endif -%}
    {%- assign fam = data.families | where: "key", best.family | first -%}
    {%- assign tkey = best.tier -%}
    {%- if best.lean and tkey contains "lean-" -%}
      {%- assign lr = site.data.lean_records | where: "name", best.lean.name | first -%}
      {%- if lr %}{% assign tkey = "lean-" | append: lr.kind %}{% else %}{% assign tkey = "" %}{% endif -%}
    {%- endif -%}
    {%- assign t = data.tiers | where: "key", tkey | first -%}
    {%- capture href %}{{ '/leaderboard.html' | relative_url }}#{{ best.id }}{% endcapture -%}
    {%- capture item -%}
      {{ best.date }}|<a href="{{ href }}">{% if g.name != "" and g.items.size > 1 %}New records: {{ fam.name }} ({{ g.items.size }} parameter points, best bound Δ ≥ {{ best.bound.decimal }}){% else %}New record: {{ fam.name }}, {{ best.parameter }}, Δ ≥ {{ best.bound.decimal }}{% endif %}</a>{% if t %} ({{ t.label }}){% endif %}
    {%- endcapture -%}
    {%- assign items = items | append: item | append: "␞" -%}
  {%- endfor -%}
{%- endfor -%}
{%- assign items = items | split: "␞" | sort | reverse -%}
{%- assign cutoff = site.time | date: "%s" | minus: 1209600 | date: "%Y-%m-%d" -%}
{%- capture list -%}
{%- for it in items -%}
  {%- assign d = it | split: "|" | first -%}
  <li><span class="date">{{ d }}</span>{{ it | remove_first: d | remove_first: "|" }}{% if d >= cutoff %} <span class="new">new</span>{% endif %}</li>␞
{%- endfor -%}
{%- endcapture -%}
{%- assign lines = list | split: "␞" %}
<ul class="news">
{%- for l in lines limit: 5 %}{{ l }}{% endfor %}
</ul>
{%- if lines.size > 5 %}
<details><summary>More news</summary>
<ul class="news">
{%- for l in lines offset: 5 %}{{ l }}{% endfor %}
</ul>
</details>
{%- endif %}

See the [open challenges]({{ '/leaderboard.html' | relative_url }}#challenges) for concrete targets and next steps.

## Documentation and resources

* **Interactive explainers:** [periodic ring proof]({{ '/explainer/' | relative_url }}) ·
  [open-chain extension]({{ '/explainer/open.html' | relative_url }})
* **Lean 4 formalization:** [blueprint]({{ '/blueprint/' | relative_url }})
  ([PDF]({{ '/blueprint.pdf' | relative_url }})) ·
  [dependency graph]({{ '/blueprint/dep_graph_document.html' | relative_url }}) ·
  [API documentation]({{ '/docs/' | relative_url }})
* **Contribute:** [leaderboard]({{ '/leaderboard.html' | relative_url }}) ·
  [contribution guide]({{ '/contribute.html' | relative_url }}) ·
  [GitHub repository]({{ repo_url }})

## References

<ol class="references">
<li id="ref-1">F. D. M. Haldane, Continuum dynamics of the 1-D Heisenberg antiferromagnet: identification with the O(3) nonlinear sigma model, Phys. Lett. A <b>93</b>, 464–468 (1983). <a href="https://doi.org/10.1016/0375-9601(83)90631-X">doi:10.1016/0375-9601(83)90631-X</a></li>
<li id="ref-2">F. D. M. Haldane, Nonlinear field theory of large-spin Heisenberg antiferromagnets: semiclassically quantized solitons of the one-dimensional easy-axis Néel state, Phys. Rev. Lett. <b>50</b>, 1153–1156 (1983). <a href="https://doi.org/10.1103/PhysRevLett.50.1153">doi:10.1103/PhysRevLett.50.1153</a></li>
<li id="ref-3">OpenAI, The periodic spin-one Haldane gap, OpenAI Math Release preprint (October 6, 2026). <a href="https://github.com/openai/math/blob/main/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/paper.pdf">PDF</a></li>
<li id="ref-4">OpenAI, A boundary-field gap for the spin-one Heisenberg chain, OpenAI Math Release preprint (October 6, 2026). <a href="https://github.com/openai/math/blob/main/preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026/paper.pdf">PDF</a></li>
<li id="ref-5">S. R. White and D. A. Huse, Numerical renormalization-group study of low-lying eigenstates of the antiferromagnetic S=1 Heisenberg chain, Phys. Rev. B <b>48</b>, 3844–3852 (1993). <a href="https://doi.org/10.1103/PhysRevB.48.3844">doi:10.1103/PhysRevB.48.3844</a></li>
<li id="ref-6">K. Osterwalder and R. Schrader, Axioms for Euclidean Green's functions, Commun. Math. Phys. <b>31</b>, 83–112 (1973). <a href="https://doi.org/10.1007/BF01645738">doi:10.1007/BF01645738</a></li>
<li id="ref-7">L. de Moura and S. Ullrich, The Lean 4 theorem prover and programming language, in Automated Deduction – CADE 28, Lecture Notes in Computer Science <b>12699</b>, 625–635 (2021). <a href="https://doi.org/10.1007/978-3-030-79876-5_37">doi:10.1007/978-3-030-79876-5_37</a></li>
</ol>
