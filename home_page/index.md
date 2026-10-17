---
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}
{% assign data = site.data.records %}

A quantum spin chain has an energy *gap* when its lowest energy level stays separated from the
next one however long the chain is. Haldane predicted such a gap for the spin-1 Heisenberg
chain [[1](#ref-1), [2](#ref-2)]. A computer can only treat short chains, so proving a gap for
every length needs an argument that turns one finite computation into a statement about all
lengths.

This project gives such an argument, the *spacetime bootstrap*. It compares the partition
functions $$Z_L(\beta) = \operatorname{Tr} e^{-\beta H_L}$$ of rings of length $$L$$ through two
ratios, the spatial purity $$Z_{2n}(\beta)/Z_n(\beta)^2$$ and the thermal purity
$$Z_n(2\beta)/Z_n(\beta)^2$$. If the interaction is reflection positive [[6](#ref-6)] and both
ratios are close to 1 for one ring $$n_0$$ at one inverse temperature $$\beta_0$$, then every
even ring from $$n_0$$ on has a unique ground state and a gap above an explicit number. We check
that finite computation by computer in exact or interval arithmetic, and we prove the step from
it to all lengths in the Lean proof assistant [[7](#ref-7)].

## Headline results

* **The spin-1 Heisenberg chain.** The gap is at least 0.0729 for every even ring of at least 36
  sites. The starting inequalities come from the finite computations of the OpenAI preprint
  [[3](#ref-3)], reproduced here in exact arithmetic and independently reviewed; Lean proves the
  rest. The preprint itself gives $$\log(20)/784 \approx 0.00382$$ for rings of at least 2304
  sites. The numerical value of the gap is about 0.41 [[5](#ref-5)]. A companion preprint treats
  the open chain with fields at its ends [[4](#ref-4)].
* **A gap proved in Lean with no assumptions.** For the Ising ring
  $$-\sum_j Z_j Z_{j+1} - \tfrac14 \sum_j X_j - \tfrac1{20} \sum_j Z_j$$ the gap is at least
  0.3922 for every even length of at least 6. Lean also checks the finite computation.
* **Lean checks almost every result.** For every proven gap except one, Lean proves the step
  from the computed starting inequalities to all ring lengths.
* **{{ data.records | where_exp: "r", "r.tier != 'target'" | size }} proven gaps** for eight model
  families, from spin ladders to chains in staggered fields. See the
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
      {{ best.date }}|<a href="{{ href }}">{% if g.name != "" and g.items.size > 1 %}New records: {{ fam.name }}, {{ g.items.size }} points, the best with a proven gap of at least {{ best.bound.decimal }}{% else %}New record: {{ fam.name }}, {{ best.parameter }}, proven gap at least {{ best.bound.decimal }}{% endif %}</a>{% if t %} ({{ t.label }}){% endif %}
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

Want to set the next record? See the [open challenges]({{ '/leaderboard.html' | relative_url }}#challenges).

## Read

* [Explainer]({{ '/explainer/' | relative_url }}): an illustrated walk through the proof for
  undergraduates, with a sequel on [open chains]({{ '/explainer/open.html' | relative_url }}).
{% comment %} The note PDFs are temporary: the paper will replace these links. {% endcomment %}
* Notes (PDF): [main note]({{ '/notes/thermal-mpo-gap-note.pdf' | relative_url }}) ·
  [scope]({{ '/notes/thermal-mpo-gap-scope.pdf' | relative_url }}) ·
  [certificate]({{ '/notes/thermal-mpo-gap-certificate.pdf' | relative_url }}) ·
  [2+1 dimensions]({{ '/notes/thermal-mpo-gap-2plus1d.pdf' | relative_url }})

## Lean

* [Blueprint]({{ '/blueprint/' | relative_url }}) ([PDF]({{ '/blueprint.pdf' | relative_url }})) ·
  [dependency graph]({{ '/blueprint/dep_graph_document.html' | relative_url }}) ·
  [documentation]({{ '/docs/' | relative_url }})

## Take part

* [Leaderboard]({{ '/leaderboard.html' | relative_url }}) ·
  [how to contribute]({{ '/contribute.html' | relative_url }}) ·
  [GitHub repository]({{ repo_url }})

## References

<ol class="references">
<li id="ref-1">F. D. M. Haldane, Continuum dynamics of the 1-D Heisenberg antiferromagnet: identification with the O(3) nonlinear sigma model, Phys. Lett. A <b>93</b>, 464–468 (1983). <a href="https://doi.org/10.1016/0375-9601(83)90631-X">doi:10.1016/0375-9601(83)90631-X</a></li>
<li id="ref-2">F. D. M. Haldane, Nonlinear field theory of large-spin Heisenberg antiferromagnets: semiclassically quantized solitons of the one-dimensional easy-axis Néel state, Phys. Rev. Lett. <b>50</b>, 1153–1156 (1983). <a href="https://doi.org/10.1103/PhysRevLett.50.1153">doi:10.1103/PhysRevLett.50.1153</a></li>
<li id="ref-3">OpenAI, The periodic spin-one Haldane gap, preprint (September 24, 2026). <a href="https://github.com/openai/math/blob/main/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/paper.pdf">PDF</a></li>
<li id="ref-4">OpenAI, A boundary-field gap for the spin-one Heisenberg chain, preprint (September 24, 2026). <a href="https://github.com/openai/math/blob/main/preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026/paper.pdf">PDF</a></li>
<li id="ref-5">S. R. White and D. A. Huse, Numerical renormalization-group study of low-lying eigenstates of the antiferromagnetic S=1 Heisenberg chain, Phys. Rev. B <b>48</b>, 3844–3852 (1993). <a href="https://doi.org/10.1103/PhysRevB.48.3844">doi:10.1103/PhysRevB.48.3844</a></li>
<li id="ref-6">K. Osterwalder and R. Schrader, Axioms for Euclidean Green's functions, Commun. Math. Phys. <b>31</b>, 83–112 (1973). <a href="https://doi.org/10.1007/BF01645738">doi:10.1007/BF01645738</a></li>
<li id="ref-7">L. de Moura and S. Ullrich, The Lean 4 theorem prover and programming language, in Automated Deduction – CADE 28, Lecture Notes in Computer Science <b>12699</b>, 625–635 (2021). <a href="https://doi.org/10.1007/978-3-030-79876-5_37">doi:10.1007/978-3-030-79876-5_37</a></li>
</ol>
