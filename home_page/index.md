---
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}

This project proves lower bounds on the spectral gap of quantum spin chains. Each bound holds for
all even ring lengths above a stated length. A computer checks one finite calculation, and Lean 4
checks the rest of the proof.

Haldane predicted in 1983 that the spin-1 Heisenberg chain has a gap that stays open as the chain
grows [[1](#ref-1), [2](#ref-2)]. OpenAI published a proof in October 2026, as
[result 268](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/CONTENTS.md?plain=1#L6747) of its math release [[3](#ref-3), [4](#ref-4)]. We call the method of this proof the
*spacetime bootstrap*.

The method is not specific to the spin-1 chain. It applies to many other spin chains, and the
[leaderboard]({{ '/leaderboard.html' | relative_url }}) lists the bounds that we have so far. It can do more.
Sharper finite calculations give larger bounds, and an extension to two-dimensional lattices is in
progress.

## Method

Let $$Z_n(\beta) = \mathrm{Tr}\, e^{-\beta H_n}$$ be the partition function of a ring of $$n$$
unit cells. Draw it as a torus: $$n$$ cells around the ring in space, and a circle of length
$$\beta$$ in imaginary time. The method compares tori through two ratios:

<div class="purities">
  <div class="purity time">
    <span class="def">\(\class{t}{P_{\mathrm T}}(n, \beta) = \dfrac{Z_n(2\beta)}{Z_n(\beta)^2} =\)</span>
    <span class="frac"><span>{% include torus.html w=2 h=2 %}</span><span>{% include torus.html w=2 h=1 %}<sup>2</sup></span></span>
    <span class="gloss">thermal purity: two tori glued in time</span>
  </div>
  <div class="purity space">
    <span class="def">\(\class{s}{P_{\mathrm S}}(n, \beta) = \dfrac{Z_{2n}(\beta)}{Z_n(\beta)^2} =\)</span>
    <span class="frac"><span>{% include torus.html w=4 h=1 %}</span><span>{% include torus.html w=2 h=1 %}<sup>2</sup></span></span>
    <span class="gloss">spatial purity: two tori glued in space</span>
  </div>
</div>

<div class="method">
<ol>
  <li>Show that the interaction between neighboring cells is reflection positive
    [<a href="#ref-6">6</a>].</li>
  <li>Pick one ring length \(n_0\) and one inverse temperature \(\beta_0\), the <em>seed</em>.
    Calculate upper bounds on the defects \(p = 1 - P_{\mathrm S}(n_0, \beta_0)\) and
    \(q = 1 - P_{\mathrm T}(2n_0, \beta_0)\). The method works if \((p, q)\) is in the green
    region, for example if both are below \(0.127\).</li>
  <li>Double \(n\) and \(\beta\) together. Each doubling moves \((p, q)\) closer to 0. The result
    is a lower bound on the gap for all even \(L \ge n_0\).</li>
</ol>
<figure>
  {% include seed-region.svg %}
  <figcaption>The spin-1 chain. The path starts at the seed, \(n = 36\). Each step doubles
    \(n\) and \(\beta\).</figcaption>
</figure>
</div>

Only step 2 uses a numerical calculation. Lean 4 [[7](#ref-7)] proves step 3.

## Main results

| Model | Gap bound | Ring lengths | Proof |
|---|---|---|---|
| Spin-1 Heisenberg chain | $$\Delta_L > 0.0729$$ | even $$L \ge 36$$ | Lean 4; seed checked by computer, from the data of [[3](#ref-3)] |
| Ising chain, $$g = 1/4$$, $$h = 1/20$$ | $$\Delta_L > 0.392$$ | even $$L \ge 6$$ | Lean 4, seed included |

The numerical value of the spin-1 gap is approximately $$0.4105$$ [[5](#ref-5)]. The Ising chain
is $$H_L = -\sum_j Z_j Z_{j+1} - g \sum_j X_j - h \sum_j Z_j$$. The
[leaderboard]({{ '/leaderboard.html' | relative_url }}) lists all bounds and the history of the
spin-1 record.

## Links

* Explainer: [periodic rings]({{ '/explainer/' | relative_url }}),
  [open chains]({{ '/explainer/open.html' | relative_url }})
* Two-page note: [a short diagrammatic explanation of the spin-1 Haldane gap proof]({{ '/notes/spin1-gap-tn-diagrams.pdf' | relative_url }}),
  by S. Lu and GPT-6-Astra
* Lean 4: [blueprint]({{ '/blueprint/' | relative_url }}),
  [blueprint PDF]({{ '/blueprint.pdf' | relative_url }}),
  [dependency graph]({{ '/blueprint/dep_graph_document.html' | relative_url }}),
  [API documentation]({{ '/docs/' | relative_url }})
* Contributions: [contribution guide]({{ '/contribute.html' | relative_url }}),
  [source on GitHub]({{ repo_url }})

## News

{%- assign news = site.data.news | sort: "date" | reverse %}
<ul class="news">
{%- for n in news limit: 3 %}{% include news-item.html item=n %}{% endfor %}
</ul>
<details><summary>Earlier news</summary>
<ul class="news">
{%- for n in news offset: 3 %}{% include news-item.html item=n %}{% endfor %}
</ul>
</details>

## Acknowledgments

We thank Xiao-Liang Qi for discussions. Sirui Lu did this work while visiting Xiao-Liang Qi at
Stanford University.

## References

<ol class="references">
<li id="ref-1">F. D. M. Haldane, Continuum dynamics of the 1-D Heisenberg antiferromagnet: identification with the O(3) nonlinear sigma model, Phys. Lett. A <b>93</b>, 464–468 (1983). <a href="https://doi.org/10.1016/0375-9601(83)90631-X">doi:10.1016/0375-9601(83)90631-X</a></li>
<li id="ref-2">F. D. M. Haldane, Nonlinear field theory of large-spin Heisenberg antiferromagnets: semiclassically quantized solitons of the one-dimensional easy-axis Néel state, Phys. Rev. Lett. <b>50</b>, 1153–1156 (1983). <a href="https://doi.org/10.1103/PhysRevLett.50.1153">doi:10.1103/PhysRevLett.50.1153</a></li>
<li id="ref-3">OpenAI, The periodic spin-one Haldane gap, OpenAI Math Release, result 268, preprint (October 6, 2026). <a href="https://github.com/openai/math/blob/main/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/paper.pdf">PDF</a></li>
<li id="ref-4">OpenAI, A boundary-field gap for the spin-one Heisenberg chain, OpenAI Math Release, result 268, preprint (October 6, 2026). <a href="https://github.com/openai/math/blob/main/preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026/paper.pdf">PDF</a></li>
<li id="ref-5">S. R. White and D. A. Huse, Numerical renormalization-group study of low-lying eigenstates of the antiferromagnetic S=1 Heisenberg chain, Phys. Rev. B <b>48</b>, 3844–3852 (1993). <a href="https://doi.org/10.1103/PhysRevB.48.3844">doi:10.1103/PhysRevB.48.3844</a></li>
<li id="ref-6">K. Osterwalder and R. Schrader, Axioms for Euclidean Green's functions, Commun. Math. Phys. <b>31</b>, 83–112 (1973). <a href="https://doi.org/10.1007/BF01645738">doi:10.1007/BF01645738</a></li>
<li id="ref-7">L. de Moura and S. Ullrich, The Lean 4 theorem prover and programming language, in Automated Deduction – CADE 28, Lecture Notes in Computer Science <b>12699</b>, 625–635 (2021). <a href="https://doi.org/10.1007/978-3-030-79876-5_37">doi:10.1007/978-3-030-79876-5_37</a></li>
</ol>
