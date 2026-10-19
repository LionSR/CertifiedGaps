---
usemathjax: true
---

{% assign repo_url = 'https://github.com/' | append: site.repository %}

<p>In 1983 Haldane predicted that the antiferromagnetic spin-1 Heisenberg chain has a
spectral gap that stays open however long the chain becomes [<a href="#ref-1">1</a>,
<a href="#ref-2">2</a>]. In October 2026 OpenAI released two preprints that prove it, for periodic
rings [<a href="#ref-3">3</a>] and for open chains with fields at the ends
[<a href="#ref-4">4</a>]. Both proofs combine one finite computation, at a fixed ring length and
temperature, with an argument that doubles the length and the inverse temperature together, step by
step, out to every length.</p>

We state that argument as a criterion for general spin chains and ladders, which we call the
*spacetime bootstrap*, prove the criterion in Lean 4 [[7](#ref-7)], and carry out the finite
computation, in exact or interval arithmetic, for several further models. Each result is a theorem
that holds for every even ring length beyond a stated one.

## The criterion

Take a ring of $$L$$ unit cells with Hamiltonian $$H_L$$ and partition function
$$Z_L(\beta) = \mathrm{Tr}\, e^{-\beta H_L}$$. The spacetime bootstrap compares a ring with one
twice as long, and a temperature with one half as high, through two ratios:

$$
P_{\mathrm{s}}(n, \beta) = \frac{Z_{2n}(\beta)}{Z_n(\beta)^2}, \qquad
P_\tau(n, \beta) = \frac{Z_n(2\beta)}{Z_n(\beta)^2}.
$$

When the interaction between neighboring cells is reflection positive [[6](#ref-6)], the defects
$$1 - P_{\mathrm{s}}$$ and $$1 - P_\tau$$ shrink quadratically each time $$n$$ and $$\beta$$ are
doubled together. It is therefore enough to know them at a single seed $$(n_0, \beta_0)$$: if
$$1 - P_{\mathrm{s}}(n_0, \beta_0)$$ and $$1 - P_\tau(2n_0, \beta_0)$$ are both below $$0.127$$,
every even ring with $$L \ge n_0$$ has a unique ground state and a gap bounded below by an explicit
number. The seed is the only numerical input to the proof.

## Results

For the spin-1 Heisenberg chain we prove $$\Delta_L > 0.0729$$ for every even $$L \ge 36$$; the
numerical value of the gap is about $$0.4105$$ [[5](#ref-5)]. The finite inputs are those of the
periodic preprint [[3](#ref-3)]. We reproduce them in exact integer arithmetic, combine them more
efficiently, and prove in Lean 4 that the bound follows from them; an independent audit in exact
rational and in ball arithmetic confirms it. The main bounds stated in the preprints,
$$\Delta_L > (\log 20)/784$$ for even $$L \ge 2304$$ on rings and $$\Delta_n > (\log 10)/392$$ for
odd $$n \ge 1921$$ on open chains [[4](#ref-4)], are likewise proved in Lean 4 from their finite
inputs.

For the Ising chain $$H_L = -\sum_j Z_j Z_{j+1} - \tfrac14 \sum_j X_j - \tfrac1{20} \sum_j Z_j$$,
Lean 4 proves $$\Delta_L > 0.392$$ for every even $$L \ge 6$$ with no hypotheses at all: the finite
computation at the seed $$(n_0, \beta_0) = (6, 97/4)$$ is itself checked inside the Lean kernel.

The [leaderboard]({{ '/leaderboard.html' | relative_url }}) lists the other models, each at
several couplings: spin-1/2 chains with a staggered field or with alternating bonds, the dimerized
spin-1 chain, the spin-1 chain with single-ion anisotropy, and two-leg Heisenberg and Ising ladders.
It also lists the [open problems]({{ '/leaderboard.html#challenges' | relative_url }}).

## Reading

The [explainer]({{ '/explainer/' | relative_url }}) walks through the proof for periodic rings
with pictures and short films, and a [second part]({{ '/explainer/open.html' | relative_url }})
treats open chains. A [two-page note]({{ '/notes/spin1-gap-tn-diagrams.pdf' | relative_url }}) by
S. Lu and GPT-6-Astra writes the spin-1 proof in tensor-network diagrams. The Lean 4 formalization
has a [blueprint]({{ '/blueprint/' | relative_url }}) (also as a
[PDF]({{ '/blueprint.pdf' | relative_url }})), a
[dependency graph]({{ '/blueprint/dep_graph_document.html' | relative_url }}), and
[API documentation]({{ '/docs/' | relative_url }}). To add or improve a bound, see the
[contribution guide]({{ '/contribute.html' | relative_url }}); code, certificates and proofs are on
[GitHub]({{ repo_url }}).

## News

{%- assign news = site.data.news | sort: "date" | reverse %}
<ul class="news">
{%- for n in news limit: 5 %}{% include news-item.html item=n %}{% endfor %}
</ul>
<details><summary>Earlier</summary>
<ul class="news">
{%- for n in news offset: 5 %}{% include news-item.html item=n %}{% endfor %}
</ul>
</details>

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
