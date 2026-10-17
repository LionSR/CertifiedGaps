---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

# layout: home
usemathjax: true
---

This project formalizes in Lean 4 the **spacetime bootstrap** (the *two-purity doubling
criterion*): a reflection-positive quantum spin ring whose partition function
Z<sub>L</sub>(β) = Tr e<sup>−βH<sub>L</sub></sup> has spatial purity
P<sub>s</sub>(n₀, β₀) = Z<sub>2n₀</sub>(β₀)/Z<sub>n₀</sub>(β₀)² and thermal purity
P<sub>τ</sub>(2n₀, β₀) = Z<sub>2n₀</sub>(2β₀)/Z<sub>2n₀</sub>(β₀)² close enough to 1 at a single
seed scale (n₀, β₀) has a unique ground state and a uniform spectral gap on every even ring of
length L ≥ n₀.

## Headline theorems

* **Unconditional:** for every even L ≥ 6, the Ising ring
  −∑ Z<sub>j</sub>Z<sub>j+1</sub> − ¼ ∑ X<sub>j</sub> − ¹⁄₂₀ ∑ Z<sub>j</sub> has a unique ground
  state and gap Δ<sub>L</sub> > log 13480 / (97/4) ≈ 0.392
  (`SeedKernel.Ising.ising_gap`; the seed is checked in the Lean kernel, no `sorry`).
* **Conditional on two seed inequalities:** for every even L ≥ 2304, the spin-1 Heisenberg ring
  ∑ **S**<sub>j</sub> · **S**<sub>j+1</sub> has a unique ground state and
  Δ<sub>L</sub> > log 20 / 784 ≈ 0.00382, given P<sub>s</sub>(2304, 784) ≥ 0.99 and
  P<sub>τ</sub>(4608, 784) ≥ 0.99
  (`CertifiedGaps.SpinOneHeisenberg.spin_one_heisenberg_gap`).
* Seven further model families (spin-1 with single-ion anisotropy, dimerized spin-1 and
  spin-½ chains, spin-½ chain in a staggered field, spin-½ and transverse-field Ising ladders,
  Ising chain in transverse and longitudinal fields), each with its two seed inequalities as the
  only open hypotheses. The seeds are certified by computer in the repository's `models/`.

Useful links:

* [Blueprint]({{ site.url }}/blueprint/)
* [Blueprint as pdf]({{ site.url }}/blueprint.pdf)
* [Dependency graph]({{ site.url }}/blueprint/dep_graph_document.html)
* [Doc pages for this repository]({{ site.url }}/docs/)
