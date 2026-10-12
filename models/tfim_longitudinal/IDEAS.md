# Candidate routes

One site is one cell. The theorem requires even seed n; odd n will be diagnostic only.

- Full dense ED: exact matrix entries are rational. Exponential storage and cubic time
  limit the largest ring; start small because a longitudinal field can make small seeds work.
- Translation and reflection: translation sectors have dimension approximately 2^L/L.
  Reflection maps momentum k to -k, so independent reflection parity only exists at
  k=0 and pi. Spin-flip parity is broken for h>0 and must never be used. Dense sectors
  at L=20--24 can still exceed available memory/time. Measure before extending.
- Rigorous eigenvalue enclosures: certify complete dense eigensystems using a residual
  and orthogonality defect, including arithmetic rounding. An exact rational computational
  basis avoids algebraic errors in symmetry basis transformations. Use this for the
  best accessible seed, with rigorously bounded dot-product errors for speed.
- Free fermions: only h=0 is solvable by the supplied formula; use g=2 as a cross-check,
  never substitute it for longitudinal-field data.
- Cluster expansions or transfer-kernel contraction could reach larger sizes, but require
  a model-specific quantitative remainder absent from this repository. Not assumed here.

Chosen route: numerical complete-spectrum ED survey, then rigorous enclosure of the best
seed per parameter. Increase sizes only when cost permits. Report actual reached sizes,
failures, and finite-size diagnostics rather than claiming the requested 20--24 by default.
Uniqueness at every finite length follows independently by Perron--Frobenius: g>0 gives
negative edges on the connected hypercube, so a scalar shift of -H is irreducible nonnegative.
This alone does not establish a uniform gap.
