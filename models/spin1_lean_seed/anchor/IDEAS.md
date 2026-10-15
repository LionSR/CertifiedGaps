# Anchor routes (written before experiments)

All floating-point spectra and optimizations here are NON-CERTIFYING. Julia exact rational
checks establish finite identities only; they are not Lean kernel certificates.

1. Two disjoint open halves: replace each cut bond by -2 I. This has a tensor-sum
   spectrum, tiny exact local PSD checks, and a likely surface-energy penalty of order two.
2. Optimize two-site SU(2)-invariant boundary corrections. For each cut use a four-site
   remainder h23 - B12 - B34 - c I >= 0. The halves receive B at both ends plus
   c per cut. All local checks fit 81 by 81, and half-chain Sz blocks have at most
   141 states for length six. Optimize bilinear and biquadratic coefficients, then
   round rationally and perform exact PSD checks. Compare ground optimization with
   the actual spectral-count objective at the thermal cutoff.
3. Three or more disjoint blocks: same construction, more boundaries and smaller
   block spectra; probably a worse count. Measure this rather than assume it.
4. Overlapping Anderson clusters: good scalar ground bounds, but overlapping
   operators do not have a tensor-sum spectrum. Test scalar bounds and distinguish
   ground-energy certificates from usable excited-eigenvalue anchors.
5. Effective edge spin halves: a projection alone gives Ritz upper bounds. A lower
   comparison needs certified discarded-space and coupling control. Do not count
   an unproved effective model as a valid anchor.
6. If a one-shot count is too large, interpolate from the block comparison to H.
   The difference is PSD with an exactly bounded norm; propagate counting anchors
   using monotonicity and refine using residual/inertia or Lehmann bounds. Test
   stage counts and retained-space sizes, explicitly distinguishing a numerical
   feasibility study from a certified continuation proof.

Chosen first: routes 1 and 2, with complete block spectra and ring low spectra in Sz
sectors, then quantitative comparison against routes 3 and 6. New files stay in this
folder. No changes to lean/ or models/spin1/results.
