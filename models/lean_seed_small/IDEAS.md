# Candidate routes: small seeds in Lean

All numerical searches are untrusted and non-certifying. The deliverable must separate
checked integer/rational facts from the analytic bridge to the original `RPRing`.

1. **Full dense eigen-certificate.** For the longitudinal Ising grid the largest matrix is
   4096 by 4096; for spin-one dimerization it is 6561 by 6561. Exact Gram and residual
   certificates avoid trusting BLAS, but full eigenvector storage is quadratic and checking
   is cubic before packing. Translation reduction needs an additional basis proof.
2. **Low-energy Ritz plus a certified anchor.** Cost depends on the number of levels below
   `(log dimension + log(1/tail_budget))/beta` above the ground energy. Inspect the existing
   spectra before choosing. An unproved anchor is not a certificate, even if the residual
   is tiny. Packed dot products can reduce kernel overhead, not the analytic proof burden.
3. **Exact characteristic polynomial / inertia on the smaller ring.** Rational symmetric
   matrices allow fraction-free exact certificates, without floating point or a generator
   correctness proof. Useful as a first trustworthy component; 4096-dimensional generic
   full characteristic polynomials remain expensive.
4. **Perturbation from independent spins or dimers.** Cheap operator bounds can eliminate
   large spectral computations, but may lose the small finite-size cancellation needed
   for the exact existing spatial seed. Test the error budget rather than assume it works.
5. **Partition enclosure interface.** Prove the exact seed inequalities from explicit
   partition-function enclosures and kernel-check the rational arithmetic. This is useful
   partial progress only: the interface hypotheses must remain visible until their matrix
   exponential semantics are proved. No theorem with these hypotheses is an unconditional seed.

Initial choice: inventory exact seed targets and low-energy counts, then use the cheapest
small-ring model for an exact certificate component and the original-RPRing interface.
Prefer a concrete counterexample to a proposed shortcut over an optimistic runtime estimate.
Keep all source, downloads, caches, and build outputs inside this directory; never edit `lean/`.

## Outcome

Route 2 (low-energy Ritz plus a classical anchor) at the longitudinal Ising point
`g = 1/4, h = 1/20, (n₀, β₀) = (6, 97/4)`, with only two levels per ring, plus route 5 (the
partition-enclosure interface). The integer Gram/residual/anchor arithmetic is kernel-checked
at 6 and 12 sites, and the exact seed conclusions for the first grid point are proved
conditional on four `LogEnclosures` hypotheses. No seed is unconditional: the matrix
identification, min–max counting and scalar `exp` enclosures remain (METHOD.md). The
43-point inventory shows that full k×k Gram checks reach hours to ~1,000 hours for the
65,536- and 531,441-dimensional seeds, so routes 1–2 do not scale there without clustering or
symmetry reductions.
