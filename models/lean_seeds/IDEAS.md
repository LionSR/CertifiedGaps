# Toward a Lean-checked seed: route memos

Five route memos, each written before its study's computations (originally in separate
directories: `spin1_lean_seed/{seed_design,spectrum,anchor,lean_prototype}` and
`lean_seed_small`). They are kept as written, apart from headings; "this directory" meant the
study's own directory, and the Lean prototypes now live in `lean/SeedKernel`. The common plan:
Julia is an untrusted generator, and the Lean kernel (`decide +kernel`, no `native_decide`, no
new axioms) checks every inequality. The outcomes of all five are collected at the end.

## seed_design: which finite inputs does the spin-1 seed need?

All floating-point calculations are non-certifying. Exact scalar checks remain
conditional on the thermal enclosures, the trial inequality, and analytic transfer
lemmas; Julia is not the Lean trusted kernel.

1. Restrict existing source-1 data to even lengths through 10 or 8. Optimize
   nonnegative filters x^(2s) q(x)^2, including even q, and preserve the threefold
   N multiplicity using J=(Z+3P)/4 and N=(Z-P)/4. This avoids complex C traces.
   Cost is dominated by 3^n times the low-state count, with smaller polynomial
   degree reducing both table count and scalar proof complexity. Risk: the cap
   on J is too far above one, and its power grows at the trial length.
2. Compare the full I/O/N split, using all three twists, with the coarser J/N
   split and a twist-free bound. It may tighten caps enough to compensate for
   the extra Hamiltonian and complex arithmetic.
3. Raise temperature before extracting purities. Scan base beta, multiplier,
   moment order, and trial length. A longer trial improves spatial isolation but
   makes an imperfect upper cap more costly in the temporal estimate.
4. Change the energy shift only within a justified trial bound. Shifts rescale
   every transfer eigenvalue and moment; their benefit comes from slack in the
   variational inequality, not from changing a purity. Retain error rescaling.
5. Compare capped versus uncapped mass estimates, lower-degree filters, and
   relaxed absolute table errors. A cheap, slightly weaker seed is preferable
   to maximizing the gap if it removes the n=12 computations.

Chosen workflow: reuse exact stored thermal centers first, generate new small
ring spectra with Float64 magnetization-block ED/Lanczos only where necessary,
search scalar bounds, then replay 2–3 useful designs with exact rationals.
Report failures and conditional assumptions explicitly, along with work proxies
and the analytic lemmas each winning design needs. Do not infer a runtime for
an unimplemented Lean spectral checker from Julia timings.

## spectrum: how many low states, how precise, how expensive?

1. Dense diagonalization in each Sz sector through n=8 gives an independent
   multiplicity check and complete thermal sums. Cost is cubic in sector size.
2. Sparse restarted Lanczos/Arnoldi through n=12 gives low spectra at cost roughly
   nnz(H) times Krylov iterations, with explicit residual checks. Repeated
   eigenvalues and missing states are risks: validate against dense spectra and
   request guard eigenvalues above every counting threshold. No Float64 count is
   a certified spectral anchor.
3. Chebyshev block filtering from individual computational-basis states tests the
   proposed small recipe. Oversampling and QR are necessary to avoid losing low
   states. Measure total vector Hv calls, not just polynomial degree. Compare
   floor-rounded integer recurrences with Float64, separating vector generation
   from numerical QR and from the eventual exact Gram checker.
4. A sharper tail bound using sector thresholds may save work; retain the proposed
   dimension-times-exponential bound as the baseline. Quantify the difference
   between absolute shifted-trace errors and relative errors, particularly odd n.
5. Packed ternary basis has cheap fixed shifts but 3^n lanes; compressed Sz basis
   has fewer lanes but irregular permutations. Do not assume both advantages at
   once. Report both layouts and treat kernel runtime as conditional on an actual
   packing implementation.

Chosen route: dense/sparse numerical spectra, followed by measured Chebyshev
filtering and exact integer recurrence probes. This is an experiment and cost
audit, not a certificate. All new files remain in this task directory.

## anchor: can a cheap comparison anchor the eigenvalue count?

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

## spin1_kernel: does packed arithmetic fit the kernel?

Task 4 tests the proposed Lean seed certificate; it does not claim a completed gap proof.
All new files remain in this directory; the original Lean project is read-only.

1. Full base-3 packing: digit-pattern masks have closed geometric-series recipes and
   nearest-neighbor transitions become fixed shifts. Cost is O(n) bignum operations per
   matvec on 3^n lanes. Chosen as the first executable prototype.
2. Compact magnetization packing: reduces n=12 from 531441 to 73789 lanes, but rank
   compression destroys constant shift distances. The advertised O(50) shifts needs an
   additional construction; audit this rather than assuming it.
3. Nonnegative shifted/Marshall-transformed operator: removes signed intermediate
   subtraction and permits cheap power-filter generation from product states. Signed
   arithmetic still needs explicit lane bounds; use separate positive/negative limbs.
4. Ritz plus residual: prove exact integer Gram inequalities first. Global low-eigenvalue
   identification additionally requires a counted anchor or a complementary-space bound.
   A residual alone is not a lower bound on the lowest eigenvalues.
5. Full spectral certification: exact LDL/inertia at rational endpoints is a fallback,
   but sparse elimination fill and kernel cost may defeat the proposed schedule.
6. Minimal honest baseline: kernel-check H6/H8 arithmetic and crude global bounds,
   recording explicitly any missing bridge or strong low-level lower certificate.

Scaling variables: vector length 3^n, lane width set by filter growth and dot-product
carry bounds, k retained vectors, filter degree, and recursive anchor count. Beta affects
required spectral precision and energy cutoff, not a single matvec directly.

Experiments using Float64 are non-certifying. Julia is only a generator; no generated
number is accepted as a theorem without a Lean kernel proof. No native_decide or sorry.

## ising: the cheapest model end to end

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

### seed_design

Route 2 (I/O/N split with the C twist) won; route 1 (real twists only) works only with the
length-12 tables. Three designs replay exactly (README); the recommended one uses nine traces
at n = 6, 8, 10 with absolute error 3/10000 and gives a conditional gap > 0.005966152 for every
even L ≥ 60. Length 8 fails in the tested family even at zero error (best maximum defect about
0.529), and length 10 without the C twist fails as well (about 0.296). Route 4 (moving the
shift) only relabels the error budget. Details and the Pareto front are in METHOD.md.

### spectrum

No certificate (none was attempted). The low-state counts are small (n = 12, β = 12.25 keeps 26
states for I and 34 for P; source 2 at n = 10 keeps 22/30/26 for I/P/C), and Q48 integer
recurrences with degree-16 cycles meet the numerical residual targets in every active sector.
Two premises of the proposed Lean seed route fail: direct compressed-Sz masked shifts need
thousands of displacement masks (3,816 at n = 12, Sz = 0; 13,620 at Sz = 2) instead of about 50,
and the elementary two-half-chain anchor leaves 242 Sz = 0 levels below the n = 12/I anchor
where the ring has 10. A kernel-checkable spectral anchor is the open problem (README, "Result").

### anchor

- Route 2 works locally: the four-site remainder is exactly PSD for the selected boundary
  parameters, and the half-chain spectra (Sz blocks ≤ 141) are cheap to enclose exactly.
  Optimising the comparison for the ground energy is the wrong objective; optimising the
  count at the tail cutoff cuts the Sz=0 count at n=12 from 191 (bare) / 212 (ground-opt) to 121.
- A binned tail using the comparison density of states lowers that to 49 (n=12, Sz=0).
- A one-shot count does not by itself give ordered lower bounds (missing-level obstruction);
  route 6 (continuation) is required. Its Float64 version completed in all 26 nonzero
  sectors for n=8,10,12 (2–73 stages per sector), as a design test only.
- Routes 3 (more blocks: larger counts), 4 (overlapping Anderson clusters: scalar only) and
  5 (edge spin-1/2: no lower comparison) were rejected in their cheap forms.

### spin1_kernel

Route 1 (full base-3 packing) works as an arithmetic primitive: dense n = 12 matvec + dot
checks by `decide +kernel` in tens of seconds, and power-filter Ritz/residual witnesses for
H6, H8, H10 are kernel-checked. Route 2 fails as stated: compact Sz = 0 packing needs 3804
distinct shifts at n = 12, not about 50. Route 4 is blocked by the anchor: with the untuned
cut comparison (two open half chains minus 4I), a Float64 scan of the full H6 and H8 spectra
finds no rank k with the separation the counting lemma needs, even for exact Ritz vectors.
Routes 5 and 6 were not completed. No seed, eigenvalue enclosure or gap is certified; the
long review and the implementation budget (8.5–15k Lean lines) are in METHOD.md.

### ising

Route 2 (low-energy Ritz plus a classical anchor) at the longitudinal Ising point
`g = 1/4, h = 1/20, (n₀, β₀) = (6, 97/4)`, with only two levels per ring, plus route 5 (the
partition-enclosure interface). The integer Gram/residual/anchor arithmetic is kernel-checked
at 6 and 12 sites, and the exact seed conclusions for the first grid point are proved
conditional on four `LogEnclosures` hypotheses. No seed is unconditional: the matrix
identification, min–max counting and scalar `exp` enclosures remain (METHOD.md). The
43-point inventory shows that full k×k Gram checks reach hours to ~1,000 hours for the
65,536- and 531,441-dimensional seeds, so routes 1–2 do not scale there without clustering or
symmetry reductions.

### Together

The plan survives in outline (Julia generates, Lean checks, packed `Nat` lanes carry the
arithmetic), but no seed is closed. Target design C of `seed_design` (real twists, n ≤ 12); the
bottleneck is the counting anchor, for which the local PSD comparison of `anchor` is exact and
cheap but leaves too many levels for a one-shot count, so a certified continuation (or a better
comparison) is the next research question. The missing Lean linear algebra (min–max, Weyl,
interlacing) and the packing/Hamiltonian bridges are the formalization cost (README, "Open
obligations").
