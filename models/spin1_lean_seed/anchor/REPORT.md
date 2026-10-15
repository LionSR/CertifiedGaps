# Task 2: cheaply checkable spectral anchors

The local comparison works, but the proposed one-shot Ritz/counting argument does not
close merely by retaining the number of comparison levels below a cutoff. The best
tested two-half comparison requires **121 vectors in Sz=0 at length 12**, versus **8
physical ring levels** below the original absolute-error cutoff. A cheaper tail bound
using the comparison density of states reduces the initial count to **49**, with half
the trace-error budget reserved for the low spectrum. Direct certification of the
four-to-six-site halves is cheap; recursion is unnecessary at these sizes.

The deliverable includes exact Julia arithmetic checks for the local comparison, half-chain
enclosures, and counts; numerical continuation experiments; and reproducible Julia scripts.
**There is no Lean kernel certificate, certified ring spectrum, or completed seed proof in
this directory.** All ring spectra, optimizations, and continuation stages are
NON-CERTIFYING Float64 experiments. Exact Julia checks remain untrusted relative to the
requested Lean end goal.

## The comparison and its exact check

Use the integer spin-one bond `h = Sz⊗Sz + (S+⊗S- + S-⊗S+)/2`. Its eigenvalues are
`-2,-1,1`. The pi twist negates only exchange entries. In a half of length `l=n/2`, set

```
B = u h + v h²
A_l = H_open,l + B_12 + B_(l-1,l)
H0 = A_l ⊗ I + I ⊗ A_l + 2 c I
R = h_23 - B_12 - B_34 - c I.
```

If `R >= 0`, the two embedded copies of R at the cuts sum exactly to `H-H0`.
They act on four sites each. For n=8,10,12 their supports are disjoint, although disjointness
is not needed for positivity. For the twisted cut, conjugate R by the diagonal pi rotation
on sites 3 and 4: the crossing bond twists while B_34 is invariant. The script checks this
integer identity. It does not assume that the twisted ring is SU(2)-invariant.

Every two-site SU(2)-invariant self-adjoint operator is a linear combination of `I,h,h²`;
the scalar part is absorbed into c. Thus this is the full **symmetric two-site boundary
family**, but not the set of all conceivable local comparisons. Four-start direct
minimization optimizes the ground bound. A separate 438-point grid/refinement per length
minimizes `(Sz=0 count, total count)` at the original tail cutoff. This finite search does
not certify global count optimality. It does not explore independent asymmetric boundaries.

The selected count parameters are particularly simple:

| n | u | v | c (exact) |
|---:|---:|---:|---:|
| 8 | -3/10 | -1/10 | -2.0307364 |
| 10 | -1/5 | 0 | -2.2154884 |
| 12 | -1/5 | 1/20 | -2.4266500 |

The four-site matrix decomposes by Sz into blocks of dimensions `19,16,10,4,1` for
nonnegative Sz. Negative sectors are related by spin flip. `certify.jl` and `check.jl`
perform exact Rational{BigInt} LDL checks on these blocks, including singular-pivot
handling, and check the twisted conjugation. The JSON stores all pivots as exact strings.
The ground-optimized parameters also pass these checks.

Optimizing the ground bound is the wrong objective for an anchor count. At n=12 it
improves the ground bound from -18.740550 (bare cuts) to -17.527216, but gives 914 states
below the original cutoff across all Sz, versus 773 for bare cuts. The count-optimized
comparison has a weaker ground bound, -17.690178, and only 475 states below that cutoff.
The corresponding Sz=0 counts are 212, 191, and 121 respectively (see generated tables
for the measured values; the total counts include negative Sz).

## What rho and k mean

Let `e_(m,i)` be the ordered spectrum of H0 in sector m, and let
`K_m(rho)=#{i:e_(m,i)<rho}`. Min-max gives

```
lambda_(K_m+1)(H_m) >= lambda_(K_m+1)(H0_m) >= rho.
```

This is an **upper bound on the number of H levels below rho**. It is not an assertion
that H actually has K_m such levels. The physical count is often much smaller. Both
counts, the anchor available at the physical count, the stronger anchor available at
the comparison count, and the first five ordered energy losses are recorded for every
`(n,g,Sz)` in [results/sectors.csv](results/sectors.csv). There are 66 nonnegative-sector
rows; negative Sz gives the identical row. `results/analysis.json` also contains losses
for bare-cut, ground-optimized, and count-optimized comparisons separately.

For the absolute shifted-trace budget epsilon=3e-5, the uniform-dimension cutoff is

```
rho_uniform = -n a + log(3^n/epsilon)/beta,
a=700741/500000, beta=49/4.
```

It equals approximately -9.644247, -12.267846, -14.891444 for n=8,10,12. These are not
exactly E0+2, and they use the *entire* error budget for the tail. In particular, n=8
untwisted has only two Sz=0 levels below the absolute-error cutoff, but five below E0+2.
The exact anchor certificate uses rational cutoffs rounded upward to six decimals and
proves the comparison counts agree on both sides of the half-spectrum enclosures.
The transcendental derivation of this baseline cutoff is only numerical here.

The failure in the proposed one-shot counting inference is elementary. Take
`H=diag(-2,-1,10)`, let V consist of the second and third coordinate vectors, and use
`rho=0 <= lambda_3(H)`. The two Ritz vectors have zero residual, yet the lowest Ritz
value -1 is not a lower bound on the true ground energy -2. A coarse count of two below
rho does not identify the missing eigenvalue. In this study, the upper Ritz levels in
a space of dimension K_m generally lie above the comparison anchor, so the simple
separation hypothesis needed to turn that count into ordered lower bounds fails.

This does not prove that every possible one-shot method fails. A certified lower bound
on the orthogonal complement, or a separate inertia/Birman-Schwinger calculation, could
replace continuation. Such a bound is additional work and does not follow from the
scalar counting anchor. For the count-only route actually proposed, keep homotopy in
the plan as a required component rather than an unspecified fallback.

## A cheaper tail using the comparison spectrum

Suppose continuation has established that all but r ring levels lie above a common
target rho. For i>r, min-max and this count give

```
lambda_i(H) >= max(rho,e_i(H0)).
sum_(i>r) exp[-beta(lambda_i(H)+na)]
  <= sum_i exp[-beta(max(rho,e_i(H0))+na)].
```

Bound the last sum with 32 bins of width 1/8 plus one final bin. Bin zero includes all
comparison energies below `rho+1/8`, charged at rho; later bins are charged at their
lower edges. The final bin charges all remaining states at `rho+4`. Multiplicities
are exact sector/tensor-product counts, not the full dimension at every energy.

The following targets reserve 1.5e-5 for the tail and leave 1.5e-5 for low-spectrum errors:

| n | rho | initial K by Sz=0,1,2,3,4 | total including negative Sz |
|---:|---:|---|---:|
| 8 | -9.925 | 28,22,11,2,0 | 98 |
| 10 | -12.687 | 38,32,15,4,0 | 140 |
| 12 | -15.453 | 49,40,20,6,1 | 183 |

These parameters were selected after the uniform-cutoff optimization; no global joint
optimization of boundary parameters and binned tail is claimed. The exact bin populations
use lower half-chain eigenvalue enclosures, so rounding can only weaken the tail bound.
Each scalar bound uses `exp(-x) <= 1/sum_(j=0)^256 x^j/j!` for positive rational x, rounded
up to denominator 10^24. Julia checks all rational arithmetic and the final budget.
The associated analytic Taylor inequality still needs a Lean lemma. No floating exponential
is accepted by the finite checker. The exact tail totals are just below 1.5e-5.

## Homotopy that was actually tested

Use `H(t)=H0+t(H-H0)`; local PSD gives monotonicity. At a stage retain K vectors and a
strict upper counting barrier `rho_bar < lambda_(K+1)(H(t))`. For trial matrix V form

```
M = V' (H(t)-rho_bar) V
N = ((H(t)-rho_bar)V)' ((H(t)-rho_bar)V).
```

Require M negative definite. This proves there are at least K levels below the barrier;
the propagated anchor proves there are at most K. Apply ordinary Ritz to
`(H-rho_bar)^(-1)` in the image space `(H-rho_bar)V`. If its K generalized Ritz values
`tau_i` are negative, `rho_bar+1/tau_i`, in reversed order, are lower bounds on the first
K eigenvalues. Levels whose lower bounds exceed the target can be removed, giving the
next counting barrier. Monotonicity propagates that barrier to the next t. If the count
reaches zero, monotonicity completes the rest of the path without another eigensolve.

The scripts implement this logic with Float64 generalized eigenproblems and adaptive
steps. Every accepted stage checks negativity; rejected steps are halved. A numerical
Hellmann-Feynman predictor chooses the next step but is not used as a proof. Dropping
levels immediately at the target can collapse the margin and stall at degeneracies;
the implemented rule removes them near the current barrier, preserving a buffer.
The final retained K can therefore exceed the number of levels strictly below the target.
Scalar Lanczos also stalled at some degenerate ring endpoints when asked for exactly K
eigenvectors. Oversampling small trial spaces resolved these stalls; the saved runs include
the stalled steps and subsequent resumes. A clean run with the revised solver can be cheaper.

These are genuine end-to-end **numerical tests of the continuation design**, not certified
continuations. Residual norms and acceptance margins are recorded. A Lean implementation
must replace every floating generalized eigenproblem and definiteness test with checked
integer/rational inequalities, account for Gram conditioning, and prove the counting
indices and strict barriers. In particular, the numerical safety padding is not a proof
of enclosure. The generated tables give stages, peak/final K, and timing per sector;
all K values at every stage are available in the corresponding JSON files.

The uniform-tail Sz=0 tests required 75 stages for n=8 and 73 for n=10 (untwisted).
The length-12 uniform pilot was stopped after two accepted stages: its initial
121-vector eigensolve alone took about 63 seconds. It is retained as an explicitly
incomplete pilot, not included among successful continuation tests.

## Other candidate routes

Three or more blocks were tested with bare cuts and the selected two-site boundary
parameters. The tested partitions include 4+2+2 and 3+3+2 at n=8; 4+3+3 and 4+4+2 at
n=10; and 4+4+4, 6+3+3, 6+4+2, and 3+3+3+3 at n=12. Additional cuts increased the counts.
The best sampled results are in the generated table and all samples are in
`results/analysis.json`. This is not a global optimization over multiblock parameters.

Anderson overlapping clusters can give a better scalar ground bound. For example the
five-site open-chain ground energy is approximately -5.830213, giving the scalar bound
`H_n >= (n/4) E_open,5 I`. Its largest Sz block is only 51. However, a scalar comparison
has every one of the 3^n levels below the thermal cutoff, so it supplies no useful low-rank
excited-level count. Keeping overlapping cluster operators instead of replacing them by
scalars produces noncommuting terms: their spectrum is not a tensor convolution. Weighted
overlaps do not remove that missing spectral calculation.

An edge-spin-1/2 compression similarly gives variational upper bounds by itself. A lower
comparison requires a certified discarded-space lower bound and coupling control. No
unproved effective edge model was admitted as an anchor. These observations reject these
routes in their cheap, directly computable forms; they do not exclude stronger constructions.

## Half-chain certification and the remaining Lean work

The maximum half-chain Sz dimensions are 19, 51, and 141 for lengths 4, 5, and 6. A
floating eigensolve proposes Q and diagonal D, then both are rounded to denominator 2^32.
The checker evaluates exactly

```
delta = ||Q Q' - I||_infinity
e = ||A - Q D Q'||_infinity
radius = e + 3 ||D|| delta.
```

Here Q denotes the scaled matrix. When delta<1/2, polar decomposition gives an orthogonal
U with `||Q-U|| <= delta` and `||Q||<2`. Thus
`||A-U D U'|| <= radius`, and Weyl encloses the entire ordered spectrum. Integer products
use Int128 only after a BigInt bound proves that every product and partial sum fits;
norm sums and final comparisons use BigInt rationals. This avoids trusting the floating
eigensolver for eigenvalue completeness or accuracy.

For the selected length-six half, the largest raw enclosure radius is about 2.78e-7.
Stored endpoints are rounded outward to a 1e-6 grid. The warm arithmetic check of the
141-dimensional block took roughly 0.04 seconds in an uncontended Julia run; timings in
later concurrent runs vary. All half-spectrum/tensor counts pass exact checks. Recursive
half-chain certification is therefore not needed for this experiment.

The current JSON stores enclosures and local pivots, but regenerates the integer Q proposals
using Julia/LAPACK. A standalone Lean checker must receive or regenerate these Q matrices
by a trusted deterministic recipe and then check them. For the three selected half-chain
families there are 53,346 Q entries: about 417 KiB if stored as Int64, before intervals and
metadata. Including the ground-optimized families doubles that. The current approximately
200 KiB JSON is not evidence that a few-KiB standalone Lean certificate has been achieved.

No Lean runtime has been measured here. Two dense products for all selected half sectors
require about 11.9 million scalar multiply-add terms. At the supplied 3e4 scalar-loop
iterations/second, a naive port already spends roughly seven minutes on this small part;
packed arithmetic remains necessary. For the continuation, the generated work table gives
`sum 3 K(K+1)/2`, the packed-dot count for a straightforward three-Gram-matrix implementation.
At n=12 a 64-bit-lane vector is about 4.72 million bits. Applying a 2.5-second *budget*
per product to those counts gives hours to days, not 15–60 minutes. This is a cost scenario,
not a lower bound or a measured Lean timing; the supplied “under 2.5 seconds” does not
establish that every product costs 2.5 seconds. Batched Gram products, reuse across stages,
or a stronger one-shot complement certificate could substantially change this budget.

The earlier 4.5–7.5k Lean LoC / 1.5–3 week estimate is not supported for the count-only
route. A planning allowance for this anchor/continuation component is 3.5–6k Lean LoC:
roughly 0.6–1k for sector/local-PSD lifting, 0.8–1.4k for exact spectral enclosures,
0.4–0.7k for tensor counts and bins, 1.2–2k for inverse-Ritz/counting, and 0.5–0.9k for
continuation bookkeeping. Depending on overlap with existing infrastructure, budget
approximately 8–13k total Lean LoC and 4–8 weeks for the full seed route. These are
engineering estimates, not code produced or a guarantee. Replace the 15–60 minute kernel
promise with an explicit benchmark milestone; even an optimistic target of 1–8 hours
requires substantial batching beyond the naive continuation implementation.

## Validation and scope

The six numerical shifted traces agree with the reproduction's tabulated centers to less
than 10^-6, providing a useful check of the bond, normalization, magnetization multiplicities,
and pi twist. Small-sector dense diagonalization, Hamiltonian identities, negative-Sz
symmetry, exact-checker rejection tests, and the missing-level counterexample are covered
by `test/runtests.jl`. Small residuals and trace agreement do not certify numerical counts.

The pinned paper source directory was absent in this checkout. The conventions and trace
centers were read from `models/spin1/src/thermal.jl`, and the pinned upstream identity from
`models/spin1/NOTICE`. No paper theorem has been silently treated as formalized. Nothing
under `lean/` or `models/spin1/results/` was changed.

Open risks are the Lean inverse-Ritz/counting development, strict barrier checks near
degeneracies, conditioning and generation of all ring trial spaces, kernel cost of repeated
Gram/residual calculations, and the lack of a global optimum certificate for boundary
parameters. The exact local positivity and small-half spectra are no longer the main
uncertainties. The useful plan change is to retain the optimized comparison, adopt the
binned tail, and make counting/continuation and its kernel benchmark explicit milestones.

<!-- GENERATED TABLES -->

## Numerical and exact-check results

Ring energies and losses in the next table are NON-CERTIFYING. K counts are for the uniform-cutoff comparison.

| n | g | E0(H) | ground-opt loss | count-opt loss | physical k (m=0) | bare K0 | ground-opt K0 | count-opt K0 |
|---|---|---|---|---|---|---|---|---|
| 8 | I | -11.336956 | 0.663044 | 1.227759 | 2 | 44 | 63 | 34 |
| 8 | P | -11.156192 | 0.843808 | 1.408523 | 3 | 44 | 63 | 34 |
| 10 | I | -14.094130 | 0.515383 | 0.743774 | 5 | 86 | 115 | 62 |
| 10 | P | -13.982142 | 0.627371 | 0.855762 | 6 | 86 | 115 | 62 |
| 12 | I | -16.869556 | 0.657659 | 0.820622 | 8 | 191 | 212 | 121 |
| 12 | P | -16.797640 | 0.729575 | 0.892538 | 8 | 191 | 212 | 121 |

Exact comparison counts and next-level lower anchors (both twists). Decimals here display the rational values in certificate.json. All unlisted sectors have K=0 at both cutoffs.

| n | Sz | uniform K | uniform next rho | binned K | binned next rho |
|---|---|---|---|---|---|
| 8 | 0 | 34 | -9.643366 | 28 | -9.799542 |
| 8 | 1 | 26 | -9.643366 | 22 | -9.799542 |
| 8 | 2 | 13 | -9.643366 | 11 | -9.799542 |
| 8 | 3 | 2 | -9.643365 | 2 | -9.643365 |
| 10 | 0 | 62 | -12.243159 | 38 | -12.657827 |
| 10 | 1 | 50 | -12.243159 | 32 | -12.657827 |
| 10 | 2 | 27 | -12.243159 | 15 | -12.657827 |
| 10 | 3 | 8 | -12.243159 | 4 | -12.473613 |
| 10 | 4 | 1 | -11.954917 | 0 | -12.473613 |
| 12 | 0 | 121 | -14.890385 | 49 | -15.440805 |
| 12 | 1 | 100 | -14.890385 | 40 | -15.440805 |
| 12 | 2 | 56 | -14.890385 | 20 | -15.373806 |
| 12 | 3 | 18 | -14.890385 | 6 | -15.373806 |
| 12 | 4 | 3 | -14.855058 | 1 | -14.961896 |

Exact binned-tail arithmetic totals; the analytic Taylor and spectral inequalities remain to be formalized in Lean.

| n | target rho | tail upper | budget |
|---|---|---|---|
| 8 | -9.925 | 1.48995173385e-05 | 0.000015 |
| 10 | -12.687 | 1.49126012899e-05 | 0.000015 |
| 12 | -15.453 | 1.49383528391e-05 | 0.000015 |

First five ordered losses in Sz=0 for the count-optimized comparison (NON-CERTIFYING). Full per-sector losses and both definitions of k are in [results/sectors.csv](results/sectors.csv).

| n | g | loss 1 | loss 2 | loss 3 | loss 4 | loss 5 |
|---|---|---|---|---|---|---|
| 8 | I | 1.227759 | 1.333551 | 2.480392 | 2.245592 | 2.245592 |
| 8 | P | 1.408523 | 1.306997 | 2.215279 | 2.235966 | 2.235966 |
| 10 | I | 0.743774 | 1.268582 | 2.240853 | 1.952080 | 1.952080 |
| 10 | P | 0.855762 | 1.252706 | 2.053984 | 1.937183 | 1.937183 |
| 12 | I | 0.820622 | 1.114069 | 1.970001 | 1.822924 | 1.822924 |
| 12 | P | 0.892538 | 1.104232 | 1.831977 | 1.811095 | 1.811095 |

Best sampled multiblock comparison at the uniform cutoff (NON-CERTIFYING selection).

| n | partition | parameters | Sz=0 K | total K |
|---|---|---|---|---|
| 8 | 3+3+2 | boundary | 104 | 434 |
| 10 | 4+3+3 | boundary | 207 | 901 |
| 12 | 4+4+4 | boundary | 407 | 1819 |

Binned-tail continuation, all initially nonzero sectors (**NON-CERTIFYING**). Final K is the retained count, which may exceed the actual number below the target. Timings include numerical eigensolves, compilation, and contention; they are not Lean timings. Negative Sz shares the same spectrum. Each linked file lists k and rho at every stage.

| n | g | Sz | stages | initial → final K | physical k | Julia seconds | min margin | data |
|---|---|---|---|---|---|---|---|---|
| 8 | I | 0 | 28 | 28 → 2 | 2 | 21.0 | 5.7e-05 | [stages](results/homotopy_n8_gI_m0_binned.json) |
| 8 | I | 1 | 24 | 22 → 1 | 1 | 30.0 | 5.3e-05 | [stages](results/homotopy_n8_gI_m1_binned.json) |
| 8 | I | 2 | 14 | 11 → 0 | 0 | 0.9 | 0.00015 | [stages](results/homotopy_n8_gI_m2_binned.json) |
| 8 | I | 3 | 3 | 2 → 0 | 0 | 0.0 | 0.012 | [stages](results/homotopy_n8_gI_m3_binned.json) |
| 8 | P | 0 | 25 | 28 → 3 | 2 | 48.4 | 0.00027 | [stages](results/homotopy_n8_gP_m0_binned.json) |
| 8 | P | 1 | 17 | 22 → 2 | 2 | 3.4 | 0.00061 | [stages](results/homotopy_n8_gP_m1_binned.json) |
| 8 | P | 2 | 13 | 11 → 0 | 0 | 1.7 | 6.2e-05 | [stages](results/homotopy_n8_gP_m2_binned.json) |
| 8 | P | 3 | 2 | 2 → 0 | 0 | 0.0 | 0.026 | [stages](results/homotopy_n8_gP_m3_binned.json) |
| 10 | I | 0 | 51 | 38 → 2 | 2 | 54.8 | 4.2e-05 | [stages](results/homotopy_n10_gI_m0_binned.json) |
| 10 | I | 1 | 48 | 32 → 1 | 1 | 23.3 | 4.2e-05 | [stages](results/homotopy_n10_gI_m1_binned.json) |
| 10 | I | 2 | 23 | 15 → 0 | 0 | 6.2 | 5.2e-05 | [stages](results/homotopy_n10_gI_m2_binned.json) |
| 10 | I | 3 | 5 | 4 → 0 | 0 | 0.3 | 0.00088 | [stages](results/homotopy_n10_gI_m3_binned.json) |
| 10 | P | 0 | 42 | 38 → 3 | 3 | 68.1 | 3.9e-05 | [stages](results/homotopy_n10_gP_m0_binned.json) |
| 10 | P | 1 | 29 | 32 → 2 | 2 | 41.2 | 5.3e-05 | [stages](results/homotopy_n10_gP_m1_binned.json) |
| 10 | P | 2 | 21 | 15 → 0 | 0 | 7.3 | 5.6e-05 | [stages](results/homotopy_n10_gP_m2_binned.json) |
| 10 | P | 3 | 4 | 4 → 0 | 0 | 0.6 | 0.0029 | [stages](results/homotopy_n10_gP_m3_binned.json) |
| 12 | I | 0 | 73 | 49 → 5 | 5 | 666.2 | 1.2e-05 | [stages](results/homotopy_n12_gI_m0_binned.json) |
| 12 | I | 1 | 63 | 40 → 4 | 4 | 361.3 | 1.6e-05 | [stages](results/homotopy_n12_gI_m1_binned.json) |
| 12 | I | 2 | 25 | 20 → 1 | 1 | 56.6 | 0.00015 | [stages](results/homotopy_n12_gI_m2_binned.json) |
| 12 | I | 3 | 9 | 6 → 0 | 0 | 16.8 | 0.00043 | [stages](results/homotopy_n12_gI_m3_binned.json) |
| 12 | I | 4 | 2 | 1 → 0 | 0 | 1.6 | 0.016 | [stages](results/homotopy_n12_gI_m4_binned.json) |
| 12 | P | 0 | 72 | 49 → 5 | 5 | 675.5 | 1.1e-05 | [stages](results/homotopy_n12_gP_m0_binned.json) |
| 12 | P | 1 | 33 | 40 → 4 | 4 | 222.1 | 3.1e-05 | [stages](results/homotopy_n12_gP_m1_binned.json) |
| 12 | P | 2 | 24 | 20 → 0 | 0 | 46.7 | 0.00016 | [stages](results/homotopy_n12_gP_m2_binned.json) |
| 12 | P | 3 | 6 | 6 → 0 | 0 | 16.2 | 0.00024 | [stages](results/homotopy_n12_gP_m3_binned.json) |
| 12 | P | 4 | 2 | 1 → 0 | 0 | 0.5 | 0.012 | [stages](results/homotopy_n12_gP_m4_binned.json) |

Naive three-Gram packed-dot work for completed binned runs. Only nonnegative Sz is computed; negative Sz is reused. The final column applies a hypothetical 2.5-second budget per dot, not a measured cost.

| n | g | sector-stages | sum K | pair products | budget hours |
|---|---|---|---|---|---|
| 8 | I | 69 | 747 | 19560 | 13.6 |
| 8 | P | 57 | 659 | 17940 | 12.5 |
| 10 | I | 127 | 1858 | 63882 | 44.4 |
| 10 | P | 96 | 1457 | 50649 | 35.2 |
| 12 | I | 172 | 3200 | 140313 | 97.4 |
| 12 | P | 137 | 2616 | 115320 | 80.1 |
