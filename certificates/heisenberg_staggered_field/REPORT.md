# Uniform gap certificates: Heisenberg chain in a staggered field

For current-main compatibility and separately recorded validation of the
CertifiedGaps package rename, see MIGRATION.md. Historical results and hashes
remain unchanged; migration outputs have a separate migration_ prefix.

## Outcome and scope

All five requested fields have rigorous certificates from 8- and 16-site rings.
The final bounds apply to **every even L >= 8 one-spin cells = sites**. The physical
Hamiltonian has a two-site field pattern, but its rotated representative has a
one-spin unit cell. The theorem only uses even n0; odd n0 are not valid seeds here.
No threshold occurs within the requested field set: the smallest certified field
is 1/4. This is a certificate/search threshold statement, not a physical critical field.

A separate numerical search reaches 20 sites using Sz and translation-by-two
momentum sectors. It finds stronger candidates but uses partial Krylov spectra
for blocks above dimension 2500; these candidates are **numerical, not certified**.
We certify the best seed on the complete dense grid through 16 sites, rather than
claim certification of the best candidate on the larger, partial-spectrum grid.
The requested fully dense certification of the best 20-site candidates remains
unperformed. No 24-site computation is claimed.
The largest block grows from 1620 to 18504: a single complex dense matrix then
occupies about 5.48 GB, and the current verifier retains several such arrays.
The naive cubic eigensolve cost grows by about 1490 times per largest block.
These are resource estimates, not measured failed dense twenty-site runs.

## Model, exact RP identity, and length convention

H = sum_j S_j dot S_(j+1) - h sum_j (-1)^j Sz_j, with periodic bonds,
spin 1/2, exchange strength 1, and h in {1,3/4,1/2,7/20,1/4}.
The two-site ring includes both periodic bonds, as does the core implementation.
Define exact rational matrices

```
X = [0  1/2;  1/2  0],
J = [0  1/2; -1/2  0] = i Sy,
Z = [1/2 0; 0 -1/2],   U = 2J.
```

Rotate odd sites by the real orthogonal U (a pi rotation about y). Both orientations
of an adjacent rotated/unrotated pair give

```
h_rot = -X tensor X - J tensor J - Z tensor Z,
-h_rot = X tensor X + J tensor J + Z tensor Z.
```

Every coefficient c is exactly 1 and every C is real. J is antisymmetric; the
hypothesis requires real C, not symmetric C. U Z U^T = -Z, so the staggered field
becomes the uniform real symmetric on-site term -h Z. `exact_rp` checks the bond
identities in Rational arithmetic in both orientations, and checks the on-site
term. Thus the stated spatial transfer self-adjointness hypothesis holds. The
positive two-purity bound then gives uniqueness as well as the gap, using the
repository criterion; no numerical assumption of uniqueness is substituted.

## Dense-grid certificates

Grid: n0 in {2,4,6,8}, beta in {0.25,0.50,...,32}. `scan_seeds` and `best_seed`
select the largest bound. Only that seed is rerun with dense eigenvectors and
outward error bounds for each requested field. All exact rational parameters
and seed defects are strings in the certificate JSON.

| h | n0 | beta | rigorous bound | numerical 16-site gap | bound/gap | certification seconds |
|---|---:|---:|---:|---:|---:|---:|
| 1//1 | 8 | 21//4 | 1.148261 | 1.87351193 | 0.612892543 | 79.2549119 |
| 3//4 | 8 | 21//4 | 0.890782 | 1.57089206 | 0.567054967 | 80.252265 |
| 1//2 | 8 | 11//2 | 0.594785 | 1.23131976 | 0.483046933 | 50.400383 |
| 7//20 | 8 | 23//4 | 0.383508 | 0.994380622 | 0.385675562 | 47.0755551 |
| 1//4 | 8 | 25//4 | 0.214486 | 0.812031433 | 0.264135226 | 44.339813 |

Displayed rigorous bounds are rounded DOWN to six decimal places from the smaller original and freshly replayed bound. Original raw outputs remain in certificate JSONs; fresh raw outputs and exact seed differences are in `results/spectral_replay.json`. The full spectral replay found small numerical differences, so the earlier nine-decimal display is superseded; see AUDIT.md.

For h=1, the exact stored upper seed defects are epsilon_s = `3124087506499//4503599627370496` and epsilon_t = `6599110342439//9007199254740992`, with beta = `21//4`.

| h=1 ring | max entry error | max residual upper | max orthogonality upper | max eigenvalue radius |
|---:|---:|---:|---:|---:|
| 8 | 1.5034199e-15 | 7.08836598e-12 | 1.83181771e-12 | 2.2315322e-11 |
| 16 | 5.27181194e-14 | 5.1456421e-08 | 9.56335631e-09 | 2.10443015e-07 |

The h=1 survey value supplied in the task is 1.147821, with reference gap 1.87351.
The current repository recurrence gives 1.1490006457432576 from unverified dense
seed values at (8,5.25). Our verified bound includes enclosure loss and is the
value above. We do not silently substitute the survey number. Reference gaps and
ratios in this report are numerical ED values, not independent uniform gap proofs.

## Complete sectors and multiplicity

In the physical frame, p counts down spins, m=L/2, and k indexes momentum under
translation by two sites. An orbit of period ell contributes iff k*ell is divisible
by m. Its normalized state is sum_r exp(-2 pi i k r/m) |T^(2r) a>/sqrt(ell).
For a transition from representative a into T^(2t)b, the matrix element is the
local amplitude times sqrt(ell_a/ell_b) exp(2 pi i k t/m). This formula is used
directly; a small-ring test compares the full merged spectrum with the independent
core Sz implementation.

Global spin flip composed with one-site translation maps p to L-p and preserves
H and momentum. Complex conjugation maps k to m-k. We diagonalize one dense
representative per such exact equivalence, with multiplicities 1,2, or 4; no
uncomputed inequivalent sector is omitted. Exact integer orbit counts and the
assertion sum(dimension*multiplicity)==2^L run for every ring. No Krylov solve is
called on the certificate path.

| ring | representative blocks | largest block | weighted dimension | required dimension |
|---:|---:|---:|---:|---:|
| 8 | 13 | 20 | 256 | 256 |
| 16 | 41 | 1620 | 65536 | 65536 |

Full sector dimensions, multiplicities, eigenvalues, radii and timings are stored per block in each certificate JSON. For L=16:

| p | dimensions for k=0,1,2,3,4 | charge multiplicity |
|---:|---|---:|
| 0 | 0:1 | 2 |
| 1 | 0:2, 1:2, 2:2, 3:2, 4:2 | 2 |
| 2 | 0:16, 1:14, 2:16, 3:14, 4:16 | 2 |
| 3 | 0:70, 1:70, 2:70, 3:70, 4:70 | 2 |
| 4 | 0:232, 1:224, 2:230, 3:224, 4:232 | 2 |
| 5 | 0:546, 1:546, 2:546, 3:546, 4:546 | 2 |
| 6 | 0:1008, 1:994, 2:1008, 3:994, 4:1008 | 2 |
| 7 | 0:1430, 1:1430, 2:1430, 3:1430, 4:1430 | 2 |
| 8 | 0:1620, 1:1600, 2:1616, 3:1600, 4:1618 | 1 |

Momentum multiplicities for k=0,1,2,3,4 are 1,2,2,2,1. Absent entries have dimension zero.

## Eigenvalue-enclosure theorem and implementation

All stored Float64 entries are interpreted as exact binary numbers. Let H be the
exact Hermitian momentum block and A the Hermitian matrix submitted to LAPACK.
`block(...; verified=true)` evaluates rational coefficients, orbit normalization
square roots, and trigonometric phases in outward-rounded IntervalArithmetic.
Pi itself is enclosed, not replaced with a single floating-point approximation.
Each exact upper-triangular entry is enclosed; the lower triangle is its conjugate.
This gives a Frobenius upper bound a >= ||H-A||_2, stored as `entry_error`.

Let Q and Lambda be the full computed eigenvector/eigenvalue arrays. The program
explicitly calculates Rhat=fl(AQ-Q Lambda) and Ehat=fl(Q* Q-I). It never assumes
that LAPACK's answer is exact. Each Frobenius norm is bounded using upward scalar
multiplication, addition, and square root (`nextfloat` after each operation,
including a subnormal allowance). Complex real and imaginary components are
accounted for separately. Write q,h,l for these upper norms of Q,A,Lambda.
With u=2^-53 and g=gamma_(32d+32)=(32d+32)u/(1-(32d+32)u), the implementation uses

```
r = ||Rhat||_F,up + g (h q + q l) + 10^-280,
e = ||Ehat||_F,up + g (q^2 + sqrt(d)) + 10^-280.
```

These follow from the standard dot-product rounding bound
|fl(sum products)-sum products| <= gamma_s sum |products|, applied to real and
imaginary components and their subtractions, and || |A| |B| ||_F <= ||A||_F||B||_F.
The factor 32d+32 conservatively exceeds the number of elementary operations per
complex dot product and subtraction. Sparse AQ has no longer dot products than
the dense dimension. The sqrt(d) term explicitly budgets subtraction of I.
The added 10^-280 exceeds the aggregate possible subnormal
absolute errors at all reached dimensions. Inputs/results are finite and the
chosen gamma denominator is positive. No fast-math transformations are requested.
This is a rigorous norm error bound on ordinary matrix multiplication, rather
than merely trusting a residual accumulated in the same precision.

**Polar/Weyl argument (with multiplicity):** assert e<1. The polar unitary V of Q
exists and ||Q-V||_2 <= dQ=e/(1+sqrt(1-e)). Therefore

```
||A - V Lambda V*||_2
 = ||AV - V Lambda||_2
 <= r + (||A||_2 + ||Lambda||_2) dQ.
```

Use ||A||_2 <= L(3/4+h/2)+a and ||Lambda||_2=max|lambda_i|. Weyl's inequality
for the two Hermitian matrices H and V Lambda V* encloses each ordered eigenvalue,
with multiplicity, in lambda_i +/- eta, where

```
eta = a + r + (L(3/4+h/2)+a+max|lambda_i|) dQ.
```

All scalar operations defining r,e,dQ,eta are outward interval operations. This
argument does not rely on simple eigenvalues or on matching individual residuals
to distinct eigenvalues. Every eigenvalue of each block is present.

For beta>0, `enclosed_logz` sums exp(-beta*[lambda_i-eta,lambda_i+eta]) with exact
symmetry multiplicities and a common exact binary shift, then takes an interval
log. Both lower and upper logZ endpoints are retained. The lower purity bounds
come from logZ(2n,beta)-2logZ(n,beta) and logZ(2n,2beta)-2logZ(2n,beta).
Upper endpoints of 1-exp(these intervals) are converted to exact Rational{BigInt}
values, with no intermediate decimal rounding. Those rational seeds and rational
beta are passed to the repository's `certified_gap_bound`, including its rigorous
tail correction. The fast test recomputes this propagation from stored enclosures.

## Twenty-site exploration and threshold

The following results are **numerical, not certified**. Full dense spectra are used
through 16 sites. At 20 sites only blocks of dimension <=2500 are dense; larger
blocks retain 80 low levels from KrylovKit. All sectors are visited, with the same
exact dimension count. A common truncation cutoff and a dimension-weighted tail
are included via core `PartialSpectrum`/`logZ_bounds`, but Krylov convergence and
multiplicities are not verified, so these are not rigorous partition bounds.
Exact twenty-site dimensions are recorded below.

| h | best n0 | beta | numerical candidate bound | numerical 20-site reference gap |
|---|---:|---:|---:|---:|
| 1//1 | 10 | 6.25 | 1.27832796 | 1.87351191 |
| 3//4 | 10 | 6.5 | 1.01075578 | 1.57089173 |
| 1//2 | 10 | 6.75 | 0.707243353 | 1.23131113 |
| 7//20 | 10 | 7 | 0.493229271 | 0.994292667 |
| 1//4 | 10 | 7.25 | 0.326952394 | 0.811528896 |
| 1//5 | 10 | 7.5 | 0.231456632 | 0.708641218 |
| 3//20 | 10 | 8 | 0.119366516 | 0.594804607 |
| 1//10 | — | — | no grid seed | 0.466432604 |

Actual L=20 largest block: 18504; weighted total: 1048576 = 2^20.

Additional fields 1/5,3/20,1/10 probe the numerical search threshold beyond the requested set. A failed grid point is not evidence of a gapless phase.

The sampled twenty-site search fails at h=0.1 and succeeds at h=0.15. This brackets an apparent finite-grid threshold between those sampled fields; no monotonicity between them or physical critical field is proved.

For every explored field, the following gives max_beta min(Ps,Pt) and delta(n)=E0(2n)-2E0(n). Partial-spectrum values remain numerical.

| h | n0 | max min purity | delta(n) |
|---|---:|---:|---:|
| 1//1 | 2 | 0.529873172 | 0.505022849 |
| 1//1 | 4 | 0.94067001 | 0.0221827722 |
| 1//1 | 6 | 0.99334286 | 0.00159103285 |
| 1//1 | 8 | 0.999275644 | 0.00013128857 |
| 1//1 | 10 | 0.999923253 | 1.14962961e-05 |
| 3//4 | 2 | 0.447640968 | 0.613999064 |
| 3//4 | 4 | 0.890324463 | 0.0389149028 |
| 3//4 | 6 | 0.983514815 | 0.00385786603 |
| 3//4 | 8 | 0.997555655 | 0.000440747681 |
| 3//4 | 10 | 0.999636029 | 5.36822693e-05 |
| 1//2 | 2 | 0.399066156 | 0.754873673 |
| 1//2 | 4 | 0.800725359 | 0.0769369751 |
| 1//2 | 6 | 0.952122862 | 0.0113905614 |
| 1//2 | 8 | 0.988846519 | 0.00192196739 |
| 1//2 | 10 | 0.997543855 | 0.000347058281 |
| 7//20 | 2 | 0.361338696 | 0.850802846 |
| 7//20 | 4 | 0.699446862 | 0.124633121 |
| 7//20 | 6 | 0.893328792 | 0.0251937055 |
| 7//20 | 8 | 0.967347486 | 0.00566434475 |
| 7//20 | 10 | 0.990324488 | 0.00135842657 |
| 1//4 | 2 | 0.31998312 | 0.913017541 |
| 1//4 | 4 | 0.602523494 | 0.177017263 |
| 1//4 | 6 | 0.818031626 | 0.0466268272 |
| 1//4 | 8 | 0.926399235 | 0.0132591633 |
| 1//4 | 10 | 0.971872848 | 0.00397943173 |
| 1//5 | 2 | 0.303223146 | 0.940972384 |
| 1//5 | 4 | 0.560054301 | 0.2117524 |
| 1//5 | 6 | 0.757748899 | 0.0653997426 |
| 1//5 | 8 | 0.881241013 | 0.0214187006 |
| 1//5 | 10 | 0.947706946 | 0.00732718927 |
| 3//20 | 2 | 0.289667018 | 0.965093873 |
| 3//20 | 4 | 0.495941276 | 0.252166166 |
| 3//20 | 6 | 0.674748062 | 0.0934440211 |
| 3//20 | 8 | 0.814781105 | 0.0361784563 |
| 3//20 | 10 | 0.902221022 | 0.0144183736 |
| 1//10 | 2 | 0.279689865 | 0.983874689 |
| 1//10 | 4 | 0.442024104 | 0.295440509 |
| 1//10 | 6 | 0.569803257 | 0.134401906 |
| 1//10 | 8 | 0.7028359 | 0.0639327614 |
| 1//10 | 10 | 0.80680422 | 0.0308278889 |

## Predicted scale, extrapolation, and limits

Using the supplied Delta proportional to h^(2/3) scaling and calibrating the
**quality of the h=1 survey seed**, the naive comparable-quality scales are
n0 approximately 8 h^(-2/3), beta approximately 5.25 h^(-2/3). Logarithmic corrections
and changes in amplitudes are omitted. This predicts a scale, not a necessary
condition for any positive bound. Weaker positive certificates can occur earlier.

| h | predicted n0 (continuous) | next even n0 | predicted beta | largest searched n0 |
|---|---:|---:|---:|---:|
| 1//1 | 8 | 8 | 5.25 | 10 |
| 3//4 | 9.69130983 | 10 | 6.35992207 | 10 |
| 1//2 | 12.6992084 | 14 | 8.33385552 | 10 |
| 7//20 | 16.1081114 | 18 | 10.5709481 | 10 |
| 1//4 | 20.1587368 | 22 | 13.229171 | 10 |

For h=1//10, fitting the last three delta(n) values gives log delta approximately 0.199436518 + (-0.368103712) n, or correlation scale 2.71662569.
A heuristic model Ps=exp[-beta delta(n)] and Pt=(1+4n exp[-2beta gap])/(1+4n exp[-beta gap])^2 first yields a predicted seed at (n0,beta)=(14,17.0), requiring 28 sites. This is an extrapolation, not an existence proof or certificate.

## Reproduction, environment, and evidence status

Julia 1.12.4 was verified. The approved juliaup installer returned HTTP 403; the
official julialang-s3 binary archive succeeded. Dependencies were instantiated and
precompiled from the unchanged julia project. OPENBLAS_NUM_THREADS=1 and
BLAS.set_num_threads(1) were used. The environment has approximately 17 GiB RAM
and five CPUs. Independent parameter runs used separate processes. No individual
run approached the 45-minute limit. Source SHA-256 hashes, package versions,
per-block times, and Julia/BLAS thread counts are recorded in JSON.

Run commands are in README.md. `test/runtests.jl` reruns the exact RP identities,
compares small complete momentum spectra with the independent core ED, checks
sector totals, recomputes enclosed partition functions, reconstructs rational
seeds, and reruns the certified recurrence. It does not reconstruct the large
Q matrices; rerun `certify.jl` to audit residual calculations from scratch.

**Certified:** exact RP at the five requested rational fields; complete weighted
sector coverage on certificate rings; eigenvalue, partition-function and purity
enclosures under IEEE arithmetic and the stated error analysis; the reported
uniform lower bounds for every even L>=8 sites, with uniqueness from the criterion.

**Numerical only:** grid optimization, all reference gaps and ratios, twenty-site
partial spectra and seed candidates, finite-grid threshold observations,
finite-size delta fits, and h^(2/3) scale predictions.

**Remaining limit:** the stronger best candidates on the twenty-site search grid
have not been rerun through complete dense verification. No physical threshold
or impossibility of certificates at weaker fields is proved. No existing core,
other certificate directory, or notes were edited; nothing was merged or published.

Actual checker result: **71 tests passed**; 22.408746091 seconds; all tests passed. See `results/tests.log`. Source-hash audit and full dependency versions are in `results/environment.json`.

The separate independent certification review and its acceptance status are recorded in [AUDIT.md](AUDIT.md).
