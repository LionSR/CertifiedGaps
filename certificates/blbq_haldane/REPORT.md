# Spin-1 bilinear–biquadratic Haldane survey

## Outcome and scope

**No eligible seed was found on the computed grid. No uniform spectral gap is certified.**
Complete dense spectra reach **12 sites**, with numerical references through **16 sites**.
The constructive SU(2) reduction passes small-spectrum tests and 14-site pilots; its
complete 14-site sweep is still being collected. No partial sector set is used as Z14.


One unit cell is one spin-1 site. The tested Hamiltonian is
`H_L = sum_j [cos(theta) X_(j,j+1) + sin(theta) X_(j,j+1)^2]`,
`X = S·S`, with periodic boundaries and angles
`theta = 0, -1/10, -3/10, -1/2, -3/5, -7/10` radians.
For L=2 both directed ring bonds are included, as in the core convention.
There is no on-site term. No energy shift or rescaling is applied.

The scope note's theorem requires **even n0**. Every reported odd n0 is only a
finite-ring diagnostic; odd-n spatial ratios need not be purities of a positive operator.
The hypothetical validity range would be **even L ≥ n0 cells = sites**. Since there is
no certified seed, **there is no asserted validity range or positive bound-to-gap ratio**.

## Exact reflection positivity

Use Cartesian generators `(G_a)_(bc)=epsilon_(abc)`; explicitly

```
G1 = [0 0 0; 0 0 1; 0 -1 0]
G2 = [0 0 -1; 0 0 0; 1 0 0]
G3 = [0 1 0; -1 0 0; 0 0 0]
```

Here `S_a=-i G_a`. Define the **symmetrized**, real symmetric quadrupoles
`Q_ab = -(G_a G_b + G_b G_a)/2` (not the unsymmetrized products).
With `A=sum_a kron(G_a,G_a)` and `B=sum_ab kron(Q_ab,Q_ab)`, the exact identities are

```
X = -A
X^2 = A/2 + B
-h = (cos(theta)-sin(theta)/2) A - sin(theta) B.
```

`model.jl:exact_rp` verifies the two coefficient identities using `Rational{BigInt}`;
this proves the identity for symbolic real c=cos(theta), s=sin(theta), without rounding
transcendental coefficients to rational numbers. For nonzero rational angles these
coefficients are not being claimed to be algebraic. The generator and quadrupole entries
are exact integers/half-integers. The on-site zero matrix is exactly real symmetric.
For -pi/4 < theta ≤ 0, sin(theta) ≤ 0 and cos(theta) > 0, so both weights are nonnegative.

The skew-symmetric and symmetric operator spaces are orthogonal. The three G's span the
former, and six independent Q's span the latter, as checked by a nonzero **exact rational
9-by-9 determinant** in the tests. Consequently the realigned bond is PSD iff the two
weights are nonnegative: a negative weight has a negative direction in its own orthogonal
space. This proves the iff condition for this Cartesian bond, rather than relying on a
floating-point RP test or the numerical converse in the scope note.
An independent exact realignment identity gives eigenvalues `2c-s` (multiplicity 3),
`-s` (multiplicity 5), and `-4s` (multiplicity 1), where c=cos(theta), s=sin(theta).

## Complete spectra and symmetry counting

All partition functions below use every eigenvalue from dense real symmetric LAPACK
solves. The reduction uses total Sz=M, cyclic translation k, spin reversal within M=0,
and reflection at k=0 and k=pi. Opposite M and conjugate k share spectra, with explicit
multiplicities 2 where appropriate. At general k, reflection followed by conjugation is
an antiunitary involution; its fixed basis makes the Hamiltonian real. It does not itself
halve the block. Sparse involution projectors supply the additional ± sectors.

Each result records every block dimension and multiplicity, and asserts both
`sum(dim * multiplicity) == 3^L` and that the expanded spectrum has exactly `3^L`
states. SU(2) outputs retain energies with exact integer multiplet/momentum weights.
The tests compare complete reduced and unreduced spectra at L=2,3,4,5 for three
angles (including -0.7). Matrix Hermiticity, real-basis consistency, symmetry invariance
and involution identities are also checked during construction. Complete spectra also
pass analytic first- and second-moment checks at infinite temperature: for L>=3,
`mean(H)=4 L sin(theta)/3` and
`var(H)=L [4 cos(theta)^2/3 - 4 cos(theta)sin(theta)/3 + 8 sin(theta)^2/9]`;
the L=2 variance has an extra factor two because both bonds coincide. These floating checks
validate the numerical implementation; they are not verified eigenvalue enclosures.

| L (sites) | Total dimension | Largest dense block | Stored blocks | Wall seconds, six angles |
|---:|---:|---:|---:|---:|
| 2 | 9 | 2 | 5 | 1.31596 |
| 3 | 27 | 2 | 9 | 0.00296807 |
| 4 | 81 | 5 | 20 | 0.00273705 |
| 5 | 243 | 9 | 22 | 0.00369382 |
| 6 | 729 | 21 | 40 | 0.00802088 |
| 7 | 2187 | 51 | 39 | 0.0163689 |
| 8 | 6561 | 127 | 61 | 0.148689 |
| 10 | 59049 | 835 | 86 | 4.07773 |
| 12 | 531441 | 5798 | 115 | 1067.58 |

## Purity scaling table

For each n the entry is the **grid maximum** of min(Ps(n,beta),Pt(2n,beta)), with its
maximizing beta. Grid: beta = 1/8, 2/8, ..., 64 (512 points), using core `scan_seeds`
and `best_seed`. Maxima are numerical, not upper bounds on the continuous-beta optimum.
All even grid points were checked with the actual asymmetric recurrence, not just the
symmetric 0.12731 sufficient threshold. No eligible positive recurrence bound appeared.

| theta | n0 | eligible even seed? | beta at max | max min(Ps,Pt) | Ps | Pt |
|---:|---:|:---:|---:|---:|---:|---:|
| 0 | 2 | yes | 1 | 0.179556 | 0.190126 | 0.179556 |
| 0 | 3 | no | 44.875 | 1 | 1.02534e+51 | 1 |
| 0 | 4 | yes | 2.25 | 0.261921 | 0.261921 | 0.285953 |
| 0 | 5 | no | 60.25 | 1 | 1.0005e+27 | 1 |
| 0 | 6 | yes | 3.25 | 0.328774 | 0.328774 | 0.344587 |
| -0.1 | 2 | yes | 0.875 | 0.174944 | 0.179517 | 0.174944 |
| -0.1 | 3 | no | 40.25 | 1 | 4.40885e+61 | 1 |
| -0.1 | 4 | yes | 1.875 | 0.243032 | 0.264827 | 0.243032 |
| -0.1 | 5 | no | 59 | 1 | 3.66317e+39 | 1 |
| -0.1 | 6 | yes | 2.875 | 0.291791 | 0.291791 | 0.302041 |
| -0.3 | 2 | yes | 0.75 | 0.144661 | 0.144661 | 0.185162 |
| -0.3 | 3 | no | 32 | 1 | 5.29978e+58 | 1 |
| -0.3 | 4 | yes | 1.625 | 0.227725 | 0.227725 | 0.254607 |
| -0.3 | 5 | no | 50.75 | 1 | 2.21151e+59 | 1 |
| -0.3 | 6 | yes | 2.375 | 0.263705 | 0.263705 | 0.266475 |
| -0.5 | 2 | yes | 0.625 | 0.157937 | 0.158698 | 0.157937 |
| -0.5 | 3 | no | 26 | 1 | 1.82802e+44 | 1 |
| -0.5 | 4 | yes | 1.375 | 0.232111 | 0.240975 | 0.232111 |
| -0.5 | 5 | no | 42.375 | 1 | 1.3138e+49 | 1 |
| -0.5 | 6 | yes | 2.125 | 0.248583 | 0.248583 | 0.265919 |
| -0.6 | 2 | yes | 0.625 | 0.134015 | 0.134015 | 0.174979 |
| -0.6 | 3 | no | 23.75 | 1 | 3.47498e+38 | 1 |
| -0.6 | 4 | yes | 1.375 | 0.217236 | 0.217236 | 0.257173 |
| -0.6 | 5 | no | 39 | 1 | 3.32114e+42 | 1 |
| -0.6 | 6 | yes | 2 | 0.256224 | 0.256224 | 0.257587 |
| -0.7 | 2 | yes | 0.625 | 0.116593 | 0.116593 | 0.189237 |
| -0.7 | 3 | no | 21.875 | 1 | 3.61373e+33 | 1 |
| -0.7 | 4 | yes | 1.25 | 0.221996 | 0.251182 | 0.221996 |
| -0.7 | 5 | no | 36.125 | 1 | 6.493e+36 | 1 |
| -0.7 | 6 | yes | 1.875 | 0.241512 | 0.27508 | 0.241512 |

No n0=7,8 grid maximum is available until complete Z14/Z16 spectra exist.
In particular,
the larger-ring Krylov results are **not** silently treated as complete spectra or used
in `logZ_bounds`. Odd-n near-unit diagnostic ratios cannot establish the theorem.

### Best eligible point per parameter (ranked by balanced purity)

These are failed candidate seeds, not certificates. Their stored rational defects are
upward decimal roundings of numerical defects, **not upper enclosures of true defects**.
The checker recomputes `certified_gap_bound` on those exact rational inputs and verifies
that `.bound` is `nothing`. No rational seed is advertised as validated.

| theta | n0 | beta | epsilon_s | epsilon_t | certified bound | bound / ED gap |
|---:|---:|---:|---|---|:---:|:---:|
| 0 | 6 | 3.25 | 167806529617609//250000000000000 | 655412896527563//1000000000000000 | none | N/A |
| -0.1 | 6 | 2.875 | 708209229326509//1000000000000000 | 697959185079679//1000000000000000 | none | N/A |
| -0.3 | 6 | 2.375 | 736294903548947//1000000000000000 | 733525378728861//1000000000000000 | none | N/A |
| -0.5 | 6 | 2.125 | 187854195229087//250000000000000 | 367040343114593//500000000000000 | none | N/A |
| -0.6 | 6 | 2 | 371887996489099//500000000000000 | 185603260026111//250000000000000 | none | N/A |
| -0.7 | 6 | 1.875 | 724919880870063//1000000000000000 | 758488383533373//1000000000000000 | none | N/A |

The closest computed eligible point is **theta=0, n0=6, beta=3.25**, with balanced purity **0.328774**. This ranking is numerical and limited to the stated grid.

## Ground-energy mismatch and ED reference gaps

`delta(n)=E0(2n)-2 E0(n)` includes the finite-size nonextensive ground-energy contribution.
Positive even-n delta makes Ps fall as exp(-beta delta) in the ground-state-dominated
regime. This competes with the larger beta required for thermal purity. Odd-n deltas are
shown separately through the parity column and must not be fit together with even n.
The reference gap column is the **finite-ring** E1-E0 at 2n, not a bulk gap assertion.

At L=14 and 16, `reference.jl` searches all momentum and spin-reversal/reflection blocks
in M=0, with six converged low levels per block. Every integer-spin SU(2) multiplet has an
M=0 state, so this covers low multiplet energies numerically. Krylov convergence and
residual norms are recorded, but multiplicities are not certified. These are numerical
ED references only, never partition-function inputs or uniform gap lower bounds.

| theta | n | parity | delta(n) | ED gap at 2n | full spectrum at 2n? |
|---:|---:|:---:|---:|---:|:---:|
| 0 | 2 | even | 2 | 1 | yes |
| 0 | 3 | odd | -2.61742 | 0.720627 | yes |
| 0 | 4 | even | 0.663044 | 0.593555 | yes |
| 0 | 5 | odd | -1.03187 | 0.524808 | yes |
| 0 | 6 | even | 0.36529 | 0.484196 | yes |
| 0 | 7 | odd | -0.510737 | 0.458965 | no, reference only |
| 0 | 8 | even | 0.227105 | 0.442796 | no, reference only |
| -0.1 | 2 | even | 2.50886 | 1.17499 | yes |
| -0.1 | 3 | odd | -3.52649 | 0.803729 | yes |
| -0.1 | 4 | even | 0.833391 | 0.63211 | yes |
| -0.1 | 5 | odd | -1.54405 | 0.535898 | yes |
| -0.1 | 6 | even | 0.482296 | 0.476218 | yes |
| -0.1 | 7 | odd | -0.835243 | 0.436861 | no, reference only |
| -0.1 | 8 | even | 0.320877 | 0.409849 | no, reference only |
| -0.3 | 2 | even | 3.38206 | 1.55259 | yes |
| -0.3 | 3 | odd | -4.29421 | 1.01185 | yes |
| -0.3 | 4 | even | 1.12472 | 0.757481 | yes |
| -0.3 | 5 | odd | -2.69257 | 0.610228 | yes |
| -0.3 | 6 | even | 0.680991 | 0.514999 | yes |
| -0.3 | 7 | odd | -1.63748 | 0.448967 | no, reference only |
| -0.3 | 8 | even | 0.48199 | 0.400934 | no, reference only |
| -0.5 | 2 | even | 4.0591 | 1.92963 | yes |
| -0.5 | 3 | odd | -4.00439 | 1.24593 | yes |
| -0.5 | 4 | even | 1.34513 | 0.921334 | yes |
| -0.5 | 5 | odd | -2.72087 | 0.730862 | yes |
| -0.5 | 6 | even | 0.827735 | 0.605661 | yes |
| -0.5 | 7 | odd | -2.09587 | 0.517228 | no, reference only |
| -0.5 | 8 | even | 0.599448 | 0.45157 | no, reference only |
| -0.6 | 2 | even | 4.32085 | 2.10766 | yes |
| -0.6 | 3 | odd | -3.8291 | 1.36407 | yes |
| -0.6 | 4 | even | 1.42743 | 1.00992 | yes |
| -0.6 | 5 | odd | -2.56682 | 0.801343 | yes |
| -0.6 | 6 | even | 0.880739 | 0.663674 | yes |
| -0.6 | 7 | odd | -1.95282 | 0.565998 | no, reference only |
| -0.6 | 8 | even | 0.640803 | 0.493127 | no, reference only |
| -0.7 | 2 | even | 4.5295 | 2.27455 | yes |
| -0.7 | 3 | odd | -3.63279 | 1.47935 | yes |
| -0.7 | 4 | even | 1.49079 | 1.09953 | yes |
| -0.7 | 5 | odd | -2.40723 | 0.875251 | yes |
| -0.7 | 6 | even | 0.920176 | 0.726833 | yes |
| -0.7 | 7 | odd | -1.813 | 0.621247 | no, reference only |
| -0.7 | 8 | even | 0.670672 | 0.542253 | no, reference only |

## Conditional required-scale forecasts

These are **model-dependent extrapolations, not evidence of an actual seed at those sizes**.
Fit the positive even-n mismatches (n≥4) to both `A exp(-n/xi)` and
`A exp(-n/xi)/sqrt(n)`. Independently extrapolate the last three even finite-size gaps with
`Delta(L)=Delta_inf+a/L` and `Delta(L)=Delta_inf+a/L^2`. The finite-size data are short,
especially near -pi/4, so none of these extrapolations has controlled asymptotic error.
A negative extrapolated gap produces no forecast, rather than being clipped positive.

For an assumed gap g, use the independent-triplet ansatz
`z_L(beta)=(1+3 exp(-beta*g))^L` after removing the ground energy. This gives
`Ps=exp(-beta*delta(n))` and
`Pt=[(1+3 exp(-2 beta*g))/(1+3 exp(-beta*g))^2]^(2n)`.
Scan even n=8:2:2000 and beta*g=2:0.125:30, and report the first n at which the
actual recurrence returns a positive conditional bound. This ansatz ignores dispersion,
interactions, and the detailed low-energy density of states; it is a transparent scale
estimate, not a fit to the full partition function. The variation across four choices is
model sensitivity, **not a confidence interval**.

| theta | gap fit | assumed bulk gap | delta fit | fitted xi | first n0 | beta0 |
|---:|---|---:|---|---:|---:|---:|
| 0 | gap intercept vs L^(-1) | 0.317551 | exponential | 3.73333 | 22 | 23.6183 |
| 0 | gap intercept vs L^(-1) | 0.317551 | exponential_sqrt | 5.51834 | 24 | 25.5864 |
| 0 | gap intercept vs L^(-2) | 0.389494 | exponential | 3.73333 | 20 | 20.5395 |
| 0 | gap intercept vs L^(-2) | 0.389494 | exponential_sqrt | 5.51834 | 24 | 19.8976 |
| -0.1 | gap intercept vs L^(-1) | 0.209794 | exponential | 4.19091 | 26 | 38.1326 |
| -0.1 | gap intercept vs L^(-1) | 0.209794 | exponential_sqrt | 6.58032 | 32 | 39.3242 |
| -0.1 | gap intercept vs L^(-2) | 0.325032 | exponential | 4.19091 | 24 | 24.6129 |
| -0.1 | gap intercept vs L^(-2) | 0.325032 | exponential_sqrt | 6.58032 | 30 | 24.6129 |
| -0.3 | gap intercept vs L^(-1) | 0.0581775 | exponential | 4.7205 | 36 | 148.253 |
| -0.3 | gap intercept vs L^(-1) | 0.0581775 | exponential_sqrt | 7.98732 | 48 | 158.996 |
| -0.3 | gap intercept vs L^(-2) | 0.256094 | exponential | 4.7205 | 30 | 31.2385 |
| -0.3 | gap intercept vs L^(-2) | 0.256094 | exponential_sqrt | 7.98732 | 38 | 33.1909 |
| -0.5 | gap intercept vs L^(-1) | -0.0109547 | exponential | — | no forecast | — |
| -0.5 | gap intercept vs L^(-1) | -0.0109547 | exponential / sqrt(n) | — | no forecast | — |
| -0.5 | gap intercept vs L^(-2) | 0.256349 | exponential | 4.94906 | 32 | 31.2075 |
| -0.5 | gap intercept vs L^(-2) | 0.256349 | exponential_sqrt | 8.66439 | 42 | 33.1579 |
| -0.6 | gap intercept vs L^(-1) | -0.0186599 | exponential | — | no forecast | — |
| -0.6 | gap intercept vs L^(-1) | -0.0186599 | exponential / sqrt(n) | — | no forecast | — |
| -0.6 | gap intercept vs L^(-2) | 0.277173 | exponential | 4.99432 | 32 | 29.3138 |
| -0.6 | gap intercept vs L^(-2) | 0.277173 | exponential_sqrt | 8.80404 | 42 | 31.5687 |
| -0.7 | gap intercept vs L^(-1) | -0.0115617 | exponential | — | no forecast | — |
| -0.7 | gap intercept vs L^(-1) | -0.0115617 | exponential / sqrt(n) | — | no forecast | — |
| -0.7 | gap intercept vs L^(-2) | 0.308603 | exponential | 5.00763 | 32 | 25.9233 |
| -0.7 | gap intercept vs L^(-2) | 0.308603 | exponential_sqrt | 8.84551 | 42 | 27.9485 |

For theta=-0.5,-0.6,-0.7, the incompatible signs of the two extrapolated gap intercepts
leave the required scale **unresolved**. The finite values shown for those angles are
only outputs of the positive-intercept ansatz. The short series cannot distinguish
the assumed exponential mismatch tail from a long preasymptotic regime.

## Resources and reachable sizes

The environment has Julia **1.12.4**, a **16 GiB** cgroup memory limit and **4 CPU** quota.
The juliaup installer returned HTTP 403. The official Julia 1.12.4 S3 tarball succeeded;
`Pkg.instantiate()` and `Pkg.precompile()` completed using the unmodified core project
and manifest. The runtime/depot are ignored and confined to this certificate directory.
All runs set `BLAS.set_num_threads(1)` and `OPENBLAS_NUM_THREADS=1`; dense parameter solves
use four Julia threads. No existing tracked repository file was modified.

At L=14 the initial M=1,k=1 block has **41,835** states before SU(2) reduction.
One real dense matrix requires **13.04 GiB**; matrix construction plus eigensolver copies
in the implemented path requires more than the available memory (two matrices alone:
**26.08 GiB**). Cubic work is approximately **376 times** the corresponding L=12 block.
At L=16 that block has **310,572** states and needs **718.65 GiB** for one matrix.
`preflight.json` supplies exact integer character counts and the resource projections.

A numerical SU(2) highest-weight construction is now implemented in `su2.jl`.
Its polynomial Casimir projector is checked exactly on small charge spaces, and weighted
complete spectra agree with unreduced dense ED through L=7. A Householder-QR pilot at
L=14, S=2, k=1 produced an 8,106-state dense block in 649.6 seconds including one angle's
eigensolve. The faster CholeskyQR2 version passed the same small tests; its first full
L=14, S=0, k=0 block took 61.8 seconds for all six angles. No incomplete sector collection
is used for Z14. The full sweep is resumable with `run_su2_grid.sh`, with a
45-minute timeout per sector process. See `review_notes.md` and the SU2 pilot JSON.
The exact character counts give largest `(S,k)` highest-weight multiplicities **8,482**
at L=14 (0.536 GiB per dense matrix) and **59,143** at L=16 (26.06 GiB). These counts use
`d(S,k)=d(M=S,k)-d(M=S+1,k)` and assert `sum_(S,k) (2S+1)d(S,k)=3^L`.
Thus this work does **not** establish that full 14-site ED is intrinsically infeasible:
the former basis-construction blocker has been resolved numerically.
The remaining full sector sweep is the current blocker to the larger-ring purity table. No
memory-overcommitted 14-site dense solve or incomplete-spectrum purity calculation was run.
With one-site cells, full Z14 adds the odd-n0=7 diagnostic; the next eligible seed
beyond n0=6 is n0=8, which needs complete Z16 under the stated two-purity workflow.

## Rigor boundary and eigenvalue-enclosure theorem

**Proved exactly:** the RP coefficient identities, real on-site term, and coefficient
signs in the stated parameter window. Integer Hilbert-space counts are asserted.
**Numerical only:** matrix symmetry checks, every computed spectrum, all purities, failed
seed defects, finite-size gaps, delta fits, rankings, and extrapolated required scales.
**Not established:** a true-spectrum eigenvalue enclosure, an admissible rational seed,
a uniform gap, uniqueness from the criterion, or any thermodynamic extrapolation error bar.

The requested best-seed rigorous rerun is not triggered because no eligible candidate
passes. Consequently no residual/orthogonality eigenvalue enclosure is claimed or
implemented. For clarity, the theorem such a future path must apply is the polar-factor
Weyl bound: for square approximate eigenvectors Q and real diagonal Lambda, if verified
`||Q'*Q-I||_2 <= eta < 1`, `||H*Q-Q*Lambda||_2 <= rho`, then there is a unitary polar factor U
with `||U-Q||_2 <= 1-sqrt(1-eta)`, and sorted eigenvalues of H lie within
`rho+(||H||_2+||Lambda||_2)*(1-sqrt(1-eta))` of sorted Lambda. This follows by bounding
`U'*H*U-Lambda` and applying Weyl. All norms, matrix assembly errors, trigonometric
coefficients and arithmetic errors must be enclosed outward. Merely accumulating in
BigFloat without a rounding bound is not sufficient.

Given eigenvalue intervals, logZ lower and upper bounds follow by summing exp(-beta E)
with directed rounding, and purity lower bounds use
`exp(logZ2n_lo-2 logZn_hi)` and `exp(logZ2n_2beta_lo-2 logZ2n_beta_hi)`.
Only then may rational upper defect bounds be passed to the certified recurrence.
The numerical `PartialSpectrum.exact=true` flag here means *complete*, not *certified*.

## Reproduction and validation

From the repository root (the wrapper selects the downloaded Julia and private depot):

```sh
certificates/blbq_haldane/run.sh -t 4 certificates/blbq_haldane/scan.jl 2 3 4 5 6 7 8 10 12
certificates/blbq_haldane/run.sh certificates/blbq_haldane/reference.jl 14
certificates/blbq_haldane/run.sh certificates/blbq_haldane/reference.jl 16
certificates/blbq_haldane/run.sh certificates/blbq_haldane/preflight.jl
certificates/blbq_haldane/run.sh certificates/blbq_haldane/analyze.jl
certificates/blbq_haldane/run.sh certificates/blbq_haldane/report.jl
certificates/blbq_haldane/run.sh certificates/blbq_haldane/test/runtests.jl
```

Each invocation runs `julia --project=julia`; each completed numerical run is below the
45-minute cap. Stored timers cover the calculation functions, excluding environment setup
and initial package imports. Result JSON uses the shared core's `Results.write_result`
(ThermalMPOGap historically, CertifiedGaps after the rename). Large complete spectra are
losslessly gzip-compressed; `archive.json` records original and compressed hashes.
and contains source hashes and run metadata. Exact rational tokens are strings. Spectrum
files retain every eigenvalue and every sector dimension, enabling independent numerical
reanalysis without repeating ED. `scaling.json` hashes its spectrum/reference inputs.
The fast test recomputes exact RP, small-system spectrum crosschecks, complete dimension
counts, and the failed recurrence outputs from stored exact rational tokens. It does not
pretend to certify the rejected numerical input defects. See `validation_current_main.json`
for current-main test timing and `historical/thermal_validation.json` for the original
snapshot. `MIGRATION.md` describes the independent source/dependency comparison. No merge
has been performed; draft publication is tracked in `PR.md` when available.

### Larger-ring reference runtimes

| L | Wall seconds |
|---:|---:|
| 14 | 95.5662 |
| 16 | 1209.74 |
