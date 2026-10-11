# Spin-1 dimerized ring: uniform-gap certificate

The smallest successfully certified parameter on the requested grid is **delta=7/10**.
The seed n=3 cells, beta=17/4 gives
`certified_gap_bound(...).bound = 0.2187941698399417`. In particular,
**Delta_L >= 0.218794169839** with a unique ground state for every even L>=3
cells: physical spin counts N=8,12,16,... . The exact rational seeds are
epsilon_s=58597799297/500000000000 and epsilon_t=48363543279/1000000000000.
See [the certificate](results/certificate_d7.json) and [parameter summary](results/summary.json).

The twelve-spin ED reference gap is 1.2659290109040136; the certificate is
approximately 17.2833% of this **numerical finite-size** reference. The full
verified rerun took 1017.27 seconds (16.95 minutes). All five parameter scans
completed through six cells. No seed was found for delta=0.6 or 0.5 on the stated
grid; that is numerical evidence at the size cap, not a no-go theorem.

## Model, lengths and scope

We use spin one, H = sum_j [1+(-1)^j delta] S_j.S_(j+1), periodic,
with two spins per cell (cell dimension 9), and delta in {9/10,4/5,7/10,3/5,1/2}.
All energies and gaps are in the normalization of this Hamiltonian. A ring of L
cells has N=2L spins. Ring sizes 2,3,4,6 cells support seed sizes n=2,3 only:
a seed n needs the ring 2n. In particular, the six-cell cap does **not** allow
seed sizes 4,5,6, which would require 16,20,24 spins respectively.

The beta grid is 1/4,1/2,...,32. The search uses `Seeds.scan_seeds` and
`Seeds.best_seed`. Those search outputs are numerical, not certified, even though
they call the outward-rounded recurrence. Only results named `certificate_*.json`
have verified spectral enclosures and rational upper bounds on the seed defects.
A finite grid failure is not a proof that the model cannot be certified.

## Exact reflection positivity, including the cell

In the Cartesian basis let

```
Gx = [0 0 0; 0 0 1; 0 -1 0]
Gy = [0 0 -1; 0 0 0; 1 0 0]
Gz = [0 1 0; -1 0 0; 0 0 0]
Sa = -i Ga, so Ga = i Sa.
```

These integer matrices are real and antisymmetric. The exact checks include
[Sx,Sy]=iSz and sum_a Sa^2=2I. For either c=1-delta or c=1+delta,

    -c S.S = c sum_a Ga tensor Ga,    c>0.

All coefficients are the stated positive rationals. This is an identity over
complex rationals, not a floating point `rp_test`.

Inside a cell (a,b), D=-(1-delta) sum_a Ga tensor Ga is real symmetric.
Let P swap the two sites and Theta=P K. Set A_a=I tensor Ga, B_a=Ga tensor I.
The intercell bond obeys

    -h_inter=(1+delta) sum_a A_a tensor Theta(A_a),
    Theta(A_a)=P conjugate(A_a) P=B_a,    Theta(D)=D.

`rp.jl` checks these identities exactly, including the 81-by-81 intercell
interaction. The cell reflection is essential: the physical intercell interaction
is **not** a sum of identical cell operators A_a tensor A_a in the naive product
cell basis. We use the internal-reflection extension in the scope note.

In fact one obtains the stronger positivity needed for an odd-cell seed. Let
F(a,b) be the one-site traced event kernel, with event histories on the two incident
bonds. Each insertion is one of the same real Ga matrices irrespective of which
bond supplied it. Therefore F is real and F(a,b)=F(b,a), preserving time order.
Let mu_- and mu_+ be the positive event measures for weak and strong bonds. Cutting
at the strong bonds and integrating the weak-bond history u gives the cell kernel

    X(a,b)=integral F(a,u) F(u,b) dmu_-(u)=K K* (a,b),

on L2(mu_+). Bounded kernels and finite event measures make K Hilbert-Schmidt.
Thus X is positive trace-class and Z_L=Tr X^L for every integer L>=1.
The only use of even n in the supplied self-adjoint-transfer proof is positivity
of X^n and the real-eigenvalue moment interpolation. Here all eigenvalues of X
are nonnegative, so these steps hold for odd n as well. The same recurrence and
its certified implementation consequently apply to n=3. We report only the
requested uniform conclusion for even L>=n, not a stronger length range.
Uniqueness in this range follows from the two-purity criterion itself.

## Complete dense spectra and exact symmetries

In the Sz product basis each allowed exchange of neighboring spin-one digits has
matrix element J_j exactly: the two square-root ladder coefficients cancel the
factor 1/2. Consequently `integer_hamiltonian` builds an exact symmetric integer
matrix A with H=A/denominator(delta). Its diagonal is sum_j J_j m_j m_(j+1).
There is no numerical threshold for dropping entries.

The numerical search uses total Sz, translation by two sites, spin flip, and
spatial reflection. Charge M>0 is paired with -M with multiplicity 2; spin flip
is split into parities at M=0. Generic momenta k and L-k carry identical spectra,
so a real reflection-even copy is diagonalized and counted twice. At k=0 and
k=L/2 both reflection parities are separately diagonalized. These multiplicities
follow from the exact dihedral symmetry, not an assumption about eigenvalue
nondegeneracy. Translation by two preserves bond parity; reversal j -> N+1-j
maps bond j to N-j and also preserves bond parity; spin flip sends m -> -m and
preserves both diagonal and exchange entries. These prove commutation with H.
Every independent block is diagonalized with dense LAPACK; no Lanczos/ARPACK,
truncated spectrum, missing-state estimate or DMRG is used.

The verified small-ring path uses complete dense Sz blocks directly. The large-ring
path uses real character projectors on individual permutation-group orbits. For
L in {2,3,4,6}, the coefficients 2 cos(2 pi k r/L) are exact integers. Each
projector is represented as Pnum/den with integer entries, and symmetry,
Pnum^2=den Pnum and integer trace/rank are asserted exactly for every orbit.
The numerical basis has precisely this exact rank. Sector dimensions with all
multiplicities are asserted to sum to 9^L, both in the search and verification.
Full block records are retained in the JSON results.

## Eigenvalue enclosure theorem and implemented rounding bounds

Let H be the exact real symmetric block, Q the square numerical eigenvector matrix,
and Lambda the diagonal matrix of sorted numerical eigenvalues. Suppose verified
upper bounds r>=||H Q-Q Lambda||_F and e>=||Q'Q-I||_F satisfy e<1.
Writing Q=U S with U orthogonal and S=(Q'Q)^(1/2), put

    sigma=sqrt(1-e),    d=e/(1+sigma).

Then ||Q^(-1)||<=1/sigma and ||S-I||<=d. The identity

    H-U Lambda U' = R Q^(-1) + U (S Lambda S^(-1)-Lambda) U'

gives

    ||H-U Lambda U'|| <= eta=(r+2||Lambda|| d)/sigma.

Weyl's variational eigenvalue inequality now implies, **with multiplicity and
in sorted order**, |lambda_i(H)-Lambda_i|<=eta. This is not an unverified
residual-per-vector claim near degenerate eigenvalues.

`eigen_enclosure` computes R and E=Q'Q-I in Float64, then adds rigorous
floating point accumulation error bounds before applying this theorem.
Let u=2^-53 and gamma_k=k u/(1-k u). For a d-dimensional block it uses
gamma_(2d+8), the Frobenius product inequalities, and

    r <= (||Rcomputed||_F + gamma*(||A||_F ||Q||_F
                    + den ||Q||_F ||Lambda||) + d^3 floatmin)/den,
    e <= ||Ecomputed||_F + gamma*(||Q||_F^2+sqrt(d)) + d^3 floatmin.

These conservative dot-product bounds allow separate multiplies/additions or FMA.
For each computed norm, the Float64 sum of squares with m entries is corrected
by gamma_(2m+2), plus an absolute underflow allowance, before a directed-rounded
MPFR square root. Specifically sqrt((sum_computed+m*floatmin)/(1-gamma)) is
an upper bound. All scalar calculations use 128-bit BigFloat directed rounding;
denominators needing lower bounds are separately rounded down. Nonfinite sums
and excessive orthogonality defects abort. No `@fastmath` is used, and BLAS is
single-threaded. This is the rigorous-norm-bound alternative in requirement 3(c),
not a claim that unverified Float64 residual norms suffice. Radii are rounded up
to exact rational multiples of 10^-12.

For large-ring symmetry bases an additional construction error is included.
For exact orbit projector P and floating basis B, verified bounds on
||PB-B||=q and ||B'B-I||=e imply

    ||(PB)'(PB)-I|| <= w=e+2||B||q+q^2,
    ||B-U|| <= q+w/(1+sqrt(1-w)),

where U is an exactly orthonormal basis of the exact projector range. These
small-matrix residuals include IEEE dot-product and exact-projector conversion
errors using interval arithmetic. On disjoint orbits the maximum such distance
bounds the whole basis distance. If that bound is d_B, replacing B by U changes
B' H B by at most ||H||(2d_B+d_B^2), with the rigorous bound ||H||<=4L.
Sparse projection multiplication and rational division errors are bounded using
gamma_(2 maxnnz(A row)+2 maxnnz(B column)+12) times
2 ||A||_F ||B||_F^2/den. The factor two also covers copying the upper triangle
into the lower triangle. This construction allowance is added to every block's
Weyl radius. Thus numerical orbit eigenvectors and momentum transformations
are not silently assumed exact.

## Partition functions, seeds and bound

For all enclosed eigenvalues, `verified_logz` sums exp(-beta E) with
IntervalArithmetic, after a common exactly represented energy shift, including
every multiplicity. It stores rational lower and upper endpoints for log Z_n(beta),
log Z_2n(beta), and log Z_2n(2 beta). The seed upper bounds are

    epsilon_s >= 1-exp(lower(log Z_2n(beta))-2 upper(log Z_n(beta))),
    epsilon_t >= 1-exp(lower(log Z_2n(2 beta))-2 upper(log Z_2n(beta))).

The actual implementation performs the equivalent interval operations and takes
the upper endpoint, rounded up to a rational with denominator 10^12.
These rational defects and exact rational beta are passed to
`ThermalMPOGap.certified_gap_bound`. Its `.bound` is the reported uniform lower
bound. The reference ED gap and bound/reference ratio are numerical; neither is
used as an input to the certificate.

## Reproduction and what the fast checker verifies

Run from the repository root, with Julia 1.12.4 and the core environment:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/spin1_dimerized/scan.jl 7 2,3,4,6
julia --project=julia certificates/spin1_dimerized/certify.jl 7
julia --project=julia certificates/spin1_dimerized/test/runtests.jl
julia --project=julia certificates/spin1_dimerized/test/verified_blocks.jl
```

Use parameter integers 9,8,7,6,5 for delta=integer/10. Existing spectrum files are
reused by the scan; remove only that parameter's spectrum files to rerun ED.
The optional certificate arguments select a stored-grid seed, e.g. `9 2 7//2`.
The fast checker reruns exact RP, checks symmetry spectra against direct dense ED
at four spins and the decoupled-dimer spectrum, checks exact projector identities,
checks stored dimensions, and recomputes the bound from stored rational seeds.
It does **not** redo twelve-spin diagonalization or independently recertify the
stored partition sums; that requires the full `certify.jl` run.

Final validation: `test/runtests.jl` passed all 115 assertions in 9.37 seconds,
including test compilation, and `test/verified_blocks.jl` passed all 16 assertions
(6.4 seconds in its test set). [The test result](results/tests.json) records
provenance and timing. The main checker verifies both stored certificates and
their source hashes. Scan/verification timings in the tables exclude Julia
startup and package precompilation. Every individual calculation was below the
requested 45-minute run limit. There are no remaining execution blockers.

Source SHA-256 hashes, Julia/package versions, thread counts and run configurations
are stored through `ThermalMPOGap.Results.write_result`/`run_metadata`. Exact
rationals are strings. All repository edits are confined to this certificate
folder. Subsequent draft-PR preparation, upstream package-name compatibility,
independent review and replay are recorded in [PUBLICATION.md](PUBLICATION.md).
No changes to the shared core are proposed.

Environment: the approved juliaup installer URL returned HTTP 403. The official
Julia 1.12.4 Linux x86-64 tarball was downloaded from julialang-s3.julialang.org
and unpacked under /tmp; `Pkg.instantiate(); Pkg.precompile()` succeeded using
`JULIA_DEPOT_PATH=/tmp/julia-depot`. The project and manifest remained unchanged.

Certified: exact RP/cell reflection, the positive-cell-transfer argument above,
complete symmetry counts, the selected verified spectral enclosures and partition
bounds, rational seed bounds, and the resulting uniform gap/uniqueness conclusion.
Numerical: search-grid values, unverified ED reference gaps, parameter ranking
outside the certified seed, and any extrapolation beyond the six-cell ring cap.

For commands in this particular cloud environment, first set
`PATH=/tmp/julia-1.12.4/bin:$PATH` and `JULIA_DEPOT_PATH=/tmp/julia-depot`.
The route memo was committed locally as required by `AGENTS.md`. Draft publication
status is tracked separately; nothing is merged into the target branch.

## Independent audit and fresh replay

[The independent audit](AUDIT.md) found no substantive mathematical blocker.
A fresh full spectral computation against `CertifiedGaps` completed in 949.35
seconds and reproduced every mathematical output exactly, including all sector
error bounds, log-partition endpoints, rational seeds and the final bound.
The replay comparison passed 40 assertions; current-core fast and supplemental
checks passed 115 and 16 assertions respectively. Original outputs are preserved.
See [PUBLICATION.md](PUBLICATION.md) and the separate
[comparison record](results/replay_comparison.json) for provenance and limitations.

## Results from completed runs

Best numerical seed per parameter, over n=2,3 and beta=0.25:0.25:32:

| delta | n | beta | epsilon_s | epsilon_t | numerical recurrence bound | ED gap (12 spins) |
|---|---:|---:|---:|---:|---:|---:|
| 9//10 | 3 | 4.75 | 0.003735505 | 0.004751553 | 0.891357627 | 1.762670627 |
| 4//5 | 3 | 4.25 | 0.031253772 | 0.022772241 | 0.544718970 | 1.517667420 |
| 7//10 | 3 | 4.25 | 0.117194900 | 0.048360599 | 0.218800397 | 1.265929011 |
| 3//5 | — | — | — | — | no seed | 1.010968921 |
| 1//2 | — | — | — | — | no seed | 0.764485184 |

Verified certificates (the 0.9, n=2 result was verified before the larger-ring scan):

| delta | n | exact beta | certified lower bound | ratio to 12-spin ED gap | verification seconds | validity |
|---|---:|---:|---:|---:|---:|---|
| 9//10 | 2 | 7//2 | 0.528801949327 | 0.300000 | 1.01 | even L >= 2 cells; N divisible by 4, N >= 4 spins |
| 7//10 | 3 | 17//4 | 0.218794169839 | 0.172833 | 1017.27 | even L >= 3 cells; N divisible by 4, N >= 8 spins |

The exact Float64 lower bound is stored also as an exact binary rational in each certificate JSON; these table bounds are rounded down to twelve decimal places.

| delta | max min purity, n=2 | max min purity, n=3 | E0(4)-2E0(2) | E0(6)-2E0(3) | extrapolated (n,beta), numerical only |
|---|---:|---:|---:|---:|---|
| 9//10 | 0.950284354 | 0.996076163 | 0.0144862082 | 0.0007826475 | not needed |
| 4//5 | 0.829648851 | 0.968746228 | 0.0621074387 | 0.0073561969 | not needed |
| 7//10 | 0.673733633 | 0.895128915 | 0.1465020299 | 0.0289358279 | not needed |
| 3//5 | 0.539672632 | 0.762071657 | 0.2639962827 | 0.0779594100 | (5, 7.25) |
| 1//2 | 0.409109875 | 0.585111913 | 0.3984628798 | 0.1636639315 | (7, 9.75) |

The extrapolation fits delta(n)=A exp(-kappa n) through n=2,3 and holds the numerical twelve-spin gap fixed. It estimates epsilon_t=12n exp(-beta gap) and epsilon_s=1-exp(-beta delta(n)), targeting both below 0.05. Two sizes cannot validate this scaling or an infinite-volume gap; these are planning estimates, not demonstrated seeds.

| cells | spins | full dimension | largest numerical block | largest full Sz block |
|---:|---:|---:|---:|---:|
| 2 | 4 | 81 | 7 | 19 |
| 3 | 6 | 729 | 42 | 141 |
| 4 | 8 | 6561 | 254 | 1107 |
| 6 | 12 | 531441 | 11596 | 73789 |

Full sector dimension lists and multiplicities are in `results/spectrum_d*_c*.json` and `results/certificate_d*.json`.

| delta | dense search seconds (sum across sizes) | twelve-spin seconds |
|---|---:|---:|
| 9//10 | 549.40 | 549.20 |
| 4//5 | 538.11 | 537.91 |
| 7//10 | 536.69 | 536.48 |
| 3//5 | 546.82 | 546.55 |
| 1//2 | 537.86 | 537.65 |
