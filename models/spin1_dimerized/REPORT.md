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

`src/rp.jl` checks these identities exactly, including the 81-by-81 intercell
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
`CertifiedGaps.certified_gap_bound`. Its `.bound` is the reported uniform lower
bound. The reference ED gap and bound/reference ratio are numerical; neither is
used as an input to the certificate.

## Validation and provenance

Commands and run times are in README.md. The fast checker (`check.jl` = `test/runtests.jl`,
123 tests, about 55 s including compilation) reruns exact RP, checks the symmetry spectra
against direct dense ED at four spins and the decoupled-dimer spectrum, checks the exact
projector identities and the verified block path on a small ring, checks stored dimensions,
recomputes the rational seeds from the stored log-partition enclosures and the bound from the
seeds, verifies source provenance, and recertifies delta = 9/10 from scratch (rings of 4 and 8
spins). For delta = 7/10 it does **not** redo the twelve-spin diagonalization; `check.jl
--full` does (17–36 minutes). A fresh recertification is not expected to be bit-identical
on another machine, since LAPACK output differs in the last bits; the checks therefore
require the fresh certificate to be valid on its own and its log Z enclosures to overlap the
stored ones (both enclose the same exact numbers). On an Apple M-series machine the fresh
delta = 9/10 bound is 0.5288019471885558 against the stored 0.5288019493276747, and the
fresh delta = 7/10 bound (`check.jl --full`, 36 min under load) is 0.21879416840454358 against
the stored 0.2187941698399417.

Source SHA-256 hashes, Julia/package versions, thread counts and run configurations are
stored through `CertifiedGaps.Results.write_result`/`run_metadata`. Exact rationals are
strings. The stored certificates were produced when the core was named `ThermalMPOGap` and
this directory was `certificates/spin1_dimerized/`; their source hashes refer to archived
commit `6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2` and are checked against it when that commit
is in the local object store. Scan and verification timings exclude Julia startup and package
precompilation; every individual calculation was below 45 minutes.

Certified: exact RP/cell reflection, the positive-cell-transfer argument above, complete
symmetry counts, the selected verified spectral enclosures and partition bounds, rational
seed bounds, and the resulting uniform gap/uniqueness conclusion. Numerical: search-grid
values, unverified ED reference gaps, parameter ranking outside the certified seed, and any
extrapolation beyond the six-cell ring cap.

## Independent audit (2026-10-11)

Reviewer: an independent delegated review agent that did not author the calculation and did
not run a duplicate twelve-spin job. It inspected `rp.jl`, `common.jl`, `verified_blocks.jl`
(now in `src/`), the certificate code, both test files, the stored delta = 7/10 certificate,
and the transfer, coupling, interpolation and gap-recovery arguments in the repository notes.

**Finding:** no substantive mathematical obstruction to the odd three-cell seed. The positive
transfer proof below supplies the extra hypothesis that the generic self-adjoint-transfer
theorem does not supply. The enclosure algorithm has the required complete-spectrum and
outward-bound structure.

### 1. Exact model and transfer orientation

Write J-=1-delta and J+=1+delta. The integer Cartesian matrices in `rp.jl`
are real antisymmetric, with Sa=-i Ga, the spin-one commutation relation and
sum Sa^2=2I. Thus -J Sa tensor Sa=J Ga tensor Ga, with positive rational J
for every requested parameter. This is an exact identity, not a tolerance test.
For a cell with weak internal bond, D=-J- sum Ga tensor Ga is real symmetric.
The internal site swap P maps A=I tensor Ga to B=Ga tensor I and fixes D.
Consequently -h_inter=J+ sum A tensor (P conjugate(A) P).
The naive claim that the intercell operators are identical on both cell factors
would be false; the implemented internal-reflection identity has the correct
orientation.

Here is a direct transfer construction that avoids relying on that claim.
For J>0 let mu_J be the **unnormalised** measure on finite histories of labelled
events in [0,beta], with weight J^m on each ordered m-event simplex, and a unit
atom at the empty history. There are three labels. Its mass is exp(3 beta J).
For two histories a,u, define F(a,u) as the trace of the time-ordered product
of the Ga insertions from their union. Coincident event times have measure zero
and may be assigned a consistent order. Since both incident bonds insert exactly
the same Ga on this site, exchanging their history names changes no ordered
product: F(a,u)=F(u,a). All Ga are real, so F is real. Since ||Ga||=1,
|F|<=3, including the empty product.

Set H+=L2(mu_J+) and H-=L2(mu_J-). The integral operator K:H- -> H+ with
kernel F is Hilbert-Schmidt since both measures are finite and F is bounded.
Integrating the weak history inside a cell gives the strong-boundary kernel

    X(a,b) = integral F(a,u) F(u,b) dmu_J-(u) = (K K*)(a,b).

Therefore X is positive and trace class. Expanding the bounded finite-ring
Hamiltonian exponential into bond events and taking the physical trace gives
exactly Z_L(beta)=Tr X^L. Independent bond histories are equivalent to the
global ordered Dyson expansion; there is no Poisson normalisation factor
because the measures above are unnormalised. Absolute convergence follows from
the bounded insertions and finite history measures. This identity holds for
every positive integer number of cells (with both periodic bonds retained for
the one-cell ring); only L>=3 is needed here. In particular all transfer
eigenvalues lambda_i(beta) are nonnegative at every beta>0.

### 2. Why n0=3 has the same recurrence and uniform conclusion

For every integer n>=1, X^n/Z_n is a density operator. Its purity is
P_s(n,beta), so P_s<=1 and the spatial squaring lemma applies at odd n.
Thermal density operators give P_t<=1 without any parity assumption.
Partition-function cancellation gives the two coupling identities in the note
for every such n. Define

    z(s,t)=1-(1-s)^2(1-t),    f(z)=z^2/[2(1-z)^2].

Starting from defects bounded by (s,t), the first coupling gives a defect at
(n,2 beta) bounded by z(s,t). Spatial squaring gives
s'=f(z(s,t)) at (2n,2 beta). The second coupling gives a thermal defect at
(4n,beta) bounded by z(t,s'), and temporal squaring gives
t'=f(z(t,s')). This is precisely the unequal-defect recurrence used by the
certified implementation. Its derivation has no remaining parity restriction.

For n_j=2^j n0, beta_j=2^j beta0, and n_j<=L<=2n_j, nonnegative transfer
eigenvalues and monotonicity of sequence p-norms imply

    Z_L(beta_j) <= Z_nj(beta_j)^(L/n_j),
    Z_L(2 beta_j) >= Z_2nj(2 beta_j)^(L/(2 n_j)).

Hence, writing alpha=L/(2n_j)<=1,

    P_t(L,beta_j)
      >= [P_t(2n_j,beta_j) P_s(n_j,beta_j)^2]^alpha
      >= (1-t_j)(1-s_j)^2 = 1-z(s_j,t_j).

The certified orbit requires z(s_j,t_j)<1/2 at every scale. Since the largest
Gibbs eigenvalue w0 is at least its purity, w0>1/2. A degenerate ground state
would have w0<=1/2, proving uniqueness without an external uniqueness
assumption. Moreover

    exp(-beta_j Delta_L) <= (1-w0)/w0
                        <= z(s_j,t_j)/(1-z(s_j,t_j)).

Thus Delta_L>=b_j, where b_j is the rate bounded by the core recurrence's
admissibility and tail theorem. Its certified lower bound on inf_j b_j works
for every L>=n0 by selecting its enclosing dyadic interval. Restricting this
conclusion to the requested even cell lengths gives L=4,6,8,... for n0=3,
or N=8,12,16,... physical spins. It does not accidentally claim N=6 or all
even physical spin counts. No extension of the generic theorem has been
assumed without supplying the additional positivity and interpolation argument.

### 3. Complete spectra and numerical verification

The Sz Hamiltonian uses integer digits m=-1,0,1. Each allowed nearest-neighbour
exchange has coefficient J: the spin-one square roots cancel the 1/2 in the
ladder expression. Diagonal entries are the exact weighted m_j m_(j+1) sums.
For the specified sizes and rational parameters the integer entries fit exactly
in Float64. Charge reversal, two-site translation, and reversal of site order
preserve this Hamiltonian. Reversal preserves the alternating bond parity.

The character projectors use exact integer twice-cosines at cell periods
2,3,4,6; symmetry, idempotence after scaling, and trace/rank are asserted per
orbit. Positive and negative generic momenta have equivalent spectra. Restricting
to reflection + takes one copy, with multiplicity two: on a generic momentum
pair, T-T^(-1) gives an invertible intertwiner between the two reflection
parities. Both reflection parities are retained at k=0 and k=L/2. Spin flip
supplies multiplicity two for M>0 and is split explicitly at M=0. These are
exact representation facts, not numerical near-degeneracy assumptions.

An independent JSON read found 7 small-ring Sz blocks with weighted dimension
729 and maximum dimension 141, and 79 large-ring blocks with weighted dimension
531441 and maximum dimension 11596. These agree with 9^3 and 9^6. Dimension
counts alone would not prove completeness, but together with these disjoint
charge/character projectors and the exact rank assertions they do.

The polar/Weyl argument in the report is valid. If ||Q'Q-I||<=e<1 and
||HQ-Q Lambda||<=r, then Q=US has invertible S, sigma_min(S)>=sqrt(1-e),
and ||S-I||<=e/(1+sqrt(1-e)). Comparing H with U Lambda U' gives the stated
radius (r+2||Lambda|| ||S-I||)/sqrt(1-e). Weyl bounds every ordered eigenvalue
with multiplicity, including clusters. It is not merely a residual test for
isolated eigenvectors.

For symmetry bases, the independently checked projector residual controls PB-B;
the Gram defect of PB is at most e+2||B||q+q^2. The polar basis of PB spans the
exact projector range because its dimension equals the exact rank and this
defect is below one. Disjoint orbits make the full operator-norm distance the
maximum orbit distance. The bound 4L on ||H|| is valid: ||Sa.Sa||=2 for a
spin-one bond and the sum of alternating positive couplings is 2L.
The projected-matrix construction allowance and eigenvector allowance therefore
refer to the exact model, rather than silently treating numerical bases as exact.

Residuals and Gram defects are accumulated in Float64 and enlarged with normwise
IEEE error bounds; this is the allowed rigorous-norm-bound route, not high
precision matrix multiplication. The gamma factors conservatively include dot
products, multiplication, subtraction and scaling. Norm sums include an absolute
underflow allowance. Scalar bounds use directed MPFR or interval arithmetic,
including lower-rounded polar denominators. The projection bounds use substantial
gamma-factor slack. This audit assumes standard IEEE arithmetic, no fast-math
reassociation invalidating the error model, and the stated single-thread BLAS
configuration; it is not a formal verification of BLAS or MPFR. No substantive
direction error was found in these bounds.

### 4. Partition directions, stored seed, and audit limits

Each exact spectrum lies within the computed centres plus/minus rational radii.
Exponentials, summation, logarithms, and the common energy shift are enclosed by
IntervalArithmetic. For the spatial purity the safe exponent is lower(logZ_2n)
minus twice upper(logZ_n); for thermal purity it is lower(logZ_2n_2beta) minus
twice upper(logZ_2n_beta). The implementation uses these correct directions
and rounds the resulting defect upper endpoint upward to a rational.

The independently read original stored inputs are

    epsilon_s = 58597799297//500000000000
    epsilon_t = 48363543279//1000000000000
    beta0     = 17//4

with stored bound 0.2187941698399417. The decimal lower claim
0.218794169839 is weaker and correctly rounded downward. Numerical reference
gaps, scan rankings, and extrapolated seeds remain numerical observations.

The fast checker verifies exact RP, rational-seed calculations, counts and source
provenance, plus small numerical comparisons. It does not reconstruct twelve-spin
eigenvalue radii from model inputs. Source hashes attest provenance, not correctness;
the stored `positive_transfer_gram_proof` boolean is metadata, not a proof checker.
The mathematical argument above supplies that proof. No Lean formalisation of this
model, transfer construction, or numerical enclosure algorithm is claimed.

### 5. Computational replay

A fresh complete `certify.jl 7` run against the renamed core (949 s) was compared
with the preserved original `results/certificate_d7.json`, recursively excluding only
metadata and timings: **all mathematical fields were exactly identical**, including
all sector error bounds, partition enclosures, rational defects and the bound
0.2187941698399417. The replay output is not stored in the repository (see README).
No mathematical bound drift or remaining substantive audit blocker was found.

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

Full sector dimension lists and multiplicities are in `results/certificate_d*.json` and in the regenerated `results/spectrum_d*_c*.json` (`scripts/scan.jl`).

| delta | dense search seconds (sum across sizes) | twelve-spin seconds |
|---|---:|---:|
| 9//10 | 549.40 | 549.20 |
| 4//5 | 538.11 | 537.91 |
| 7//10 | 536.69 | 536.48 |
| 3//5 | 546.82 | 546.55 |
| 1//2 | 537.86 | 537.65 |
