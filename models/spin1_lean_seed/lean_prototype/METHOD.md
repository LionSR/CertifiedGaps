# Lean seed feasibility review: method, audit and budget

The long companion to [README.md](README.md): the Mathlib audit, the packed arithmetic and its
trust boundary, the packing and anchor experiments, a mathematical review of the proposed
Peierls/Weyl/counting route, benchmarks, and the remaining implementation budget. Lean files
named below live in `src/SeedPrototype/` (library) and `src/Bench/` (benchmark-only files).

The paper excerpts were fetched from the repository-pinned OpenAI commit
`adc7f1241b42e322a6451854ab7e4b4c146bf78a`, files
`preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/build/`
`{initialization,transfer,thermal72,trials}.tex`
(`scripts/fetch_references.sh` downloads them into the git-ignored `references/` and checks
`results/paper-source-sha256.txt`). Their inputs confirm 18 thermal tables, including **odd
lengths**, and the π twist about z.

## 1. Mathlib audit

Paths below are relative to `Mathlib/` in the Mathlib checkout pinned by the repository's
`lean/` project (commit `c55e6e78`, the same pin this prototype was written against).
Declaration names are fully qualified. `results/mathlib-audit.txt` records matching source
lines; `scripts/audit_mathlib.sh` reproduces the search (re-run in the cleanup with identical
output; the negative search, `results/mathlib-negative-search.txt`, has zero hits). A wider
manual search for interlacing, minimax and eigenvalue monotonicity also found nothing; the
only related extras are `Matrix.eigenvalue_mem_ball` (Gershgorin, `LinearAlgebra/Matrix/
Gershgorin.lean`) and `LinearMap.IsSymmetric.card_filter_eigenvalues_eq` (multiplicities).
Negative search is evidence of an API gap, not a logical proof that no equivalent statement
can be derived elsewhere.

| Need | Exact declarations already present | Missing work and incremental Lean LoC estimate |
|---|---|---|
| Hermitian diagonalization and eigenvectors | `Matrix.IsHermitian.spectral_theorem`, `Matrix.IsHermitian.eigenvectorBasis`, `Matrix.IsHermitian.mulVec_eigenvectorBasis`, `Matrix.IsHermitian.eigenvectorUnitary` in `Analysis/Matrix/Spectrum.lean` | Thin matrix/operator wrappers, 80–150 |
| Eigenvalue ordering | `Matrix.IsHermitian.eigenvalues₀`, `Matrix.IsHermitian.eigenvalues₀_antitone`; `LinearMap.IsSymmetric.eigenvalues`, `LinearMap.IsSymmetric.eigenvalues_antitone` in `Analysis/InnerProductSpace/Spectrum.lean` | Ascending-index wrapper and reindexing, 50–100. Existing order is **descending**. `Matrix.IsHermitian.eigenvalues` is reindexed onto the original index type and must not be treated as the ascending list. |
| Trace and spectral invariance | `Matrix.IsHermitian.trace_eq_sum_eigenvalues`, `Matrix.IsHermitian.eigenvalues_eq_eigenvalues_iff`, `Matrix.IsHermitian.roots_charpoly_eq_eigenvalues₀` | Unitary/basis transport wrappers, 100–200; repository already has useful wrappers listed below |
| Rayleigh quotient; extremal eigenvalue existence | `ContinuousLinearMap.rayleighQuotient`, `ContinuousLinearMap.rayleighQuotient_le_norm`, `ContinuousLinearMap.norm_eq_iSup_rayleighQuotient`; `LinearMap.IsSymmetric.hasEigenvalue_iInf_of_finiteDimensional`, `LinearMap.IsSymmetric.hasEigenvalue_iSup_of_finiteDimensional` in `Analysis/InnerProductSpace/Rayleigh.lean` | Ground-state Rayleigh upper-bound wrapper, 100–200. These do not give the full indexed Ritz interlacing theorem. |
| Courant–Fischer / indexed min–max | No ready-made declaration found in the spectral/inner-product/matrix API or wider text search | Dimension-intersection argument and both extrema on k-dimensional subspaces, 600–1000 |
| Weyl inequalities | No Hermitian ordered-eigenvalue perturbation/addition theorem found; Weyl-group hits are unrelated | Operator-order monotonicity and norm perturbation, 250–450 after min–max; full two-index Weyl addition is unnecessary for this task |
| Cauchy interlacing / indexed Rayleigh–Ritz | No ready-made declaration found | Orthogonal compression, ascending indices, 200–350 after min–max |
| PSD / Loewner order | `Matrix.le_iff`, `Matrix.nonneg_iff_posSemidef` with scoped `MatrixOrder` in `Analysis/Matrix/Order.lean`; `Matrix.posSemidef_iff_dotProduct_mulVec`, `Matrix.PosSemidef.dotProduct_mulVec_nonneg` in `LinearAlgebra/Matrix/PosDef.lean`; `Matrix.IsHermitian.posSemidef_iff_eigenvalues_nonneg`, `Matrix.PosSemidef.eigenvalues_nonneg` in `Analysis/Matrix/PosDef.lean` | Exact rational LDL/congruence checker still needed, 350–650; positivity of a matrix and efficient decidable certification are distinct |
| Trace monotonicity | `Matrix.PosSemidef.trace_nonneg`, `Matrix.tracePositiveLinearMap` | Trace of an **antitone scalar spectral function** for A >= B follows from ordered eigenvalues, 100–200. Not available merely from trace positivity, and false for arbitrary f. |
| Jensen | `ConvexOn.map_sum_le`, `ConvexOn.map_centerMass_le` in `Analysis/Convex/Jensen.lean`; `convexOn_exp` in `Analysis/Convex/SpecificFunctions/Basic.lean` | Apply scalar Jensen to spectral weights, normalize rows/columns and sum: Peierls 200–350 |
| Nonorthogonal Peierls / Ritz data | Above spectral, PSD, Jensen building blocks | Gram whitening or frame correction; generalized residual and small-matrix enclosure interface, 500–900 |

Relevant existing repository declarations (read-only):

- `CertifiedGaps.trace_exp_smul_eq_sum_of_isHermitian` in `lean/CertifiedGaps/Ring/Model.lean`.
- `CertifiedGaps.eigenvalues_eq_of_conj`, `CertifiedGaps.piTensor_conj`,
  `CertifiedGaps.conj_bondOp_kronecker`, `CertifiedGaps.piTensor_conjTranspose_mul_self`
  in `lean/CertifiedGaps/Ring/Conjugation.lean`.
- `SpinOneHeisenberg.hamiltonian_isHermitian`, `SpinOneHeisenberg.trace_exp_hamiltonian`,
  `SpinOneHeisenberg.partitionFunction_eq_sum` in `lean/CertifiedGaps/Heisenberg/Hamiltonian.lean`.

A focused spectral foundation is approximately **2.8–4.9k Lean lines**, including overlap
between the rows above. Do not add these audit estimates to the module budget below;
they describe the same work from an API perspective.

## 2. Executable prototype and trust boundary

`Packed.lean` implements:

- little-endian base-3 configurations in a full `3^n`-lane Nat, radix `B = 2^w`;
- digit-pattern masks via finite geometric series, avoiding per-configuration mask generation;
- standard-Sz `H v` with two nonnegative limbs representing signed lanes, including a
  pi-twisted seam;
- coordinatewise signed floor using an aligned bias, with masking after the global shift;
- butterfly reversal of vector lanes in O(log d) big operations;
- one big multiplication for each nonnegative dot product (four for a signed pair);
- Gram matrices and product-state power-filter recipes.

`Optimized.lean` hoists masks out of the filter and hoists reversal out of the k-by-k Gram
loop. It also provides reversal masks built by repeated doubling (shifts and ORs),
avoiding the expensive big divisions in the initial reversal. `FastReverse8.lean` checks
this against the original on a non-palindromic regenerated vector. General correctness
is specified by `FastReverseCorrect` and remains open. This removes measured division work.
It does not solve compact-sector addressing.
`PreparedH8.lean` kernel-checks equality of all three optimized
and original H8 moment integers.

The tested trial is

`v = (n I - H_M)^q e_Neel`,

where `H_M` is the Marshall conjugate of the untwisted even-ring Hamiltonian. This is an
exact integer power filter, not a Chebyshev or tuned thermal certificate. Set
`F = nI-H_M`, `G=vᵀv`, `A=vᵀFv`, `C=(Fv)ᵀFv`. Then the intended energy and residual are

`theta = n - A/G`, `residual² = C/G - (A/G)²`.

The integer checker proves `G>0`, `den*A >= (den*n+num)*G`,
`epsDen²*(C*G-A²) <= epsNum²*G²`, and `A²<=C*G`.
That final inequality prevents a Nat subtraction from silently truncating a negative
residual numerator. The raw exact G/A/C values are in `results/exact_witnesses.tsv`.
The Julia generator independently computes them using uncompressed sparse rows and BigInts;
its agreement is a regression check, not part of Lean's trust base.

| n | Degree q | Lane bits | Checked intended Ritz upper endpoint | Checked intended residual upper endpoint |
|---|---:|---:|---:|---:|
| 6 | 24 | 262 | -8.495447 | 0.274967 |
| 8 | 32 | 345 | -11.179504 | 0.279904 |
| 10 | 16 | 222 | -13.548710 | 0.632355 |

**These table entries are arithmetic witnesses, not proved physical eigenvalue bounds.**
Only one trial per size is tested; multiple low levels and a full spectral counting
certificate are not delivered. The residuals are far too loose for the thermal tables.

`PackingProofs.lean` proves `radixPack_head`, `radixPack_tail`, `radixPack_extract`,
`add_noCarry`, and `pack_eq_radixPack`. `Obligations.lean` gives precise propositions
`DotCorrect`, `ReverseCorrect`, `FastReverseCorrect`, `MaskCorrect`, `FloorCorrect`,
`HApplyCorrect`, and `FilterCorrect`.
These are `def ... : Prop`, never `axiom` or theorem-by-placeholder.
They are not assumed by the checkers. The bitwise/matrix interpretations remain unproved.

Important explicit bounds:

- Inputs must have no lanes above d, and each limb must have lanes at most M.
- H accumulation is safe if `3*n*M < B`; the positive filter uses `4*n*M < B`.
- Dot products require **every convolution coefficient**, not just the center, below B.
  `d*M² < B` suffices. The trial modules kernel-check
  `3^n*(4*n)^(2*(q+1)) < 2^w` as a conservative arithmetic bound.
- Signed floor requires `s <= biasBits < w`, negative lanes `<=2^biasBits`, and
  positive lanes plus the bias `< B`. The bias must be divisible by `2^s`.
- With d=531441 and lanes of size `2^55`, a dot needs at least **130 bits**, before
  allowing H-vector growth. The 128-bit dense microbenchmark uses tiny lanes and is safe;
  it is not evidence that 128 bits suffices for production Q=2^55.

`Tests.lean` checks signed H against a separate row formula for n=2,3,4 and both twists;
Marshall filtering for n=2,4; negative rounding `[3,-3,17,-17]/2 -> [1,-2,8,-9]`;
reversal; and a Gram matrix. The n=2 case checks the double bond. A deliberately undersized
radix gives a packed dot of 0 instead of the true 18, exhibiting why carry bounds matter.
`Gram6.lean` additionally checks a 4-by-4 Gram matrix on 729-lane vectors.

`Counting.lean` proves the one-Ritz-value **conditional** counting step in scaled integer
coordinates, with Weyl, separation and anchor as explicit premises. It checks a concrete
counterexample to omitting separation. It does not instantiate a physical spectral theorem.

All proof computations use `decide +kernel`. Printed axiom sets are subsets of
`propext`, `Classical.choice`, `Quot.sound`; no `sorryAx` or native-decide axiom is used.

## 3. Packing and anchor experiments

### Compact-sector addressing

Counts below are exact Julia enumeration of the Sz=0 sector, not kernel-proved counts.
The packing order is increasing base-3 index restricted to that sector. Each entry counts
**distinct** target-rank minus source-rank differences.

| n | Full dimension | Sz=0 dimension | Distinct compact shifts | Closing-bond shifts | Full vector MiB at 128 bits |
|---|---:|---:|---:|---:|---:|
| 6 | 729 | 141 | 32 | 18 | 0.0111 |
| 8 | 6561 | 1107 | 120 | 94 | 0.1001 |
| 10 | 59049 | 8953 | 614 | 574 | 0.9010 |
| 12 | 531441 | 73789 | 3804 | 3744 | 8.1091 |

At n=12 the per-bond counts are `2,4,6,8,10,12,10,8,6,4,2,3744`.
This refutes the straightforward compact-sector implementation of “about 50 masked shifts.”
It is not a lower bound on every possible sector permutation circuit: a fast cyclic rotation,
bit-permutation network, or a different representation could recover a cheap seam.
Such a construction has to be supplied and benchmarked. Full packing costs 7.202 times as
many lanes at n=12; it was chosen here because its transitions really are fixed shifts.

### Untuned cut anchor (all floating-point results NON-CERTIFYING)

The always-valid elementary comparison is

`H >= H_open(n/2) tensor I + I tensor H_open(n/2) - 4I`,

because each of the two removed exchange bonds has minimum eigenvalue -2. Removing bonds
without this compensation is not a lower comparison. The eigenvalues of a bond are
-2,-1,1, with multiplicities 1,3,5. Twisting a bond by a unitary preserves this fact.

For `a=700741/500000`, `beta=49/4`, `tau=3e-5`, use the conservative full-dimension cutoff
`rho = -n*a + log(3^n/tau)/beta`. This spends the entire tolerance on the tail;
a real error budget must reserve some for retained levels and scalar enclosures.

| n | H0 ground energy | rho | # H0 levels below rho | # below -na |
|---|---:|---:|---:|---:|
| 6 | -10.000000000000 | -7.020648170750 | 64 | 15 |
| 8 | -13.291502622129 | -9.644246899131 | 162 | 32 |
| 10 | -15.660425045542 | -12.267845627512 | 322 | 64 |
| 12 | -18.740549938849 | -14.891444355892 | 773 | 86 |

A large overcount is not merely a request for more Ritz vectors. The counted-anchor lemma
below requires that the kth Ritz value be below the anchor. We diagonalized all Sz blocks
of H6/H8 in Float64, sorted their full spectra with multiplicities, and compared every k:

| n | max over 1<=k<d of lambda_(k+1)(H0) - lambda_k(H) | Best k | Positive margins |
|---|---:|---:|---:|
| 6 | -1.382576818186 | 1 | 0 |
| 8 | -1.445376754516 | 1 | 0 |

Thus **even ideal Ritz vectors cannot satisfy the separation condition with this H0** in
these two examples. This is a robust numerical diagnostic with margins over one energy unit,
not a certified impossibility theorem. Tuned boundary terms were not optimized or excluded.
A proposed tuning should first demonstrate a positive separated margin, with its exact local
PSD certificate, before a recursive Lean implementation is budgeted.

## 4. Mathematical review of proposed steps 1–4

### Step 1: Peierls with integer trial families

For orthonormal vectors q_i, scalar spectral Jensen gives

`exp(-beta <q_i,H q_i>) <= <q_i, exp(-beta H) q_i>`.

Summing and using the orthogonal projection <= I proves the proposed lower bound. For a
nonorthogonal family the uncorrected sum is false: H=0 on C² and three copies of one unit
vector give 3 on the proposed lower side but trace 2.

Two valid repairs:

1. If V has full column rank, G=V*V>0 and Q=V G^(-1/2) is orthonormal. Diagonalize
   T=G^(-1/2)(V*HV)G^(-1/2); then `Z >= sum exp(-beta theta_i)` for its Ritz values.
   An integer family costs no Gram penalty when handled through this generalized problem.
2. Normalize columns to w_i and prove `WW* <= c I` (equivalently the nonzero spectrum of
   W*W is <=c). Then `Z >= (1/c) sum exp(-beta <w_i,H w_i>)`. The Gram here is the
   **normalized** Gram, not raw V*V. The frame penalty c can destroy the target accuracy.

A near-identity Gram estimate alone does not justify dividing signed Rayleigh energies by
one common endpoint. Use generalized PSD comparisons or a proven shift to positivity.
Rank and normalization checks are mandatory, even though the generator needs no proof.

### Step 2: the Weyl/counting lemma, and what the anchor must prove

Write the ordered physical eigenvalues **ascending** as lambda_1,...,lambda_d.
For orthonormal Q with k columns, let T=Q*HQ have Ritz eigenvalues theta_1<=...<=theta_k,
R=HQ-QT, and `||R|| <= eps`. Decompose into trial and complement, and delete the
cross block to get `B = T direct_sum C`. Then `||H-B||=||R||`.

A sufficient certificate is:

- `rho <= lambda_(k+1)(H)` with multiplicity, for k<d;
- **strict separation** `theta_k + eps < rho`;
- correctly ordered enclosures for theta_i and a spectral/operator norm residual bound.

Weyl implies `lambda_(k+1)(B) >= rho-eps`. Since all k Ritz eigenvalues are below
`rho-eps`, C cannot contribute an eigenvalue below `rho-eps`. Therefore the lowest k
values of B are exactly theta_i, and Weyl gives

`theta_i-eps <= lambda_i(H) <= theta_i`, for 1<=i<=k,

where the upper bound is Rayleigh–Ritz. Repeated eigenvalues and threshold equality must be
handled explicitly; strict separation is a convenient safe formulation.

Without separation, `H=diag(-100,0)` and Q=e₂ have theta=0, residual=0 and the valid anchor
`rho=0<=lambda_2(H)`, but lambda_1=-100. Without counting, zero residual can similarly
select any excited invariant subspace. Bounding each column residual separately is not a
bound of eps on the whole block: a Frobenius bound, or an operator-norm Gram bound, is needed.

For integer V, put `G=V*V`, `A=V*HV`, `C=(HV)^*(HV)`. The correct generalized residual test is

`C - A G^(-1) A <= eps² G`, with `G>0`.

Whitening turns this into `R*R<=eps²I`. If approximate Ritz coordinates or approximate
normalization are used, include those errors too. A small-matrix LDL/inertia certificate
can avoid computing square roots inside Lean. None of this supplies the missing anchor.
Lehmann–Goerisch is also not anchor-free; it needs a separately justified spectral count or
complement estimate. Calling homotopy a fallback does not supply a crossing-count proof.

### Tail and trace monotonicity

With lower enclosures l_i for the retained eigenvalues and
`rho<=lambda_(k+1)(H)`, beta>0,

`Z_shifted <= sum_(i<=k) exp(-beta*(l_i+n*a)) + (d-k)*exp(-beta*(rho+n*a))`.

The shift n*a is essential. Use sector dimensions and multiplicities consistently; a
sectorwise k and rho cannot be inserted into a full-space tail without summing all sectors.
The trace monotonicity requested in the task holds for antitone f, including
`f(t)=exp(-beta*t)`. It is false for unrestricted f (take f(t)=t). No operator-monotonicity
of the exponential is needed: ordered eigenvalue monotonicity suffices after taking trace.

At n=12 the crude tail threshold is -14.891444355892, or shifted energy 1.926339644108.
A ground-energy upper estimate cannot by itself replace an absolute certified rho.
Avoid counting “states near a numerical E0” unless all endpoints and multiplicities are
certified. Scalar exp values may be supplied by rational Taylor/interval certificates, but
the analytic exp enclosure theorem still has to be proved; exp does not disappear from the
mathematical trace bound.

Precision is more severe than the energy window suggests. The n=4 untwisted shifted trace
is approximately 124.885379. A uniform eigenvalue error delta changes its trace bound by
roughly beta*Z*delta, giving delta about `1.96e-8` for a 3e-5 budget, before allocating
other errors. For n=12, Z approximately 1.900547 gives about `1.29e-6`. The observed
0.27–0.63 residuals do not test the needed precision. Filter degree, k, and Gram conditioning
must be chosen from a complete error budget, not just a cutoff near E0+2.

### Local comparison and recursion

A <=81-by-81 local PSD check is useful only after an exact identity decomposes H-H0 into
embeddings of those PSD operators. Boundary terms assigned to the two halves must cancel
in that identity; overlapping supports and scalar allocations cannot be counted twice.
A half-chain eigenvalue list must include certified multiplicities, all relevant Sz sectors,
and a certified remainder count before pair-sum ordering gives an anchor. Recursive anchors
can accumulate losses and rank inflation. The measured untuned comparison is already too
weak for the required separation at n=6 and 8.

A cheaper development order is to search for a usable tuned comparison numerically, emit
small rational local PSD certificates, and test the separation inequality first. If no such
comparison exists, investigate a verified inertia/Schur-complement scheme using the low
subspace of H0, a complement coercivity bound, or a genuinely counted homotopy. These are
new routes requiring their own cost measurements, not implemented fallbacks here.

### Step 3: untrusted regeneration

It is correct that V may be generated by any terminating recipe, including arbitrary
rounding: Rayleigh–Ritz does not require proving the recipe approximates eigenvectors.
The checker must still establish G>0, valid arithmetic/carry bounds for its interpretation
of V/HV, its Gram/residual inequalities, and the anchor. Floor errors need not be bounded
for V generation, but cannot be ignored if HV or the reported Gram entries are rounded.
This distinction can simplify the generator substantially.

### Step 4 and connection to the existing Hamiltonian

The existing `SpinOneHeisenberg.hamiltonian` is Cartesian, complex and **untwisted**. A
base-3 matrix is not definitionally this matrix. Required bridges:

1. Let m_j=s_j-1 and identify index `sum s_j*3^j` with `Fin n -> Fin 3`.
2. Check a local Cartesian-to-weight unitary. Compatible columns for m=-1,0,+1 are
   `(1,-i,0)/sqrt(2)`, `(0,0,1)`, `(-1,-i,0)/sqrt(2)`. With the repository's
   `S_alpha=-i epsilon_(alpha s s')`, these give S^z=diag(-1,0,1) and positive sqrt(2)
   ladder entries. Their product cancels the 1/2 in the exchange to give integer hops +1.
3. Lift this identity sitewise using the existing tensor conjugation lemmas. Restrict to
   n>=2; preserve the n=2 double bond. A blanket n=1 formula is unsafe.
4. Define a twisted Hamiltonian: conjugate one endpoint of the closing exchange by
   `diag(-1,1,-1)`. Only its transverse hops change sign. Prove trace/spectrum invariance
   under the global rotation identifying the paper's P with this z-axis pi twist.
5. Prove Sz conservation, direct-sum decomposition, and multiplicity-preserving spectral
   and trace recombination. Sz=0 alone is not enough for a full thermal trace. Avoiding
   SU(2) theory is possible, but then certify every Sz sector.
6. The positive Marshall filter is valid only for **even untwisted** rings. Source 1
   includes n=5,7,9,11 and pi twists; signed generation or a frustrated seam is needed.
   Source 2's third twist has complex phases; an integer pair representation in 1,omega
   is possible, but its conjugation, Gram and norm form `a²-ab+b²` require extra proofs.

The hApply tests exercise twisted signs in the weight representation. They do not prove
any of these global bridges. The paper's symmetry-resolved spatial-transfer I/O/N moments
are also additional structure beyond the already-proved two-purity transfer inequality.

## 5. Benchmark table and extrapolation

Two serial runs, one fresh Lean process per check: the original PR run (Linux container, Intel
Xeon Platinum 8573C, 16 GiB, `results/benchmarks_original_linux.tsv`) and the cleanup re-run
(Apple M5, 32 GiB, `results/benchmarks.tsv`, `scripts/benchmark.jl`). Seconds / peak RSS GiB.

| Check | file | Linux s | Linux GiB | M5 s | M5 GiB |
|---|---|---:|---:|---:|---:|
| Imports / trivial check | Baseline | 4.205 | 0.798 | 3.24 | 0.60 |
| n=6 dense matvec + dot | Bench6 | 2.265 | 0.800 | 0.52 | 0.61 |
| n=8 dense matvec + dot | Bench8 | 2.500 | 0.800 | 0.52 | 0.64 |
| n=10 dense matvec + dot | Bench10 | 6.457 | 1.118 | 1.33 | 1.03 |
| n=12 dense matvec + dot | Bench12 | 72.996 | 5.548 | 26.38 | 3.70 |
| H6, degree 24, 3 moments | H6 | 3.596 | 0.799 | — | — |
| H8, degree 32, 3 moments | H8 | 6.584 | 1.557 | — | — |
| H10, degree 16, 3 moments | H10 | 21.589 | 4.045 | — | — |
| n=12 matvec + checksum | Matvec12 | 27.960 | 3.761 | 11.34 | 2.51 |
| n=12 reversal + dot | Dot12 | 56.997 | 2.617 | 20.64 | 2.28 |
| n=12 prepared masks, matvec + dot | Prepared12 | 60.279 | 4.694 | 28.34 | 3.41 |
| H8 optimized/original equality | PreparedH8 | 8.111 | 2.352 | 1.43 | 2.29 |
| n=6, four-vector Gram | Gram6 | 2.568 | 0.800 | 1.36 | 0.61 |
| Signed/twisted/floor/Gram cases | Tests | 3.852 | 0.810 | — | — |
| n=12 dot, reversed operand supplied | PlainDot12 | 5.307 | 0.823 | 2.46 | 0.77 |
| n=12 dot with division-free reversal | FastDot12 | 9.957 | 2.856 | 2.96 | 2.57 |
| n=12 prepared matvec + fast dot | Fast12 | 30.369 | 5.105 | 9.63 | 3.56 |
| n=8 fast/original reversal equality | FastReverse8 | 2.785 | 0.800 | 1.24 | 0.67 |
| n=12 prepared matvec + checksum | PreparedMatvec12 | 18.042 | 2.961 | 6.98 | 2.99 |

H6, H8, H10 and Tests are now library modules built by `lake build`; on the M5 a clean build
of the whole library (all of them) takes 8.5 s, with H10 the largest module at 3.2 s.

The dense microbenchmark computes one filter matvec on a full repunit vector, reverses a
vector, multiplies, and checks the closed sum `17*n*3^(n-2)`. Every theorem is elaborated
and kernel-checked afresh. The Hn tests regenerate the entire trial and three dot products.
Times include Lean/Lake startup and imports, elaboration, proof construction and checking;
there is no claim that the table isolates kernel normalization time. RSS is the maximum
reported by Linux child rusage; startup contributes a large floor. Small-n timing variance
is dominated by startup. The final serial run avoids concurrent eigensolvers or bignum jobs.

One 4-million-bit multiplication does not characterize a whole matvec. At n=12 this
prototype uses 68,024,448-bit full vectors at 128 bits per lane, repeated divisions for
masks, many shifts/ands, and reversal. Peak resident memory is orders of magnitude larger
than one packed vector. Hoisting masks and reversals is necessary. The optimized
n=12 dense matvec-plus-dot takes **30.369 s** versus **72.996 s** initially
(2.40x faster), while peak RSS is **5.105 GiB** versus **5.548 GiB**. A dot with a
pre-reversed operand takes **5.307 s**; the division-free reversal-plus-dot takes **9.957 s**,
versus **56.997 s** for the first reversal. Hoisting reversed vectors out of a Gram loop
is consequently a material optimization, not a cosmetic one.

Single-check peak memory
cannot be extrapolated linearly to thousands of theorems; checkpointed modules may be needed.

For illustrative full-packing workloads with k=50 and q=200 (Linux timings; the M5 re-run
is 2.8–3.2x faster at n=12, which divides the hours below but does not change the conclusion):

- `q*k = 10000` matvecs. At roughly 30 seconds per n=12 matvec, this is about **83 hours**
  for one n/twist, before Gram work. This is a conditional workload calculation, not a
  measured production runtime; q and k are not established by these experiments.
- The prepared-mask matvec plus checksum measures **18.042 seconds**, including mask
  preparation. Repeating it 10000 times would take **50.1 hours**. Reusing masks across
  all filter iterations can reduce that cost, but was not measured at n=12 over q=200.
- Three symmetric Gram families need up to `3*k*(k+1)/2 = 3825` dot products. At the
  observed 5.307 seconds per dot with an already-reversed operand, a naive individual-check
  extrapolation is **5.64 hours**, before reversal/preparation. Startup amortization and
  shared expressions can change this; it is not a measured batched Gram runtime.
- Using the optimized full matvec-plus-dot as a 30.369-second budget gives **84.4 hours**
  for 10000 such operations. A complete filter should reuse prepared masks and omit dots
  at intermediate steps, so this is not a lower bound. There is no measured 200-step,
  50-vector n=12 kernel run proving that it fits memory or meets this estimate.
- Compact sector packing, reusable masks and reversals, narrower certified lanes, and
  carefully selected k/q could reduce this substantially. The sector seam circuit and
  anchor are currently unknown, so a defensible optimized total is **unknown**, not a
  confident time interval. Redesigning n_max=10 or 8 is strongly favored by the measured
  size growth, but must revalidate the full analytic seed and precision budget.

## 6. Module-by-module remaining implementation budget

These are engineering estimates for reviewed Lean source, not measured coding throughput.
Ranges account for missing reusable APIs and proof representation choices. They exclude
this small prototype and do not include a general-purpose full two-index Weyl theorem.

| Proposed module | New Lean LoC | Why it is needed |
|---|---:|---|
| PackedRadix / Carry / SignedFloor | 350–650 | Connect proved reference extraction to shift/mask operations, convolution carry bounds |
| Masks / Reverse / Gram | 350–650 | Geometric masks, reversal, signed products, reuse and Gram interpretation |
| SpinOnePacked / SectorDecomposition | 600–1000 | H action, seams, bounded accumulation, all Sz sectors |
| CartesianWeightBridge | 400–700 | Local sqrt(2) unitary and global indexed matrix transport |
| OrderedMinMax | 700–1100 | Ascending eigenvalues, subspace dimension arguments, Ritz |
| EigenvalueOrder / Weyl / Interlacing | 450–750 | Loewner monotonicity, operator perturbation, compression |
| RationalPSD / GeneralizedRitz | 500–850 | LDL/congruence certificates, G>0, rational spectral intervals |
| ResidualCounting | 300–500 | Generalized block residual and separated counted-anchor theorem |
| Peierls / TraceTail / ScalarExp | 500–900 | Nonorthogonal treatment, finite trace/tail and exp enclosures |
| LocalAnchor / RecursiveCount | 900–1800 | Local comparison identity, pair sums, all-sector recursive certification |
| Recipe / CertificateInterface | 250–450 | Bound checks, failure cases, replay and theorem-facing interface |
| ThermalTableAssembly | 300–600 | 18 enclosures, shifts, multiplicities, error budgeting |
| **Thermal route subtotal** | **5600–9950** | Does not solve the numerical anchor-design problem |
| PeriodicMPS / ExactTransfer | 1000–1800 | Norm and insertion contractions, integer matrix/block checker |
| SymmetryTransfer / CappedMass / DominantIsolation | 1600–3000 | I/O/N moment identities and finite-seed analytic bridge |
| SeedAssembly / ExistingChainGlue | 250–450 | Feed (36,49/4) seed into existing doubling chain |
| **Whole seed-route total** | **8450–15200** | Source-1 scope; extra third-twist work for source 2 |

No 1.5–3-week completion promise is justified. Conditional on finding a usable anchor and
representation, a planning range is **8–16 expert weeks** for this broader proof work,
plus potentially unbounded mathematical/certificate search. The estimate corresponds to
roughly 200 reviewed Lean lines per working day across the 8450–15200-line range;
that throughput is an explicit planning assumption, not a measured productivity result. This is
a staffing estimate, not evidence that an anchor exists.
A first milestone should be a fully bridged H6
counted enclosure, before expanding to H8 or implementing all transfer/trial modules.

## 7. Plan changes and open risks

1. Put the anchor/separation test first. The untuned comparator fails before precision
   becomes relevant. Require a concrete tuned local certificate and counted margin.
2. Replace the assumed compact O(50)-shift matvec with an actual seam circuit, or budget
   the measured full-base-3 costs. Cache mask construction and vector reversal.
3. Use generalized Ritz eigenvalues for the lower trace; avoid an uncontrolled frame
   penalty. Budget residual precision against beta times the shifted trace.
4. Prove a reusable ordered min–max/Weyl layer and the Cartesian/twist bridge. Do not
   present arithmetic witnesses as physical eigenvalue theorems.
5. Make the first complete acceptance test H6: all sectors, a counted lower enclosure,
   a Ritz upper enclosure, a shifted trace interval, and `#print axioms` on that theorem.
6. Explore n_max<=10 or source 2 only after charging its complex twist, trial length,
   analytic transfer identities and precision. The existing optimized scalar note is
   useful evidence for redesign, not a Lean-certified shortcut by itself.

Remaining risks: tuned anchor existence, rank inflation, norm-vs-column residual errors,
Gram conditioning, odd rings and twists, degeneracies and threshold equalities, all-sector
multiplicities, carry/borrow proofs, expression-retention memory, Mathlib API integration,
scalar exp bounds, and the unimplemented analytic/MPS seed bridge. The project deliberately
keeps all of these visible; none is hidden behind an extra axiom or an unchecked generator.

## 8. Validation record

Re-run in the cleanup (Apple M5, Julia 1.12, Lean v4.35.0-rc3):

- `lake build` in `src/` succeeds from clean in 8.5 s with no Mathlib dependency (the PR's
  declared but unused Mathlib requirement was dropped; nothing imports it).
- `#print axioms` for all 11 audited declarations (`src/AxiomAudit.lean`, also
  `results/axioms.txt`) shows only subsets of `propext`, `Classical.choice`, `Quot.sound`;
  no `sorryAx`, no `Lean.ofReduceBool` (no `native_decide`).
- Every `src/Bench` file checks (exit 0), the slowest being Bench12 (26 s) and Prepared12
  (28 s); none needed the 5-minute limit.
- `scripts/generate.jl` regenerates `src/SeedPrototype/H{6,8,10}.lean` and
  `results/exact_witnesses.tsv` byte-identically (105 s). The witnesses are tight: raising
  the H6 Ritz numerator by 1 (8495447 to 8495448) makes `decide +kernel` report the
  proposition false.
- `scripts/experiments.jl` reproduces the structure and anchor tables exactly and the
  Float64 level table to about 1e-13.
- `test/runtests.jl` (24 exact tests, 2 s) checks sparse-row/packed agreement, digit masks,
  reversal and dots for n = 2..6, the carry-overflow counterexample, and the H6 recipe.
- `check.jl` verifies the source hashes in `results/witnesses.json`, regenerates the three
  witnesses and their Lean files, rebuilds the library from clean, and checks the axiom
  reports; `--full` adds the axiom audit and every benchmark file.

The float eigenvalue/anchor tables remain non-certifying.
