# Dimerized Heisenberg: verified gaps under general reflected RP

**The physical model satisfies the repository's general reflection-positivity
hypothesis.** What fails is only the original task's stronger requirement of
real channels with literally identical factors `C tensor C` in a common cell
basis. An independent mathematical review and exact Julia checks confirm an
eligible **internally reflected** formulation `C tensor Theta(C)`, with real
channels and `Theta = F_cell K`. No change to the Hamiltonian or core is needed.

Uniform-gap certificates are now produced for delta=0.5, 0.3, and 0.25 under
that **general reflected RP theorem**, using complete dense spectra with
verified eigenvalue and purity enclosures. The original stronger real
identical-factor condition remains false and is explicitly marked false in
each certificate. These are not certificates of that stronger hypothesis.
The enclosure and asymmetric recurrence proof is in
[ENCLOSURE_PROOF.md](ENCLOSURE_PROOF.md).

## Verified uniform bounds (general reflected RP)

All three use n0=4 cells and apply to **every even L>=4 cells**, equivalently
8,12,16,... spins. The exact binary rational value of each returned `.bound`
is stored in its JSON; the table shows the round-trip Float64 output.

| delta | beta0 | `certified_gap_bound(...).bound` | reference 16-spin ED gap | bound/reference | verified run seconds |
|---|---|---|---|---|---|
| 1/2 | 6 | 0.5326577398571627 | 1.188651652 | 0.448119 | 47.38 |
| 3/10 | 13/2 | 0.17291646480187198 | 0.811144396 | 0.213176 | 43.23 |
| 1/4 | 29/4 | 0.057966591686137485 | 0.709347628 | 0.081718 | 43.94 |

The exact outward rational seeds are:

| delta | epsilon_s upper | epsilon_t upper |
|---|---|---|
| 1/2 | 50996126991899/4503599627370496 | 13777527984355/1125899906842624 |
| 3/10 | 857117450063167/9007199254740992 | 466505354307053/9007199254740992 |
| 1/4 | 373168159679129/2251799813685248 | 115354947892529/2251799813685248 |

Thus 0.25 is the smallest delta certified **under general reflected RP on the
tested even-seed grid**. The original stronger identical-factor assumption
is not satisfied at any of these deltas. The unsuccessful smaller deltas and
their heuristic extrapolations remain numerical only.

`verified_spectra.jl` constructs exact-algebraic block intervals, diagonalizes
every block densely, and implements a rigorous norm upper bound for residual
and orthogonality errors. The polar/Weyl enclosure, outward partition sums,
and exact rational seed conversion are detailed in ENCLOSURE_PROOF.md.
The computation uses classical binary64 OpenBLAS GEMM with an explicit
rounding/underflow bound, not an assumption that LAPACK eigenvalues are exact.
The largest 16-spin eigenvalue enclosure radii are, respectively,
2.724e-5, 1.821e-5, and 1.715e-5 (rounded upward here for display). Each
sector has its own exact rational radius; the 8-spin radii are below 4.1e-11.
The reference gaps and bound/reference ratios remain numerical comparison
data, not claims about the infinite-chain gap.

## Model, cell orientation, and exact reflection

Take cells `(2j,2j+1)` with the strong bond inside each cell. Translating the
physical ring by one spin exchanges origin conventions without changing its
partition functions. Write `s=1+delta`, `w=1-delta`,
`A_a=S_a tensor I`, and `B_a=I tensor S_a`, where `S_a` are Pauli matrices/2.
Before rotation,

```
D = s sum_a A_a B_a,       h = w sum_a B_a tensor A_a.
```

Rotate the second spin of every cell by `J=[0 1;-1 0]`, so `R=I tensor J`.
In the rotated product basis,

```
Drot = s (-A_x B_x + A_y B_y - A_z B_z),
-hrot = w (B_x tensor A_x - B_y tensor A_y + B_z tensor A_z).
```

Let `F_cell` exchange the two spins **inside** a cell. This is distinct from
`F_outer`, which exchanges two complete cells. Define

```
Theta = F_cell K,
(C_x,C_y,C_z) = (B_x, i B_y, B_z).
```

All `C_a` and `Drot` are real; `Drot` is symmetric. Exactly,

```
Theta^2 = 1,   Theta(Drot) = Drot,
Theta(C_x,C_y,C_z) = (A_x, i A_y, A_z),
-hrot = w sum_a C_a tensor Theta(C_a).
```

`rp.jl` verifies these identities with `Complex{Rational{BigInt}}` at all six
parameters. The two-site cell is palindromic (one internal bond). This is
precisely the general antiunitary hypothesis in Appendix E, Eq. (rpform), of
`outputs/thermal-mpo-gap-scope.tex`. Its mirror lemma permits `Theta=VK` and
`Theta^2=+/-1`; it does not require `Theta(C)=C`. Here the transfer kernel is
also real: the channels, their reflected counterparts, and the on-site
propagator are real. Internal reflection exchanges left/right insertions and
leaves the on-site term invariant, making that real kernel symmetric under
the repository proof's conventions.

### Singlet-i basis: the equivalent conjugate form

In product basis `(uu,ud,du,dd)`, let

```
M = [1  1  0   0;
     0  0  1   i;
     0  0  1  -i;
     1 -1  0   0],     W=M/sqrt(2),     U=R W.
```

The first three columns of `W` are real triplets; the fourth is the singlet
multiplied by `i`. `M' M=2I` exactly. `W` converts `Theta` to ordinary complex
conjugation. In the original physical basis transformed by `U`,

```
U' D U = s diag(-3/4,1/4,1/4,1/4),
C_x = (1/2) [0 0 1  0; 0 0 0 i; 1 0 0 0; 0 -i 0 0],
C_y = (1/2) [0 0 0  i; 0 0 1 0; 0 -1 0 0; i 0 0 0],
C_z = (1/2) [0 1 0  0; 1 0 0 0; 0 0 0 -i; 0 0 i 0],
-h' = w sum_a C_a tensor conj(C_a).
```

These channels are **complex**, not real. All conjugations are evaluated with
`M` and rational denominators, so the checks use no approximation to `sqrt(2)`.
Both exact representations satisfy the broader repository RP condition.

## Precise scope of the identical-factor obstruction

Any `sum c_k C_k tensor C_k` commutes with `F_outer`, regardless of positivity
or reality of the matrices. But for the oriented physical cell bond,

```
||h - F_outer h F_outer||_F^2 = 6 (1-delta)^2 > 0.
```

This follows from `Tr(A_a' A_b)=Tr(B_a' B_b)=delta_ab` and
`Tr(A_a' B_b)=0`; the six tensor products in the difference are orthogonal.
It is also checked by exact construction of the 16 by 16 matrix. A uniform
basis change `U tensor U` commutes with `F_outer` and preserves the norm.
Thus the obstruction excludes the **uniform literal identical-factor form**,
not an internally reflected or complex-conjugate form. Scalar shifts cannot
fix it; redistributing on-site terms cannot cancel the two-cell traceless
component either.

| delta | w | exact squared outer-swap defect |
|---|---|---|
| 1/2 | 1/2 | 3/2 |
| 3/10 | 7/10 | 147/50 |
| 1/4 | 3/4 | 27/8 |
| 1/5 | 4/5 | 96/25 |
| 3/20 | 17/20 | 867/200 |
| 1/10 | 9/10 | 243/50 |

### Other orientations and larger blocks

`reflection_conventions.jl` checks alternatives exactly:

- Reversing every cell or shifting the origin exchanges endpoint labels or
  strong/weak boundary weights; it does not erase a nonzero outer-swap defect.
- Reversing alternating cells produces alternating `B tensor B` and
  `A tensor A` bonds. Each has real identical channels, but these are **two
  different bond operators**. This does not prove the homogeneous identical-
  channel assumption. We do not assume an alternating product of symmetric
  transfer operators is itself symmetric.
- For a contiguous block of `c>=2` spins, dimension `d=2^c`, the boundary
  endpoint operators give squared outer-swap defect `3 d^2 w^2 / 8`.
  Larger uniform blocking cannot restore the literal identical-factor form.
- Even blocks retain a palindromic alternating internal Hamiltonian. Internal
  reversal composed with conjugation supplies the general reflected RP form.
  The code checks two- and four-spin blocks, including on-site reflection
  invariance at all six deltas. The general defect formula follows from the
  same endpoint-support trace argument. Arbitrary nonlocal encodings or
  unrelated transfer gauges have not been exhaustively excluded.

An independent subagent read the repository criterion and rederived the
`Theta=F_cell K` identity separately, checking it in exact Julia arithmetic.
Its conclusion agrees: general repository RP is available; the stronger
original task requirement is incompatible with this common-cell bond.

## Dense numerical spectra and sizes reached

All following spectral and purity values are **numerical, not certified**.
`numerical_ed.jl` constructs every Sz / two-spin-translation momentum /
spin-flip block and uses dense LAPACK eigenvalues only. Nonzero Sz sectors
are counted with multiplicity two using exact physical spin-flip symmetry;
zero Sz is split into spin-flip parity sectors. Group stabilizer compatibility
is checked with integer modular arithmetic. The script asserts each fixed-Sz
dimension sum and the complete weighted sum `2^L`. Complete dimension does
not make Float64 eigenvalues rigorous enclosures.

| spins | cells | full dimension | dense blocks solved | largest block |
|---|---|---|---|---|
| 4 | 2 | 16 | 6 | 3 |
| 6 | 3 | 64 | 13 | 5 |
| 8 | 4 | 256 | 21 | 14 |
| 12 | 6 | 4096 | 43 | 132 |
| 16 | 8 | 65536 | 73 | 1430 |

Every sector's labels, dimension, multiplicity, Hermiticity discrepancy and
runtime are stored in the per-parameter JSON files. Independent comparison
against the core's dense charge-sector spectra passes at 4, 6, 8 and 10 spins;
trace and squared-trace identities also pass. No Lanczos/ARPACK was used.

Seeds are scanned at `n0=2,3,4` cells and `beta0=0.25:0.25:32` using core
`scan_seeds` and `best_seed`. Odd n0=3 is diagnostic only, excluded from
criterion-eligible best seeds. The largest ring is 16 spins, not 20. A seed
at n0 cells needs 4 n0 spins; the original 20-spin limit permits n0<=5, and
only n0=2,4 are even. Thus the completed scan includes every even seed size
within a literal 20-spin ceiling, although spectra at 18/20 spins and the
odd n0=5 diagnostic have not been computed. If "up to 10 cells" instead means
n0=10, the largest required ring is 40 spins, not 20.

### Best even seeds and reference-gap ratios

The recurrence outputs below use floating-point seeds. Calling the outward-
rounded recurrence does **not** certify their input spectra or purities.
Reference gaps are ED gaps of the 16-spin ring, not extrapolated infinite-
volume gaps or certified reference bounds.

| delta | best (n0 cells,beta0) | numerical recurrence output | 16-spin ED gap | output/gap | spectral runtime s |
|---|---|---|---|---|---|
| 0.50 | (4,6.00) | 0.534795335 | 1.188651652 | 0.449918 | 25.99 |
| 0.30 | (4,6.50) | 0.173204312 | 0.811144396 | 0.213531 | 20.06 |
| 0.25 | (4,7.25) | 0.058232184 | 0.709347628 | 0.082093 | 20.00 |
| 0.20 | none | — | 0.604614582 | — | 18.94 |
| 0.15 | none | — | 0.498326272 | — | 19.23 |
| 0.10 | none | — | 0.395028708 | — | 20.69 |

Smallest successful **numerical** delta is 0.25 on this tested grid. Every
successful even seed first appears at n0=4 (8 spins for the smaller ring,
16 spins for the doubled ring). This is not a smallest certified delta.
The unequal-defect recurrence can accept delta=0.25 even though its optimized
minimum purity below is less than the simple equal-defect sufficient cutoff.
The supplied survey values and the current core recurrence agree in scale;
we report actual current outputs rather than replacing them with the survey.

### Purity-grid maxima and ground-energy differences

Each purity entry is `max_beta min(Ps(n,beta),Pt(2n,beta))`; each energy entry
is `D(n)=E0(2n)-2E0(n)` with n in cells.

| delta | purity n=2 | purity n=3 | purity n=4 | D(2) | D(3) | D(4) |
|---|---|---|---|---|---|---|
| 0.50 | 0.793458 | 0.946829 | 0.988292 | 0.0765253 | 0.0115238 | 0.00186628 |
| 0.30 | 0.597391 | 0.794418 | 0.912550 | 0.181631 | 0.0519727 | 0.0155181 |
| 0.25 | 0.535246 | 0.734203 | 0.864767 | 0.216162 | 0.0726371 | 0.0253836 |
| 0.20 | 0.488038 | 0.660558 | 0.789594 | 0.252362 | 0.0996171 | 0.0410052 |
| 0.15 | 0.454453 | 0.579736 | 0.701708 | 0.287892 | 0.133086 | 0.0650049 |
| 0.10 | 0.429492 | 0.492905 | 0.583603 | 0.319045 | 0.170599 | 0.0992359 |

### Heuristic extrapolation and delta^(-2/3) comparison

`summarize_numerics.jl` fits `D(n)=A exp(-n/ell)` to n=3,4. This is a
short-size effective energy-defect decay scale, **not a measured correlation
length**. It fits an independent-excitation prefactor `a` to the actual
16-spin thermal purity at beta=8, using the 16-spin gap `g`. It then models

```
Ps(n,beta) ~= exp(-beta A exp(-n/ell)),
Pt(2n,beta) ~= [(1+a exp(-2 beta g))/(1+a exp(-beta g))^2]^(2n).
```

The model assumes extensive excitation contributions cancel from spatial
purity, retains a finite-ring gap as a proxy, and neglects interactions and
finite-size thermal corrections. Scanning it over even n<=40 and beta<=64
gives the following **unvalidated predictions**, not guaranteed seed scales.

| delta | ell in cells | ell * delta^(2/3) | first predicted even seed (n0,beta0) |
|---|---|---|---|
| 0.50 | 0.54931 | 0.34604 | (4,5.75) |
| 0.30 | 0.82733 | 0.37076 | (4,6.25) |
| 0.25 | 0.95114 | 0.37746 | (4,7.00) |
| 0.20 | 1.12659 | 0.38529 | (6,9.25) |
| 0.15 | 1.39561 | 0.39399 | (6,10.50) |
| 0.10 | 1.84565 | 0.39763 | (8,14.00) |

The roughly stable scaled decay estimate is qualitatively compatible with
delta^(-2/3), but neither validates that asymptotic law nor identifies ell
with xi. The observed successful n0 values alone are all 4 and cannot fit a
scaling exponent. Predicted n0=6,8 require 24,32-spin doubled rings, beyond the
20-spin ceiling. No extrapolation repairs the literal identical-factor
requirement; these predictions concern numerical seeds for the general RP
formulation only.

## Rigor, tests, environment, and scope

**Exactly checked:** all displayed local reflected/conjugate identities,
real symmetric on-site terms, the narrower obstruction, and explicit
alternating/two-/four-spin-block convention checks. Exact rational audit
values are strings in `results/rp_audit.json`.

**Verified under general reflected RP:** complete spectral enclosures for the
8- and 16-spin seed rings at delta=0.5, 0.3, 0.25; outward partition sums;
exact rational purity-defect upper bounds; and the uniform recurrence bounds
listed below. Their validity range is **even L>=4 cells**, i.e. physical
site counts **8,12,16,...**. Uniqueness follows from the same purity argument.

**Numerical only:** the exploratory `numerical_*.json` files, reference gaps,
ratios using those numerical reference gaps, and extrapolations. The original
real identical-factor condition remains exactly disproved. The general
certificate does not silently assert it or a different site-length range.

All JSON files were written with the core `Results.write_result` (then named
`ThermalMPOGap`) and contain source hashes and metadata. `numerical_summary.json`
also hashes its input results. BLAS threads are one. All six full-spectrum runs
took about 19–26 seconds each after startup; the JSON files record runtimes.

The earlier 127-test audit is superseded by the expanded enclosure tests
reported below. The fast checker reconstructs scalar spectral error witnesses,
all three partition-function intervals and rational seeds from stored
spectral centers, then recomputes `certified_gap_bound`. Independent tests
also evaluate small-block residuals in exact rational arithmetic and compare
with the complete analytic four-spin spectrum. Reproducing the eigenvector
solves themselves uses `run_verified.jl`; no eigenvectors are needed by the
fast checker.

### Final validation run

The complete `test/runtests.jl` passed **2425/2425 tests** in about a minute
wall time (current `CertifiedGaps`, Apple M-series, one BLAS thread): 106 exact audit/provenance checks, 20 numerical sector cross-checks, 2
alternative-convention checks (plus their internal exact assertions), 108
independent exact-residual/analytic-spectrum checks, and 2189 stored
spectral-witness/partition/seed/recurrence/summary checks. It reproduced all three
returned bounds exactly.

Each verified result records its source SHA-256 values, package versions,
loaded BLAS configuration, all sector dimensions and multiplicities,
per-block residual and orthogonality witnesses, exact eigenvalue centers and
radii, exact log-partition endpoints, rational seeds, and the exact binary
rational representation of `.bound`. The main mathematical proof is
`ENCLOSURE_PROOF.md`, also hashed in each verified result.

There is no remaining enclosure or general-theorem blocker for these three
parameters. The outstanding qualification is the original stronger
identical-factor requirement, which is false and explicitly not claimed.
Smaller-delta extrapolations remain unverified; no enlarged seed grid or
unsupported 40-spin calculation is claimed. No file outside this directory is
modified.

### Independent acceptance

The independent review passed **4012/4012 checks in 145.5 seconds**. It reran
all six dense spectra, proved fresh eigenvalue intervals contained in every
stored interval, and checked exact rational residual/Gram norms in all 63
eight-spin sectors. Some delta=0.5 eigensolver results differed bitwise;
interval containment resolves this rigorously. See
[`review/independent_review.md`](review/independent_review.md); its 10 MB output
`review/acceptance.json` is not committed and is regenerated by
`review/independent_review.jl`. The verdict is PASS for general reflected RP
and FAIL for the original strict identical-real-factor hypothesis. The
classical binary64 OpenBLAS arithmetic model remains an explicit trusted
assumption, not a formally verified BLAS implementation.

### Provenance after the core rename and the move to `models/`

The stored results were produced under the core's former name `ThermalMPOGap`
(core commit `ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b`) with sources under
`certificates/dimerized_heisenberg/` (commit
`fd3c3d5592a91c0d37f9f4dfb85048ca5d9cfba2`). Their `source_sha256` records are
left untouched. The checker compares them with the git objects of those two
commits (`provenance.jl`) and recomputes everything else with the current
`CertifiedGaps`: all scalar witnesses, partition sums, rational seeds and
the three bounds agree exactly.
