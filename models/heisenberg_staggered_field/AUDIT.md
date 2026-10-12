# Independent certification audit

**PASS: accept the five conservative bounds below at h = 1, 3/4, 1/2, 7/20,
and 1/4, computed from complete 8- and 16-site spectra.** Their uniform gap bounds apply
to every even L >= 8 one-spin cells/sites. This acceptance does not extend to
the numerical twenty-site search candidates.

The separate reviewer inspected AGENTS.md, REPORT.md, common.jl, certify.jl,
test/runtests.jl, the core certified recurrence, and the scope note's transfer
operator and two-purity criterion. No defect requiring a mathematical source
or certificate-data change was found.

## Mathematical and numerical review

- **RP and physical normalization:** the exact rational matrices X, J=i Sy,
  and Z are spin-1/2 matrices. The real orthogonal rotation U=2J on odd sites
  transforms either bond orientation to -X tensor X-J tensor J-Z tensor Z.
  Thus the negative rotated bond is a sum of three real tensor squares with
  exact coefficient 1. The antisymmetry of J is allowed. U Z U^T=-Z converts
  the staggered field to the uniform real symmetric term -h Z. The L=2
  implementation correctly includes both periodic bonds.
- **Length and uniqueness:** on an even ring the alternating rotation closes
  consistently and gives a translation-invariant one-spin-cell Hamiltonian.
  Consequently the even n0=8 seed has the claimed even-L>=8-site scope.
  The scope note's criterion derives uniqueness along with the gap; numerical
  evidence of uniqueness is not substituted for this conclusion.
- **Complete sectors:** the translation-by-two orbit compatibility condition,
  normalized momentum matrix-element formula, p versus L-p equivalence under
  spin flip followed by one-site translation, and k versus -k equivalence under
  complex conjugation are consistent. Their combined multiplicities account
  for every state. Exact integer weighted dimension checks give 2^L, and every
  inequivalent block in the certificate path is diagonalized densely.
- **Exact matrix enclosures:** interval evaluation includes rational entries,
  orbit-normalization square roots, and trigonometric phases with an interval
  enclosure of pi. Upper-triangle enclosures and Hermitian completion correctly
  bound the submitted matrix's distance from the exact block in Frobenius norm.
- **Residual and orthogonality rounding:** the implemented upward Frobenius
  reductions and gamma_(32d+32) allowances give valid norm bounds. A complex-dot
  component involves at most 2d real products; component norm combination,
  sparse accumulation, diagonal multiplication and subtraction fit comfortably
  within this conservative budget. The sqrt(d) contribution explicitly covers
  subtraction of the identity in the orthogonality defect. At the largest
  certified dimension d=1620, even O(d^3) elementary underflow contributions are
  far below the added absolute allowance 10^-280.
- **Eigenvalue theorem:** e<1 permits polar orthogonalization of Q, with
  ||Q-V|| <= e/(1+sqrt(1-e)). The residual, polar correction, and exact-entry
  error therefore bound ||H-V Lambda V*||. Hermitian Weyl inequalities then
  enclose every ordered eigenvalue, including multiplicities and degeneracies.
- **Partition functions and seeds:** energy shifts are exact stored binary
  numbers; exponentiation and summation include every level and its exact
  symmetry multiplicity. Both log-partition endpoints are retained. Subtracting
  denominator log intervals gives the correct lower purity bounds; the upper
  endpoint of each defect is converted exactly to a rational binary fraction.
  Rational seeds and beta are passed to the outward-rounded recurrence with
  its tail correction. No decimal truncation is used to obtain the seeds.

The floating-point rounding proof assumes IEEE binary64 round-to-nearest
arithmetic, ordinary correctly rounded scalar operations and square root,
standard finite BLAS arithmetic without unsafe fast-math transformations, and
the documented outward arithmetic of IntervalArithmetic. It does not accept
the eigensolver output as exact or assume that residual accumulation in ordinary
precision is intrinsically rigorous: explicit norm-error allowances are used.

## Independent execution and supporting checks

The reviewer reran the saved checker without overwriting artifacts:

- `check.jl`: **71/71 tests passed**, approximately 21 seconds of reported test
  execution. These cover RP, comparison with independent small-ring ED, sector
  counts, stored enclosure propagation, exact rational seeds, and recurrence.
- Independently recomputed and verified **235 stored source-hash entries**.
- Verified Julia **1.12.4** and checked that `iv(pi)` contains higher-precision pi.

The coordinating agent also ran the separate `audit_bounds.jl` audit:
**17/17 exact-normalization and direct scalar interval residual checks passed**.
Its evidence is in `results/audit_bounds.json`.
These checks supplement the mathematical rounding analysis; finite test cases
alone would not establish a general rounding theorem.

## Limits of acceptance

The initial fast checker propagates saved large-block eigenvalue enclosures;
it does not reconstruct their eigenvectors or original residuals. The separate
17-test audit uses direct scalar interval accumulation on small blocks. Neither
initial check alone constituted a fresh large-block spectral execution.

## Completed full spectral replay

The reviewer subsequently executed `spectral_replay.jl` on all five fields.
For each field this rebuilt every exact 8/16-site momentum block, computed full
dense Q and Lambda, recomputed R and E and their rounding allowances, and
propagated fresh eigenvalue radii through logZ, rational seeds and recurrence.
Original eigenvalues, radii and seeds were comparison data, not inputs to these
fresh spectral stages. This is fresh execution of the same reviewed verifier,
**not an independently implemented second spectral verifier**. No search ran.

A first replay asserted bitwise equality with the original certificates and stopped at
h=3/4, where the fresh seeds differ in the last bits; a second stopped because the fresh
h=3/4 bound did not support the original nine-decimal display. The replay was then completed
with the acceptance rule below (positive fresh bound, overlapping enclosures). Those attempts
are in branch history (commit `cb635a1`), not in the tree.

| h | fresh raw bound | full spectral stage seconds | accepted bound |
|---|---:|---:|---:|
| 1 | 1.1482614900945929 | 47.904 | 1.148261 |
| 3/4 | 0.8907821437909684 | 31.552 | 0.890782 |
| 1/2 | 0.5947852422865794 | 41.783 | 0.594785 |
| 7/20 | 0.38350830551700044 | 33.863 | 0.383508 |
| 1/4 | 0.21448610675192437 | 33.450 | 0.214486 |

Accepted values are rounded **down to six decimals from min(original,fresh)**,
with exact rational values stored in the replay JSON. Every fresh bound is
positive, and every matched original/fresh eigenvalue enclosure overlaps.
Overlap is only a consistency check; the fresh certificate itself supplies the
proof. All original certificate files remain unchanged.

Only h=1 reproduced all original spectral fields, seeds and raw bound exactly.
At h=3/4 the original raw bound was 0.8907822037452694; at h=7/20 it was
0.3835083070577823. Their fresh bounds are smaller and do **not** reproduce the
original nine-decimal claims 0.890782203 and 0.383508307. The conservative values
above supersede those display claims. The other two changed fresh bounds are
slightly stronger. Exact fresh-minus-original seed differences are:

| h | spatial defect difference | thermal defect difference |
|---|---|---|
| 1 | 0//1 | 0//1 |
| 3/4 | 2991185//9007199254740992 | 5980245//4503599627370496 |
| 1/2 | -2009921//9007199254740992 | -4027265//4503599627370496 |
| 7/20 | 907975//9007199254740992 | 3666675//9007199254740992 |
| 1/4 | -431123//9007199254740992 | -1822071//9007199254740992 |

The observed h=3/4 maximum eigenvalue difference is approximately
1.95399e-14. Its largest radius increase occurs at L=16,p=8,k=2:
1.9753391920702282e-10. In that block the entry-error bound is unchanged,
the residual upper bound changes by 3.6492607308249224e-15, and the
orthogonality upper bound increases by 1.471519586164419e-11. Thus the
observed radius variation there is dominated by the orthogonality contribution,
not by shifts of the eigenvalues. These reported differences are numerical
diagnostics; no unverified explanation about a library-level cause is asserted.

The five conservative 8/16-site certificate bounds above are accepted. Reference gaps, optimization,
twenty-site partial-spectrum candidates, apparent search thresholds, and scaling
extrapolations remain numerical. In particular the stronger best twenty-site
candidates have not undergone complete dense certification, and no physical
threshold is established.

## CertifiedGaps migration review

After the core was renamed from `ThermalMPOGap` to `CertifiedGaps`, the reviewer checked that
every hashed source changed only by the package name, and a fresh 270-block replay on the
renamed core (189.418 seconds, no resumed witnesses) reproduced the preceding replay's raw
bounds exactly. Its compact summary is `results/spectral_replay_summary.json`.
