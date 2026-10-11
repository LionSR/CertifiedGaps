# Independent certification audit

**PASS: accept the five certificates at h = 1, 3/4, 1/2, 7/20, and 1/4,
computed from complete 8- and 16-site spectra.** Their uniform gap bounds apply
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
Its evidence is in `results/audit_bounds.json` and `results/audit_bounds.log`.
These checks supplement the mathematical rounding analysis; finite test cases
alone would not establish a general rounding theorem.

## Limits of acceptance

The fast checker propagates the saved large-block eigenvalue enclosures; it
does not reconstruct their large eigenvector matrices or recompute all original
residuals. Source hashes establish which implementation produced the stored
records. Re-running `certify.jl` is the available full reconstruction path.

The five 8/16-site certificates are accepted. Reference gaps, optimization,
twenty-site partial-spectrum candidates, apparent search thresholds, and scaling
extrapolations remain numerical. In particular the stronger best twenty-site
candidates have not undergone complete dense certification, and no physical
threshold is established.
