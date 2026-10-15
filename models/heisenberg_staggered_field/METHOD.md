# Method: uniform gap certificates for the Heisenberg chain in a staggered field

This document describes the reflection-positivity transformation, the symmetry-block construction,
the eigenvalue-enclosure theorem, the 20-site numerical exploration, and the independent audit for
the spin-1/2 Heisenberg chain in a staggered longitudinal field.

## Outcome and scope

All five target staggered fields $h \in \{1, 3/4, 1/2, 7/20, 1/4\}$ are certified from complete
verified 8- and 16-site ring spectra at seed length $n_0 = 8$. Because the alternating sublattice
rotation produces a homogeneous reflection-positive ring with a single-spin unit cell on any even
ring length $L$, the resulting uniform spectral gap bounds apply to **every even chain length
$L \ge 8$ sites**, establishing both ground-state uniqueness and $\Delta_L \ge \Delta_{\mathrm{cert}}$.

In addition to the certified $n_0 = 8$ seeds, an exploratory numerical scan (`scripts/scan.jl`)
reaches $2n_0 = 20$ sites ($n_0 = 10$) by combining dense diagonalization for symmetry blocks of
dimension $\le 2{,}500$ with 80-state `KrylovKit` Lanczos runs for larger blocks (the largest
$L = 20$ block has dimension $18{,}504$, where a single dense `ComplexF64` matrix occupies
$\approx 5.48\text{ GB}$ and cubic eigensolve work increases by $\approx 1{,}490\times$ relative to
$L = 16$). Those 20-site candidate bounds are reported below as non-certified numerical benchmarks.

## Model, exact RP identity, and length convention

We consider the periodic spin-1/2 chain
$$H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} - h \sum_{j=1}^L (-1)^j S_j^z, \qquad h \in \{1, 3/4, 1/2, 7/20, 1/4\},$$
where the $L = 2$ ring includes both periodic bonds. Define the exact rational spin-1/2 matrices
$$X = S^x = \begin{pmatrix} 0 & 1/2 \\ 1/2 & 0 \end{pmatrix}, \qquad
G_y = i S^y = \begin{pmatrix} 0 & 1/2 \\ -1/2 & 0 \end{pmatrix}, \qquad
Z = S^z = \begin{pmatrix} 1/2 & 0 \\ 0 & -1/2 \end{pmatrix}, \qquad U = 2G_y.$$
(In `src/common.jl`, $G_y$ is denoted `J`.) Rotating every odd site by the real orthogonal matrix
$U$ (a $\pi$-rotation about the $y$-axis) transforms both odd-to-even and even-to-odd nearest-neighbor
bonds into the identical homogeneous form
$$h_{\mathrm{rot}} = -(X \otimes X + G_y \otimes G_y + Z \otimes Z), \qquad
-h_{\mathrm{rot}} = X \otimes X + G_y \otimes G_y + Z \otimes Z = \sum_{\alpha \in \{x,y,z\}} C_\alpha \otimes \overline{C_\alpha}.$$
Every bond factor $C_\alpha \in \{X, G_y, Z\}$ has exact rational entries (note that reflection
positivity requires real $C_\alpha$ when $F = I$, not symmetric $C_\alpha$, so the antisymmetry of
$G_y$ is completely compatible). Since $U Z U^\top = -Z$, the staggered field becomes the uniform
real symmetric on-site term $o = -h Z$. The function `exact_rp` in `src/common.jl` verifies the
bond identity in both orientations and the on-site identity in exact `Rational{BigInt}` arithmetic.

## Dense-grid certificates ($n_0 = 8$)

Candidate seeds were scanned over $n_0 \in \{2, 4, 6, 8\}$ and $\beta_0 \in \{0.25, 0.50, \dots, 32.0\}$
via `scan_seeds` and `best_seed`, and the optimal $(n_0, \beta_0)$ for each field was certified
with full dense eigenvector residual and orthogonality bounds. All exact rational parameters and
seed defects are stored as `"num//den"` strings in `results/certificate_<h>.json`.

| $h$ | $n_0$ | $\beta_0$ | Certified gap $\Delta_L \ge$ | 16-site ED gap | Bound / ED gap | Certification time (s) |
|---|---:|---:|---:|---:|---:|---:|
| `1//1` | 8 | `21//4` | 1.148261 | 1.87351193 | 0.612892543 | 79.25 |
| `3//4` | 8 | `21//4` | 0.890782 | 1.57089206 | 0.567054967 | 80.25 |
| `1//2` | 8 | `11//2` | 0.594785 | 1.23131976 | 0.483046933 | 50.40 |
| `7//20` | 8 | `23//4` | 0.383508 | 0.994380622 | 0.385675562 | 47.08 |
| `1//4` | 8 | `25//4` | 0.214486 | 0.812031433 | 0.264135226 | 44.34 |

The displayed six-decimal bounds are rounded **down** from the minimum of the original certificate
and an independent full spectral replay (`results/spectral_replay_summary.json`); see
"Independent audit" below. For $h = 1$, the stored upper seed defects at $\beta_0 = 21/4$ are
$\varepsilon_{\mathrm{s}} = \text{3124087506499}//\text{4503599627370496}$ and
$\varepsilon_\tau = \text{6599110342439}//\text{9007199254740992}$, with the following worst-case
block error enclosures:

| $h = 1$ ring size $L$ | Max entry error $a$ | Max residual upper $r$ | Max orthogonality upper $e$ | Max eigenvalue radius $\eta$ |
|---:|---:|---:|---:|---:|
| 8 | $1.5034199 \times 10^{-15}$ | $7.08836598 \times 10^{-12}$ | $1.83181771 \times 10^{-12}$ | $2.2315322 \times 10^{-11}$ |
| 16 | $5.27181194 \times 10^{-14}$ | $5.1456421 \times 10^{-8}$ | $9.56335631 \times 10^{-9}$ | $2.10443015 \times 10^{-7}$ |

## Complete symmetry sectors and multiplicities

In the physical basis, let $p \in \{0, \dots, L\}$ denote the number of down spins, $m = L/2$ the
number of two-site unit cells of the physical field pattern, and $k \in \{0, \dots, m-1\}$ the
momentum index under translation by two sites ($T^2$). A $T^2$-orbit of representative configuration
$a$ with period $\ell$ supports momentum $k$ if and only if $k\ell$ is a multiple of $m$, with
normalized basis state $\frac{1}{\sqrt{\ell}} \sum_{r=0}^{\ell-1} e^{-2\pi i k r / m} |T^{2r} a\rangle$.
For a Hamiltonian transition from representative $a$ into $T^{2t} b$, the projected matrix element
equals the local transition amplitude multiplied by $\sqrt{\ell_a / \ell_b}\,e^{2\pi i k t / m}$.
Small-ring unit tests verify that the merged spectrum across all $(p, k)$ blocks matches full dense
diagonalization in the shared core.

Global spin flip composed with a one-site translation commutes with $H_L$ and $T^2$ while mapping
$p \mapsto L - p$, and complex conjugation maps $k \mapsto m - k \bmod m$. We diagonalize one dense
representative block for each equivalence class, assigning exact multiplicity $1, 2,$ or $4$, and
assert $\sum_{\mathrm{blocks}} \mathrm{multiplicity} \times d_{\mathrm{block}} = 2^L$ on every ring.

| Ring size $L$ | Representative blocks | Largest block dimension | Weighted dimension sum | Full Hilbert dimension $2^L$ |
|---:|---:|---:|---:|---:|
| 8 | 13 | 20 | 256 | 256 |
| 16 | 41 | 1620 | 65536 | 65536 |

For $L = 16$ ($m = 8$), the representative block dimensions by down-spin count $p \in \{0, \dots, 8\}$
and momentum $k \in \{0, 1, 2, 3, 4\}$ (which carry momentum multiplicities $1, 2, 2, 2, 1$) are:

| Down spins $p$ | Block dimensions for $k = 0, 1, 2, 3, 4$ | Charge multiplicity ($p \leftrightarrow L - p$) |
|---:|---|---:|
| 0 | `0:1` | 2 |
| 1 | `0:2, 1:2, 2:2, 3:2, 4:2` | 2 |
| 2 | `0:16, 1:14, 2:16, 3:14, 4:16` | 2 |
| 3 | `0:70, 1:70, 2:70, 3:70, 4:70` | 2 |
| 4 | `0:232, 1:224, 2:230, 3:224, 4:232` | 2 |
| 5 | `0:546, 1:546, 2:546, 3:546, 4:546` | 2 |
| 6 | `0:1008, 1:994, 2:1008, 3:994, 4:1008` | 2 |
| 7 | `0:1430, 1:1430, 2:1430, 3:1430, 4:1430` | 2 |
| 8 | `0:1620, 1:1600, 2:1616, 3:1600, 4:1618` | 1 |

## Eigenvalue-enclosure theorem and implementation

All `Float64` arrays returned by LAPACK are treated as untrusted approximate proposals whose exact
binary64 entries are subsequently certified. Let $H$ denote the true mathematical Hermitian momentum
block and $A$ the `ComplexF64` Hermitian matrix passed to LAPACK. The function `block(...; verified=true)`
in `src/common.jl` evaluates all rational couplings, orbit-normalization square roots
$\sqrt{\ell_a / \ell_b}$, and trigonometric phases $\exp(2\pi i k t / m)$ in outward-rounded
`IntervalArithmetic` (using an interval enclosure of $\pi$) to enclose each upper-triangular entry
of $H$ and bound the entry-construction error $a \ge \|H - A\|_2 \le \|H - A\|_F$ (stored as
`entry_error`).

Given the approximate eigenvector matrix $Q$ and real eigenvalue vector $\Lambda$ from LAPACK, the
verifier computes the floating-point residual $\hat{R} = \operatorname{fl}(AQ - Q\Lambda)$ and
orthogonality defect $\hat{E} = \operatorname{fl}(Q^* Q - I)$. Upper bounds $q \ge \|Q\|_F$,
$h_{\mathrm{norm}} \ge \|A\|_F$, $\ell \ge \|\Lambda\|_F$, $\|\hat{R}\|_{F,\mathrm{up}}$, and
$\|\hat{E}\|_{F,\mathrm{up}}$ are accumulated using directed upward rounding (`nextfloat` after
each scalar operation, separating real and imaginary components). With unit roundoff $u = 2^{-53}$
and $\gamma_{32d+32} = (32d+32)u / (1 - (32d+32)u)$, the true residual and orthogonality norms
satisfy
$$\|AQ - Q\Lambda\|_F \le r := \|\hat{R}\|_{F,\mathrm{up}} + \gamma_{32d+32}(h_{\mathrm{norm}} q + q\ell) + 10^{-280},$$
$$\|Q^* Q - I\|_F \le e := \|\hat{E}\|_{F,\mathrm{up}} + \gamma_{32d+32}(q^2 + \sqrt{d}) + 10^{-280}.$$
This follows from the standard IEEE binary64 dot-product rounding bound
$|\operatorname{fl}(\sum x_i y_i) - \sum x_i y_i| \le \gamma_s \sum |x_i y_i|$ applied to the real
and imaginary parts and the submultiplicativity $\||A|\,|B|\|_F \le \|A\|_F \|B\|_F$, where $32d+32$
conservatively bounds the operation count per complex inner product and subtraction, $\sqrt{d}$
accounts for subtracting the identity matrix $I$, and $10^{-280}$ covers gradual underflow.

**Polar-factor and Weyl inequality (with multiplicity):** Asserting $e < 1$, the unitary polar
factor $V = Q(Q^* Q)^{-1/2}$ exists and satisfies $\|Q - V\|_2 \le d_Q := e / (1 + \sqrt{1 - e})$.
Consequently,
$$\|A - V \Lambda V^*\|_2 = \|AV - V\Lambda\|_2 \le r + (\|A\|_2 + \|\Lambda\|_2) d_Q.$$
Using $\|A\|_2 \le L(3/4 + h/2) + a$ and $\|\Lambda\|_2 = \max_i |\lambda_i|$, Weyl's inequality
for the Hermitian matrices $H$ and $V \Lambda V^*$ guarantees that the $i$th ordered eigenvalue of
$H$ (counting multiplicity) lies in $[\lambda_i - \eta, \lambda_i + \eta]$ with radius
$$\eta = a + r + \left(L(3/4 + h/2) + a + \max_i |\lambda_i|\right) d_Q.$$
All scalar operations defining $r, e, d_Q,$ and $\eta$ are performed with outward interval rounding.

For $\beta > 0$, `enclosed_logz` sums $\exp(-\beta [\lambda_i - \eta, \lambda_i + \eta])$ weighted
by each block's exact symmetry multiplicity (after factoring out a common exact binary64 energy
shift) to obtain rigorous intervals for $\log Z_{n_0}(\beta_0)$, $\log Z_{2n_0}(\beta_0)$, and
$\log Z_{2n_0}(2\beta_0)$. The upper endpoints of the resulting defect intervals
$1 - P_{\mathrm{s}}(n_0, \beta_0)$ and $1 - P_\tau(2n_0, \beta_0)$ are converted directly to exact
`Rational{BigInt}` fractions (`epsilon_s` and `epsilon_t`) without decimal truncation and passed to
`CertifiedGaps.certified_gap_bound`.

## Twenty-site numerical exploration and scaling

The 20-site scan (`scripts/scan.jl 20`) combines full dense diagonalization for blocks of dimension
$\le 2{,}500$ with 80 low-lying Ritz values from `KrylovKit` for larger blocks (up to dimension
$18{,}504$, with total weighted dimension $\sum \mathrm{mult} \times d = 1{,}048{,}576 = 2^{20}$).
Because Krylov convergence and eigenvalue completeness are not certified, the $n_0 = 10$ rows below
are exploratory numerical estimates:

| $h$ | Best $n_0$ | $\beta_0$ | Numerical candidate bound | Numerical 20-site ED gap |
|---|---:|---:|---:|---:|
| `1//1` | 10 | 6.25 | 1.27832796 | 1.87351191 |
| `3//4` | 10 | 6.5 | 1.01075578 | 1.57089173 |
| `1//2` | 10 | 6.75 | 0.707243353 | 1.23131113 |
| `7//20` | 10 | 7 | 0.493229271 | 0.994292667 |
| `1//4` | 10 | 7.25 | 0.326952394 | 0.811528896 |
| `1//5` | 10 | 7.5 | 0.231456632 | 0.708641218 |
| `3//20` | 10 | 8 | 0.119366516 | 0.594804607 |
| `1//10` | — | — | no grid seed at $n_0 \le 10$ | 0.466432604 |

At $n_0 \le 10$, the sampled $\beta_0$ grid finds an admissible seed at $h = 3/20 = 0.15$ and none
at $h = 1/10 = 0.10$. The table below records $\max_\beta \min(P_{\mathrm{s}}(n_0, \beta), P_\tau(2n_0, \beta))$
and the ground-state spatial energy defect $\delta(n_0) = E_0(2n_0) - 2E_0(n_0)$ across all scanned
fields and seed lengths $n_0 \in \{2, 4, 6, 8, 10\}$:

| $h$ | $n_0$ | $\max_\beta \min(P_{\mathrm{s}}, P_\tau)$ | $\delta(n_0) = E_0(2n_0) - 2E_0(n_0)$ |
|---|---:|---:|---:|
| `1//1` | 2 | 0.529873172 | 0.505022849 |
| `1//1` | 4 | 0.94067001 | 0.0221827722 |
| `1//1` | 6 | 0.99334286 | 0.00159103285 |
| `1//1` | 8 | 0.999275644 | 0.00013128857 |
| `1//1` | 10 | 0.999923253 | 1.14962961e-05 |
| `3//4` | 2 | 0.447640968 | 0.613999064 |
| `3//4` | 4 | 0.890324463 | 0.0389149028 |
| `3//4` | 6 | 0.983514815 | 0.00385786603 |
| `3//4` | 8 | 0.997555655 | 0.000440747681 |
| `3//4` | 10 | 0.999636029 | 5.36822693e-05 |
| `1//2` | 2 | 0.399066156 | 0.754873673 |
| `1//2` | 4 | 0.800725359 | 0.0769369751 |
| `1//2` | 6 | 0.952122862 | 0.0113905614 |
| `1//2` | 8 | 0.988846519 | 0.00192196739 |
| `1//2` | 10 | 0.997543855 | 0.000347058281 |
| `7//20` | 2 | 0.361338696 | 0.850802846 |
| `7//20` | 4 | 0.699446862 | 0.124633121 |
| `7//20` | 6 | 0.893328792 | 0.0251937055 |
| `7//20` | 8 | 0.967347486 | 0.00566434475 |
| `7//20` | 10 | 0.990324488 | 0.00135842657 |
| `1//4` | 2 | 0.31998312 | 0.913017541 |
| `1//4` | 4 | 0.602523494 | 0.177017263 |
| `1//4` | 6 | 0.818031626 | 0.0466268272 |
| `1//4` | 8 | 0.926399235 | 0.0132591633 |
| `1//4` | 10 | 0.971872848 | 0.00397943173 |
| `1//5` | 2 | 0.303223146 | 0.940972384 |
| `1//5` | 4 | 0.560054301 | 0.2117524 |
| `1//5` | 6 | 0.757748899 | 0.0653997426 |
| `1//5` | 8 | 0.881241013 | 0.0214187006 |
| `1//5` | 10 | 0.947706946 | 0.00732718927 |
| `3//20` | 2 | 0.289667018 | 0.965093873 |
| `3//20` | 4 | 0.495941276 | 0.252166166 |
| `3//20` | 6 | 0.674748062 | 0.0934440211 |
| `3//20` | 8 | 0.814781105 | 0.0361784563 |
| `3//20` | 10 | 0.902221022 | 0.0144183736 |
| `1//10` | 2 | 0.279689865 | 0.983874689 |
| `1//10` | 4 | 0.442024104 | 0.295440509 |
| `1//10` | 6 | 0.569803257 | 0.134401906 |
| `1//10` | 8 | 0.7028359 | 0.0639327614 |
| `1//10` | 10 | 0.80680422 | 0.0308278889 |

### Asymptotic $h^{2/3}$ scaling comparison

In the continuum sine-Gordon description at weak staggered field $h \ll 1$, the spectral gap and
inverse correlation length scale as $\Delta \propto h^{2/3}$ (up to logarithmic corrections).
Calibrating to the $h = 1$ seed $(n_0, \beta_0) = (8, 5.25)$ gives the constant-defect-quality
scaling locus $n_0(h) \approx 8 h^{-2/3}$ and $\beta_0(h) \approx 5.25 h^{-2/3}$:

| $h$ | Scaled $n_0 \approx 8 h^{-2/3}$ | Next even $n_0$ | Scaled $\beta_0 \approx 5.25 h^{-2/3}$ | Certified $n_0$ |
|---|---:|---:|---:|---:|
| `1//1` | 8 | 8 | 5.25 | 8 |
| `3//4` | 9.69130983 | 10 | 6.35992207 | 8 |
| `1//2` | 12.6992084 | 14 | 8.33385552 | 8 |
| `7//20` | 16.1081114 | 18 | 10.5709481 | 8 |
| `1//4` | 20.1587368 | 22 | 13.229171 | 8 |

Because the admissible region allows seed defects up to $\approx 0.12731$ (whereas the $h = 1$ seed
has defects $\sim 7 \times 10^{-4}$), $n_0 = 8$ remains sufficient all the way down to $h = 1/4$,
with the certified bound-to-ED-gap ratio decreasing from $0.613$ at $h = 1$ to $0.264$ at $h = 1/4$.
At $h = 1/10$, fitting $\log \delta(n_0) \approx 0.19944 - 0.36810 n_0$ over $n_0 \in \{6, 8, 10\}$
gives a decay scale of $2.72$ sites and predicts a first admissible seed near
$(n_0, \beta_0) \approx (14, 17.0)$ (requiring a 28-site doubled ring).

## Independent audit

An independent audit reviewed `src/common.jl`, `certify.jl`, `test/runtests.jl`, `audit_bounds.jl`,
and the core recurrence, and executed a complete fresh spectral replay (`scripts/spectral_replay.jl`)
across all five fields:

1. **Mathematical and numerical verification:** Confirmed the exact rational reflection-positivity
   identity under the odd-site rotation $U = 2G_y$, the even-$L$ closure condition, the completeness
   of the $(p, k)$ symmetry-sector decomposition ($\sum \mathrm{mult} \times d = 2^L$), the interval
   enclosure of the momentum-block entries, the $\gamma_{32d+32}$ plus $10^{-280}$ rounding bounds
   for residual and orthogonality Frobenius norms, the polar-factor/Weyl eigenvalue enclosure
   theorem, and the outward-rounded $\log Z$ and rational seed propagation.
2. **Direct interval residual cross-check (`scripts/audit_bounds.jl`):** Passed 17/17 exact
   normalization and direct scalar `IntervalArithmetic` residual checks on small blocks
   (`results/audit_bounds.json`).
3. **Full spectral replay (`scripts/spectral_replay.jl`):** Rebuilt every 8- and 16-site momentum
   block from scratch across all five fields, re-solved the dense eigenproblems, recomputed the
   residual and orthogonality error bounds, and propagated the fresh eigenvalue enclosures through
   $\log Z$, rational seed extraction, and `certified_gap_bound`:

| $h$ | Fresh replay raw bound | Replay wall time (s) | Accepted 6-decimal bound |
|---|---:|---:|---:|
| $1$ | 1.1482614900945929 | 47.904 | 1.148261 |
| $3/4$ | 0.8907821437909684 | 31.552 | 0.890782 |
| $1/2$ | 0.5947852422865794 | 41.783 | 0.594785 |
| $7/20$ | 0.38350830551700044 | 33.863 | 0.383508 |
| $1/4$ | 0.21448610675192437 | 33.450 | 0.214486 |

The fresh replay reproduced the $h = 1$ seed defects and raw bound bit-for-bit, while at the other
four fields LAPACK eigenvector differences in the last bits ($\sim 10^{-14}$) shifted the rational
seed defects by a few parts in $10^{10}$:

| $h$ | Spatial defect difference (fresh $-$ original) | Thermal defect difference (fresh $-$ original) |
|---|---|---|
| $1$ | `0//1` | `0//1` |
| $3/4$ | `2991185//9007199254740992` | `5980245//4503599627370496` |
| $1/2$ | `-2009921//9007199254740992` | `-4027265//4503599627370496` |
| $7/20$ | `907975//9007199254740992` | `3666675//9007199254740992` |
| $1/4$ | `-431123//9007199254740992` | `-1822071//9007199254740992` |

Every matched original and fresh eigenvalue interval overlaps, and rounding down to six decimal
places from $\min(\text{original}, \text{fresh})$ gives the accepted headline bounds verified by
`check.jl`. Following the rename of the core package from `ThermalMPOGap` to `CertifiedGaps`, a
second 270-block replay (189.4 s) reproduced the replay bounds identically
(`results/spectral_replay_summary.json`).

## Cross-check with the shared core (`CertifiedGaps.GapCertificates`)

The script `scripts/core_certificate.jl` independently certifies the $h = 1$ seed $(8, 21/4)$ using
`CertifiedGaps.GapCertificates`: it constructs the rotated single-spin-cell `RingModel`, verifies
reflection positivity via `RingModels.verify_rp`, encloses the 8- and 16-site spectra using the
shared core's $S^z$-parity and single-site momentum blocks (`VerifiedSpectra`), and runs the exact
rational doubling chain (`rational_gap_chain`), yielding $\Delta_L > 1.147699$ for every even
$L \ge 8$ (`results/core_certificate_1_1.json`). The fast checker verifies its RP identity, model
equivalence, and rational chain, and `check.jl --full` re-solves and re-verifies all of its symmetry
blocks from scratch.
