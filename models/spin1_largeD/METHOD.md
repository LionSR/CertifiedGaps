# Method: model, enclosures, transfer squaring, and Weyl transport

This document presents the mathematical derivations for both certification stages of the spin-1
large-$D$ chain. Stage 1 (`certify.jl`) runs the shared core pipeline
`CertifiedGaps.GapCertificates.certify_gap` on the exact model of `src/model.jl` at $n_0 = 6$ and
stores the verified spectra of the rings $L = 6, 12$. Stage 2 (`src/TransferCertificate.jl`,
`scripts/transfer_certify.jl`) reuses those stored spectra to certify lifted seeds at
$n_0 = 12, 24, 48, 96$ and continuous $D$ intervals without diagonalizing rings larger than
$L = 12$.

## 1. Exact model, frame, and reflection positivity

The core `RingModel` is written in the $S^z$ eigenbasis $(|{+1}\rangle, |0\rangle, |{-1}\rangle)$,
where $`S^\pm = \sqrt{2}\,S_\pm`$ with the integer matrices
$S_+ = |{+1}\rangle\langle 0| + |0\rangle\langle {-1}|$ (`SP`) and $S_- = S_+^T$ (`SM`). Then
```math
h = \mathbf{S} \cdot \mathbf{S} = S^z \otimes S^z + S_+ \otimes S_- + S_- \otimes S_+, \qquad
o = D (S^z)^2 = D \mathrm{diag}(1, 0, 1),
```
both with exact integer (rational) entries, and the charge $`m \in \{1, 0, -1\}`$ is conserved.

**Frame.** The Lean model uses the Cartesian basis $(|x\rangle, |y\rangle, |z\rangle)$ with
$(S^\alpha)_{ab} = -i\varepsilon_{\alpha ab}$. With
```math
W = \begin{pmatrix} -i & 0 & i \\ 1 & 0 & 1 \\ 0 & i & 0 \end{pmatrix}, \qquad
W^\dagger W = N = \mathrm{diag}(2, 1, 2),
```
the matrix $U = W N^{-1/2}$ is unitary and maps the Cartesian frame to the $S^z$ frame:
`cartesian_frame` checks in exact Gaussian-integer arithmetic that $W^\dagger S^z W = N S^z$,
$W^\dagger S^\pm W = 2 S_\pm$ and $W^\dagger (S^z)^2 W = N \mathrm{diag}(1, 0, 1)$, i.e.
$U^\dagger S^z U = S^z$, $U^\dagger S^\pm U = \sqrt{2} S_\pm$ and
$U^\dagger (S^z)^2 U = \mathrm{diag}(1,0,1)$.
Hence $U^{\otimes L}$ conjugates the Cartesian ring Hamiltonian into the core one, and the two have
the same partition functions. `test/frame.jl` also compares the verified spectra with dense
diagonalization of the Cartesian Hamiltonian for $L \le 6$.

**Reflection positivity.** Let $R = \mathrm{antidiag}(1, -1, 1)$, a real symmetric
involution commuting with $o$. Then $R S^z R = -S^z$ and $R S_\pm R = -S_\mp$, so
```math
-h = S^z \otimes R S^z R + S_+ \otimes R S_+ R + S_- \otimes R S_- R
= \sum_{C \in \{S^z, S_+, S_-\}} C \otimes R\,\overline{C}\,R,
```
the "reflected" RP form of the core (`rp_data`, checked exactly by
`CertifiedGaps.RingModels.verify_rp`; the exact LDL factorization of the realignment finds the
same form, while without $R$ the bond is not RP in this basis). In the Cartesian frame this is
the familiar $-h = \sum_\alpha G_\alpha \otimes G_\alpha$ with $G_\alpha = i S^\alpha$.

## 2. Symmetry-sector decomposition

The core `ring_blocks` splits the ring exactly by total charge $q$, by lattice momentum
$k = 2\pi m/L$ (for this real model only $0 \le m \le L/2$, with multiplicity $2$ for
$0 < m < L/2$), and uses the spin flip $F = \mathrm{antidiag}(1, 1, 1)$ (`flip`, a charge-
reversing permutation commuting with $h$ and $o$): the sectors $q$ and $-q$ are built once, and the
self-conjugate sector $q = 0$ is split by flip parity. At $L = 12$ this gives 92 blocks of
dimension at most $5{,}798$, the same decomposition as the earlier model-specific code; every run
checks that the blocks cover all $3^L$ states.

## 3. Verified eigenvalue enclosures (Stage 1)

Each block is enclosed by `CertifiedGaps.VerifiedSpectra.verified_eigen`: its interval entries
(exact couplings times enclosed phases and orbit-norm ratios) give an exactly Hermitian Float64
center $M$ (real for $m = 0, L/2$) and an exact bound on $`\|A - M\|_2`$ by the Frobenius norm of
the entrywise radii. With the LAPACK eigenpairs $(\Lambda, X)$ of $M$, the residual
$`\|MX - X\Lambda\|_2`$ and orthogonality $`\|X^\dagger X - I\|_2`$ are bounded by their computed
Frobenius norms plus a priori IEEE binary64 rounding bounds ($\gamma_{N+8}$ for real blocks,
$\gamma_{2N+8}$ for complex blocks evaluated with four real matrix products, valid for any
summation order and FMA, with an additive underflow term). The radius is the smaller of the
Ostrowski and polar-factor Weyl bounds plus the matrix error; every sorted eigenvalue of the
true block lies within it of the sorted centers. See the docstring of
`julia/src/VerifiedSpectra.jl` for the inequalities.

## 4. Partition functions and Stage-1 seeds

`seed_bounds` encloses $\log Z_6(\beta_0)$, $\log Z_{12}(\beta_0)$ and $\log Z_{12}(2\beta_0)$ by
`log_partition_enclosure` (outward-rounded interval sums over every block, each with its radius
and multiplicity) and rounds the defect enclosures
$1 - P_{\mathrm{s}}(6, \beta_0) = 1 - e^{\log Z_{12}(\beta_0) - 2\log Z_6(\beta_0)}$ and
$1 - P_\tau(12, \beta_0) = 1 - e^{\log Z_{12}(2\beta_0) - 2\log Z_{12}(\beta_0)}$ up to the grid
$10^{-15}$. `certify_gap` then runs the exact rational chain (`rational_gap_chain`) and the
interval recurrence (`certified_gap_bound`). Both bounds are rigorous; the tables report the
interval bound and list the chain bound. `certify.jl` writes both spectra with `write_spectrum`
and pins their SHA-256 in the certificate (`spectrum_sha256`).

## 5. Spatial transfer squaring and purity coupling without $L = 24$ (Stage 2)

For even $n$, $\sigma_n = X_\beta^n / Z_n(\beta)$ is positive semidefinite with unit trace. Writing
its eigenvalues as $p_i \ge 0$ ($\sum_i p_i = 1$), $P = \sum_i p_i^2 = P_{\mathrm{s}}(n, \beta)$
and $P' = \sum_i p_i^4 / P^2 = P_{\mathrm{s}}(2n, \beta)$ satisfy
```math
1 - P' = \frac{2\sum_{i < j} p_i^2 p_j^2}{P^2} \le \frac{(1 - P)^2}{2P^2},
```
because $(1 - P)^2 = (2\sum_{i < j} p_i p_j)^2 \ge 4\sum_{i < j} p_i^2 p_j^2$. A spatial defect
bound $s \ge 1 - P_{\mathrm{s}}(n, \beta)$ therefore lifts to
$f(s) = \min(1, s^2 / [2(1 - s)^2]) \ge 1 - P_{\mathrm{s}}(2n, \beta)$ (`square_defect`).

From the four finite log-partition enclosures $\log Z_6(\beta_0)$, $\log Z_{12}(\beta_0)$,
$\log Z_6(2\beta_0)$, and $\log Z_{12}(2\beta_0)$, we obtain upper bounds on
$s = 1 - P_{\mathrm{s}}(6, \beta_0)$, $s_2 = 1 - P_{\mathrm{s}}(6, 2\beta_0)$, and
$t = 1 - P_\tau(12, \beta_0)$. Using the exact partition-function identity
```math
P_\tau(4n, \beta_0) = \frac{P_{\mathrm{s}}(2n, 2\beta_0)\,P_\tau(2n, \beta_0)^2}{P_{\mathrm{s}}(2n, \beta_0)^2}
```
together with $P_{\mathrm{s}}(2n, \beta_0) \le 1$ yields the coupled spatial-lift map at fixed
inverse temperature $\beta_0$:
```math
s' = f(s) \ge 1 - P_{\mathrm{s}}(12, \beta_0), \qquad
s_2' = f(s_2) \ge 1 - P_{\mathrm{s}}(12, 2\beta_0), \qquad
t' = 1 - (1 - s_2')(1 - t)^2 \ge 1 - P_\tau(24, \beta_0).
```
Iterating this map (`seed_ladder`) generates certified seed pairs $(s^{(k)}, t^{(k)})$ at
$n_0 = 12, 24, 48, 96$ from the $L = 6, 12$ spectra alone, which are then fed into the standard
space-and-temperature doubling recurrence `certified_gap_bound`.

## 6. Centered Weyl parameter transport in $D$ (Stage 2)

Writing $H_L(D) = H_L(D^*) + \delta D \cdot V$ with $V = \sum_{j=1}^L (S_j^z)^2$ and
$0 \le V \le L I$, centering $V$ by $\frac{L}{2}I$ gives
$`\|\delta D (V - \frac{L}{2}I)\|_2 \le \frac{1}{2}L|\delta D|`$. By Weyl's inequality, every ordered
eigenvalue of $`H_L(D) - \frac{L\,\delta D}{2}I`$ lies within $\frac{1}{2}L|\delta D|$ of the
corresponding eigenvalue of $H_L(D^*)$, so $\log Z_L(\beta)$ of the shifted Hamiltonian lies within
$\frac{1}{2}\beta L |\delta D|$ of the reference enclosure at $D^*$ (and the constant energy shift
cancels identically in $P_{\mathrm{s}}$ and $P_\tau$). For a compact interval $[D_{\mathrm{lo}}, D_{\mathrm{hi}}]$,
we replace $|\delta D|$ by $\max_{D \in [D_{\mathrm{lo}}, D_{\mathrm{hi}}]} |D - D^*|$
(`transported_logz`).


## 7. Finite-ring gap enclosures (Stage 2)

Stage 2 sums $\log Z$ over the complete stored spectra (no truncation). The finite-ring gap
$\Delta_L = E_1(L) - E_0(L)$ is enclosed by the core `spectral_gap_enclosure` (second-smallest
lower endpoint minus smallest upper endpoint, with multiplicities), widened by $L|\delta D|$ when
transported in $D$. At $D = 4$, verified finite-ring gaps at $L = 6, 8, 10$ (`bridge_certificate`)
exceed the $n_0 = 12$ uniform bound, extending that certificate to all even $L \ge 6$.
