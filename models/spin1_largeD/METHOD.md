# Method: model, enclosures, transfer squaring, and Weyl transport

This document presents the mathematical derivations for both certification stages of the spin-1
large-$D$ chain: Stage 1 (`certify.jl`, `seed.jl`) constructs complete verified spectra at
$L = 6, 12$ and certifies the direct $n_0 = 6$ seeds, while Stage 2 (`src/TransferCertificate.jl`,
`scripts/transfer_certify.jl`) reuses those verified spectra to certify lifted seeds at
$n_0 = 12, 24, 48, 96$ and continuous $D$ intervals without diagonalizing rings larger than
$L = 12$.

## 1. Exact model and reflection positivity

In the Cartesian spin-1 basis $\{|x\rangle, |y\rangle, |z\rangle\}$, define the real, integer
antisymmetric matrices $G_\alpha = i S^\alpha$ ($S^\alpha = -i G_\alpha$):
$$G_x = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & -1 & 0 \end{pmatrix}, \qquad
G_y = \begin{pmatrix} 0 & 0 & -1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \qquad
G_z = \begin{pmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.$$
In `model.jl`, `rp_identity` verifies in exact integer arithmetic that $[S^x, S^y] = i S^z$,
$\sum_\alpha (S^\alpha)^2 = 2I$, that the two-site bond operator
$h = P - |\Omega\rangle\langle\Omega|$ (where $P$ is the swap and
$|\Omega\rangle = |xx\rangle + |yy\rangle + |zz\rangle$) equals $\sum_\alpha S^\alpha \otimes S^\alpha$,
and that
$$-h = G_x \otimes G_x + G_y \otimes G_y + G_z \otimes G_z = \sum_{\alpha \in \{x,y,z\}} G_\alpha \otimes \overline{G_\alpha}.$$
For any rational $D$, the single-ion anisotropy is the real symmetric on-site operator
$o = D(S^z)^2 = D\operatorname{diag}(1, 1, 0)$.

In the normalized $S^z$ eigenbasis $m \in \{-1, 0, 1\}$, each allowed nearest-neighbor spin-exchange
transition $(m_j, m_{j+1}) \mapsto (m_j \pm 1, m_{j+1} \mp 1)$ has unit matrix element $1$, and the
diagonal is $\sum_{j=1}^L m_j m_{j+1} + D \sum_{j=1}^L m_j^2$. The unit test suite (`test/reviewer.jl`)
verifies the exact unitary equivalence between the Cartesian and $S^z$ representations, including
the $L = 2$ convention of retaining both periodic bonds.

## 2. Symmetry-sector decomposition

Basis product states in the $S^z$ basis are partitioned by total magnetization $q = \sum_j m_j$;
sectors $q$ and $-q$ are exchanged by global spin flip $m_j \mapsto -m_j$ and have identical
spectra. Within the $q = 0$ sector, blocks are further labeled by single-site translation momentum
$k$ and spin-flip parity $p = \pm 1$ (with character $e^{2\pi i k t / L} p^f$), whereas for $q \ne 0$
translation momentum $k$ alone is used. A group orbit of representative configuration $a$
contributes to sector $(k, p)$ if and only if the character is trivial on the stabilizer of $a$.
In the normalized orbit basis, a Hamiltonian transition from representative $a$ to $T^u F^f b$ has
matrix element
$$h_{ab} \sqrt{\frac{|O_a|}{|O_b|}}\,e^{-2\pi i k u / L}\,p^f.$$
Complex conjugation pairs momentum $k$ with $L - k \bmod L$, so each independent block carries
charge multiplicity $1$ or $2$ and momentum multiplicity $1$ or $2$. Every run asserts that
$\sum_{\mathrm{blocks}} \mathrm{multiplicity} \times d_{\mathrm{block}} = 3^L$; at $L = 12$, this
partitions the $531{,}441$ states into $92$ independent blocks of dimension at most $5{,}798$.
`test/small.jl` verifies all block eigenvalues against full Hilbert-space diagonalization through
$L = 6$, and `test/symmetry.jl` checks the orbit counts and symmetry covariances through $L = 12$.

## 3. Verified eigenvalue enclosures (Stage 1)

Every rational transition amplitude, trigonometric phase, and orbit-ratio square root in a symmetry
block $H_{\mathrm{block}}$ is evaluated in outward-rounded `IntervalArithmetic`. The `ComplexF64`
Hermitian midpoint matrix $A$ is formed from the upper triangle, and the maximum row sum of
entrywise interval radii $|\operatorname{Re}| + |\operatorname{Im}|$ yields an upper bound
$c \ge \|H_{\mathrm{block}} - A\|_2$. Since each Heisenberg bond has eigenvalues in $\{-2, -1, 1\}$
and $0 \le (S_j^z)^2 \le I$, $\|H_{\mathrm{block}}\|_2 \le L(2 + |D|)$ and hence
$\|A\|_2 \le M := L(2 + |D|) + c$.

Given the approximate dense eigensystem $AQ \approx Q\Lambda$ of an $m \times m$ block, with unit
roundoff $u = 2^{-53}$ and $\gamma_{8m+32} = (8m + 32)u / (1 - (8m + 32)u)$, the standard IEEE
binary64 / single-threaded BLAS dot-product rounding model bounds the residual and orthogonality
errors by
$$\|AQ - Q\Lambda\|_2 \le r := \|\operatorname{fl}(AQ - Q\Lambda)\|_F + \gamma_{8m+32}(\sqrt{m}\,M \|Q\|_F + \|Q\|_F \|\Lambda\|_2) + 10^{-280},$$
$$\|Q^* Q - I\|_2 \le e := \|\operatorname{fl}(Q^* Q - I)\|_F + \gamma_{8m+32}(\|Q\|_F^2 + \sqrt{m}) + 10^{-280},$$
where each floating-point sum of $N$ squared magnitudes is enlarged by $\gamma_{4N+16}$ before
taking an outward-rounded square root, and $10^{-280}$ covers gradual underflow.

**Polar-factor and Weyl theorem:** Asserting $e < 1$, $Q$ is invertible and its unitary polar
factor $U = Q(Q^* Q)^{-1/2}$ satisfies $\|Q - U\|_2 \le d_Q := e / (1 + \sqrt{1 - e})$. Thus
$\|A - U\Lambda U^*\|_2 = \|AU - U\Lambda\|_2 \le r + (M + \|\Lambda\|_2)d_Q$. By Weyl's
inequality, each sorted eigenvalue of $H_{\mathrm{block}}$ (counting multiplicity) lies within
$$\eta = c + r + (M + \|\Lambda\|_2)\frac{e}{1 + \sqrt{1 - e}}$$
of the corresponding sorted diagonal entry of $\Lambda$. Merging and sorting all block centers with
their multiplicities yields the complete verified spectrum of $H_L$.

## 4. Partition functions and Stage-1 ED seeds

For a complete verified spectrum, $\sum_i \mathrm{mult}_i \cdot \exp(-\beta[\lambda_i - \eta_i, \lambda_i + \eta_i])$
is summed in outward-rounded interval arithmetic after factoring out an exact binary64 energy shift,
yielding rigorous intervals for $\log Z_n(\beta)$. The upper endpoints of
$$\varepsilon_{\mathrm{s}} = \operatorname{upper}\!\big(1 - \exp(\log Z_{2n}(\beta) - 2\log Z_n(\beta))\big), \qquad
\varepsilon_\tau = \operatorname{upper}\!\big(1 - \exp(\log Z_{2n}(2\beta) - 2\log Z_{2n}(\beta))\big)$$
are converted exactly to `Rational{BigInt}` and passed with the rational $\beta_0$ to
`CertifiedGaps.certified_gap_bound`.

## 5. Spatial transfer squaring and purity coupling without $L = 24$ (Stage 2)

For even $n$, $\sigma_n = X_\beta^n / Z_n(\beta)$ is positive semidefinite with unit trace. Writing
its eigenvalues as $p_i \ge 0$ ($\sum_i p_i = 1$), $P = \sum_i p_i^2 = P_{\mathrm{s}}(n, \beta)$
and $P' = \sum_i p_i^4 / P^2 = P_{\mathrm{s}}(2n, \beta)$ satisfy
$$1 - P' = \frac{2\sum_{i < j} p_i^2 p_j^2}{P^2} \le \frac{(1 - P)^2}{2P^2},$$
because $(1 - P)^2 = (2\sum_{i < j} p_i p_j)^2 \ge 4\sum_{i < j} p_i^2 p_j^2$. A spatial defect
bound $s \ge 1 - P_{\mathrm{s}}(n, \beta)$ therefore lifts to
$f(s) = \min(1, s^2 / [2(1 - s)^2]) \ge 1 - P_{\mathrm{s}}(2n, \beta)$ (`square_defect`).

From the four finite log-partition enclosures $\log Z_6(\beta_0)$, $\log Z_{12}(\beta_0)$,
$\log Z_6(2\beta_0)$, and $\log Z_{12}(2\beta_0)$, we obtain upper bounds on
$s = 1 - P_{\mathrm{s}}(6, \beta_0)$, $s_2 = 1 - P_{\mathrm{s}}(6, 2\beta_0)$, and
$t = 1 - P_\tau(12, \beta_0)$. Using the exact partition-function identity
$$P_\tau(4n, \beta_0) = \frac{P_{\mathrm{s}}(2n, 2\beta_0)\,P_\tau(2n, \beta_0)^2}{P_{\mathrm{s}}(2n, \beta_0)^2}$$
together with $P_{\mathrm{s}}(2n, \beta_0) \le 1$ yields the coupled spatial-lift map at fixed
inverse temperature $\beta_0$:
$$s' = f(s) \ge 1 - P_{\mathrm{s}}(12, \beta_0), \qquad
s_2' = f(s_2) \ge 1 - P_{\mathrm{s}}(12, 2\beta_0), \qquad
t' = 1 - (1 - s_2')(1 - t)^2 \ge 1 - P_\tau(24, \beta_0).$$
Iterating this map (`seed_ladder`) generates certified seed pairs $(s^{(k)}, t^{(k)})$ at
$n_0 = 12, 24, 48, 96$ from the $L = 6, 12$ spectra alone, which are then fed into the standard
space-and-temperature doubling recurrence `certified_gap_bound`.

## 6. Centered Weyl parameter transport in $D$ (Stage 2)

Writing $H_L(D) = H_L(D^*) + \delta D \cdot V$ with $V = \sum_{j=1}^L (S_j^z)^2$ and
$0 \le V \le L I$, centering $V$ by $\frac{L}{2}I$ gives
$\|\delta D (V - \frac{L}{2}I)\|_2 \le \frac{1}{2}L|\delta D|$. By Weyl's inequality, every ordered
eigenvalue of $H_L(D) - \frac{L\,\delta D}{2}I$ lies within $\frac{1}{2}L|\delta D|$ of the
corresponding eigenvalue of $H_L(D^*)$, so $\log Z_L(\beta)$ of the shifted Hamiltonian lies within
$\frac{1}{2}\beta L |\delta D|$ of the reference enclosure at $D^*$ (and the constant energy shift
cancels identically in $P_{\mathrm{s}}$ and $P_\tau$). For a compact interval $[D_{\mathrm{lo}}, D_{\mathrm{hi}}]$,
we replace $|\delta D|$ by $\max_{D \in [D_{\mathrm{lo}}, D_{\mathrm{hi}}]} |D - D^*|$
(`transported_logz`).

## 7. Partition tails and finite-ring gap enclosures (Stage 2)

To accelerate repeated $\log Z_{12}$ evaluations in Stage 2, the verifier sums the lowest
$K = 8{,}192$ sorted levels explicitly in interval arithmetic and bounds the remaining
$3^{12} - K = 523{,}249$ levels using their certified minimum energy: since every level $i > K$
satisfies $\lambda_i \ge E_{K+1} - \max_{i > K} r_i$, their shifted Boltzmann sum is enclosed in
$[0,\; (3^{12} - K)\exp(-\beta(E_{K+1} - r_{\max} - E_{\mathrm{shift}}))]$. Because the entire
531,441-state spectrum is enclosed, this tail bound is completely rigorous.

Similarly, to enclose a finite-ring spectral gap $\Delta_L = E_1(L) - E_0(L)$, the verifier applies
the global maximum radius across all blocks to the two lowest sorted eigenvalue centers (plus
$L|\delta D|$ when transported in $D$). At $D = 4$, verified finite-ring gaps at $L = 6, 8, 10$
(`bridge_certificate`) exceed the $n_0 = 12$ uniform bound $0.899143399$, extending the certificate
to all even $L \ge 6$.
