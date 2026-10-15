# Spin-1 dimerized ring: uniform-gap certificate

The smallest certified dimerization on the five-point grid $\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}$
within a 6-cell ($12$-spin) ring ceiling is **$\delta = 7/10$**. At seed $(n_0, \beta_0) = (3\text{ cells}, 17/4)$,
`certified_gap_bound` gives `0.2187941698399417`, establishing a unique ground state and uniform
spectral gap **$\Delta_L \ge 0.218794169839$** for every even chain length $L \ge 3$ cells
($L = 4, 6, 8, \dots$ cells, i.e., physical spin counts $N = 8, 12, 16, \dots$), from the exact
rational seed defects $\varepsilon_{\mathrm{s}} = \text{58597799297}//\text{500000000000}$ and
$\varepsilon_\tau = \text{48363543279}//\text{1000000000000}$ (`results/certificate_d7.json`). At
$\delta = 9/10$, the $n_0 = 2$ cell seed $(2, 7/2)$ certifies **$\Delta_L \ge 0.528801949327$** for
every even $L \ge 2$ cells (`results/certificate_d9.json`).

At $\delta = 7/10$, the certified bound equals $17.28\%$ of the 12-spin numerical reference gap
($1.2659290109040136$), and the full 12-spin verified eigensolve completed in $1017.27\text{ s}$
($16.95\text{ min}$, or $949\text{ s}$ in the computational replay). Across all five grid points,
numerical scans were completed through $L = 6$ cells ($12$ spins); at $\delta = 3/5$ and $1/2$,
admissible seeds require larger rings than $6$ cells.

## Model, lengths, and scope

We consider the periodic spin-1 chain
$$H_L = \sum_{j=1}^{2L} \left[1 + (-1)^j \delta\right] \mathbf{S}_j \cdot \mathbf{S}_{j+1}, \qquad \mathbf{S}_j^2 = 2I,$$
with two spins per unit cell (local cell dimension $3^2 = 9$) and
$\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}$. A ring of $L$ cells contains $N = 2L$ spins. Computing
ring spectra at $L \in \{2, 3, 4, 6\}$ cells supports seed sizes $n_0 = 2$ (using $L = 2, 4$) and
$n_0 = 3$ (using $L = 3, 6$), since certifying $P_\tau(2n_0, \beta_0)$ requires the doubled ring
$2n_0$.

Candidate inverse temperatures $\beta_0 \in \{1/4, 1/2, \dots, 32\}$ are screened numerically using
`Seeds.scan_seeds` and `Seeds.best_seed`, and the selected $(n_0, \beta_0)$ seeds at $\delta = 9/10$
and $\delta = 7/10$ are then certified with verified full-spectrum eigenvalue enclosures and exact
rational defect bounds in `results/certificate_d9.json` and `results/certificate_d7.json`.

## Exact reflection positivity and positive cell-transfer operator

In the Cartesian spin-1 basis, define the real, integer antisymmetric matrices $G_\alpha = i S^\alpha$
($S^\alpha = -i G_\alpha$):
$$G_x = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & -1 & 0 \end{pmatrix}, \qquad
G_y = \begin{pmatrix} 0 & 0 & -1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \qquad
G_z = \begin{pmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.$$
The module `src/rp.jl` verifies in exact integer arithmetic that $[S^x, S^y] = i S^z$,
$\sum_\alpha (S^\alpha)^2 = 2I$, and that for any coupling $c \in \{1 - \delta, 1 + \delta\} > 0$,
$$-c\,\mathbf{S}_i \cdot \mathbf{S}_j = c \sum_{\alpha \in \{x,y,z\}} G_\alpha \otimes G_\alpha.$$
Inside a two-spin cell $(a, b)$ joined by the weak bond $1 - \delta$, the on-site cell operator
$o = -(1 - \delta)\sum_\alpha G_\alpha \otimes G_\alpha$ is real symmetric. Let $F = P$ be the
two-site swap operator inside the cell and $\Theta(\cdot) = F \overline{(\cdot)} F$. Setting
$A_\alpha = I \otimes G_\alpha$ and $B_\alpha = G_\alpha \otimes I$, the strong inter-cell bond
satisfies the exact $81 \times 81$ identity
$$-h_{\mathrm{inter}} = (1 + \delta) \sum_{\alpha \in \{x,y,z\}} A_\alpha \otimes \Theta(A_\alpha), \qquad
\Theta(A_\alpha) = F \overline{A_\alpha} F = B_\alpha, \qquad \Theta(o) = o.$$

### Why the odd seed $n_0 = 3$ cells is valid ($X_\beta = K K^* \ge 0$)

For a general reflection-positive ring, the spatial transfer operator $X_\beta$ is self-adjoint (so
its eigenvalues are real, and $X_\beta^{n_0} \ge 0$ requires $n_0$ to be even). For this dimerized
spin-1 chain in the Cartesian basis, $X_\beta$ is in fact **positive semidefinite** ($X_\beta = K K^* \ge 0$):
- Let $F_{\mathrm{site}}(a, u)$ be the single-site traced Dyson event kernel at a spin-1 site, where
  $a$ and $u$ denote the Poisson event histories on the two incident bonds (one strong, one weak).
  Because both incident bonds insert the exact same real matrices $G_\alpha$ on the site,
  $F_{\mathrm{site}}(a, u)$ is real-valued and symmetric under exchanging the two incident bond
  histories: $F_{\mathrm{site}}(a, u) = F_{\mathrm{site}}(u, a)$.
- Let $\mu_-$ and $\mu_+$ denote the positive Poisson event measures on $[0, \beta]$ for the weak
  ($1 - \delta > 0$) and strong ($1 + \delta > 0$) bonds, respectively. Cutting the ring at the
  strong bonds and integrating out the internal weak-bond history $u$ inside each two-spin cell
  yields the cell transfer kernel on $L^2(\mu_+)$:
  $$X_\beta(a, b) = \int F_{\mathrm{site}}(a, u)\,F_{\mathrm{site}}(u, b)\,d\mu_-(u) = (K K^*)(a, b),$$
  where $K : L^2(\mu_-) \to L^2(\mu_+)$ is the Hilbert–Schmidt integral operator with kernel
  $F_{\mathrm{site}}$.
- Consequently, $X_\beta = K K^*$ is positive trace-class and $Z_L(\beta) = \operatorname{Tr}(X_\beta^L)$
  for every integer $L \ge 1$, with nonnegative transfer eigenvalues $\lambda_i(\beta) \ge 0$.
  Because the only steps in the two-purity doubling proof that use even $n_0$ are the positivity of
  $X_\beta^{n_0}$ (so that $X_\beta^{n_0} / Z_{n_0}(\beta)$ is a density operator) and the
  $\ell^p$-norm moment interpolation $Z_L \le Z_n^{L/n}$, nonnegativity of the eigenvalues of
  $X_\beta$ makes the exact same recurrence and uniform gap bound hold for the odd seed $n_0 = 3$
  cells.

## Complete dense spectra and exact integer dihedral projectors

In the $S^z$ product basis ($m_j \in \{-1, 0, 1\}$), each allowed nearest-neighbor spin exchange
has matrix element equal to the bond coupling $J_j$ (since the two spin-1 ladder factors of
$\sqrt{2}$ cancel the $1/2$ in $\frac{1}{2}(S^+ S^- + S^- S^+)$), while the diagonal entry is
$\sum_j J_j m_j m_{j+1}$. Thus `integer_hamiltonian` in `src/common.jl` constructs an exact
symmetric integer matrix $A$ such that $H_L = A / \mathrm{den}(\delta)$.

The Hamiltonian commutes with total magnetization $M = \sum_j m_j$, translation $T_{\mathrm{cell}}$
by one two-spin cell (two sites), site-order reversal $j \mapsto N + 1 - j$ (which preserves bond
parity), and global spin flip $m_j \mapsto -m_j$:
- Sectors with $M > 0$ are paired with $-M$ with exact multiplicity $2$; at $M = 0$, global
  spin-flip parity is resolved into $\pm 1$ sectors.
- Under cell translation and site reversal (forming the dihedral group $D_L$), generic momentum
  pairs $\{k, L - k\}$ ($0 < k < L/2$) form two-dimensional irreducible representations where the
  two reflection parities are intertwined by $T_{\mathrm{cell}} - T_{\mathrm{cell}}^{-1}$; we
  diagonalize the real reflection-even subspace and assign multiplicity $2$. At $k = 0$ and (for
  even $L$) $k = L/2$, both reflection parities are diagonalized separately.

On small rings ($L \le 4$ cells), the verified certificate path block-diagonalizes $H_L$ in complete
$S^z$ sectors directly. On the 6-cell ($12$-spin) ring, `src/verified_blocks.jl` uses exact
dihedral character projectors on each permutation orbit: because $2\cos(2\pi k r / L)$ is an exact
integer for every $L \in \{2, 3, 4, 6\}$, each orbit projector takes the form $P_{\mathrm{num}} / \mathrm{den}$
with exact integer entries. The verifier checks exact symmetry, $P_{\mathrm{num}}^2 = \mathrm{den}\cdot P_{\mathrm{num}}$,
and integer trace/rank on every orbit, and confirms that the multiplicity-weighted block dimensions
sum to $9^L$ ($531{,}441$ at $L = 6$, across 79 blocks of maximum dimension $11{,}596$).

## Eigenvalue-enclosure theorem and rounding error bounds

Let $H$ be an exact real symmetric block of dimension $d$, $Q$ its approximate `Float64`
eigenvector matrix, and $\Lambda = \operatorname{diag}(\Lambda_1, \dots, \Lambda_d)$ the sorted
approximate eigenvalues. Given verified upper bounds $r \ge \|H Q - Q\Lambda\|_F$ and
$e \ge \|Q^\top Q - I\|_F$ with $e < 1$, write the polar decomposition $Q = U S$ with $U$
orthogonal and $S = (Q^\top Q)^{1/2}$. Setting $\sigma = \sqrt{1 - e}$ and $d_S = e / (1 + \sigma)$,
we have $\|Q^{-1}\|_2 \le 1/\sigma$ and $\|S - I\|_2 \le d_S$, so the identity
$$H - U \Lambda U^\top = (H Q - Q\Lambda) Q^{-1} + U (S \Lambda S^{-1} - \Lambda) U^\top$$
implies
$$\|H - U \Lambda U^\top\|_2 \le \eta := \frac{r + 2\|\Lambda\|_2 d_S}{\sigma}.$$
Weyl's inequality then guarantees that every ordered eigenvalue of $H$ (with multiplicity) satisfies
$|\lambda_i(H) - \Lambda_i| \le \eta$.

In `src/enclosure.jl`, `eigen_enclosure` computes $R_{\mathrm{comp}} = \operatorname{fl}(A Q - \mathrm{den}\cdot Q\Lambda)$
and $E_{\mathrm{comp}} = \operatorname{fl}(Q^\top Q - I)$ in single-threaded `Float64` BLAS and
adds explicit IEEE binary64 rounding bounds ($u = 2^{-53}$, $\gamma_k = k u / (1 - k u)$):
$$r \le \frac{\|R_{\mathrm{comp}}\|_F + \gamma_{2d+8}(\|A\|_F \|Q\|_F + \mathrm{den}\cdot \|Q\|_F \|\Lambda\|_2) + d^3 \mathrm{floatmin}}{\mathrm{den}},$$
$$e \le \|E_{\mathrm{comp}}\|_F + \gamma_{2d+8}(\|Q\|_F^2 + \sqrt{d}) + d^3 \mathrm{floatmin},$$
where each Frobenius norm sum of $m$ squares is enlarged by $(1 - \gamma_{2m+2})^{-1}$ plus
$m\,\mathrm{floatmin}$ before taking a directed-rounded 128-bit `BigFloat` square root.

For the 12-spin symmetry-projected blocks (`src/verified_blocks.jl`), an additional basis-construction
allowance is added to $\eta$: for each orbit projector $P$ and floating-point basis $B$, verified
bounds $q \ge \|PB - B\|_F$ and $e_B \ge \|B^\top B - I\|_F$ bound the distance from $B$ to an
exact orthonormal basis $U_{\mathrm{proj}}$ of $\operatorname{ran}(P)$ by
$$d_B = q + \frac{w}{1 + \sqrt{1 - w}}, \qquad w = e_B + 2\|B\|_F q + q^2 < 1.$$
Taking the maximum $d_B$ over disjoint orbits, replacing $B$ by $U_{\mathrm{proj}}$ shifts
$B^\top H_L B$ in operator norm by at most $\|H_L\|_2 (2d_B + d_B^2) \le 4L (2d_B + d_B^2)$
(since $\|\mathbf{S}_i\cdot\mathbf{S}_{i+1}\|_2 = 2$ and the couplings sum to $2L$), and the sparse
projection multiplication $B^\top A B / \mathrm{den}$ is bounded using
$\gamma_{2\,\mathrm{maxnnz}(A) + 2\,\mathrm{maxnnz}(B) + 12}$. All final eigenvalue radii $\eta$ are
rounded up to exact rational multiples of $10^{-12}$.

## Partition functions, rational seeds, and gap bounds

For each ring, `verified_logz` sums $\exp(-\beta E)$ over all enclosed eigenvalue intervals in
outward-rounded `IntervalArithmetic` (after factoring out an exact binary64 energy shift) weighted
by each block's exact multiplicity, storing rational lower and upper bounds on $\log Z_{n_0}(\beta_0)$,
$\log Z_{2n_0}(\beta_0)$, and $\log Z_{2n_0}(2\beta_0)$. The upper endpoints of the resulting
defect intervals
$$\varepsilon_{\mathrm{s}} \ge 1 - \exp\!\big(\underline{\log Z_{2n_0}(\beta_0)} - 2\,\overline{\log Z_{n_0}(\beta_0)}\big), \qquad
\varepsilon_\tau \ge 1 - \exp\!\big(\underline{\log Z_{2n_0}(2\beta_0)} - 2\,\overline{\log Z_{2n_0}(\beta_0)}\big)$$
are rounded up to exact rationals with denominator $10^{12}$ and passed to
`CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)`.

## Independent audit (2026-10-11)

An independent audit reviewed `src/rp.jl`, `src/common.jl`, `src/enclosure.jl`,
`src/verified_blocks.jl`, `certify.jl`, `check.jl`, and `results/certificate_d7.json`:
1. **Positive transfer operator and odd seed ($n_0 = 3$):** Confirmed that in the Cartesian basis
   $G_\alpha = i S^\alpha$, the single-site Dyson trace kernel $F_{\mathrm{site}}(a, u)$ is real
   and symmetric in its two incident bond histories, so integrating out the weak intra-cell bond
   history expresses the cell transfer operator as $X_\beta = K K^* \ge 0$ on $L^2(\mu_+)$. Because
   all transfer eigenvalues are nonnegative, $X_\beta^{n_0} / Z_{n_0}(\beta)$ is a positive density
   operator and $\ell^p$-norm moment interpolation holds at $n_0 = 3$, justifying the unequal-defect
   recurrence and uniform gap bound for all even $L \ge 3$ cells ($N = 8, 12, 16, \dots$ spins).
2. **Symmetry completeness and enclosure bounds:** Verified that the 7 small-ring $S^z$ blocks
   (weighted dimension $729 = 9^3$, max block $141$) and 79 large-ring dihedral blocks (weighted
   dimension $531{,}441 = 9^6$, max block $11{,}596$) partition the full Hilbert space via disjoint
   exact integer character projectors, and confirmed all polar-factor/Weyl, $\gamma_k$, and
   projector-basis error bounds.
3. **Full computational replay:** A fresh execution of `certify.jl 7` on the renamed `CertifiedGaps`
   core (949 s) reproduced every mathematical field of `results/certificate_d7.json`—all 79 block
   error bounds, $\log Z$ intervals, rational seed defects, and the bound `0.2187941698399417`—bit
   for bit.

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
