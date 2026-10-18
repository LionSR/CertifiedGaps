# Spin-1 dimerized ring: uniform-gap certificate

The smallest certified dimerization on the five-point grid
$`\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}`$ within a 6-cell ($12$-spin) ring ceiling is
**$\delta = 7/10$**. At the odd seed $(n_0, \beta_0) = (3\text{ cells}, 17/4)$ the core pipeline
gives the interval bound **$\Delta_L \ge 0.218793401$** (exact rational chain bound
$\Delta_L > 0.218607065$) and a unique ground state for every even $L \ge 3$ cells, from
$\varepsilon_{\mathrm{s}} = 117195681829365/10^{15}$ and $\varepsilon_\tau = 48363910725858/10^{15}$
(`results/certificate_d7.json`). At $\delta = 9/10$ the seed $(2, 7/2)$ gives
**$\Delta_L \ge 0.528802288$** (chain $0.528349759$) for every even $L \ge 2$ cells
(`results/certificate_d9.json`).

At $\delta = 7/10$ the bound is $`17.28\%`$ of the 12-spin numerical reference gap
($1.2659290109040136$); the core enclosure of the 6-cell ring took $10{,}339$ s on a heavily
shared machine.
Across all five grid points, numerical scans were completed through $L = 6$ cells ($12$
spins); at $\delta = 3/5$ and $1/2$, admissible seeds require larger rings than $6$ cells.

## Model, lengths, and scope

We consider the periodic spin-1 chain
```math
H_L = \sum_{j=1}^{2L} \left[1 + (-1)^j \delta\right] \mathbf{S}_j \cdot \mathbf{S}_{j+1}, \qquad \mathbf{S}_j^2 = 2I,
```
with two spins per unit cell (local cell dimension $3^2 = 9$) and
$`\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}`$. A ring of $L$ cells contains $N = 2L$ spins. Computing
ring spectra at $`L \in \{2, 3, 4, 6\}`$ cells supports seed sizes $n_0 = 2$ (using $L = 2, 4$) and
$n_0 = 3$ (using $L = 3, 6$), since certifying $P_\tau(2n_0, \beta_0)$ requires the doubled ring
$2n_0$.

Candidate inverse temperatures $`\beta_0 \in \{1/4, 1/2, \dots, 32\}`$ are screened numerically using
`Seeds.scan_seeds` and `Seeds.best_seed`, and the selected $(n_0, \beta_0)$ seeds at $\delta = 9/10$
and $\delta = 7/10$ are then certified with verified full-spectrum eigenvalue enclosures and exact
rational defect bounds in `results/certificate_d9.json` and `results/certificate_d7.json`.

## Exact reflection positivity and positive cell-transfer operator

In the Cartesian spin-1 basis, define the real, integer antisymmetric matrices $G_\alpha = i S^\alpha$
($S^\alpha = -i G_\alpha$):
```math
G_x = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & -1 & 0 \end{pmatrix}, \qquad
G_y = \begin{pmatrix} 0 & 0 & -1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \qquad
G_z = \begin{pmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
```
`test/runtests.jl` verifies in exact integer arithmetic that $[S^x, S^y] = i S^z$,
$\sum_\alpha (S^\alpha)^2 = 2I$, and that for any coupling $`c \in \{1 - \delta, 1 + \delta\} > 0`$,
```math
-c\,\mathbf{S}_i \cdot \mathbf{S}_j = c \sum_{\alpha \in \{x,y,z\}} G_\alpha \otimes G_\alpha.
```
Inside a two-spin cell $(a, b)$ joined by the weak bond $1 - \delta$, the on-site cell operator
$o = -(1 - \delta)\sum_\alpha G_\alpha \otimes G_\alpha$ is real symmetric. Let $F = P$ be the
two-site swap operator inside the cell and $\Theta(\cdot) = F \overline{(\cdot)} F$. Setting
$A_\alpha = I \otimes G_\alpha$ and $B_\alpha = G_\alpha \otimes I$, the strong inter-cell bond
satisfies the exact $81 \times 81$ identity
```math
-h_{\mathrm{inter}} = (1 + \delta) \sum_{\alpha \in \{x,y,z\}} A_\alpha \otimes \Theta(A_\alpha), \qquad
\Theta(A_\alpha) = F \overline{A_\alpha} F = B_\alpha, \qquad \Theta(o) = o.
```

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
  ```math
  X_\beta(a, b) = \int F_{\mathrm{site}}(a, u)\,F_{\mathrm{site}}(u, b)\,d\mu_-(u) = (K K^*)(a, b),
  ```
  where $K : L^2(\mu_-) \to L^2(\mu_+)$ is the Hilbert–Schmidt integral operator with kernel
  $F_{\mathrm{site}}$.
- Consequently, $X_\beta = K K^*$ is positive trace-class and $Z_L(\beta) = \mathrm{Tr}(X_\beta^L)$
  for every integer $L \ge 1$, with nonnegative transfer eigenvalues $\lambda_i(\beta) \ge 0$.
  Because the only steps in the two-purity doubling proof that use even $n_0$ are the positivity of
  $X_\beta^{n_0}$ (so that $X_\beta^{n_0} / Z_{n_0}(\beta)$ is a density operator) and the
  $\ell^p$-norm moment interpolation $Z_L \le Z_n^{L/n}$, nonnegativity of the eigenvalues of
  $X_\beta$ makes the exact same recurrence and uniform gap bound hold for the odd seed $n_0 = 3$
  cells.

## Core pipeline

The certificates `results/certificate_d9.json` and `results/certificate_d7.json` are produced by
`GapCertificates.certify_point` (which runs `certify_gap`) on the `RingModel` of
`src/DimerizedCertificate.jl` (`dimerized(δ)`), in the rotated $S^z$ frame
$U = \prod_{\text{cells}} (1 \otimes V)$, $V = e^{-i\pi S^y}$, described in `README.md`:

- **Exact model.** With $X' = \sqrt2 S^x$, $`A' = \sqrt2\,iS^y`$ and $S^z$ integer matrices,

  $`U (J\,\mathbf{S}_i\cdot\mathbf{S}_j) U^\dagger = -J(\tfrac12 X'\otimes X' + \tfrac12 A'\otimes A' + S^z\otimes S^z)`$
  for every bond, since $V$ reverses $S^x$ and $S^z$ and keeps $S^y$. The on-site term is this
  bond with $J = 1 - \delta$ on $(a, b)$, the ring bond is it with $J = 1 + \delta$ on
  $(b, a')$. The core checks exactly that the cell charge $z_a - z_b$ is conserved, that the
  flip $P \otimes P$ reverses it and commutes with $h$ and $o$, and that the lattice reflection
  with on-site part $\mathrm{SWAP}\,(P \otimes P)$ (site reversal of the chain composed with
  the flip, which keeps the charge) maps $h$ to its mirror image and commutes with $o$; the
  tests check the frame change exactly (`is_frame_change` from the physical cell model).
- **RP.** $-h = \sum_\alpha w_\alpha C_\alpha \otimes R C_\alpha R$ with
  $C_\alpha = 1 \otimes (X', A', S^z)$, $w = ((1+\delta)/2, (1+\delta)/2, 1+\delta)$ and $R$ the
  cell swap, verified exactly by `verify_rp`.
- **Spectra.** `verified_ring_spectrum` splits the rings of $n_0$ and $2n_0$ cells exactly into
  charge sectors (only $q \ge 0$, by the flip; $q = 0$ split by flip parity) and cell momenta
  $0 \le k \le L/2$ (the reflection pairs $\pm k$ and makes every block real symmetric) and
  encloses each block with `verified_eigen` (residual and orthogonality in binary64 with a
  priori $\gamma_{N+8}$ / $\gamma_{2N+8}$ bounds and Frobenius norms, then Ostrowski or
  polar-factor and Weyl). The 4-cell ring has 28 blocks (largest 254, radius
  $\le 2.5 \times 10^{-10}$); the 6-cell ring has 53 real blocks, the largest $11{,}596$, with
  radii up to $6.6 \times 10^{-7}$.
- **Seed and chain.** $\log Z$ enclosures (`log_partition_enclosure`), the seed defects rounded
  up to the grid $10^{-15}$, the interval bound `certified_gap_bound` and the exact rational
  doubling chain. For the odd seed $n_0 = 3$ the core requires the `positive_transfer` reference
  to the $X_\beta = K K^* \ge 0$ argument above and stores it in the certificate. The checker's
  `check_point` compares the stored model, RP data, seed and `params` with the point, and
  `check_certificate` reruns the exact parts (fast) and the spectra (`--full`).

## Former local pipeline

Until 2026-10 both certificates came from a local pipeline (removed; in the git history of this
directory). It used the physical $S^z$ basis, where $`\mathrm{den}(\delta)\,H_L`$ is an integer
matrix, and split the 6-cell ring by $S^z$, spin-flip parity, cell momentum and site reversal
into 79 real blocks (maximum dimension $11{,}596$) with exact integer dihedral character
projectors on every orbit. Each block was enclosed by a local polar-factor/Weyl bound with
binary64 rounding terms plus a bound on the floating-point projector basis. Its $\delta = 7/10$
seed was $\varepsilon_{\mathrm{s}} = 58597799297/(5 \cdot 10^{11})$,
$\varepsilon_\tau = 48363543279/10^{12}$ (interval bound $0.2187941698$). The first core
certificate, before the core had lattice reflections, had complex momentum blocks $k = 1, 2$
with radii up to $2.0 \times 10^{-6}$ and the interval bound $0.218793170$; with the
reflection the core blocks are real again and the bound is $0.218793401$.

## Independent audit (2026-10-11)

An independent audit reviewed the pre-core pipeline (`src/rp.jl`, `src/common.jl`,
`src/enclosure.jl`, `src/verified_blocks.jl`, since removed), `certify.jl`, `check.jl`, and
`results/certificate_d7.json`:
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

Verified certificates (the 0.9, n=2 seed was chosen before the larger-ring scan):

| delta | n | exact beta | certified lower bound | ratio to 12-spin ED gap | verification seconds | validity |
|---|---:|---:|---:|---:|---:|---|
| 9//10 | 2 | 7//2 | 0.528349759000 | 0.299744 | 21.46 | even L >= 2 cells; N divisible by 4, N >= 4 spins |
| 7//10 | 3 | 17//4 | 0.218606834000 | 0.172685 | 6520.24 | even L >= 3 cells; N divisible by 4, N >= 8 spins |

The table gives each certificate's exact rational chain bound `gap_lower_bound`; its interval bound `certified_gap_bound` is `interval_bound_decimal` in the certificate (and in README.md). Table bounds are rounded down to twelve decimal places.

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
