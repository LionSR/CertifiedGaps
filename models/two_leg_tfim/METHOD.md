# Two-leg transverse-field Ising ladder: verification method and results

**Headline result:** $g = 3$ is the smallest certified transverse field on the tested grid $`g \in \{2, 5/2, 3, 4\}`$. The seed $(n_0, \beta_0) = (4\text{ rungs}, 3/2)$ certifies $\Delta_L \ge 0.636386863$ and a unique ground state for every even $L \ge 4$ rungs ($N = 2L$ spins, $N$ divisible by 4, $N \ge 8$). Fields $g = 2$ and $5/2$ have no admissible seed with $n_0 \le 4$ rungs, and $g = 4$ passes the numerical scan but was not run through full verification since $g = 3$ already succeeds.

## Model and scope

Each unit cell is a two-spin rung ($d = 4$). In Pauli normalization on a periodic ring of $L$ rungs,

```math
H_L = -\sum_{j=1}^L (Z_{j,a} Z_{j+1,a} + Z_{j,b} Z_{j+1,b} + Z_{j,a} Z_{j,b}) - g \sum_{j=1}^L (X_{j,a} + X_{j,b}).
```

The $L = 2$ ring includes both periodic leg bonds, following the core package convention. For every $g > 0$, the off-diagonal terms in the $\sigma^z$ basis form a connected hypercube graph with negative edge weights $-g$, so Perron–Frobenius applied to $cI - H_L$ guarantees a unique ground state on every finite ring (while the two-purity doubling criterion establishes the uniform infinite-volume gap and ground-state uniqueness).

In this Pauli normalization, DMRG finite-size scaling by Xavier, Pereira, Nunes, and Plascak ([arXiv:2111.11404](https://arxiv.org/abs/2111.11404), Eq. (1) and Table I with $J_\perp^x = 0$, $J_\perp^z = 1$, $J_4 = 0$) locates the quantum critical point at $g_c = 1.83214(5)$, placing all four tested fields $`g \in \{2, 5/2, 3, 4\}`$ in the gapped paramagnetic phase. This literature value is used only for physical context and never enters the certificate.

## Exact reflection positivity and symmetry decomposition

In the real $\sigma^z$ cell basis $(\uparrow\uparrow, \uparrow\downarrow, \downarrow\uparrow, \downarrow\downarrow)$ with $Z = \mathrm{diag}(1, -1)$ and $`X = \bigl(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\bigr)`$, setting $C_1 = Z \otimes I$ and $C_2 = I \otimes Z$ gives

```math
-h_{\text{bond}} = C_1 \otimes C_1 + C_2 \otimes C_2,
```

with real symmetric on-site rung operator $o = -Z \otimes Z - g(X \otimes I + I \otimes X)$. The function `rp_identity` in [`src/Ladder.jl`](src/Ladder.jl) checks this four-spin matrix identity over `Rational{BigInt}`, constructing the inter-rung bond independently from four $2 \times 2$ Kronecker factors. Because a single rung cell is invariant under spatial reflection along the leg direction, $F = I$.

To block-diagonalize $H_L$ for exact diagonalization, a real Hadamard rotation on every spin diagonalizes the transverse field so that $`\mathbb{Z}_2`$ spin-flip parity (the parity of the number of down spins in the rotated basis) is conserved, while the $`ZZ`$ couplings act as pair flips. Leg exchange $`S`$, rung reflection $`R: j \mapsto -j`$, and (for even $`L`$) half-ring translation $`T: j \mapsto j + L/2`$ are commuting involutions that commute with $`H_L`$. For each spin-flip parity and real character $`\chi`$ of the involution group, each permutation orbit $`O`$ with a $`\chi`$-trivial stabilizer contributes one normalized basis state $`|O, \chi\rangle = |O|^{-1/2} \sum_{s \in O} \chi(s) |s\rangle`$, and orbits with non-trivial stabilizer characters vanish. Because distinct orbits have disjoint support and characters are orthogonal, this produces an exact orthonormal decomposition into real symmetric blocks whose dimensions sum to $4^L$.

Exact block entries have the algebraic form $A_{ij} / (2\sqrt{m_i / m_j})$, where $A$ is an integer matrix and $m_i = |O_i|$ are orbit sizes; exact symmetry is verified via the integer identity $A_{ij} m_j = A_{ji} m_i$. Small-ring tests check the block spectra against the full $4^L \times 4^L$ Hamiltonian.

## Implemented enclosure theorem

Let $B$ be an exact real symmetric block of dimension $m$, $H$ its `Float64` approximation, $Q$ the computed eigenvector matrix, and $\Lambda = \mathrm{diag}(\Lambda_1, \dots, \Lambda_m)$ the computed eigenvalues. Suppose we have certified upper bounds

```math
\|B - H\|_2 \le d, \qquad \|HQ - Q\Lambda\|_2 \le r, \qquad \|Q^\top Q - I\|_2 \le e < 1.
```

Then the orthogonal polar factor $U = Q(Q^\top Q)^{-1/2}$ satisfies $`\|U - Q\|_2 \le a = e / (1 + \sqrt{1 - e})`$, which implies

```math
\|U^\top H U - \Lambda\|_2 \le r + (\|H\|_2 + \|\Lambda\|_2)\, a.
```

By Weyl's inequality and the triangle inequality, every ordered eigenvalue of $B$ (counted with multiplicity) satisfies

```math
|\lambda_i(B) - \Lambda_i| \le \eta = d + r + (\|H\|_F + \|\Lambda\|_2)\, a.
```

In the implementation:
- Each algebraic entry of $B$ is enclosed using 256-bit directed-rounding MPFR arithmetic, and $`d \ge \|B - H\|_2`$ is bounded by the maximum absolute row sum of $B - H$.
- Under the IEEE binary64 floating-point model (unit roundoff $u = 2^{-53}$, $\gamma_k = ku / (1 - ku)$ with $k = 2m + 8$, plus an additive $16m^3 \cdot 2^{-1074}$ allowance for gradual underflow), the rounding error of the floating-point matrix products $HQ - Q\Lambda$ and $Q^\top Q - I$ is bounded by $`\gamma_k(\|H\|_F \|Q\|_F + \|Q\|_F \|\Lambda\|_2)`$ and $`\gamma_k(\|Q\|_F^2 + \sqrt{m})`$, respectively. All Frobenius norms and combined upper bounds are accumulated in `BigFloat` with upward rounding.
- Directed-rounding MPFR exponentials, sums, and logarithms enclose $\log Z_n(\beta_0)$, $\log Z_{2n}(\beta_0)$, and $\log Z_{2n}(2\beta_0)$ from the eigenvalue intervals $[\Lambda_i - \eta, \Lambda_i + \eta]$. From the lower bounds on the log-purities $l_{\mathrm{s}} = \mathrm{lo}(\log Z_{2n}(\beta_0)) - 2\mathrm{hi}(\log Z_n(\beta_0))$ and $l_\tau = \mathrm{lo}(\log Z_{2n}(2\beta_0)) - 2\mathrm{hi}(\log Z_{2n}(\beta_0))$, upward-rounded values of $-\mathrm{expm1}(l_{\mathrm{s}})$ and $-\mathrm{expm1}(l_\tau)$ are converted exactly to `Rational{BigInt}` and passed to `CertifiedGaps.certified_gap_bound`.

## Validation and independent audit

The certificate [`results/certificate_g3_1.json.gz`](results/certificate_g3_1.json.gz) (uncompressed SHA-256 `aba98d9865f1b496a96e7a7e17dec0239ec07b34e012b73714d373c1e8a74928`) was generated on Julia 1.12.4 in 645.83 s, with maximum eigenvalue enclosure radii of $1.50 \times 10^{-10}$ ($L = 4$ rungs) and $3.94 \times 10^{-5}$ ($L = 8$ rungs), yielding $\varepsilon_{\mathrm{s}} \approx 0.11240836635$, $\varepsilon_\tau \approx 0.05161523886$, and $\Delta_L \ge 0.6363868632722607$.

- **Fast checker (`check.jl` / `test/runtests.jl`):** runs 70 tests in ~35 s (10 exact RP and symmetry checks, 58 stored-seed, recurrence, and provenance checks, and 2 enclosure checks on an analytic test matrix), recomputing $\log Z$ and the rational seed defects from the stored eigenvalues and radii.
- **Full spectral replay (`check.jl --full` and `scripts/replay_audit.jl`):** rebuilds all 32 dense symmetry blocks at $4$ and $8$ rungs ($256$ and $65{,}536$ states) from the Hamiltonian, recomputing the residual and orthogonality enclosures, partition intervals, rational seed defects, and recurrence bound from scratch (~445 s). On the original Linux OpenBLAS environment, the replay matched every stored floating-point and rational number bit for bit (`scripts/check_replay.jl`, 262 tests in 5.4 s). On an Apple M5 (Julia 1.12.4), platform differences in LAPACK eigenvector tie-breaking yield slightly tighter enclosure radii and a valid gap bound of $0.6364738604216351$, while all seed, bound, and source-hash checks pass.

## Results

All scan values below are **numerical, not certified**. The grid is $\beta = 1/8, 2/8, \dots, 32$; $n_0 = 2, 3, 4$.

| g | best eligible n0 | β | numerical bound | ED gap (8 rungs) | spectral runtime (s) |
|---|---:|---:|---:|---:|---:|
| 2//1 | — | — | no seed | 0.528639460 | 89.62 |
| 5//2 | — | — | no seed | 1.579207409 | 134.27 |
| 3//1 | 4 | 1.500 | 0.636872819 | 2.661781588 | 90.40 |
| 4//1 | 4 | 1.125 | 2.142199085 | 4.763362747 | 136.61 |

Spectral runtimes sum the measured dense-spectrum calls after orbit construction; they exclude Julia startup, the $\beta$-grid scan, and JSON I/O. Verification runtime includes spectra, residual validation, and seed construction, excluding startup and final JSON I/O.

### Purity and ground-energy finite-size diagnostics

| g | n | max_β min(Ps,Pτ) | δ(n)=E0(2n)−2E0(n) |
|---|---:|---:|---:|
| 2//1 | 2 | 0.454349841 | 0.753150185 |
| 2//1 | 3 | 0.526080292 | 0.393739787 |
| 2//1 | 4 | 0.564897524 | 0.262460357 |
| 5//2 | 2 | 0.532084570 | 0.717774565 |
| 5//2 | 3 | 0.683199040 | 0.292564344 |
| 5//2 | 4 | 0.780193496 | 0.142802875 |
| 3//1 | 2 | 0.602242471 | 0.634477163 |
| 3//1 | 3 | 0.797921242 | 0.200727919 |
| 3//1 | 4 | 0.893759824 | 0.073955799 |
| 4//1 | 2 | 0.722816094 | 0.492692762 |
| 4//1 | 3 | 0.905693809 | 0.106397922 |
| 4//1 | 4 | 0.969496099 | 0.026210081 |

For the fields without an $n_0 \le 4$ seed, a two-point fit $\delta(n) = \delta(4) \exp(-(n-4)/\xi_\delta)$ to $n = 3, 4$, combined with a dilute-excitation extrapolation of the 8-rung partition function, provides heuristic estimates of the first admissible even seed (`scripts/extrapolate.jl`):

- `g = 2//1`: $\xi_\delta = 2.465543$ rungs; estimated $(n_0, \beta_0) = (10, 8.250)$, requiring a 20-rung doubled ring.
- `g = 5//2`: $\xi_\delta = 1.394273$ rungs; estimated $(n_0, \beta_0) = (6, 2.500)$, requiring a 12-rung doubled ring.

### Dense sector dimensions

- `L = 2`: 5, 1, 1, 1, 2, 2, 2, 2; sum = $16 = 4^2$. Zero-dimensional sectors are omitted.
- `L = 3`: 14, 6, 6, 6, 10, 10, 6, 6; sum = $64 = 4^3$. Zero-dimensional sectors are omitted.
- `L = 4`: 34, 18, 14, 14, 10, 10, 14, 14, 24, 24, 16, 16, 8, 8, 16, 16; sum = $256 = 4^4$. Zero-dimensional sectors are omitted.
- `L = 6`: 320, 256, 256, 256, 240, 240, 240, 240, 272, 272, 272, 272, 240, 240, 240, 240; sum = $4096 = 4^6$. Zero-dimensional sectors are omitted.
- `L = 8`: 4384, 4128, 4064, 4064, 4000, 4000, 4064, 4064, 4224, 4224, 4096, 4096, 3968, 3968, 4096, 4096; sum = $65536 = 4^8$. Zero-dimensional sectors are omitted.

### Certified results

- **`g = 3//1`: $\Delta_L \ge 0.636386863$** for all even $L \ge 4$ rungs ($2L$ spin sites, multiples of 4, $\ge 8$ sites). Seed $(n_0, \beta_0) = (4, 3/2)$; verification time $645.83\text{ s}$; ratio of certified bound to 8-rung ED gap $= 0.239083$.
  The returned `.bound` is `0.6363868632722607` (exact binary rational `179126977518527//281474976710656`), rounded downward in the headline decimal.
  Exact rational seed defects and per-sector eigenvalue enclosure radii are stored in `results/certificate_g3_1.json.gz`.

