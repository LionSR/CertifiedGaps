# Two-leg transverse-field Ising ladder: route analysis

## Candidate routes

1. **Closed forms / free fermions.** Unlike the single-chain transverse-field Ising model, the two-leg ladder is interacting after a Jordan–Wigner transformation because the rung $ZZ$ couplings become four-fermion terms, so no closed-form free-fermion partition function exists.
2. **Large-field cluster expansion.** Perturbing around the polarized product state at large $g$ could bound $\log Z_{2n} - 2\log Z_n$ for arbitrary $n$, but requires a proved Kotecký–Preiss remainder bound tailored to the two-leg lattice.
3. **Symmetry-resolved complete spectra (chosen).** We decompose the Hilbert space using global $\mathbb{Z}_2$ spin-flip parity and the commuting spatial involutions (leg exchange, rung reflection, and half-ring translation for even $L$). This reduces each sector to a real symmetric matrix whose exact entries are integer multiples of square roots of rationals, with block dimensions summing to $4^L$. Candidate seeds are screened numerically over $\beta$, and the winning even seed is certified using a polar-decomposition and Weyl eigenvalue enclosure with explicit binary64 BLAS rounding bounds and outward-rounded MPFR partition sums. Odd $n_0$ values are recorded only as diagnostics because the single-cell reflection theorem requires even $n_0$.
4. **Half-length seed from a partial spectrum of the doubled ring (follow-up).** The criterion
   needs only $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ and
   $\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$. By the coupling identity,

   ```math
   P_\tau(2n_0, \beta) = \frac{Z_{2n_0}(2\beta)}{Z_{2n_0}(\beta)^2}
     = \frac{Z_{2n_0}(2\beta)/Z_{n_0}(\beta)^4}{P_{\mathrm{s}}(n_0, \beta)^2}
     \ge \frac{Z_{2n_0}(2\beta)}{Z_{n_0}(\beta)^4}
   ```

   since $P_{\mathrm{s}}(n_0, \beta) \le 1$ (Cauchy–Schwarz for the self-adjoint transfer operator
   at even $n_0$). This equals $`P_{\mathrm{s}}(n_0, 2\beta)\,P_\tau(n_0, \beta)^2`$, so

   ```math
   \varepsilon_{\mathrm{s}} = 1 - \frac{Z_{2n_0}(\beta)}{Z_{n_0}(\beta)^2}, \qquad
   \varepsilon_\tau = 1 - \frac{Z_{2n_0}(2\beta)}{Z_{n_0}(\beta)^4}
   ```

   need the complete spectrum of the $n_0$-rung ring only, and only **lower** bounds on
   $Z_{2n_0}$ at $\beta$ and $2\beta$. Rayleigh–Ritz values are rigorous upper bounds on the
   lowest eigenvalues of each symmetry block (Courant–Fischer), so
   $`Z_{2n_0}(\beta) \ge \sum_B m_B \sum_{i \le k_B} e^{-\beta \theta_i^{\mathrm{up}}(B)}`$, and
   omitting states or blocks only loses weight. In the paramagnet the Boltzmann weight at the
   seed temperature sits in the ground state and the one- and two-quasiparticle states, so a few
   low states per block should suffice. Cost: $n_0 = 6$ needs the complete 6-rung ring
   ($4^6 = 4096$ states) and the low states of the 12-rung ring (24 spins,
   $4^{12} \approx 1.7 \times 10^7$ states, sectors of about $3.5 \times 10^5$ states under
   $\mathbb{Z}_2$ parity, momentum and leg exchange; sparse Lanczos, seconds per matrix-vector
   product). $n_0 = 8$ needs the 16-rung ring (32 spins, sectors of about $6.7 \times 10^7$
   states, so every Lanczos vector takes 0.5–1 GB and every matrix-vector product
   $\sim 2 \times 10^9$ symmetry lookups), which is beyond a 4 GB, 2-thread budget. The rigorous
   Ritz bounds (explicit rounding bounds on the Rayleigh-quotient and Gram matrices) belong in
   the shared core; `scripts/ritz_scan.jl` first screens the route in floating point.

## Outcome

Symmetry-resolved exact diagonalization certified $g = 3$ at $(n_0, \beta_0) = (4, 3/2)$ using the $4$- and $8$-rung rings ($8$ and $16$ spins), without needing the $10$-rung ring ($4^{10} = 1{,}048{,}576$ states). On the tested grid, $g = 2$ and $g = 5/2$ have no admissible seed with $n_0 \le 4$; the finite-size extrapolation in [`METHOD.md`](METHOD.md) estimates their first admissible seeds at $n_0 \approx 10$ and $n_0 \approx 6$ rungs ($20$- and $12$-rung doubled rings), which lie beyond dense exact diagonalization and would require cluster-expansion or transfer-operator methods.

Later, the local symmetry blocks (Z₂ parity × leg exchange × reflection × half translation),
the local eigenvalue enclosure and the local log Z code were replaced by the shared core
pipeline: the ladder is a `RingModel` of rungs in the Hadamard-rotated frame (RP with
$C = X \otimes 1, 1 \otimes X$; conserved $\mathbb{Z}_2$ charge), certified by
`GapCertificates.certify_point`. The core splits by charge, momentum and leg-exchange parity
(not by reflection parity), so the 8-rung blocks have about half the old dimension, and the
declared lattice reflection makes every block real. The seed is unchanged and the tighter core
enclosure lowers both defects.

### Half-length seeds at $n_0 = 6$ (route 4)

Floating-point screening (`scripts/ritz_scan.jl`, not certified): complete spectrum of the
6-rung ring, Lanczos values of the 28 sectors (charge × momentum $0 \le m \le 6$ × leg
exchange) of the 12-rung ring, coupled defects
$`\varepsilon_{\mathrm{s}} = 1 - Z_{12}(\beta)/Z_6(\beta)^2`$ and
$`\varepsilon_\tau = 1 - Z_{12}(2\beta)/Z_6(\beta)^4`$, and `certified_gap_bound` on a
$\beta$ grid of $1/8$ (refined to $1/32$ near the optimum):

| $g$ | best $\beta_0$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | float bound | previous |
|---|---:|---:|---:|---:|---:|
| 3 | 29/16 | 0.02418 | 0.07593 | 0.8494 | 0.6369 certified ($n_0 = 4$) |
| 4 | 11/8 | 0.00283 | 0.01286 | 2.5466 | 2.142 float ($n_0 = 4$), no record |
| 11/4 | 69/32 | 0.04986 | 0.12648 | 0.3712 | no record |
| 5/2 | — | 0.106 at $\beta = 5/2$ | 0.243 (minimum over $\beta$) | no seed | no seed |

- **Few states suffice.** Keeping only the lowest $k$ levels of every sector changes the bound
  by less than $10^{-4}$ already at $k = 1$ and by less than $10^{-8}$ at $k = 4$ (even
  charge) and $k = 2$ (odd charge): at the seed temperature the weight of $Z_{12}$ sits in the
  ground state and the bottom of the one-quasiparticle band, one level per sector. The
  certificates bound 6 and 3 levels per block.
- **$g = 5/2$ fails because of the coupling penalty, not the ring length.** At
  $\beta = 5/2$, $\varepsilon_{\mathrm{s}} = 0.106$, and the coupled thermal defect is
  $`1 - (1 - \varepsilon_\tau^{\mathrm{direct}})(1 - \varepsilon_{\mathrm{s}})^2 = 0.243`$, so
  the direct $`1 - P_\tau(12, 5/2) \approx 0.054`$ would be admissible. The direct defect needs an
  upper bound on $Z_{12}(\beta)$, i.e. the complete 24-spin spectrum, which Rayleigh–Ritz cannot
  give; the bound $Z_{12} \le Z_6^2$ used by the coupling is exactly what costs the factor
  $`P_{\mathrm{s}}(6, \beta)^2`$. $g = 2$ (larger $\xi$) is further away and was not rescanned.
- **$n_0 = 8$** (16-rung doubled ring, 32 spins) is out of reach of this route on the
  available budget: sectors of $\approx 6.7 \times 10^7$ states, so each Lanczos vector takes
  0.5 GB (real) and the core's dictionary-based block construction, which already peaks near
  3.7 GB for the 12-rung ring, would need about 200 times more.

**Certified (route 4).** With the core `RitzBounds` (rigorous Rayleigh–Ritz upper bounds on
6 levels per charge-0 block and 3 per charge-1 block of the 12-rung ring, from the sparse
interval-enclosed core blocks) and `coupled_seed_bounds`, plus the verified gap of the 4-rung
ring for $L = 4$: $\Delta_L \ge 0.371180737$ at $g = 11/4$, $0.849386526$ at $g = 3$ (was
$0.636870849$) and $2.546626505$ at $g = 4$, every even $L \ge 4$. Each point takes about
35–41 min and 4.2 GB, dominated by building the 28 sparse blocks of $3.5 \times 10^5$ states.
