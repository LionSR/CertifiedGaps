# Two-leg transverse-field Ising ladder: route analysis

## Candidate routes

1. **Closed forms / free fermions.** Unlike the single-chain transverse-field Ising model, the two-leg ladder is interacting after a Jordan–Wigner transformation because the rung $ZZ$ couplings become four-fermion terms, so no closed-form free-fermion partition function exists.
2. **Large-field cluster expansion.** Perturbing around the polarized product state at large $g$ could bound $\log Z_{2n} - 2\log Z_n$ for arbitrary $n$, but requires a proved Kotecký–Preiss remainder bound tailored to the two-leg lattice.
3. **Symmetry-resolved complete spectra (chosen).** We decompose the Hilbert space using global $\mathbb{Z}_2$ spin-flip parity and the commuting spatial involutions (leg exchange, rung reflection, and half-ring translation for even $L$). This reduces each sector to a real symmetric matrix whose exact entries are integer multiples of square roots of rationals, with block dimensions summing to $4^L$. Candidate seeds are screened numerically over $\beta$, and the winning even seed is certified using a polar-decomposition and Weyl eigenvalue enclosure with explicit binary64 BLAS rounding bounds and outward-rounded MPFR partition sums. Odd $n_0$ values are recorded only as diagnostics because the single-cell reflection theorem requires even $n_0$.

## Outcome

Symmetry-resolved exact diagonalization certified $g = 3$ at $(n_0, \beta_0) = (4, 3/2)$ using the $4$- and $8$-rung rings ($8$ and $16$ spins), without needing the $10$-rung ring ($4^{10} = 1{,}048{,}576$ states). On the tested grid, $g = 2$ and $g = 5/2$ have no admissible seed with $n_0 \le 4$; the finite-size extrapolation in [`METHOD.md`](METHOD.md) estimates their first admissible seeds at $n_0 \approx 10$ and $n_0 \approx 6$ rungs ($20$- and $12$-rung doubled rings), which lie beyond dense exact diagonalization and would require cluster-expansion or transfer-operator methods.

Later, the local symmetry blocks (Z₂ parity × leg exchange × reflection × half translation),
the local eigenvalue enclosure and the local log Z code were replaced by the shared core
pipeline: the ladder is a `RingModel` of rungs in the Hadamard-rotated frame (RP with
$C = X \otimes 1, 1 \otimes X$; conserved $\mathbb{Z}_2$ charge), certified by
`GapCertificates.certify_gap`. The core splits by charge, momentum and leg-exchange parity
(rung reflection is not available), so the 8-rung blocks have about half the old dimension but
are mostly complex; the seed is unchanged and the tighter core enclosure lowers both defects.
