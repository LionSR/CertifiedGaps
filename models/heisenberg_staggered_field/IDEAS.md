# Candidate routes: spin-1/2 Heisenberg chain in a staggered field

Although the staggered longitudinal field $-h \sum_j (-1)^j S_j^z$ has a two-site spatial period in
the physical basis, rotating every odd spin by $\pi$ about the $y$-axis ($U = 2i S^y$) makes the
Hamiltonian translation-invariant with a single-spin unit cell on any even ring $L$. Consequently,
the seed length $n_0$ must be an even integer and counts individual spins. In the physical basis,
total magnetization $S^z$ and translation by two sites commute with $H_L$, decomposing the
$2^L$-dimensional Hilbert space into $(S^z, k)$ symmetry blocks.

1. **Verified dense diagonalization in $(S^z, k)$ symmetry blocks:** The largest symmetry block
   at $S^z = 0$ has dimension $\approx \binom{L}{L/2}/(L/2)$ ($20$ at $L = 8$ and $1{,}620$ at
   $L = 16$), with $O(d^3)$ eigensolve cost and $O(d^2)$ storage per block. This easily covers
   $n_0 = 8$ ($2n_0 = 16$ sites) across all five staggered fields $`h \in \{1, 3/4, 1/2, 7/20, 1/4\}`$,
   using interval enclosures for the momentum-block matrix entries and a priori IEEE binary64
   rounding bounds on the eigenvector residual and orthogonality defects. At $n_0 = 10$ ($2n_0 = 20$
   sites), the largest block reaches dimension $18{,}504$, making full dense verification
   memory-intensive. **Chosen route for $n_0 = 8$.**
2. **Large-field perturbation around the Néel product state:** Perturbation theory around $h = \infty$
   scales independently of chain length $L$, but $h \le 1$ is not in the weak-coupling regime
   relative to the unit Heisenberg exchange.
3. **Convergent cluster expansion:** A low-temperature polymer expansion around the staggered Néel
   state could bound $\log Z_L(\beta)$ uniformly in $L$ at large $h$, but deriving explicit
   rigorous convergence thresholds down to $h = 1/4$ is unnecessary when $n_0 = 8$ already lies
   well inside the admissible region.
4. **Partial Krylov spectra (`KrylovKit`):** Useful for numerical seed exploration at $L = 20$
   ($n_0 = 10$) and weaker fields ($`h \in \{1/5, 3/20, 1/10\}`$), but cannot certify partition
   functions without a rigorous lower bound on all omitted eigenvalues.

## Outcome

Route 1 certified all five target fields $`h \in \{1, 3/4, 1/2, 7/20, 1/4\}`$ at seed length
$n_0 = 8$ from complete verified 8- and 16-site spectra. Exploratory 20-site scans ($n_0 = 10$)
using partial Krylov spectra for blocks of dimension $> 2{,}500$ identified stronger numerical
candidate seeds down to $h = 3/20$, while an independent shared-core cross-check
(`scripts/core_certificate.jl`) reproduced the $h = 1$ certificate using a different symmetry-block
reduction and the exact rational doubling chain.
