# Transverse-field Ising model with longitudinal field: route analysis

Each unit cell is a single spin-1/2 site. The single-cell reflection theorem requires an even seed length $n_0$; odd $n_0$ values serve only as diagnostics.

## Candidate routes

1. **Free fermions.** Only the $h = 0$ limit is solvable by Jordan–Wigner transformation into free fermions. For $h > 0$, the longitudinal field breaks $\mathbb{Z}_2$ fermion parity and renders the chain non-integrable (giving $E_8$ meson confinement in the ferromagnetic phase), so the $h = 0$ formula is used solely as a numerical cross-check.
2. **Cluster expansions or transfer-operator contraction.** High-field or strong-coupling expansions could reach large $n_0$, but require proved remainder bounds not yet implemented in the core library. Since a longitudinal field opens a substantial gap across the phase diagram, small rings are expected to suffice.
3. **Translation-resolved complete spectra with verified enclosures (chosen).** For $h > 0$, $\mathbb{Z}_2$ spin-flip parity is broken, and spatial reflection exchanges momentum $k$ with $-k$ (giving an independent parity only at $k = 0, \pi$). We therefore block-diagonalize $H_L$ into translation momentum sectors of dimension $\approx 2^L / L$ and enclose the complete spectrum of each block using a polar-decomposition and Weyl perturbation bound with explicit IEEE binary64 rounding-error allowances.

## Outcome

The even seed length $n_0 = 6$ (requiring complete spectra of $6$- and $12$-site rings, with maximum sector dimension $352$ at $L = 12$) produced admissible seeds at all 30 grid points $(g, h)$, certifying uniform spectral gaps across the entire parameter grid without needing larger rings (see [`README.md`](README.md)).
