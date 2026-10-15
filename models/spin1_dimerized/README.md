# Dimerized spin-1 Heisenberg ring

**Certified at $\delta = 9/10$ and $\delta = 7/10$:** for every even chain length $L \ge 2$ cells at
$\delta = 9/10$ ($\Delta_L \ge 0.528801949327$) and every even chain length $L \ge 3$ cells
($L = 4, 6, 8, \dots$ cells, i.e., $N = 2L = 8, 12, 16, \dots$ spins) at $\delta = 7/10$
($\Delta_L \ge 0.218794169839$), the ground state is unique and the uniform spectral gap is bounded
below as shown. Because the two-spin unit cell in the Cartesian spin-1 basis factors the cell
transfer operator into the positive semidefinite product $X_\beta = K K^*$, the doubling recurrence
and moment interpolation apply directly to the odd seed length $n_0 = 3$ cells ($6$ and $12$
spins). Full derivations, the positive cell-transfer theorem, the exact orbit projector
construction, numerical scan tables, and the independent audit are in [`METHOD.md`](METHOD.md); the
route memo is in [`IDEAS.md`](IDEAS.md).

| Dimerization $\delta$ | Seed $(n_0\text{ cells}, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $9/10$ | $(2, 7/2)$ | $0.528801949327$ | every even $L \ge 2$ cells ($N \ge 4$ spins, $4 \mid N$) |
| $7/10$ | $(3, 17/4)$ | $0.218794169839$ | every even $L \ge 3$ cells ($N \ge 8$ spins, $4 \mid N$) |

## Model

We consider the periodic dimerized spin-1 Heisenberg chain
$$H_L = \sum_{j=1}^{2L} \left[1 + (-1)^j \delta\right] \mathbf{S}_j \cdot \mathbf{S}_{j+1}, \qquad \mathbf{S}_j^2 = 2I,$$
on a ring of $L$ two-spin unit cells ($N = 2L$ spins, local cell dimension $3^2 = 9$). Each unit
cell consists of two adjacent spin-1 sites joined by the weak intra-cell bond $1 - \delta$, while
neighboring cells are linked by the strong inter-cell bond $1 + \delta$. Lengths $L$ and $n_0$
count two-spin unit cells. The parameter grid studied is $\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}$
on rings of up to $6$ cells ($12$ spins, Hilbert-space dimension $3^{12} = 531{,}441$).

## Reflection positivity

In the Cartesian spin-1 basis, the spin operators take the form $S^\alpha = -i G_\alpha$ with real,
integer antisymmetric matrices $G_\alpha = i S^\alpha$, so that every bond satisfies
$-\mathbf{S}_i \cdot \mathbf{S}_j = \sum_{\alpha \in \{x,y,z\}} G_\alpha \otimes G_\alpha$ with unit
positive weights. Grouping sites into two-spin cells $(a, b)$:
- the intra-cell on-site operator $o = (1 - \delta)\,\mathbf{S}_a \cdot \mathbf{S}_b = -(1 - \delta)\sum_\alpha G_\alpha \otimes G_\alpha$
  is real symmetric and invariant under the internal two-site swap $F = P$;
- the inter-cell bond linking site $b$ of cell $j$ to site $a$ of cell $j+1$ satisfies
  $$-h = (1 + \delta) \sum_{\alpha \in \{x,y,z\}} C_\alpha \otimes F \overline{C_\alpha} F, \qquad C_\alpha = I \otimes G_\alpha,$$
  because $F (I \otimes G_\alpha) F = G_\alpha \otimes I$.

The module `src/rp.jl` verifies all spin-1 algebraic relations and this $81 \times 81$ reflected
unit-cell identity in exact integer arithmetic. Furthermore, because both bonds incident on each
single spin-1 site insert the exact same real matrices $G_\alpha$, cutting the ring at the strong
inter-cell bonds and integrating out the weak intra-cell bond history factors the cell transfer
operator as $X_\beta = K K^* \ge 0$ (proved in [`METHOD.md`](METHOD.md)). Positivity of $X_\beta$
ensures that all transfer eigenvalues are nonnegative, removing the usual even-$n_0$ restriction and
allowing the $n_0 = 3$ cell seed at $\delta = 7/10$.

## Seeds and certified bounds

| $\delta$ | Seed $(n_0\text{ cells}, \beta_0)$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ | $\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | 12-spin ED gap |
|---|---|---|---|---:|---:|
| $9/10$ | $(2, 7/2)$ | `12428916293//250000000000` | `16027244793//500000000000` | $0.528801949327$ | $1.7627$ |
| $7/10$ | $(3, 17/4)$ | `58597799297//500000000000` | `48363543279//1000000000000` | $0.218794169839$ | $1.2659$ |

Each bound holds for every even chain length $L \ge n_0$ cells ($N = 2L$ divisible by $4$). The
certified bounds are obtained from `CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀).bound`
($0.5288019493276747$ and $0.2187941698399417$) and displayed rounded down to 12 decimal places;
the 12-spin ED gaps are numerical reference values. The exact rational seed defects and
$\log Z$ enclosures are stored in `results/certificate_d9.json` and `results/certificate_d7.json`.
Within the 6-cell ($12$-spin) ring cap, $\delta = 7/10$ is the smallest certified parameter on the
grid (at $\delta = 4/5$, the numerical scan also finds an admissible seed at $n_0 = 3, \beta_0 = 17/4$
giving $\approx 0.5447$, while $\delta = 3/5$ and $1/2$ require estimated seed sizes $n_0 \approx 5$
and $7$ cells, i.e., $20\text{--}28$ spins).

## Reproduce

From the repository root, using Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=models/spin1_dimerized -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1_dimerized models/spin1_dimerized/check.jl          # fast checker (~55 s)
julia --project=models/spin1_dimerized models/spin1_dimerized/check.jl --full   # + δ = 7/10 12-spin replay
# Recertify from scratch (rewrites results/certificate_d<d>.json unless SPIN1_DIMERIZED_REPLAY_PATH is set):
julia --project=models/spin1_dimerized models/spin1_dimerized/certify.jl 9 2 7//2
julia --project=models/spin1_dimerized models/spin1_dimerized/certify.jl 7
# Exploratory numerical scan per δ = d/10, then regenerate summary.json and METHOD.md tables:
for d in 9 8 7 6 5; do
    julia --project=models/spin1_dimerized models/spin1_dimerized/scripts/scan.jl $d 2,3,4,6
done
julia --project=models/spin1_dimerized models/spin1_dimerized/scripts/report.jl
```

Calling `certify.jl d` without explicit seed arguments reads the best candidate seed from
`results/scan_d<d>.json`. The scan caches full symmetry-sector spectra in gitignored
`results/spectrum_d<d>_c<cells>.json` files (~15 MB each for 6 cells); `scripts/report.jl` reads
those cached spectra to update `results/summary.json` and the final section of `METHOD.md`. All
commands also work with `--project=julia`.

Measured run times (Julia 1.12, 1 BLAS thread):
- `check.jl` (`test/runtests.jl`): 123 tests, ~55 s including compilation (of which ~12 s
  re-solves and recertifies $\delta = 9/10$ from scratch);
- `certify.jl 9 2 7//2`: ~1 s; `certify.jl 7`: 949–1017 s (~16 min);
- `check.jl --full`: adds the complete 12-spin $\delta = 7/10$ recertification (~16–36 min
  depending on system load);
- `scripts/scan.jl`: ~9 min per $\delta$ value (dominated by the 6-cell ring).

## What is rigorous

- **Rigorous:** the exact integer reflection-positivity and cell-reflection identities; the
  positive cell-transfer factorization $X_\beta = K K^* \ge 0$; complete symmetry-sector coverage
  with exact integer character projectors on the 12-spin ring (79 blocks of dimension at most
  $11{,}596$, weighted sum $9^6 = 531{,}441$); verified eigenvalue enclosures via polar-factor/Weyl
  bounds with explicit IEEE binary64 dot-product ($\gamma_{2d+8}$) and projector-basis construction
  error bounds; outward-rounded interval partition sums; exact rational seed defects; and the
  resulting uniform gap and ground-state uniqueness theorems.
- **Fast vs. full verification:** `check.jl` recomputes the rational seeds and gap bounds from the
  stored $\log Z$ enclosures of both certificates and recertifies $\delta = 9/10$ ($4$ and $8$
  spins) from scratch, while `check.jl --full` also rebuilds and re-verifies all 79 blocks of the
  12-spin ring for $\delta = 7/10$. Because LAPACK eigenvectors can vary in their last bits across
  platforms, `check.jl --full` verifies that the fresh $\log Z$ intervals overlap the stored ones
  and independently certify the gap bound.
- **Numerical only:** the $\beta_0$ grid scans, 12-spin reference ED gaps, and the seed-size
  extrapolations for $\delta \in \{3/5, 1/2\}$.

## Layout

- `certify.jl`, `check.jl`: certificate generator and verifier (fast mode by default; `--full`
  re-verifies the 12-spin $\delta = 7/10$ spectrum).
- `src/common.jl` (model, exact integer Hamiltonian in the $S^z$ basis, symmetry orbits),
  `src/rp.jl` (exact integer RP verification), `src/enclosure.jl` (polar/Weyl eigenvalue enclosures,
  $\log Z$ intervals, certificate writer), and `src/verified_blocks.jl` (exact integer dihedral
  character projectors and basis-error bounds for the 12-spin ring).
- `scripts/`: `scan.jl` (numerical seed scan across $L \in \{2, 3, 4, 6\}$ cells) and `report.jl`
  (updates `results/summary.json` and the `"Results from completed runs"` section of `METHOD.md`).
- `test/runtests.jl`: test suite executed by `check.jl`.
- `results/`: `certificate_d7.json` and `certificate_d9.json` (the gap certificates),
  `scan_d*.json` (numerical scan summaries), and `summary.json` (cross-parameter comparison).
- `METHOD.md`: detailed mathematical proofs, error analysis, audit report, and numerical tables.
- `IDEAS.md`: route evaluation and outcome.

## History

The certificates were generated when this directory was located at `certificates/spin1_dimerized`
and the core package was named `ThermalMPOGap`; their `source_sha256` entries reference archived
commit `6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2` and are verified against that commit when present
in the local git object store.
