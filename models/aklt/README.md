# Hexagonal spin-3/2 AKLT model

This package (`AKLTPilot`) is a numerical pilot of the $2+1$-dimensional thermal-purity recursion
developed in [`notes/thermal-mpo-gap-2plus1d.tex`](../../notes/thermal-mpo-gap-2plus1d.tex) for the
periodic spin-$3/2$ Affleck–Kennedy–Lieb–Tasaki (AKLT) model on the honeycomb lattice. It is a
floating-point implementation check, not a certificate: standard MPS thermal purification on
snake-ordered two-dimensional torus quotients suffers from severe bond-compression bias and does not
reproduce the smallest exact reference, so no seed or gluing diagnostic is resolved.

| parameter | seed (n₀, β₀) | certified gap ≥ | lengths covered |
|---|---|---|---|
| honeycomb AKLT | none (8-spin candidate $R(64) \approx 0.0083$, floating point) | no certificate | — |

| geometry | spins | method | status |
|---|---:|---|---|
| $2\times 2$ honeycomb torus | $8$ | dense exact diagonalization (`Float64`) | $\Delta \approx 0.09203$, thermal excess $R(\beta) \le 0.01$ at $\beta \approx 62.0$ (unverified `Float64`) |
| $2\times 2, 4\times 2, 4\times 4$ tori | $8, 16, 32$ | ITensors thermal MPO purification ($\chi \le 128$) | negative result: compression bias misses the 8-spin $\log Z(\beta=2)$ reference by $-0.333$ at $\chi = 128$ |

## Model

Each primitive unit cell at $(x, y) \in \mathbb{Z}_{L_x} \times \mathbb{Z}_{L_y}$ contains two spin-$3/2$
sites, $A(x, y)$ and $B(x, y)$ (local Hilbert-space dimension $d = 4$ per spin, $16$ per cell). Each site
$A(x, y)$ is coupled to its three nearest neighbors $B(x, y)$, $B(x-1, y)$, and $B(x, y-1)$ with periodic
boundary conditions in both directions. The Hamiltonian is a sum of spin-$3$ bond projectors:

$$H = \sum_{\langle i, j \rangle} P_{ij}^{(S=3)}$$

where $P_{ij}^{(S=3)}$ is the orthogonal projector onto total spin $S = 3$ on bond $\langle i, j \rangle$,
normalized with eigenvalues $0$ and $1$ so that the frustration-free AKLT ground state has exact energy
$E_0 = 0$. The pilot studies $2\times 2$, $4\times 2$, and $4\times 4$ primitive-cell tori ($8$, $16$,
and $32$ spins). Because these small periodic quotients have short non-contractible cycles, they serve as
implementation benchmarks and finite-size diagnostics rather than bulk thermodynamic-limit evidence.

### Diagnostics of the 2+1D recursion

In $2+1$ dimensions, the note [`notes/thermal-mpo-gap-2plus1d.tex`](../../notes/thermal-mpo-gap-2plus1d.tex)
formulates a recursion on the **partition-function excess** $R(n, m, \beta) = Z(n, m, \beta) - 1$ for
frustration-free Hamiltonians with $E_0 = 0$, tracking:
- the excess $R(n, m, \beta) = Z(n, m, \beta) - 1$ (total Boltzmann weight outside a single zero-energy ground state),
- the doubling coefficient $C = [Z(2n, 2m, 2\beta) - 1]/R(n, m, \beta)^2$, and
- the spatial gluing defect $G = \log Z(2n, 2m, 2\beta) - 4 \log Z(n, m, 2\beta) = 2G_x + G_y$.

To close the recursion, one needs $C R < 1$ at a seed scale together with a uniform spatial gluing bound
across larger scales. For the one-step projection-filter transfer operator at width $2$ (transfer matrix
dimension $256$), an auxiliary check documented in the note shows that neither $T$ nor $T^2$ admits an
intertwining positive metric, ruling out a naive self-adjoint spatial transfer argument for that filter,
although direct characteristic-polynomial bounds at width $2$ yield $d_{2L} \le d_L^2/4$ for even $L \ge 4$.

## Reflection positivity

Not available. The 2+1D note proves that the standard time-reversal construction fails for one
AKLT bond: with \(\theta(A)=U\bar AU^\dagger\), \(U=e^{-i\pi S^y}\), and the quadrupole
\(Q=(S^z)^2-\tfrac54\), \(\operatorname{Tr}[e^{-tP^{(3)}}(Q\otimes\theta Q)]=\tfrac75(e^{-t}-1)<0\)
for every \(t>0\). `test/runtests.jl` checks the ingredients in floating point
(\(\operatorname{Tr}[P^{(3)}(Q\otimes Q)]=7/5\) and the projector polynomial). This excludes only
that local construction; the scalar gluing inequalities the 2+1D recursion needs remain open.

## Result

- **No certified gap or partition-function enclosure.** All outputs in this directory are exploratory
  `Float64` calculations stored with `certified = false`.
- **Negative result for 2D snake-ordered thermal MPO purification:** On the 8-spin ($2\times 2$-cell)
  torus at $\beta = 2$, dense diagonalization gives $\log Z = 5.31205522615$. Evolving the purified
  amplitude $e^{-\beta H/2}$ with step $\delta\beta = 1/8$ at maximum bond dimensions
  $\chi \in \{16, 32, 64, 128\}$ underestimates $\log Z$ by $-1.773$, $-1.424$, $-0.587$, and $-0.333$,
  respectively. Halving the Trotter step to $\delta\beta = 1/16$ does not remove the error, confirming
  that the discrepancy comes from SVD bond-compression truncation across long-range snake bonds rather
  than Trotter splitting error. Avoiding site swap restoration between nonlocal gates (`thermal-swap-refinement`)
  gives $\log Z = 4.91941051889$ versus $4.97859393485$ at $\chi = 128$.
- **Gluing diagnostics at $\chi \le 128$ are dominated by truncation bias:** Because the purified $\log Z$
  values are not converged even on 8 spins at $\beta = 2$, observed negative values of $G$ or small values
  of $C$ at fixed bond dimension cannot be interpreted as physical gluing bounds.
- **Eight-spin floating-point reference:** Full dense diagonalization of the 8-spin Hilbert space ($65{,}536$
  states split across conserved total-$S^z$ sectors) finds a unique ground state with numerical energy
  $1.35 \times 10^{-14}$, a finite-size gap $\Delta \approx 0.0920298633241$, and thermal excess
  $R(64) = 0.00831113960$ ($R(\beta) = 0.01$ is crossed near $\beta \approx 61.9932$).

## Reproduce

`AKLTPilot` depends on the shared core package `CertifiedGaps` (`../../julia`) for spin matrices, sector
diagonalization, thermal purification, and JSON output. From the repository root:

```sh
julia --startup-file=no --project=models/aklt -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/aklt models/aklt/test/runtests.jl   # 38 tests, ≈3.5 min
julia --startup-file=no --project=models/aklt models/aklt/scripts/run.jl \
  models/aklt/scripts/configs/exact-small.toml models/aklt/results/exact-small.json
```

There is no `certify.jl` (nothing is certified) and no `check.jl` (there are no tracked results to
check); the test suite is the verification. Additional configuration files in `scripts/configs/`
are:
- `thermal-pilot.toml`: $2\times 2$ base geometry, $\chi \in \{16, 32\}$,
  $\delta\beta \in \{1/8, 1/16\}$;
- `thermal-reference-refinement.toml`: 8-spin torus, $\chi \in \{64, 128\}$;
- `thermal-swap-refinement.toml`: 8-spin torus without swapping moved sites back between gates.

Each run records source and environment hashes, parameter settings, package versions, thread counts,
geometry, a `complete` flag, and `certified = false`, and rewrites the output JSON after each
completed sector or temperature step, so an interrupted run is visible; use a fresh filename to keep
an earlier run.

## What is rigorous

Nothing in this directory is a rigorous certificate. The exact identities checked by the tests
(projector polynomial, $\operatorname{Tr}[P^{(3)}(Q\otimes Q)]=7/5$, and the single-bond partition
function) hold in exact arithmetic, but the tests evaluate them in `Float64`.
- **Reference diagonalization (`src/exact.jl`):** Builds $P_{ij}^{(S=3)}$ from its exact degree-3
  polynomial in $\mathbf{S}_i \cdot \mathbf{S}_j$ (independently of the Clebsch–Gordan gate
  decomposition) and diagonalizes every total-$S^z$ sector in `Float64` using spin-flip symmetry for
  $S^z < 0$. There is no Hilbert-space truncation on 8 spins, but the eigenvalues are unverified
  dense `Float64` eigensolves.
- **Thermal MPO purification (`src/thermal.jl`):** Wraps
  `CertifiedGaps.Purification.purified_log_partition`. Each two-site bond gate uses the exact
  projector identity $e^{-t P} = I + (e^{-t} - 1) P$, but the product over noncommuting honeycomb
  bonds incurs both uncontrolled second-order Trotter splitting error and uncertified SVD
  bond-truncation error.

## Layout

| file | role |
|---|---|
| `src/AKLTPilot.jl` | Top-level module re-exporting model, exact diagonalization, and thermal purification routines |
| `src/model.jl` | Periodic honeycomb torus geometry, perfect-matching bond layers, and snake ordering |
| `src/exact.jl` | Polynomial spin-3 bond projector and sector-by-sector `Float64` reference diagonalization |
| `src/thermal.jl` | Thermal MPO purification $e^{-\beta H/2}$ and excess/gluing diagnostics ($R, C, G_x, G_y$) |
| `scripts/run.jl` | Command-line runner executing a `.toml` configuration from `scripts/configs/` and writing JSON results |
| `scripts/configs/*.toml` | The four run configurations listed under "Reproduce" |
| `test/runtests.jl` | Unit tests verifying projector identities, lattice geometry, and small-system consistency |
| `.JuliaFormatter.toml` | Formatter settings for this directory (editor tooling, not data) |
| `.gitignore` | Keeps `results/` local |

`results/` is gitignored and kept locally. The stored files were produced before the move to
`models/`, with the `spin2` environment, so their `environment/` hashes refer to
`spin2/Project.toml` and `spin2/Manifest.toml`; `pilot-summary.json` and the two projection files
there came from local Python helpers that are not part of the repository. Rerun the configurations
to regenerate the Julia outputs.

## History

- Added with the notes and the Julia core; moved under `models/` when the core was renamed
  `CertifiedGaps`.
- #15: README rewritten to state that nothing is certified and to record the negative results;
  an unused dependency removed.
- #19: aligned with the common directory format.
- Run configurations moved from `configs/` to `scripts/configs/`, and `IDEAS.md` reconstructed,
  to match the directory format.
