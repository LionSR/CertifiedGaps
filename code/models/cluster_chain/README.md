# Cluster chain in a transverse field (free fermions)

**Certified at $g = 1/2$ (cluster phase).** Every even periodic ring of $L \ge 128$ two-site
cells has a unique ground state and a spectral gap $\Delta_L \ge 0.95156$, against the exact gap
$2(1 - g) = 1$ and the survey estimate $0.932$. Reflection positivity holds in an exact rotated
cell frame. The partition functions of the 256- and 512-site seed rings come from the
free-fermion closed form, enclosed in 512-bit BigFloat interval arithmetic.

| parameter | seed $(n_0, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| $g = 1/2$ | $(128 \text{ cells},\ 175/2)$ | $0.95156$ (exact gap $1$) | every even $L \ge 128$ cells |

## Model

On a periodic ring of $N = 2L$ spin-1/2 sites, grouped into $L$ two-site cells $(2c, 2c+1)$
(lengths count cells),

```math
H = -\sum_{j=1}^{N} Z_{j-1} X_j Z_{j+1} - g \sum_{j=1}^{N} X_j .
```

In cell form (core `RingModel` with $d = 4$, `physical_model` in `src/common.jl`) the bond between
neighboring cells is $`h = -[(1 \otimes Z) \otimes (X \otimes Z) + (Z \otimes X) \otimes (Z \otimes 1)]`$
and the cell term is $o = -g(X \otimes 1 + 1 \otimes X)$. Jordan–Wigner with $X_j$ as the
occupation parity makes $Z_{j-1} X_j Z_{j+1}$ a next-nearest-neighbor Majorana coupling, and as
for the transverse-field Ising chain ([`tfim_free_fermions`](../tfim_free_fermions/), whose
`METHOD.md` gives the derivation)

```math
Z_N(\beta) = \tfrac12\Big[\prod_{\mathrm{NS}} 2\cosh\tfrac{\beta\varepsilon_k}{2}
  + \prod_{\mathrm{NS}} 2\sinh\tfrac{\beta\varepsilon_k}{2}
  + \prod_{\mathrm{R}} 2\cosh\tfrac{\beta\varepsilon_k}{2}
  - \prod_{\mathrm{R}} 2\sinh\tfrac{\beta\varepsilon_k}{2}\Big],
\qquad \varepsilon_k = 2\sqrt{1 + g^2 - 2g\cos 2k},
```

with $k = 2\pi(m + 1/2)/N$ (NS) and $k = 2\pi m/N$ (R). On even rings the Ramond modes
$k = 0, \pi$ share the signed energy $2(g - 1)$, so their signs cancel in the $\sinh$ product, and
the unsigned roots serve for every $g \ne 1$. The tests fix these conventions against dense
exact diagonalization.

## Reflection positivity

The product basis is not reflection positive (the core `rp_certificate` rejects it, and a test
checks this). In the cell frame

```math
U = (H \otimes H)\,\mathrm{CS}, \qquad H \otimes H = \tfrac12\begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & -1 & 1 & -1 \\ 1 & 1 & -1 & -1 \\ 1 & -1 & -1 & 1 \end{pmatrix},
\qquad \mathrm{CS} = \mathrm{diag}(1, 1, 1, i),
```

an exact Gaussian-rational unitary, the rotated cell term $U^\dagger o U$ is real symmetric. The
core's exact rational LDL† (`rp_certificate`) factors the rotated bond as
$`-h' = C_1 \otimes \overline{C_1} + C_2 \otimes \overline{C_2}`$ with weights $1, 1$ and $4 \times 4$
factors with entries in $\{0, \pm 1, \pm i\}$. `rp_record` checks that $U$ is an exact frame change
(core `is_frame_change`) and that the factorization is exact (core `verify_rp`). Each certificate
stores the record, and `check.jl` recomputes and compares it.

## Seeds and certified bounds

$n_0 = 128$ cells is the survey's seed length, and $\beta_0$ maximizes the interval bound on a grid
of step $1/4$ (`scripts/scan.jl`). The defects are the upper ends of the interval enclosures,
rounded up to 8 significant digits (exact values in [`results/summary.json`](results/summary.json)):

| $g$ | $n_0$ (cells) | $\beta_0$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | interval bound $\Delta_L \ge$ | exact chain $\Delta_L >$ | survey estimate |
|---|---:|---:|---:|---:|---:|---:|---:|
| $1/2$ | 128 | $175/2$ | $4.4636947 \times 10^{-38}$ | $3.0980687 \times 10^{-37}$ | $0.951567$ | $0.951422$ | $0.932$ |

The interval bound is `certified_gap_bound`. The exact chain is `rational_gap_chain` with 6
doublings on a rounding grid that resolves $\varepsilon^{2^6}$. Smaller seeds give $0.684$
(16 cells), $0.833$ (32) and $0.910$ (64) in the scan.

## Reproduce

From the repository root, in the core environment, with Julia 1.12 and
`export OPENBLAS_NUM_THREADS=1`. Times are wall-clock on a 2-core ARM server, mostly package
loading and compilation:

```sh
julia --project=code/julia code/models/cluster_chain/certify.jl          # 12 s
julia --project=code/julia code/models/cluster_chain/check.jl            # 12 s
julia --project=code/julia code/models/cluster_chain/check.jl --full     # 12 s
julia --project=code/julia code/models/cluster_chain/test/runtests.jl    # 19 s
julia --project=code/julia code/models/cluster_chain/scripts/scan.jl 1//2 "[128]" "86:1//4:90"
```

## What is rigorous

- **Checked by `check.jl`:** the recorded model-local source hashes (`src/`, `certify.jl`,
  `check.jl`; enforced), the parameter and seed, the RP record recomputed in exact arithmetic,
  a rerun of `certified_gap_bound` and `rational_gap_chain` from the stored
  $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ reproducing the stored bounds, and the bound below
  the exact gap. `results/summary.json` must match. `--full` also recomputes the three $\log Z$
  enclosures, requires them to equal the stored endpoints, and requires the fresh defects not to
  exceed the stored ones.
- **Trusted base:** the free-fermion closed form for $Z_N$ (derivation in
  `../tfim_free_fermions/METHOD.md`, with the next-nearest-neighbor coupling in place of the
  nearest-neighbor one, and zero-mode sign conventions fixed by ED); `IntervalArithmetic` 1.0.12
  on `BigFloat` (MPFR, directed rounding) at 512 bits, with the defects computed and converted
  at 512 global bits (the BigFloat precision traps of the TFIM model, asserted here too); exact
  rational arithmetic for the RP data, the rounding and the chain; and the core interval
  recurrence. No eigensolver is involved in the certificate.
- **Cross-checks (tests):** the cell `RingModel` against the dense site Hamiltonian, and the
  rotated model's spectrum against it ($L = 2$ to $4$ cells); the closed form against dense ED
  ($L = 2$ to $5$ cells, $g \in \{1/4, 1/2, 2, 3\}$, three $\beta$; relative agreement $10^{-12}$,
  with a unique ground state and gap $\ge 2|1 - g|$); the interval enclosure against Arb at
  1024 bits; tampering with the defects and the RP record.
- **Floating-point only:** none of the certificate. `scripts/scan.jl` uses 256-bit intervals but
  is only a seed search.

## Layout

Stored results live in [`data/models/cluster_chain/`](../../../data/models/cluster_chain/). This
directory's `results/` is a relative symlink to that location.

| path | role |
|---|---|
| `certify.jl` | certifies every point of `POINTS` and writes the certificate and summary |
| `check.jl` | checker (`--full` recomputes the $\log Z$ enclosures) |
| `src/FreeFermions.jl` | interval closed form: `mode_energies`, `logZ`, `exact` |
| `src/common.jl` | point and seed, physical and rotated `RingModel`, RP record, certify/check |
| `scripts/scan.jl` | seed scan (256-bit intervals; diagnostics only) |
| `test/runtests.jl` | cell model, ED and Arb cross-checks, RP record, stored certificate, tampering |
| `IDEAS.md` | route analysis |
| `results/certificate_g1_2.json.gz` | certificate (gzip JSON) |
| `results/summary.json` | exact $g, n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau$ and the bounds |

No file is hash-pinned apart from the model-local sources recorded in the certificate's
`metadata.source_sha256`. The core is recorded under `core_sha256` as provenance only. The
interval helpers duplicate those of `../tfim_free_fermions/src/FreeFermions.jl` (see the note in
`../xy_field/README.md`).
