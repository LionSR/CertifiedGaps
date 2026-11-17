# XY chain in a field (free fermions)

**Certified at three points of the paramagnetic phase with $J_x J_y < 0$.** For $J_x = -1$,
$h = 1$ and $J_y \in \{1/4, 1/2, 1\}$, every even periodic ring of $L \ge 128$ sites has a unique
ground state and a spectral gap at least the interval bound below, within $`6.8\text{–}12.5\%`$ of the
exact gap $2(h - |J_x + J_y|)$. The partition functions of the 128- and 256-site seed rings come
from the free-fermion closed form, enclosed in 512-bit BigFloat interval arithmetic.

| parameter $(J_x, J_y, h)$ | seed $(n_0, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| $(-1,\ 1/4,\ 1)$ | $(128,\ 189/4)$ | $0.43756$ (exact gap $1/2$) | every even $L \ge 128$ |
| $(-1,\ 1/2,\ 1)$ | $(128,\ 157/4)$ | $0.91890$ (exact gap $1$) | every even $L \ge 128$ |
| $(-1,\ 1,\ 1)$ | $(128,\ 121/4)$ | $1.86496$ (exact gap $2$) | every even $L \ge 128$ |

## Model

In Pauli normalization on a periodic ring of $L$ spin-1/2 sites (one site per unit cell, so
lengths count sites),

```math
H_L = \sum_{j=1}^L \big(J_x X_j X_{j+1} + J_y Y_j Y_{j+1}\big) - h \sum_{j=1}^L Z_j,
\qquad J_x J_y < 0,\quad h > |J_x + J_y|.
```

The Jordan–Wigner transformation, with $Z_j$ as the occupation parity, gives free fermions with
hopping $\propto J_x + J_y$ and pairing $\propto J_x - J_y$. As for the transverse-field Ising
chain ([`tfim_free_fermions`](../tfim_free_fermions/), whose `METHOD.md` gives the derivation),

```math
Z_L(\beta) = \tfrac12\Big[\prod_{\mathrm{NS}} 2\cosh\tfrac{\beta\varepsilon_k}{2}
  + \prod_{\mathrm{NS}} 2\sinh\tfrac{\beta\varepsilon_k}{2}
  + \prod_{\mathrm{R}} 2\cosh\tfrac{\beta\varepsilon_k}{2}
  - \prod_{\mathrm{R}} 2\sinh\tfrac{\beta\varepsilon_k}{2}\Big],
\qquad
\varepsilon_k = 2\sqrt{(h - (J_x + J_y)\cos k)^2 + (J_x - J_y)^2\sin^2 k},
```

with Neveu–Schwarz momenta $k = 2\pi(m + 1/2)/L$ and Ramond momenta $k = 2\pi m/L$. The unpaired
Ramond modes $k = 0, \pi$ carry the signed energies $2(h \mp (J_x + J_y))$, which are positive for
$h > |J_x + J_y|$. In that regime every $\varepsilon_k > 0$. Since $|J_x + J_y| < |J_x - J_y|$,
$\varepsilon_k^2$ is concave in $\cos k$, so the gap is $2(h - |J_x + J_y|)$. The code refuses
parameters outside $J_x J_y < 0$, $h > |J_x + J_y|$.

## Reflection positivity

Standard ($\sigma^z$) basis. Since $\overline{Y} = -Y$, the bond satisfies

```math
-h_{\text{bond}} = (-J_x)\, X \otimes \overline{X} + J_y\, Y \otimes \overline{Y}
\qquad (C_1 = X,\ C_2 = Y,\ \text{weights } -J_x,\ J_y \ge 0 \text{ for } J_x \le 0 \le J_y),
```

and the on-site term $o = -hZ$ is real symmetric. `rp_record` in `src/common.jl` checks both
identities in exact complex rational arithmetic. Each certificate stores the record, and
`check.jl` recomputes and compares it. The other sign pattern, $J_x > 0 > J_y$, maps to this one
under the on-site rotation $X \to Y$, $Y \to -X$ (which fixes $Z$). It is not certified here.

## Seeds and certified bounds

All points use $n_0 = 128$, with $\beta_0$ maximizing the interval bound on a grid of step $1/4$
(`scripts/scan.jl`). The defects are the upper ends of the interval enclosures, rounded up to 8
significant digits (exact values in [`results/summary.json`](results/summary.json)):

| $J_y$ | $\beta_0$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | interval bound $\Delta_L \ge$ | exact chain $\Delta_L >$ | exact gap |
|---|---:|---:|---:|---:|---:|---:|
| $1/4$ | $189/4$ | $1.3475289 \times 10^{-10}$ | $4.4330057 \times 10^{-10}$ | $0.437562$ | $0.437294$ | $1/2$ |
| $1/2$ | $157/4$ | $2.3271068 \times 10^{-17}$ | $9.3710855 \times 10^{-17}$ | $0.918905$ | $0.918583$ | $1$ |
| $1$ | $121/4$ | $2.2848453 \times 10^{-26}$ | $1.4063191 \times 10^{-25}$ | $1.864960$ | $1.864542$ | $2$ |

The interval bound is `certified_gap_bound`. The exact chain is `rational_gap_chain` with 6
doublings on a rounding grid that resolves $\varepsilon^{2^6}$. Smaller seeds give weaker bounds.
For example, at $J_y = 1/2$ the scan gives $0.736$ at $n_0 = 32$ and $0.852$ at $n_0 = 64$.

## Reproduce

From the repository root, in the core environment, with Julia 1.12 and
`export OPENBLAS_NUM_THREADS=1`. Times are wall-clock on a 2-core ARM server, mostly package
loading:

```sh
julia --project=code/julia code/models/xy_field/certify.jl          # 11 s
julia --project=code/julia code/models/xy_field/check.jl            # 10 s
julia --project=code/julia code/models/xy_field/check.jl --full     # 11 s
julia --project=code/julia code/models/xy_field/test/runtests.jl    # 16 s
julia --project=code/julia code/models/xy_field/scripts/scan.jl -1 1//2 1 "[128]" "36:1//4:42"
```

## What is rigorous

- **Checked by `check.jl`:** per point, the recorded model-local source hashes (`src/`,
  `certify.jl`, `check.jl`; enforced), the parameters, seed and exact RP record against
  `src/common.jl`, a rerun of `certified_gap_bound` and `rational_gap_chain` from the stored
  $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ reproducing the stored bounds, and the bound below
  the exact gap. `results/summary.json` must match the certificates. `--full` also recomputes the
  three $\log Z$ enclosures, requires them to equal the stored endpoints, and requires the fresh
  defects not to exceed the stored ones.
- **Trusted base:** the free-fermion closed form for $Z_L$ (classical; derivation in
  `../tfim_free_fermions/METHOD.md`); `IntervalArithmetic` 1.0.12 on `BigFloat` (MPFR, directed
  rounding) at 512 bits for the mode energies, $\log Z$ and the defects, which are computed and
  converted at 512 global bits (both BigFloat precision traps of the TFIM model are handled the
  same way and asserted); exact rational arithmetic for the RP identity, the rounding and the
  chain; and the core interval recurrence. No eigensolver is involved.
- **Cross-checks (tests):** the dense Hamiltonian against a Kronecker-product construction; the
  closed form against dense exact diagonalization (rings of 4 to 10 sites, five parameter sets
  including $J_x > 0 > J_y$ and $h = 2$, three $\beta$; relative agreement $10^{-12}$); the interval
  enclosure against an independent Arb (Nemo) evaluation at 1024 bits ($L = 16, 64, 128$).
- **Floating-point only:** none of the certificate. `scripts/scan.jl` uses 256-bit intervals but
  is only a seed search.

## Layout

Stored results live in [`data/models/xy_field/`](../../../data/models/xy_field/). This
directory's `results/` is a relative symlink to that location.

| path | role |
|---|---|
| `certify.jl` | certifies every point of `POINTS` and writes the certificates and summary |
| `check.jl` | checker (`--full` recomputes the $\log Z$ enclosures) |
| `src/FreeFermions.jl` | interval closed form: `mode_energies`, `logZ`, `exact`, `in_regime` |
| `src/common.jl` | points and seeds, RP record, `certify_point`, `check_point` |
| `scripts/scan.jl` | seed scan (256-bit intervals; diagnostics only) |
| `test/runtests.jl` | Hamiltonian, ED and Arb cross-checks, RP record, stored certificates, tampering |
| `IDEAS.md` | route analysis |
| `results/certificate_Jx<…>_Jy<…>_h<…>.json.gz` | certificate of one point (gzip JSON; `m` marks a minus sign) |
| `results/summary.json` | exact parameters, $n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau$ and the bounds |

No file is hash-pinned apart from the model-local sources recorded in each certificate's
`metadata.source_sha256`. The core is recorded under `core_sha256` as provenance only.
`src/FreeFermions.jl` shares its interval helpers (`ball`, `exact`, the log-cosh/log-sinh forms)
with `../tfim_free_fermions/src/FreeFermions.jl`. Moving them to the core would need every model
checker rerun, so they are duplicated for now.
