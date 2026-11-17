# Transverse-field Ising chain at zero longitudinal field (free fermions)

**Certified at three points of the paramagnetic phase.** For $g = 2$, $5/4$ and $11/10$, every
even periodic ring of $L \ge n_0$ sites has a unique ground state and a spectral gap at least the
interval bound below, within $`4\text{–}10\%`$ of the exact gap $2(g - 1)$. The seed rings of
$128$ to $768$ sites are far beyond exact diagonalization. The partition functions come from the
free-fermion closed form, enclosed in 512-bit BigFloat interval arithmetic.

| parameter | seed $(n_0, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| $g = 2$ | $(128,\ 175/4)$ | $1.91897$ (exact gap $2$) | every even $L \ge 128$ |
| $g = 5/4$ | $(144,\ 63)$ | $0.45128$ (exact gap $1/2$) | every even $L \ge 144$ |
| $g = 11/10$ | $(384,\ 359/2)$ | $0.18259$ (exact gap $1/5$) | every even $L \ge 384$ |

## Model

In Pauli normalization on a periodic ring of $L$ spin-1/2 sites (one site per unit cell, so
lengths count sites),

```math
H_L = -\sum_{j=1}^L Z_j Z_{j+1} - g \sum_{j=1}^L X_j, \qquad g > 1,
```

the $h = 0$ line of [`tfim_longitudinal`](../tfim_longitudinal/). The Jordan–Wigner transformation
with $X_j$ as the occupation parity gives free fermions with antiperiodic (Neveu–Schwarz,
$k = 2\pi(m + 1/2)/L$) boundary conditions on the even-parity sector and periodic (Ramond,
$k = 2\pi m/L$) ones on the odd sector, with mode energies
$`\varepsilon_k = 2\sqrt{1 + g^2 - 2g\cos k}`$. Projecting each sector onto its parity gives

```math
Z_L(\beta) = \tfrac12\Big[\prod_{\mathrm{NS}} 2\cosh\tfrac{\beta\varepsilon_k}{2}
  + \prod_{\mathrm{NS}} 2\sinh\tfrac{\beta\varepsilon_k}{2}
  + \prod_{\mathrm{R}} 2\cosh\tfrac{\beta\varepsilon_k}{2}
  - \prod_{\mathrm{R}} 2\sinh\tfrac{\beta\varepsilon_k}{2}\Big].
```

For $g > 1$ every $\varepsilon_k \ge 2(g - 1) > 0$ (for $g < 1$ the Ramond zero mode would carry
the signed energy $2(g - 1)$). The code refuses $g \le 1$. `METHOD.md` sketches the derivation.

## Reflection positivity

The certificate frame is the physical $\sigma^z$ frame. The bond satisfies

```math
-h_{\text{bond}} = Z \otimes Z = Z \otimes \overline{Z} \qquad (C_1 = Z),
```

and the on-site term $o = -gX$ is real symmetric. `rp_record` in `src/common.jl` checks both
identities in exact rational arithmetic. It stores the record in each certificate, and
`check.jl` recomputes and compares it.

## Seeds and certified bounds

The seed lengths are those of the floating-point survey that listed these points as open targets.
At each $n_0$, $\beta_0$ maximizes the interval bound on a grid of step $1/4$ ($1/2$ for
$g = 11/10$) (`scripts/scan.jl`). The defects are the upper ends of the interval enclosures,
rounded up to 8 significant digits (exact values in [`results/summary.json`](results/summary.json)):

| $g$ | $n_0$ | $\beta_0$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | interval bound $\Delta_L \ge$ | exact chain $\Delta_L >$ | survey estimate |
|---|---:|---:|---:|---:|---:|---:|---:|
| $2$ | 128 | $175/4$ | $2.2318474 \times 10^{-38}$ | $1.5490343 \times 10^{-37}$ | $1.918977$ | $1.918687$ | $1.885$ |
| $5/4$ | 144 | $63$ | $4.9855656 \times 10^{-14}$ | $1.9344913 \times 10^{-13}$ | $0.451282$ | $0.451081$ | $0.451$ |
| $11/10$ | 384 | $359/2$ | $6.0964519 \times 10^{-16}$ | $2.5256604 \times 10^{-15}$ | $0.182594$ | $0.182523$ | $0.178$ |

The interval bound is `certified_gap_bound` (48 interval doublings). The exact chain is
`rational_gap_chain` with 6 doublings, on a rounding grid fine enough to resolve
$\varepsilon^{2^6}$. Larger seeds give larger bounds at the cost of a larger $L_0$. The scan gives,
for example, $1.9498$ at $(256, 78)$ for $g = 2$, $0.4704$ at $(256, 112)$ for $g = 5/4$, and
$0.1864$ at $(512, 240)$ for $g = 11/10$. These were not certified here.

## Reproduce

From the repository root, in the core environment (this model has no environment of its own),
with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`. Times are wall-clock on a 2-core ARM server
and are dominated by package loading:

```sh
julia --project=code/julia code/models/tfim_free_fermions/certify.jl          # 12 s
julia --project=code/julia code/models/tfim_free_fermions/check.jl            # 10 s
julia --project=code/julia code/models/tfim_free_fermions/check.jl --full     # 10 s
julia --project=code/julia code/models/tfim_free_fermions/test/runtests.jl    # 14 s
julia --project=code/julia code/models/tfim_free_fermions/scripts/scan.jl 2 "[128]" "40:1//4:48"
```

## What is rigorous

- **Checked by `check.jl`:** per point, the recorded model-local source hashes (`src/`,
  `certify.jl`, `check.jl`; enforced), the point, seed and exact RP record against
  `src/common.jl`, a rerun of `certified_gap_bound` and `rational_gap_chain` from the stored
  $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ reproducing the stored bounds, and the bound below
  $2(g - 1)$. `results/summary.json` must match the certificates. `--full` also recomputes the
  three $\log Z$ enclosures from the closed form, requires them to equal the stored endpoints,
  and requires the fresh defects not to exceed the stored ones.
- **Trusted base:** the free-fermion closed form for $Z_L$ (a classical theorem; see `METHOD.md`);
  `IntervalArithmetic` 1.0.12 on `BigFloat` (MPFR with directed rounding) at 512 bits for the
  mode energies, $\log Z$ and the defects; exact rational arithmetic for the RP identity, the
  rounding of the defects and the chain; and the core interval recurrence. No floating-point
  eigensolver is involved.
- **Precision traps handled:** Julia's `Rational{BigInt}(x::BigFloat)` rounds to the *global*
  BigFloat precision, and `IntervalArithmetic`'s `expm1` on a `BigFloat` interval with more bits
  than the global precision returns a point interval that need not contain the true value. The
  defects are therefore computed and converted at 512 global bits (asserted exact in
  `seed_enclosures`), and stored endpoints use `exact`. See `METHOD.md`.
- **Cross-checks (tests):** the closed form against dense exact diagonalization (rings of
  4 to 10 sites, four values of $g$, three $\beta$; relative agreement $10^{-12}$), and the
  interval enclosure against an independent Arb (Nemo) evaluation at 1024 bits for
  $L = 16, 64, 128$.
- **Floating-point only:** none of the certificate. `scripts/scan.jl` uses 256-bit intervals but
  is only a seed search.

## Layout

Stored results live in [`data/models/tfim_free_fermions/`](../../../data/models/tfim_free_fermions/).
This directory's `results/` is a relative symlink to that location.

| path | role |
|---|---|
| `certify.jl` | certifies every point of `POINTS` and writes the certificates and summary |
| `check.jl` | checker (`--full` recomputes the $\log Z$ enclosures) |
| `src/FreeFermions.jl` | interval closed form: `mode_energies`, `logZ`, `exact` |
| `src/common.jl` | points and seeds, RP record, `certify_point`, `check_point` |
| `scripts/scan.jl` | seed scan (256-bit intervals; diagnostics only) |
| `test/runtests.jl` | ED and Arb cross-checks, precision traps, stored certificates, tampering |
| `IDEAS.md`, `METHOD.md` | route analysis; derivation and audit notes |
| `results/certificate_g<p>_<q>.json.gz` | certificate of $g = p/q$ (gzip JSON) |
| `results/summary.json` | exact $g, n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau$ and the bounds |

No file is hash-pinned apart from the model-local sources recorded in each certificate's
`metadata.source_sha256`. The core is recorded under `core_sha256` as provenance only.
