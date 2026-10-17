# Dimerized spin-1 Heisenberg ring

**Certified at $\delta = 9/10$ and $\delta = 7/10$** by the shared core pipeline
(`GapCertificates.certify_gap`): a unique ground state and $\Delta_L \ge 0.528802288$ for every
even ring of $L \ge 2$ two-spin cells at $\delta = 9/10$, and $\Delta_L \ge 0.218793170$ for every
even $L \ge 3$ cells at $\delta = 7/10$ (interval bounds; the exact rational chain bounds are
$0.528349759$ and $0.218606834$). The odd seed $n_0 = 3$ at $\delta = 7/10$ is valid because the
cell transfer operator factors as $X_\beta = K K^* \ge 0$. Derivations and the audit are in
[`METHOD.md`](METHOD.md); the route memo is in [`IDEAS.md`](IDEAS.md).

| parameter | seed $(n_0\text{ cells}, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| $\delta = 9/10$ | $(2, 7/2)$ | $0.528802288$ (chain: $0.528349759$) | every even $L \ge 2$ cells ($N = 2L$ spins) |
| $\delta = 7/10$ | $(3, 17/4)$ | $0.218793170$ (chain: $0.218606834$) | every even $L \ge 3$ cells |

## Model

$$H_L = \sum_{j=1}^{2L} \left[1 + (-1)^j \delta\right] \mathbf{S}_j \cdot \mathbf{S}_{j+1}, \qquad
\mathbf{S}_j^2 = 2I,$$
on a periodic ring of $N = 2L$ spin-1 sites. A unit cell is a pair $(a, b) = (2x - 1, 2x)$ joined
by the weak bond $1 - \delta$ (on-site term $o$, cell dimension $9$); neighboring cells are joined
by the strong bond $1 + \delta$ from $b$ to the next $a$ (bond term $h$). Lengths $L$ and $n_0$
count cells. The scan grid is $\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}$ on rings of up to 6
cells.

**Certificate frame.** The core works in the $S^z$ basis ($m = 1, 0, -1$) with the spin $b$ of
every cell rotated by $\pi$ about $y$:
$U = \prod_{\text{cells}} (1 \otimes V)$ with
$V = e^{-i\pi S^y} = \begin{pmatrix} 0&0&1\\0&-1&0\\1&0&0 \end{pmatrix}$.
With the integer matrices $X' = \sqrt2 S^x$, $A' = \sqrt2\, i S^y$ (real antisymmetric) and
$S^z = \operatorname{diag}(1, 0, -1)$, every bond becomes
$$U (J\,\mathbf{S}_i \cdot \mathbf{S}_j) U^\dagger
= -J\,(\tfrac12 X' \otimes X' + \tfrac12 A' \otimes A' + S^z \otimes S^z).$$
In this frame the cell charge $z_a - z_b$ is conserved (an additive U(1) charge), and the spin
flip $P \otimes P$ ($P$: $m \mapsto -m$) maps charge $q$ onto $-q$, so the core builds the two
sectors once and splits $q = 0$ by flip parity. `test/runtests.jl` checks
$V S^x V^\dagger = -S^x$, $V (iS^y) V^\dagger = iS^y$, $V S^z V^\dagger = -S^z$, both bond
identities, and $U H_2^{\mathrm{phys}} U^\dagger$ against the ring assembled from the
`RingModel`, all exactly; it also compares the core spectrum of 3 cells with floating-point ED
in the physical frame. The Lean statement uses the Cartesian frame $G_\alpha = iS^\alpha$; all
frames are unitarily equivalent, so purities and gaps agree.

## Reflection positivity

In the certificate frame, with $R$ the swap of the two spins of a cell,
$$-h = \tfrac{1+\delta}{2}\, C_x \otimes R C_x R + \tfrac{1+\delta}{2}\, C_A \otimes R C_A R
+ (1 + \delta)\, C_z \otimes R C_z R, \qquad C_\alpha = 1 \otimes (X', A', S^z),$$
with real integer $C_\alpha$ and real symmetric, swap-invariant $o$. The core `verify_rp` checks
this $81 \times 81$ identity, $R^2 = 1$ and $R o = o R$ exactly (`DimerizedCertificate.rp`); the
certificate stores it in the `"reflected"` form. Without $R$ the bond is not RP in this basis
(tested). In the Cartesian frame, both bonds at every site insert the same real $G_\alpha$, so
the cell transfer operator factors as $X_\beta = K K^* \ge 0$ (`METHOD.md`, "Why the odd seed
$n_0 = 3$ cells is valid"); the rotated frame differs by cell unitaries, so $X_\beta$ stays
positive. The $\delta = 7/10$ certificate passes this reference to the core as
`positive_transfer`, which the core requires for an odd $n_0$ and stores in the seed.

## Seeds and certified bounds

| $\delta$ | seed | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ | $\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$ | interval bound | exact chain bound | 12-spin ED gap |
|---|---|---|---|---:|---:|---:|
| $9/10$ | $(2, 7/2)$ | `49715646652051//10^15` | `32054414052989//10^15` | $0.528802288$ | `528349759//10^9` | $1.7627$ |
| $7/10$ | $(3, 17/4)$ | `117195722398365//10^15` | `48363998311096//10^15` | $0.218793170$ | `218606834//10^9` | $1.2659$ |

Each certificate `results/certificate_d<d>.json` stores both rigorous bounds: the
outward-rounded interval bound `certified_gap_bound(εₛ, εₜ, β₀)` (`interval_bound_decimal`, the
headline value, rounded down here) and the exact rational chain bound $\Delta_L >$
`gap_lower_bound`. At $\delta = 9/10$ the seed defects are below the values hardcoded in
`lean/CertifiedGaps/Models/SpinOneDimerized.lean` ($\varepsilon_{\mathrm{s}} \le 0.04972$,
$\varepsilon_\tau \le 0.03206$; checked by `check.jl`); the Lean file has no $\delta = 7/10$ point.
The 12-spin ED gaps are floating-point references. At $\delta = 4/5$ the scan also finds an
admissible seed $(3, 17/4)$ with $\approx 0.5447$; $\delta = 3/5$ and $1/2$ need estimated seeds of
$n_0 \approx 5$ and $7$ cells.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'     # once
julia --project=julia models/spin1_dimerized/certify.jl 9           # 10 s
julia --project=julia models/spin1_dimerized/certify.jl 7           # 109 min, 22 GB peak
julia --project=julia models/spin1_dimerized/check.jl               # 10 s
julia --project=julia models/spin1_dimerized/check.jl --full        # 115 min
julia --project=julia models/spin1_dimerized/test/runtests.jl       # 41 s
# Floating-point scan per δ = d/10 (about 9 min each), then summary.json and METHOD.md tables:
for d in 9 8 7 6 5; do
    julia --project=julia models/spin1_dimerized/scripts/scan.jl $d 2,3,4,6
done
julia --project=julia models/spin1_dimerized/scripts/report.jl
```

Times are wall clock on an otherwise idle machine, including compilation (Apple silicon, one
BLAS thread). `certify.jl 7` and `check.jl --full` enclose the 6-cell ring ($9^6 = 531{,}441$
states) in 53 core blocks, the largest of dimension $11{,}596$ (real and complex); they reached
12.5 GB resident and a 22 GB peak footprint, so run them alone on a 32 GB machine. The fast check
and the tests replay only the rings of 2 and 4 cells. The scan caches spectra in gitignored
`results/spectrum_d<d>_c<cells>.json` files, which `scripts/report.jl` reads.

## What is rigorous

- **Rigorous:** the exact frame unitary, charge and flip conditions and RP identity (exact
  rationals); complete spectra of the rings $n_0$ and $2n_0$ from the core
  (`VerifiedSpectra.verified_ring_spectrum`: exact charge, momentum and flip-parity blocks,
  eigenvalue enclosures by `verified_eigen` under the IEEE binary64 + BLAS rounding model with
  $\gamma_{N+8}$ / $\gamma_{2N+8}$ dot-product bounds and Frobenius norms, then the Ostrowski or
  polar-factor and Weyl bounds); outward-rounded $\log Z$ enclosures; exact rational seed
  defects; and the exact rational doubling chain. The trusted base is that arithmetic model,
  IntervalArithmetic, and the core code that `check.jl` reruns. Provenance:
  `Results.model_metadata` records `src/`, `certify.jl` and `check.jl` (enforced by
  `verify_model_metadata`) and the core with its environment (provenance only).
- **Analytic input:** the odd seed $n_0 = 3$ rests on the $X_\beta = K K^* \ge 0$ argument in
  `METHOD.md`, recorded as the certificate's `positive_transfer` reference, not on a computer
  check.
- **Floating point only:** the $\beta_0$ scans, the 12-spin reference gaps, and the seed-size
  extrapolations for $\delta \in \{3/5, 1/2\}$.

## Layout

- `certify.jl`, `check.jl`: certificate generator (per $\delta$) and checker (`--full` replays the
  core spectra).
- `src/DimerizedCertificate.jl`: the model in the certificate frame, its frame unitary, RP
  factors, certified points, the `positive_transfer` reference and the Lean seed values,
  `certify_point` and `check_point`.
- `src/scan_spectra.jl`: floating-point dihedral and spin-flip block ED for the scan.
- `scripts/scan.jl` (floating-point seed scan), `scripts/report.jl` (`results/summary.json` and
  the last section of `METHOD.md`).
- `test/runtests.jl`: frame, RP, scan and certificate tests, including tampering.
- `results/`: `certificate_d9.json`, `certificate_d7.json` (core format), `scan_d*.json` and
  `summary.json` (floating-point summaries). No hash-pinned files.
- No `Project.toml` or `Manifest.toml`: everything runs in the core environment
  (`--project=julia`).

## History

- The certificates were first made by a local pipeline (exact integer dihedral projectors and a
  local polar/Weyl enclosure) when this directory was `certificates/spin1_dimerized` and the core
  was `ThermalMPOGap`.
- 2026-10: both points moved to the core pipeline (rotated $S^z$ frame, `certify_gap`,
  `verified_eigen`, `positive_transfer` for the odd seed); the local Hamiltonian blocks,
  projectors, enclosure, $\log Z$, RP and provenance code and the model's own Julia environment
  were removed. Old to new, $\delta = 9/10$: $\varepsilon_{\mathrm{s}}$ 0.049715665172 to
  0.049715646652, $\varepsilon_\tau$ 0.032054489586 to 0.032054414053, interval bound
  0.528801949 to 0.528802288 (chain 0.528349759). $\delta = 7/10$: $\varepsilon_{\mathrm{s}}$
  0.117195598594 to 0.117195722398, $\varepsilon_\tau$ 0.048363543279 to 0.048363998311,
  interval bound 0.218794169 to 0.218793170 (chain 0.218606834): slightly looser, because the
  core's complex momentum blocks of dimension $11{,}596$ have radii up to $2.0 \times 10^{-6}$
  against $1.0 \times 10^{-6}$ for the old real dihedral blocks.
