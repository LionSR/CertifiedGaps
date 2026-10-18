# Dimerized spin-1 Heisenberg ring

**Certified at $\delta = 9/10$ and $\delta = 7/10$** by the shared core pipeline
(`GapCertificates.certify_point`): a unique ground state and $\Delta_L \ge 0.528802288$ for every
even ring of $L \ge 2$ two-spin cells at $\delta = 9/10$, and $\Delta_L \ge 0.218793401$ for every
even $L \ge 3$ cells at $\delta = 7/10$ (interval bounds; the exact rational chain bounds are
$0.528349759$ and $0.218607065$). The odd seed $n_0 = 3$ at $\delta = 7/10$ is valid because the
cell transfer operator factors as $X_\beta = K K^* \ge 0$. Derivations and the audit are in
[`METHOD.md`](METHOD.md); the route memo is in [`IDEAS.md`](IDEAS.md).

| parameter | seed $(n_0\text{ cells}, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| $\delta = 9/10$ | $(2, 7/2)$ | $0.528802288$ (chain: $0.528349759$) | every even $L \ge 2$ cells ($N = 2L$ spins) |
| $\delta = 7/10$ | $(3, 17/4)$ | $0.218793401$ (chain: $0.218607065$) | every even $L \ge 3$ cells |

## Model

```math
H_L = \sum_{j=1}^{2L} \left[1 + (-1)^j \delta\right] \mathbf{S}_j \cdot \mathbf{S}_{j+1}, \qquad
\mathbf{S}_j^2 = 2I,
```
on a periodic ring of $N = 2L$ spin-1 sites. A unit cell is a pair $(a, b) = (2x - 1, 2x)$ joined
by the weak bond $1 - \delta$ (on-site term $o$, cell dimension $9$); neighboring cells are joined
by the strong bond $1 + \delta$ from $b$ to the next $a$ (bond term $h$). Lengths $L$ and $n_0$
count cells. The scan grid is $`\delta \in \{9/10, 4/5, 7/10, 3/5, 1/2\}`$ on rings of up to 6
cells.

**Certificate frame.** The core works in the $S^z$ basis ($m = 1, 0, -1$) with the spin $b$ of
every cell rotated by $\pi$ about $y$:
$U = \prod_{\text{cells}} (1 \otimes V)$ with
$`V = e^{-i\pi S^y} = \begin{pmatrix} 0&0&1\\0&-1&0\\1&0&0 \end{pmatrix}`$.
With the integer matrices $X' = \sqrt2 S^x$, $`A' = \sqrt2\, i S^y`$ (real antisymmetric) and
$S^z = \mathrm{diag}(1, 0, -1)$, every bond becomes
```math
U (J\,\mathbf{S}_i \cdot \mathbf{S}_j) U^\dagger
= -J\,(\tfrac12 X' \otimes X' + \tfrac12 A' \otimes A' + S^z \otimes S^z).
```
In this frame the cell charge $z_a - z_b$ is conserved (an additive U(1) charge), and the spin
flip $P \otimes P$ ($P$: $m \mapsto -m$) maps charge $q$ onto $-q$, so the core builds the two
sectors once and splits $q = 0$ by flip parity. The model also carries the lattice reflection
(site reversal of the chain): it reverses the cell order and swaps $a \leftrightarrow b$, which
reverses the cell charge, so it is composed with the flip; its on-site part is
$\mathrm{SWAP}\,(P \otimes P)$, checked exactly by the core `RingModel`. With it the core builds
every momentum block as a real symmetric matrix. `test/runtests.jl` checks
$V S^x V^\dagger = -S^x$, $V (iS^y) V^\dagger = iS^y$, $V S^z V^\dagger = -S^z$, both bond
identities and the frame change from the physical cell model
(`RingModels.is_frame_change`), all exactly; it also compares the core spectrum of 3 cells with
dense floating-point ED of the physical 6-spin Hamiltonian. The Lean statement uses the
Cartesian frame $G_\alpha = iS^\alpha$; all frames are unitarily equivalent, so purities and
gaps agree.

## Reflection positivity

In the certificate frame, with $R$ the swap of the two spins of a cell,
```math
-h = \tfrac{1+\delta}{2}\, C_x \otimes R C_x R + \tfrac{1+\delta}{2}\, C_A \otimes R C_A R
+ (1 + \delta)\, C_z \otimes R C_z R, \qquad C_\alpha = 1 \otimes (X', A', S^z),
```
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
| $9/10$ | $(2, 7/2)$ | `49715646642544//10^15` | `32054414033594//10^15` | $0.528802288$ | `528349759//10^9` | $1.7627$ |
| $7/10$ | $(3, 17/4)$ | `117195681829365//10^15` | `48363910725858//10^15` | $0.218793401$ | `218607065//10^9` | $1.2659$ |

Each certificate `results/certificate_d<d>.json` stores both rigorous bounds: the
outward-rounded interval bound `certified_gap_bound(εₛ, εₜ, β₀)` (`interval_bound_decimal`, the
headline value, rounded down here) and the exact rational chain bound $\Delta_L >$
`gap_lower_bound`. At $\delta = 9/10$ the seed defects are below the values hardcoded in
`lean/CertifiedGaps/Models/SpinOneDimerized.lean` ($\varepsilon_{\mathrm{s}} \le 0.04972$,
$\varepsilon_\tau \le 0.03206$; read with `GapCertificates.lean_points` and checked by
`check.jl`); the Lean file has no $\delta = 7/10$ point. Each certificate records its point as
`params` (`"delta"`).
The 12-spin ED gaps are floating-point references. At $\delta = 4/5$ the scan also finds an
admissible seed $(3, 17/4)$ with $\approx 0.5447$; $\delta = 3/5$ and $1/2$ need estimated seeds of
$n_0 \approx 5$ and $7$ cells.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'     # once
julia --project=julia models/spin1_dimerized/certify.jl 9           # 2 min (16 s certify)
julia --project=julia models/spin1_dimerized/certify.jl 7           # 173 min, 10 GB peak
julia --project=julia models/spin1_dimerized/check.jl               # 2 min
julia --project=julia models/spin1_dimerized/check.jl --full        # as long as certify.jl 7
julia --project=julia models/spin1_dimerized/test/runtests.jl       # 11 min
# Floating-point scan per δ = d/10, then summary.json and METHOD.md tables:
for d in 9 8 7 6 5; do
    julia --project=julia models/spin1_dimerized/scripts/scan.jl $d 2,3,4,6
done
julia --project=julia models/spin1_dimerized/scripts/report.jl
```

Times are wall clock including compilation (Apple silicon, one BLAS thread), measured on a
heavily shared machine (load average about 200), so an idle machine is faster; the δ = 7/10
certify used 122 CPU minutes. `certify.jl 7` and `check.jl --full` enclose the 6-cell ring
($9^6 = 531{,}441$ states) in 53 real core blocks, the largest of dimension $11{,}596$, with
radii up to $6.6 \times 10^{-7}$; they reached 7.1 GB resident and a 10 GB peak footprint. The
fast check skips the spectra; the tests replay the rings of 2 and 4 cells. The scan caches
spectra in gitignored `results/spectrum_d<d>_c<cells>.json` files, which `scripts/report.jl`
reads.

## What is rigorous

- **Rigorous:** the exact frame unitary, charge, flip and lattice-reflection conditions and RP
  identity (exact rationals); complete spectra of the rings $n_0$ and $2n_0$ from the core
  (`VerifiedSpectra.verified_ring_spectrum`: exact charge, momentum and flip-parity blocks,
  real by the lattice reflection,
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
- **Floating point only:** the $\beta_0$ scans (core `float_ring_spectrum`), the 12-spin
  reference gaps, and the seed-size extrapolations for $`\delta \in \{3/5, 1/2\}`$.

## Layout

- `certify.jl`, `check.jl`: certificate generator (per $\delta$) and checker (`--full` replays the
  core spectra).
- `src/DimerizedCertificate.jl`: the model in the certificate frame with its lattice
  reflection, the physical cell model and frame unitary, RP factors, certified points, the
  `positive_transfer` reference, the Lean seed lookup, and `certify` / `check` (the core
  `certify_point` / `check_point` at a point).
- `scripts/scan.jl` (floating-point seed scan), `scripts/report.jl` (`results/summary.json` and
  the last section of `METHOD.md`).
- `test/runtests.jl`: frame, Cartesian-identity, reflection, RP, spectrum and certificate tests,
  including tampering.
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
- 2026-10: the model gained the lattice reflection (real momentum blocks) and `params`; the
  model-local certify and check loops, Lean seed constants, spin and swap matrices and the
  floating-point scan ED were replaced by the core `certify_point`, `check_point`,
  `lean_points`, `integer_spin_one`, `swap_matrix`, `is_frame_change` and
  `float_ring_spectrum`. Both certificates were regenerated. Old to new, $\delta = 9/10$:
  $\varepsilon_{\mathrm{s}}$ 0.049715646652 to 0.049715646643, $\varepsilon_\tau$ 0.032054414053
  to 0.032054414034, interval bound 0.5288022883 to 0.5288022884 (chain unchanged).
  $\delta = 7/10$: $\varepsilon_{\mathrm{s}}$ 0.117195722398 to 0.117195681829,
  $\varepsilon_\tau$ 0.048363998311 to 0.048363910726, interval bound 0.218793170 to
  0.218793401, chain 0.218606834 to 0.218607065 (block radii up to $6.6 \times 10^{-7}$).
