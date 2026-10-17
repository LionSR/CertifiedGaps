# Spin-1 large-$`D`$ chain: uniform-gap certificates

We certify uniform spectral gaps for the periodic spin-1 Heisenberg chain with single-ion
anisotropy $D(S^z)^2$ in two stages: first from complete verified dense spectra at $L = 6$ and
$L = 12$ sites (Stage 1), and second by combining spatial transfer-operator squaring, an exact
purity coupling identity that eliminates the 24-site ring, and centered Weyl perturbation transport
in $D$ (Stage 2). Every certificate establishes a unique ground state and uniform spectral gap
$\Delta_L \ge \Delta_{\mathrm{cert}}$ for every even chain length $L \ge 6$ sites. Reported decimal
bounds are rounded **down**; exact rational parameters and defects are stored in `results/`.

| Single-ion anisotropy $D$ | Seed $(n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $4$ | $(12, 29/10)$ + finite bridge at $L = 6, 8, 10$ | **$0.899143399$** | every even $L \ge 6$ sites |
| $[3.2059998, 3.206]$ | $(6, 381/125)$, Weyl-transported | **$0.001001626$** | every even $L \ge 6$ sites |

## Model

We consider the periodic spin-1 chain with single-ion uniaxial anisotropy $D$:
```math
H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} + D \sum_{j=1}^L (S_j^z)^2, \qquad \mathbf{S}_j^2 = 2I.
```
The unit cell is a single spin-1 site ($c = 1$, local dimension $3$), so all chain lengths $L$ and
seed lengths $n_0$ count sites.

## Reflection positivity

In the Cartesian spin-1 basis $`\{|x\rangle, |y\rangle, |z\rangle\}`$, the spin operators are
$S^\alpha = -i G_\alpha$ with real, integer antisymmetric matrices $G_\alpha = i S^\alpha$:
```math
G_x = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & -1 & 0 \end{pmatrix}, \qquad
G_y = \begin{pmatrix} 0 & 0 & -1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \qquad
G_z = \begin{pmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
```
In `model.jl`, `rp_identity` verifies in exact integer arithmetic that $[S^x, S^y] = i S^z$,
$\sum_\alpha (S^\alpha)^2 = 2I$, that the bond operator $h = P - |\Omega\rangle\langle\Omega|$
(where $P$ is the two-site swap and $|\Omega\rangle = |xx\rangle + |yy\rangle + |zz\rangle$) equals
$\sum_\alpha S^\alpha \otimes S^\alpha$, and that
```math
-h = \sum_{\alpha \in \{x,y,z\}} G_\alpha \otimes G_\alpha = \sum_{\alpha \in \{x,y,z\}} G_\alpha \otimes \overline{G_\alpha}.
```
In the same Cartesian basis, the single-ion anisotropy is the real diagonal on-site operator
$o = D(S^z)^2 = D\mathrm{diag}(1, 1, 0)$.

## Seeds and certified bounds

| Anisotropy $D$ | Certified gap $\Delta_L \ge$ | Valid for | Seed $(n_0, \beta_0)$ | Certificate file |
|---|---:|---|---|---|
| $4$ | **$0.899143399$** | every even $L \ge 6$ | $(12, 29/10)$ + bridge $L = 6, 8, 10$ | `results/certificate.json` |
| $[3.2059998, 3.206]$ | **$0.001001626$** | every even $L \ge 6$ | $(6, 381/125)$, Weyl-transported | `results/certificate.json` |
| $4$ | $0.694694192$ | every even $L \ge 6$ | $(6, 9/4)$, Stage-1 ED seed | `results/certificate_D4_1.json` |
| $1603/500 = 3.206$ | $0.002213456$ | every even $L \ge 6$ | $(6, 381/125)$, Stage-1 ED seed | `results/certificate_D1603_500.json` |

The exact rational seed defects ($\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ and
$`\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)`$) are:

| Seed | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ |
|---|---|---|
| $`D = 4,\; n_0 = 12,\; \beta_0 = 29/10`$ | `8376656523711495 // 2^61` | `283096258190373 // 2^53` |
| $`D \in [3.2059998, 3.206],\; n_0 = 6,\; \beta_0 = 381/125`$ | `863115144109137 // 2^52` | `418622220052377 // 2^53` |
| $`D = 4,\; n_0 = 6,\; \beta_0 = 9/4`$ (Stage 1) | `577405872923865 // 2^53` | `45632040552745 // 2^50` |
| $`D = 3.206,\; n_0 = 6,\; \beta_0 = 381/125`$ (Stage 1) | `863088513139635 // 2^52` | `418496566764289 // 2^53` |

In addition, `results/certificate.json` certifies the higher transfer-lifted seeds $n_0 = 24, 48, 96$
at $D = 4$ (with bounds $0.660233530, 0.402885117,$ and $0.063967113$), the prior survey seed
$(6, 3)$ at $D = 4$ with bound **$0.543812926$**, and the parameter interval $D \in [3.9979773, 4]$
with bound $0.001442884$. For comparison with finite-volume gaps, $\Delta_{12} \approx 2.26946$ at
$D = 4$ ($\Delta_8 \approx 2.27392, \Delta_{10} \approx 2.27024$, so the $D = 4$ bound reaches
$`39.6\%`$ of $\Delta_{12}$) and $\Delta_{12} \approx 1.54533$ at $D = 3.206$.

- **Stage 1 (verified exact diagonalization, `certify.jl` and `seed.jl`):** Conservation of total
  $S^z$, single-site translation momentum, spin inversion at $S^z = 0$, and conjugate-momentum
  pairing decompose the $3^{12} = 531{,}441$-dimensional Hilbert space at $L = 12$ into 92 symmetry
  blocks (largest dimension $5{,}798$). Every block is diagonalized densely, and each ordered
  eigenvalue is enclosed via interval block construction, IEEE binary64 residual and orthogonality
  norm bounds, and a polar-factor/Weyl inequality. Summing $\log Z$ in outward interval arithmetic
  at $L = 6$ and $L = 12$ certifies the $n_0 = 6$ seeds.
- **Stage 2 (transfer squaring and Weyl transport, `scripts/transfer_certify.jl`):** Reusing the
  stored $L = 6$ and $L = 12$ verified spectra at both $\beta_0$ and $2\beta_0$, transfer-operator
  squaring bounds $P_{\mathrm{s}}(12, \beta_0)$ and the exact coupling identity
  ```math
  P_\tau(4n, \beta_0) = \frac{P_{\mathrm{s}}(2n, 2\beta_0)\,P_\tau(2n, \beta_0)^2}{P_{\mathrm{s}}(2n, \beta_0)^2}
  ```
  bounds $P_\tau(24, \beta_0)$ **without constructing the 24-site ring** (matching Theorem 2.24(v)
  of the survey). This lifts the $D = 4$ seed from $n_0 = 6$ to $n_0 = 12$, raising the certified
  uniform bound from $0.694694192$ to $0.899143399$ for all even $L \ge 12$; verified dense spectra
  at $L = 8$ and $L = 10$ (together with $L = 6$) bridge the remaining even lengths down to $L \ge 6$.
  Finally, a centered Weyl perturbation bound ($0 \le \sum_j (S_j^z)^2 \le L$, so
  $`\|\delta D (\sum_j (S_j^z)^2 - \frac{L}{2}I)\|_2 \le \frac{1}{2}L|\delta D|`$) transports the
  $D = 3.206$ seed to the continuous interval $D \in [3.2059998, 3.206]$. Full proofs are given in
  [`METHOD.md`](METHOD.md), and the route evaluation is in [`IDEAS.md`](IDEAS.md).

### Exploratory numerical scan (floating point, non-certifying)

We scanned $`D \in \{1.0, 1.25, \dots, 4.0\} \cup \{3.205, 3.206\}`$ over $`n_0 \in \{2, \dots, 6\}`$
and $`\beta_0 \in \{0.125, 0.25, \dots, 20.0\}`$ (refined to $\beta_0 = 2.75:0.002:3.75$ at
$D = 3.205, 3.206$). The results are summarized in [`results/scan_summary.json`](results/scan_summary.json)
and plotted in [`results/transfer_curve.pdf`](results/transfer_curve.pdf):

| $D$ | Best even seed $(n_0, \beta_0)$ | Numerical bound | $\Delta_{12}$ (ED) | $\max_\beta \min(P_{\mathrm{s}}, P_\tau)$ at $n_0 = 6$ | $\delta(6) = E_0(12) - 2E_0(6)$ |
|---:|---|---:|---:|---:|---:|
| $1.0$ | none at $n_0 \le 6$ | — | $0.2554$ | $0.335$ | $0.3528$ |
| $1.5$ | none at $n_0 \le 6$ | — | $0.3418$ | $0.375$ | $0.3229$ |
| $2.0$ | none at $n_0 \le 6$ | — | $0.5840$ | $0.494$ | $0.2422$ |
| $2.5$ | none at $n_0 \le 6$ | — | $0.9452$ | $0.655$ | $0.1505$ |
| $3.0$ | none at $n_0 \le 6$ | — | $1.3641$ | $0.799$ | $0.0857$ |
| $3.205$ | none at $n_0 \le 6$ | — | $1.5444$ | $0.839$ | $0.0677$ |
| $3.206$ | $(6, 3.048)$ | $0.00283$ | $1.5453$ | $0.839$ | $0.0676$ |
| $3.25$ | $(6, 3.0)$ | $0.0931$ | $1.5845$ | $0.847$ | $0.0642$ |
| $3.5$ | $(6, 2.625)$ | $0.3232$ | $1.8095$ | $0.888$ | $0.0483$ |
| $3.75$ | $(6, 2.375)$ | $0.5121$ | $2.0381$ | $0.916$ | $0.0365$ |
| $4.0$ | $(6, 2.25)$ | $0.6948$ | $2.2695$ | $0.938$ | $0.0278$ |

At $n_0 = 6$ ($2n_0 = 12$ sites), $D = 3.206$ sits right at the boundary of the admissible region
for direct ED seeds; reaching smaller $D$ requires larger seed lengths $n_0$ (a single $L = 14$,
$S^z = 1$ momentum block has dimension $41{,}835$, requiring $28.0\text{ GB}$ for one dense complex
matrix; see [`results/feasibility.json`](results/feasibility.json)).

## Reproduce

Using Julia 1.12 from the repository root:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=models/spin1_largeD -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1_largeD models/spin1_largeD/check.jl            # fast check (~20 s)
julia --project=models/spin1_largeD models/spin1_largeD/check.jl --full     # + L=12 dense replay (~70 min)
julia --project=models/spin1_largeD models/spin1_largeD/test/runtests.jl    # unit tests (~25 s)
julia --project=models/spin1_largeD models/spin1_largeD/scripts/transfer_certify.jl  # Stage 2 (~15 s)
```

Running `scripts/transfer_certify.jl` (Stage 2) regenerates `results/certificate.json` and
`results/finite_bridge.json`. `results/transfer_curve.json` (from `scripts/report.jl`, below)
records the certificate's SHA-256 hash as provenance; a test requires its certified points to equal
the certificate's cases, so rerun `scripts/report.jl` only if a certified number changes. Running
`check.jl --full` adds a complete fresh dense eigensolve and verification of both $L = 12$ spectra
(~35 min per $D$ value with 1 BLAS thread); the default fast `check.jl` already replays both
$L = 6$ spectra and the $L = 8, 10$ bridge spectra from scratch.

To rerun the exploratory floating-point scan and plot (`pdflatex` required for PDF plot):

```sh
julia --project=models/spin1_largeD models/spin1_largeD/scripts/scan.jl              # ~4 min
julia --project=models/spin1_largeD models/spin1_largeD/scripts/scan.jl 1603//500 641//200
julia --project=models/spin1_largeD models/spin1_largeD/scripts/refine_beta.jl 1603//500 641//200
julia --project=models/spin1_largeD models/spin1_largeD/scripts/report.jl
```

To regenerate the Stage-1 enclosed spectra from scratch (~35 min per $L = 12$ run on 1 thread):

```sh
julia --project=models/spin1_largeD models/spin1_largeD/certify.jl 4 12
julia -t 4 --project=models/spin1_largeD models/spin1_largeD/certify_parallel.jl 1603//500 12
julia --project=models/spin1_largeD models/spin1_largeD/seed.jl 4 6 9//4
```

## What is rigorous

- **Rigorous (relative to the trusted base):** the exact integer reflection-positivity identity;
  symmetry-sector completeness ($\sum \mathrm{mult} \times d = 3^L$); every eigenvalue enclosure at
  $L = 6, 8, 10, 12$ via polar-factor/Weyl bounds with explicit $\gamma_{8m+32}$ IEEE binary64
  rounding bounds; outward-rounded $\log Z$ intervals (including rigorous full-spectrum tail bounds
  above the lowest $K = 8{,}192$ sorted levels in Stage 2); exact rational seed defects
  $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$; spatial transfer squaring, purity coupling, and
  centered Weyl $D$-transport; and the doubling recurrence bounds from
  `CertifiedGaps.certified_gap_bound`.
- **Fast vs. full verification:** `check.jl` verifies the SHA-256 hashes of the stored $L = 6, 12$
  enclosures, re-solves and re-verifies the $L = 6$ and $L = 8, 10$ bridge spectra from scratch, and
  recomputes all Stage-1 and Stage-2 $\log Z$ intervals, rational seeds, and gap bounds. Running
  `check.jl --full` additionally re-solves and re-verifies all 92 blocks of both $L = 12$ spectra,
  checking that every fresh eigenvalue center agrees with the stored one within the sum of their
  verified radii.
- **Numerical only:** the floating-point $(D, n_0, \beta_0)$ grid scan, the $\Delta_{12}$ values
  used as plot denominators in `results/transfer_curve.pdf`, and the seed-size extrapolations for
  $D \le 3.2$.

## Layout

| Path | Role |
|---|---|
| `model.jl` | Model definition, exact RP identity, symmetry blocks, and verified dense enclosures (**hash-pinned**) |
| `certify.jl`, `certify_parallel.jl` | Stage 1: complete enclosed spectra at $L = 6, 12$ (**hash-pinned**) |
| `seed.jl` | Stage 1: ED seed certificates `results/certificate_D*.json` (**hash-pinned**) |
| `check.jl` | Verifier for both stages (`--full` adds dense $L = 12$ spectral replay) |
| `src/TransferCertificate.jl` | Stage-2 transfer squaring, purity coupling, Weyl transport, bridge, and checker |
| `src/spectra.jl` | Byte-exact compact loader/reconstructor for the $L = 12$ spectra |
| `scripts/transfer_certify.jl` | Stage-2 generator writing `results/certificate.json` and `results/finite_bridge.json` |
| `scripts/` (others) | Exploratory floating-point scan, $\beta$ refinement, $L = 14$ feasibility, and report/plot |
| `test/` | Unit tests; `test/ed_seeds.jl` implements the Stage-1 replay invoked by `check.jl` |
| `results/` | Compact enclosed spectra, Stage-1 and Stage-2 certificates, finite bridge, scan summary, and plot |

**Hash-pinned layout exceptions:** `model.jl`, `seed.jl`, `certify.jl`, and `certify_parallel.jl`
must remain byte-identical at the top level of this directory because the Stage-1 certificates and
$L = 12$ spectrum files record their SHA-256 hashes under those exact paths (and regenerating the
$L = 12$ spectra takes ~35 minutes per $`D`$ value). The $L = 12$ spectra are stored compactly in
`results/` as `*.skeleton.json` plus run-length-encoded `*.levels.bin`; `src/spectra.jl`
reconstructs the exact original JSON bytes in memory so `check.jl` can verify their SHA-256 hash.

**Stage-2 source hashes:** `results/certificate.json` and `results/finite_bridge.json` record the
SHA-256 hashes of the model-local proof sources (`src/`, `model.jl`, `scripts/transfer_certify.jl`,
`check.jl`, `Project.toml`, `Manifest.toml`; keyed relative to this directory) under
`source_hashes`, and the checker requires exactly these keys and digests. The shared core
(`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`) is recorded under
`run_metadata.core_sha256`, keyed by repository path, as provenance only: the checker reruns the
core arithmetic and does not compare those hashes with the checkout. Notes and `METHOD.md` are not
hashed. (The Stage-1 files keep their own pinned hashes, listed above.)

## History

`results/feasibility.json` was produced on the original Linux host by the script now at
`scripts/feasibility.jl`; its recorded script hashes refer to the earlier path and content. In the
repository directory-format standardization, the Stage-2 driver moved from `transfer_certify.jl` to
`scripts/transfer_certify.jl` and its output from `certificate.json` to `results/certificate.json`,
while the hash-pinned Stage-1 top-level `.jl` files were left untouched. Both Stage-2 files were
regenerated (every case, bound, interval, and bridge gap identical; only source hashes and timings
changed), and `check.jl` gained the $L = 6$ replay and `--full`. `scripts/report.jl` was rerun on
fresh scan caches: its certified entries were identical, but on an Apple M5 the floating-point scan
reproduced the stored screening numbers only to the last bits (LAPACK), so
`results/scan_summary.json`, `results/transfer_curve.pdf`, and the numbers in
`results/transfer_curve.json` were kept and only the latter's `certificate_sha256` was updated to
the regenerated certificate.

The source-hash policy change (core hashes become provenance; the checker enforces only model-local
sources) regenerated both Stage-2 files with `scripts/transfer_certify.jl`; a deep comparison of old
and new JSON found every case, $\log Z$ enclosure, defect, bound, interval, and bridge gap
identical, with only hashes, timestamps, and timings changed. The curve test now compares the
certified numbers instead of the certificate's bytes, and the `certificate_sha256` in
`results/transfer_curve.json` was again updated by hand.
