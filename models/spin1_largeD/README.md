# Spin-1 large-D gap certificate

Model: H = Σⱼ Sⱼ·Sⱼ₊₁ + D Σⱼ (Sᶻⱼ)² on periodic rings, one site per cell.

## Results

| D (exact) | Seed (n₀, β₀) | ε_s ≥ 1−Pₛ(n₀,β₀) | ε_t ≥ 1−P_τ(2n₀,β₀) | Uniform gap bound | ED gap, L=12 |
|---|---|---|---|---:|---:|
| 4 | (6, 9/4) | 577405872923865/2⁵³ ≈ 0.0641 | 45632040552745/2⁵⁰ ≈ 0.0405 | 0.694694192 | 2.2695 |
| 1603/500 = 3.206 | (6, 381/125) | 863088513139635/2⁵² ≈ 0.1916 | 418496566764289/2⁵³ ≈ 0.0465 | 0.002213456 | 1.5453 |

Both bounds hold for **every even L ≥ 6** (sites = cells), with a unique ground state,
under the repository's two-purity theorem. Bounds are rounded down; the exact binary values
and rational seeds are in `results/certificate_D*.json`. D = 3.206 is the smallest sampled D
with an admissible seed at n₀ ≤ 6; the neighbouring sample D = 3.205 has none on the refined
β grid. That is a search threshold, not a proven minimum over real D or a phase boundary.
At D = 3.206 the initial spatial defect exceeds 0.12731; the core's asymmetric recurrence
still enters the admissible square (first admissible iteration j = 7).

[REPORT.md](REPORT.md) has the scan, ED reference gaps, negative-grid diagnostics and
extrapolations (all numerical). [METHOD.md](METHOD.md) gives the RP identity, symmetry
accounting and the eigenvalue-enclosure proof. [audit/](audit/README.md) records an
independent review and a fresh dense replay that reproduced both certificates bit for bit.

## Reflection positivity

Cartesian basis (x, y, z), real antisymmetric integer matrices

    Gx = [0 0 0; 0 0 1; 0 -1 0],  Gy = [0 0 -1; 0 0 0; 1 0 0],  Gz = [0 1 0; -1 0 0; 0 0 0],

with Sₐ = −i Gₐ. The checker verifies exactly that −h = Σₐ Gₐ ⊗ Gₐ (so C_α = G_α and
conj(C_α) = C_α), h = P − |Ω⟩⟨Ω| with Ω = Σₐ|aa⟩, and that the on-site term
o = D·diag(1,1,0) = D(Sᶻ)² is real symmetric.

## What is rigorous

- **Certified:** exact RP; complete spectra of every symmetry block at L = 6 and L = 12
  (multiplicities summing to 3^L); interval block construction plus floating-point residual
  and orthogonality norm bounds with explicit rounding-error terms, fed to a polar/Weyl
  eigenvalue enclosure; outward-rounded log Z intervals; exact rational defect upper bounds;
  `CertifiedGaps.certified_gap_bound`.
- **Trusted base:** Julia, IntervalArithmetic outward rounding, the standard IEEE binary64
  dot-product error model for the BLAS products that form the residual and Q*Q−I (bounded a
  priori by γ = (8m+32)u/(1−(8m+32)u), not by interval matrix products), and the repository
  theorem.
- **Numerical only:** the parameter scan, unsuccessful grid points, ED reference gaps,
  bound/gap ratios, the plotted curve and size extrapolations.

## Checking

From the repository root (Julia 1.12, instantiated `julia/`):

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia models/spin1_largeD/test/runtests.jl     # = check.jl, about 1 min
julia --project=julia models/spin1_largeD/audit/reviewer_checks.jl
```

`test/runtests.jl` (3353 checks) verifies the exact RP identity; that the producing scripts
match the source hashes recorded in the certificates (up to the literal core rename
`ThermalMPOGap` → `CertifiedGaps`); that the stored spectra are byte-identical to the dense
runs' output (SHA-256 against `enclosure_hashes`); full multiplicity accounting; the log Z
intervals, both defect inequalities and the recurrence bound. It then runs `test/small.jl`
(full dense comparison through L = 6, enclosures through L = 4) and `test/symmetry.jl`
(exact symmetry covariance and orbit counts through L = 12). It does not redo the L = 12
dense eigensolves.

## Storage of the L = 12 spectra

`certify.jl` writes `results/enclosed_D*_L12.json` (about 30 MB each). Their sorted level
lists repeat every block eigenvalue by its symmetry multiplicity, so `spectra.jl` stores them
as `*.skeleton.json` (the JSON text with the two large arrays replaced by markers) plus
`*.levels.bin` (run-length encoded energies and radius indices, little-endian). The loader
rebuilds the original bytes exactly; the checker compares their SHA-256 with the certificate.
If a full `.json` is present it is used instead. To compact a regenerated file:

```sh
julia --project=julia -e 'include("models/spin1_largeD/spectra.jl");
    compact_enclosure("models/spin1_largeD/results", "enclosed_D4_1_L12")'
```

## Regenerating

```sh
export OPENBLAS_NUM_THREADS=1
cd models/spin1_largeD
julia --project=../../julia certify.jl 4 6
julia --project=../../julia certify.jl 4 12                   # 36 min
julia --project=../../julia seed.jl 4 6 9//4                  # needs the full L=12 .json
julia --project=../../julia certify.jl 1603//500 6
julia --threads=2 --project=../../julia certify_parallel.jl 1603//500 12   # 35 min
julia --project=../../julia seed.jl 1603//500 6 381//125
julia --project=../../julia scan.jl                           # numerical scan, D=4:-1/4:1
julia --project=../../julia scan.jl 25//8 1603//500           # extra D samples
julia --project=../../julia refine_beta.jl 1603//500          # fine beta grid at n0=6
julia --project=../../julia report.jl                         # REPORT.md, curve.{csv,json}
```

The dense runs used Julia 1.12.4, IntervalArithmetic 1.0.12, one BLAS thread, and the core
at `ef0a3e6` (then named `ThermalMPOGap`); the only later change to the producing scripts is
the package rename. `report.jl` needs the scan outputs (`scan_*`, `summary_*`,
`beta_refinement_*`), which are not stored; the sampled D values are listed in REPORT.md.
L = 14 dense ED was infeasible on the 16 GiB host (one charge-1 momentum block needs 28 GB
as a complex matrix).

## Removed from the original PR

To keep the directory small (87 MB → 3.9 MB) the following were deleted. All remain in the
branch history at commit `cf7bdd7` and can be regenerated with the commands above:

- the full `enclosed_D{4_1,1603_500}_L12.json` (replaced by the byte-exact compact form),
  and the unused `enclosed_D13_4_L4.json`, `enclosed_D3207_1000_L6.json`;
- numerical scan caches `scan_D*_L*.json`, `summary_D*.json`, `beta_refinement_D*.json`
  (regenerate with `scan.jl`, `refine_beta.jl`);
- plotting artefacts `curve.png`, `curve.vega.json` (`report.jl` rewrites the Vega spec);
- package-rename provenance: `provenance/legacy/`, `migrate_package.jl`, `check_pinned.jl`,
  `audit_current.jl`, `results/package-migration.json`, `results/validation-current.json`;
- the per-block fresh-replay evidence `audit/results/`, its driver, assembler and checker
  (`audit/replay_worker.jl`, `assemble_replay.jl`, `check.jl`), and all logs;
- `PROMPT.txt` and `PR_BODY.md`.
