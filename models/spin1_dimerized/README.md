# Dimerized spin-1 ring

**Certified at δ = 7/10: Δ_L ≥ 0.218794169839 and a unique ground state for every even
L ≥ 3 cells** (N = 2L spins = 8, 12, 16, ...), and at δ = 9/10 with bound 0.528801949327.
The method, the enclosure theorem, the odd-seed positivity argument, the numerical tables and
the independent audit are in [REPORT.md](REPORT.md); the route memo is [IDEAS.md](IDEAS.md).

## Model

Spin one, `H = Σ_j [1 + (−1)^j δ] S_j·S_{j+1}` on a periodic ring. A cell is two sites (cell
dimension 9) joined by the weak bond `1 − δ`; the strong bond `1 + δ` joins neighbouring
cells. Lengths count cells (N = 2L spins). The grid is δ ∈ {9/10, 4/5, 7/10, 3/5, 1/2} with
rings up to six cells.

## RP basis

Cartesian spin-1 basis, `S_a = −i G_a` with integer, real antisymmetric `G_a`, so
`−S·S = Σ_a G_a ⊗ G_a` exactly. Inside the cell the on-site term `D = (1 − δ) S·S` is real
symmetric and invariant under the internal swap `P`. The intercell bond is
`−h = (1 + δ) Σ_a A_a ⊗ P conj(A_a) P` with `A_a = 1 ⊗ G_a`: the cell is internally reflected
(palindromic cell), not `A ⊗ A`. `src/rp.jl` checks all identities over exact integers. The
n₀ = 3 seed is odd; REPORT.md gives the positive cell-transfer argument (`X = K K*`) that
removes the even-n₀ restriction.

## Seed and certified bound

| δ | seed (n₀ cells, β₀) | εₛ ≥ 1 − Pₛ(n₀, β₀) | εₜ ≥ 1 − P_τ(2n₀, β₀) | certified bound | 12-spin ED gap |
|---|---|---|---|---:|---:|
| 9/10 | (2, 7/2) | `12428916293//250000000000` | `16027244793//500000000000` | 0.528801949327 | 1.7627 |
| 7/10 | (3, 17/4) | `58597799297//500000000000` | `48363543279//1000000000000` | 0.218794169839 | 1.2659 |

Each bound holds for every even L ≥ n₀ cells. The bounds are
`CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀).bound` (0.5288019493276747 and
0.2187941698399417) rounded down to twelve decimals; the ED gaps are numerical. The exact
seeds and log-partition enclosures are in `results/certificate_d9.json` and
`results/certificate_d7.json`. δ = 7/10 is the smallest certified value on the grid; the
failures at 3/5 and 1/2 (six-cell ring cap) are numerical evidence, not impossibility claims.

## How to reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=models/spin1_dimerized -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1_dimerized models/spin1_dimerized/check.jl          # fast checker
julia --project=models/spin1_dimerized models/spin1_dimerized/check.jl --full   # + δ = 7/10
# Recertify (rewrites results/certificate_d<d>.json unless SPIN1_DIMERIZED_REPLAY_PATH is set):
julia --project=models/spin1_dimerized models/spin1_dimerized/certify.jl 9 2 7//2
julia --project=models/spin1_dimerized models/spin1_dimerized/certify.jl 7
# Numerical scan per δ = d/10, then the summary and the REPORT.md tables:
for d in 9 8 7 6 5; do
    julia --project=models/spin1_dimerized models/spin1_dimerized/scripts/scan.jl $d 2,3,4,6
done
julia --project=models/spin1_dimerized models/spin1_dimerized/scripts/report.jl
```

`certify.jl d` without a seed takes the best one in `results/scan_d<d>.json`. The scan caches
its spectra in `results/spectrum_d<d>_c<cells>.json` (gitignored; the six-cell files are
15 MB each) and skips files that exist; `report.jl` needs all of them. `--project=julia` (the
core environment) also works for every command.

## Run time

Julia 1.12, one BLAS thread:

- fast checker (`check.jl`, = `test/runtests.jl`): 123 tests, about 55 s including compilation
  (of which about 12 s recertify δ = 9/10 from scratch);
- `certify.jl 7`: 1017 s in the original run, 949 s in the replay; `certify.jl 9 2 7//2`: 1 s;
- `check.jl --full`: adds the δ = 7/10 recertification (36 min on an Apple M5 under load);
- `scripts/scan.jl`: about 9 minutes per δ, almost all in the six-cell ring.

## What is rigorous

Rigorous: exact RP and cell reflection, the positive cell-transfer argument, complete
symmetry-sector counts with exact integer projectors on the twelve-spin ring, the verified
eigenvalue enclosures (polar/Weyl with IEEE rounding-error bounds and basis/construction
error, single-threaded BLAS), outward interval partition sums, rational seed bounds and the
resulting uniform gap and uniqueness conclusion. No Lanczos, ARPACK or truncated spectrum
enters.

The fast checker recomputes the seeds from the stored log-partition enclosures of both
certificates and recertifies δ = 9/10 from the model; for δ = 7/10 it treats the stored
twelve-spin enclosures as inputs, and `check.jl --full` rebuilds them. Fresh runs on another
machine differ from the stored enclosures in the last bits (LAPACK); the checks require
overlapping log Z enclosures and a valid fresh bound, not bit equality (REPORT.md,
"Validation and provenance"). A full replay of `certify.jl 7` on the renamed core
reproduced every mathematical field of `certificate_d7.json` exactly.

Numerical only: the scan grids, unverified ED reference gaps, parameter ranking and
extrapolations. No Lean formalisation of the numerical inputs is claimed.

Provenance: the certificates were made when this directory was `certificates/spin1_dimerized`
and the core was named `ThermalMPOGap`; their `source_sha256` entries are checked against the
archived commit `6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2` when it is reachable and skipped
otherwise.

## Layout

- `certify.jl`, `check.jl`: certificate entry point and checker.
- `src/common.jl` (model, exact integer Hamiltonian, symmetry bases), `src/rp.jl` (exact RP),
  `src/enclosure.jl` (eigenvalue enclosures, log Z, certificate writer),
  `src/verified_blocks.jl` (exact projectors and verified blocks for the twelve-spin ring).
- `scripts/`: `scan.jl` (numerical seed scan), `report.jl` (`results/summary.json` and the
  REPORT.md tables).
- `test/runtests.jl`: the fast checker.
- `results/certificate_d{7,9}.json`: the certificates. `results/scan_d*.json` (numerical
  scans) and `results/summary.json` (parameter comparison and extrapolation).
