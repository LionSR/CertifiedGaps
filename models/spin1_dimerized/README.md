# Spin-1 dimerized ring: uniform-gap certificates

Model: spin one, `H = sum_j [1 + (-1)^j delta] S_j . S_(j+1)` on a periodic ring, with
two-site cells (cell dimension 9) containing the weak bond `1 - delta`; the strong bond
`1 + delta` joins neighbouring cells. Ring lengths count cells (N = 2L spins).

| delta | seed (n0 cells, beta0) | certified bound | lengths | verification run |
|---|---|---|---|---|
| 9/10 | (2, 7/2) | 0.528801949327 | even L >= 2 cells (N = 4, 8, 12, ...) | 1.0 s |
| 7/10 | (3, 17/4) | 0.218794169839 | even L >= 3 cells (N = 8, 12, 16, ...) | 1017 s |

The exact rational seeds and log-partition enclosures are in
`results/certificate_d9.json` and `results/certificate_d7.json`. delta = 7/10 is the
smallest certified value on the grid {0.9, 0.8, 0.7, 0.6, 0.5}; the grid failures at
0.6 and 0.5 (six-cell ring cap) are numerical evidence, not impossibility claims.

**RP basis.** Cartesian spin-1 basis, `S_a = -i G_a` with integer, real antisymmetric
`G_a`, so `-S.S = sum_a G_a (x) G_a` exactly. Within the cell the on-site term
`D = (1 - delta) S.S` is real symmetric and invariant under the internal swap `P`. The
intercell bond is `-h = (1 + delta) sum_a A_a (x) P conj(A_a) P` with `A_a = 1 (x) G_a`:
the cell is internally reflected (palindromic cell), not `A (x) A`. `rp.jl` checks all
identities over exact integers. The n0 = 3 seed is odd; `REPORT.md` and `AUDIT.md` give
the positive cell-transfer argument (`X = K K*`) that removes the even-n0 restriction.

See [REPORT.md](REPORT.md) for the eigenvalue-enclosure proof, symmetry reduction, and
tables, [AUDIT.md](AUDIT.md) for the independent review, [IDEAS.md](IDEAS.md) for the
route memo.

## Running

From the repository root (Julia 1.12):

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia models/spin1_dimerized/test/runtests.jl        # fast checker, ~30 s
julia --project=julia models/spin1_dimerized/test/verified_blocks.jl # block tests, ~20 s
```

`check.jl` is an alias for `test/runtests.jl`. The checker verifies exact RP, the
symmetry decomposition against direct ED on small rings, exact larger-period projectors,
the stored sector dimension totals, the conversion of the stored log-partition
enclosures into the rational seeds, and recomputes the bound with the current
`CertifiedGaps` core. It treats the stored twelve-spin enclosures as inputs and does not
redo the dense eigensolves. The stored certificates were computed with the core under its
old name `ThermalMPOGap`; their `source_sha256` entries are checked against the archived
commit `6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2` (path `certificates/spin1_dimerized/`)
when that commit is reachable, and skipped otherwise.

## Removed data and how to regenerate it

To keep the directory small (about 1.7 MB instead of 76 MB) the following were removed:

- `results/spectrum_d{5,6,7,8,9}_c6.json` (15 MB each): full numerical twelve-spin dense
  spectra used by the seed scan and by `report.jl`. The smaller spectra (2, 3, 4 cells),
  the scans `results/scan_d*.json` and `results/summary.json` are kept.
- `results/replay_certificate_d7.json`, `results/replay_comparison.json`,
  `replay_check.jl`: a fresh full replay of `certify.jl 7` against the renamed core
  (949 s), which reproduced every mathematical field of `certificate_d7.json` exactly.
- `PUBLICATION.md`, `results/current_main_validation.json`, `results/tests.json`:
  migration and test-run records.

Regenerate (each twelve-spin run takes about 9 minutes; the delta = 7/10 certificate
about 17 minutes):

```sh
export OPENBLAS_NUM_THREADS=1
for d in 9 8 7 6 5; do julia --project=julia models/spin1_dimerized/scan.jl $d 2,3,4,6; done
julia --project=julia models/spin1_dimerized/report.jl     # needs all c6 spectra
julia --project=julia models/spin1_dimerized/certify.jl 7  # rewrites certificate_d7.json
julia --project=julia models/spin1_dimerized/certify.jl 9 2 7//2
```

`scan.jl` skips spectrum files that already exist. To replay a certificate without
overwriting it, set `SPIN1_DIMERIZED_REPLAY_PATH=<output.json>` for `certify.jl`.

## What is rigorous

Certified: exact RP and cell reflection, the positive cell-transfer argument, complete
symmetry-sector counts, the verified eigenvalue enclosures (polar/Weyl with IEEE
rounding-error bounds, single-threaded BLAS), outward interval partition sums, rational
seed bounds and the resulting uniform gap and uniqueness conclusion. Numerical only:
the scan grids, unverified ED reference gaps, parameter ranking and extrapolations.
No Lean formalisation of this model or its numerical inputs is claimed.
