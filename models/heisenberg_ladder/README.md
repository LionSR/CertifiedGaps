# Heisenberg ladder certificate

Spin-1/2 two-leg Heisenberg ladder, J_parallel = 1, rung coupling J_perp, periodic in the
rung direction, one rung (dimension 4) per cell. See [REPORT.md](REPORT.md) for the exact RP
identity, the enclosure theorem, the numerical threshold search and timings, and
[IDEAS.md](IDEAS.md) for the routes considered.

**Certified (J_perp = 4):** Δ_L ≥ 1.10240309174412 for every even L ≥ 4 rungs
(N = 2L ≥ 8 spins, N divisible by four).

- RP basis: checkerboard rotation exp(iπS^y) on one spin of every bond; in the rotated frame
  −h_bond = Σ C⊗C over C ∈ {X⊗I, A⊗I, Z⊗I, I⊗X, I⊗A, I⊗Z} with X = S^x, A = iS^y,
  Z = S^z (all real), and real symmetric rung term J_perp(−X⊗X − A⊗A − Z⊗Z). Checked over
  `Rational{BigInt}` by `rp.jl`.
- Seed: (n0, β0) = (4, 2), εₛ = 3508321//100000000, εₜ = 123121//5000000
  (`results/certificate.json`), from complete verified dense spectra of 4 and 8 rungs.
- Bound: `certified_gap_bound(εₛ, εₜ, β0)`; certification took 12 s.

All other couplings, the threshold near J_perp ≈ 2.488 and the finite-ring gap references
are **numerical, not certified**.

## Commands

From the repository root, with `OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=julia models/heisenberg_ladder/test/runtests.jl   # fast checker, ~50 s
julia --project=julia models/heisenberg_ladder/smoke.jl           # independent full matrices, ~20 s
julia --project=julia models/heisenberg_ladder/audit.jl           # MPFR replay, writes results/audit.json
julia --project=julia models/heisenberg_ladder/check_scans.jl     # replays the numerical grids
julia --project=julia models/heisenberg_ladder/certify.jl         # regenerates the certificate, ~15 s
julia --project=julia models/heisenberg_ladder/scan.jl '4//1'     # numerical scan (10 rungs: ~30 min)
```

`check.jl` runs the test file. `scan_lowmem.jl` is the lower-memory scan driver (sparse
real basis, `sparse_scan_basis.jl`), checked by `test_lowmem.jl`; `threshold.jl` bisects the
coupling; `summarize.jl` writes the results section of REPORT.md; `compact_results.jl`
replaces numerical eigenvalue lists by logZ grids and low levels.

The fast checker verifies exact RP at all six couplings, recomputes the bound from the stored
rationals, checks sector dimension accounting of the stored enclosures, recomputes both seeds
from all stored eigenvalue enclosures, and then reruns fresh verified 4- and 8-rung eigensolves
and checks that the stored rational seeds still dominate.

## Storage

`results/verified_L8.json.gz` (65,536 exact eigenvalue enclosures, 6.5 MB uncompressed) is
gzipped; `read_result_text` in `ladder.jl` decompresses it with the system `gzip`.
`certify.jl` writes it in this form.

Removed when the directory moved from `certificates/` to `models/` (directory now 2.2 MB,
was 15 MB):

- `legacy/`: a byte-for-byte snapshot of the pre-rename (`ThermalMPOGap`) sources and
  results, including a second 6.5 MB `verified_L8.json`. It is in git history at
  `e7e364c`; the current certificate was regenerated on `CertifiedGaps` and supersedes it.
- `PR_DESCRIPTION.md` (copy of the pull-request text).
- Earlier, `compact_results.jl` dropped the full numerical spectra of the scans (millions of
  unverified levels); rerun `scan.jl` (or `scan_lowmem.jl`) with the coupling to regenerate.

The certificate's `source_sha256` was recorded at the old path
`certificates/heisenberg_ladder` and is kept as historical provenance; the checker compares
only the shared-core (`julia/`) hashes and re-establishes the seeds by fresh verified
eigensolves. Rerunning `certify.jl` refreshes all hashes.

## What is and is not rigorous

Rigorous: the exact rational RP identity; exact symmetry-sector counting; interval-enclosed
matrix entries; residual/orthogonality (polar-factor + Weyl) eigenvalue enclosures with
outward-rounded norms and explicit IEEE dot-product error allowances; interval partition
sums; the outward-rounded recurrence. Assumptions: the abstract criterion of the repository,
and the standard IEEE round-to-nearest model for the OpenBLAS products (single BLAS thread).
The fresh audit shares the symmetry builder and the residual bound with the producer. Not
claimed: any certified gap at other couplings, a threshold, monotonicity in J_perp, or the
absence of seeds beyond the scanned grid.
