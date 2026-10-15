# Two-leg transverse-field Ising ladder

**Certified at g = 3: Δ_L ≥ 0.636386863 and a unique ground state for every even L ≥ 4 rungs**
(N = 2L spins, N divisible by four, N ≥ 8). The method, the enclosure theorem, the numerical
tables and the independent audit are in [METHOD.md](METHOD.md); the route memo is
[IDEAS.md](IDEAS.md).

| g | seed (n₀ rungs, β₀) | certified gap ≥ | lengths covered |
|---|---|---:|---|
| 3 | (4, 3/2) | 0.636386863 | every even L ≥ 4 rungs |

## Model

One cell is a rung (two spins, cell dimension 4). In Pauli normalization, on a ring of L rungs,

    H_L = −Σ_j (Z_{j,a}Z_{j+1,a} + Z_{j,b}Z_{j+1,b} + Z_{j,a}Z_{j,b}) − g Σ_j (X_{j,a} + X_{j,b}).

Lengths count rungs. The two-rung ring keeps both periodic bonds, as in the core package.

## Reflection positivity

Real σᶻ cell basis (↑↑, ↑↓, ↓↑, ↓↓). With C₁ = Z⊗I and C₂ = I⊗Z (real, diagonal),
−h_bond = C₁⊗C₁ + C₂⊗C₂ exactly, and the on-site term o = −Z⊗Z − g(X⊗I + I⊗X) is real
symmetric. `rp_identity` in `src/Ladder.jl` checks this with `Rational{BigInt}` matrices,
building the bond independently from four 2×2 Kronecker factors. A single rung is its own
mirror image, so no palindromic-cell condition arises.

## Seeds and certified bounds

| g | seed (n₀ rungs, β₀) | εₛ ≥ 1 − Pₛ(4, 3/2) | εₜ ≥ 1 − P_τ(8, 3/2) | certified bound | 8-rung ED gap |
|---|---|---:|---:|---:|---:|
| 3 | (4, 3/2) | ≈ 0.11240836635 | ≈ 0.05161523886 | 0.636386863 | 2.6618 (numerical) |

The exact rational defects are `epsilon_s` and `epsilon_t` in
`results/certificate_g3_1.json.gz`. The bound is `CertifiedGaps.certified_gap_bound(εₛ, εₜ,
3//2).bound` = `0.6363868632722607` = `179126977518527//281474976710656`; the decimal above is
rounded down. It is about 24% of the 8-rung gap. On the grid g ∈ {2, 5/2, 3, 4} the fields 2
and 5/2 have no numerical seed with n₀ ≤ 4, and g = 4 (numerical bound 2.142) was not
certified.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=models/two_leg_tfim -e 'using Pkg; Pkg.instantiate()'
julia --project=models/two_leg_tfim models/two_leg_tfim/check.jl           # fast checker
julia --project=models/two_leg_tfim models/two_leg_tfim/check.jl --full    # + fresh recertification
# Fixed-seed replay from model inputs and its checker (output paths must be new):
julia --project=models/two_leg_tfim models/two_leg_tfim/scripts/replay_audit.jl results/replay.json
julia --project=models/two_leg_tfim models/two_leg_tfim/scripts/check_replay.jl \
    results/replay.json results/replay_check.json
# Numerical seed scan for one field, then the extrapolation and the METHOD.md tables:
julia --project=models/two_leg_tfim models/two_leg_tfim/scripts/scan.jl 3//1
julia --project=models/two_leg_tfim models/two_leg_tfim/scripts/extrapolate.jl
julia --project=models/two_leg_tfim models/two_leg_tfim/scripts/summarize.jl
```

`certify.jl [g] [OUTPUT]` (default `3//1`) reads the best seed from `results/scan_g<g>.json` and
writes `OUTPUT` (default `results/certificate_g<g>_certifiedgaps.json`, which the fast checker
also verifies). `check.jl --full` runs it at g = 3 into a temporary file, never into
`results/`, checks that fresh certificate like the stored ones, and requires the same seed
(g, n₀, β₀) and a bound at least that of `results/certificate_g3_1.json.gz`. The enclosure
radii depend on the LAPACK run, so a fresh run need not be bit-identical: on a heavily loaded
machine (7040 s) it gave the valid bound 0.636473… against the stored 0.636386…. The
scan writes its spectra to `results/spectra_g<g>.json` (2.5 MB each, gitignored);
`extrapolate.jl` and `summarize.jl` need them. `--project=julia` (the core environment) also
works for every command.

Run times (Julia 1.12.4, one BLAS thread):

- fast checker (`check.jl`, = `test/runtests.jl`): 70 tests, 20–35 s including compilation;
- `certify.jl 3//1`: 646 s for the verified spectra of 4 and 8 rungs;
- `scripts/replay_audit.jl`: about 450 s; `scripts/check_replay.jl`: about 5 s;
- `scripts/scan.jl`: 90–140 s of dense spectra per field.

## What is rigorous

Rigorous: the exact RP identity; complete dense spectra of all 32 symmetry sectors at 4 and
8 rungs (dimension sums 256 and 65,536), each with a polar/Weyl whole-spectrum enclosure that
includes the algebraic matrix-construction error and explicit binary64 BLAS error terms;
outward MPFR log-partition intervals; exact rational seeds; the recurrence. No Lanczos or
ARPACK spectrum enters. The enclosures assume IEEE binary64 round-to-nearest with gradual
underflow and conventional BLAS dot products (METHOD.md, "Implemented enclosure theorem").

The fast checker trusts the stored per-sector radii and recomputes log Z, the rational seeds
and the bound from them. `check.jl --full` and `scripts/replay_audit.jl` rebuild every
eigenpair and radius from the model; the replay recorded in METHOD.md (Julia 1.12.4, Linux
OpenBLAS) reproduced the certificate bit for bit. `check_replay.jl` demands that bit identity,
so on other platforms it reports mismatches: on an Apple M5 the fresh replay gave a different
but valid bound, 0.6364738604216351 (stronger than the stored one), and every seed
recomputation and source hash in the replay checked out.

Numerical only: the seed search (g ∈ {2, 5/2, 3, 4}, n₀ ∈ {2, 3, 4}, β = 1/8:1/8:32), the
g = 4 candidate, the failures at g = 2 and 5/2, the extrapolated larger seeds, ED reference
gaps, and the literature critical field g_c ≈ 1.83214(5) (not a premise).

## Layout

- `certify.jl`, `check.jl`: certificate entry point and checker.
- `src/Ladder.jl`: RP identity, symmetry blocks, block matrices, verified eigen-enclosures,
  log Z intervals, and `seed_defects` (the one implementation of the step from two spectra to
  the rational seeds, shared by `certify.jl`, the tests and the replay scripts). `src/Provenance.jl`: gzip certificate loader and historical source hashes.
- `scripts/`: `scan.jl` (numerical seed scan), `extrapolate.jl`, `summarize.jl` (the
  "Results" section of METHOD.md), `replay_audit.jl` and `check_replay.jl` (full fixed-seed
  replay and its checker).
- `test/runtests.jl`: the fast checker.
- `results/certificate_g3_1.json.gz`: the certificate (gzip of the original 5.3 MB JSON).
  `results/scan_g*.json` (numerical scans), `extrapolation.json`, `environment.json` (Julia,
  BLAS and package versions of the certifying run).
- `METHOD.md`: method, enclosure theorem, numerical tables and the independent audit.

## History

The certificate was made when this directory was `certificates/two_leg_tfim` and the core was
named `ThermalMPOGap`. Its `source_sha256` entries refer to the audited snapshot
`8a1aad63e1a285392c97bf5cf17e8da03f0b8314`; the fast checker verifies them against the Git
objects of that commit when it is available (and says so when it is not), and checks the
SHA-256 `aba98d98…a74928` of the uncompressed certificate bytes. The rename changed only the
namespace, and a replay on `CertifiedGaps` 0.2.0 reproduced every number exactly.
`METHOD.md` was named `REPORT.md` until the format tidy-up. Later the seed step that
`certify.jl`, `test/runtests.jl`, `scripts/replay_audit.jl` and `scripts/check_replay.jl` each
carried moved into `Ladder.seed_defects` (same operations, same rounding), `certify.jl` gained
the optional output path, and `check.jl --full` stopped writing into `results/`. The stored
certificate's provenance hashes refer to the historical snapshot and are unaffected.
