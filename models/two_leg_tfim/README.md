# Two-leg transverse-field Ising ladder

This module certifies a uniform spectral gap for the two-leg transverse-field Ising ladder at transverse field $g = 3$: **$\Delta_L \ge 0.636386863$ and a unique ground state for every even $L \ge 4$ rungs** ($N = 2L$ spins, $N$ divisible by 4, $N \ge 8$). The verification method, eigenvalue enclosure theorem, numerical scan tables, and independent replay audit are documented in [`METHOD.md`](METHOD.md); the route analysis is in [`IDEAS.md`](IDEAS.md).

| $g$ | seed $(n_0 \text{ rungs}, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| 3 | $(4, 3/2)$ | 0.636386863 | every even $L \ge 4$ rungs |

## Model

Each unit cell is a two-spin rung ($d = 4$). In Pauli normalization on a periodic ring of $L$ rungs,

$$H_L = -\sum_{j=1}^L (Z_{j,a} Z_{j+1,a} + Z_{j,b} Z_{j+1,b} + Z_{j,a} Z_{j,b}) - g \sum_{j=1}^L (X_{j,a} + X_{j,b}).$$

Lengths count rungs, and the $L = 2$ ring retains both periodic leg bonds, following the core package convention.

## Reflection positivity

In the real $\sigma^z$ rung basis $(\uparrow\uparrow, \uparrow\downarrow, \downarrow\uparrow, \downarrow\downarrow)$, setting the diagonal real bond factors $C_1 = Z \otimes I$ and $C_2 = I \otimes Z$ gives

$$-h_{\text{bond}} = C_1 \otimes C_1 + C_2 \otimes C_2,$$

and the on-site rung operator $o = -Z \otimes Z - g(X \otimes I + I \otimes X)$ is real symmetric. The function `rp_identity` in [`src/Ladder.jl`](src/Ladder.jl) checks this identity over `Rational{BigInt}`, constructing the inter-rung bond independently from four $2 \times 2$ Kronecker factors. Because a single rung is invariant under spatial reflection along the leg direction, no internal cell reflection $F$ is required ($F = I$).

## Seeds and certified bounds

| $g$ | seed $(n_0 \text{ rungs}, \beta_0)$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(4, 3/2)$ | $\varepsilon_\tau \ge 1 - P_\tau(8, 3/2)$ | certified bound | 8-rung ED gap |
|---|---|---:|---:|---:|---:|
| 3 | $(4, 3/2)$ | $\approx 0.11240836635$ | $\approx 0.05161523886$ | 0.636386863 | 2.6618 (numerical) |

The exact rational defect upper bounds are stored as `epsilon_s` and `epsilon_t` in [`results/certificate_g3_1.json.gz`](results/certificate_g3_1.json.gz). Passing them to `CertifiedGaps.certified_gap_bound(εₛ, εₜ, 3//2).bound` yields the exact rational `179126977518527//281474976710656` ($\approx 0.6363868632722607$, rounded down in the table above), which is about $24\%$ of the 8-rung numerical gap. On the scanned field grid $g \in \{2, 5/2, 3, 4\}$, the weaker fields $g = 2$ and $5/2$ have no admissible seed with $n_0 \le 4$ rungs, and $g = 4$ (numerical bound $2.142$) was not run through full certification because $g = 3$ is already certified.

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

Running `certify.jl [g] [OUTPUT]` (default `3//1`) reads the best seed from `results/scan_g<g>.json`
and writes `OUTPUT` (default `results/certificate_g<g>_certifiedgaps.json`, which the fast checker
also verifies). `check.jl --full` runs it at $g = 3$ into a temporary file, never into `results/`,
checks that fresh certificate like the stored ones, and requires the same seed $(g, n_0, \beta_0)$
and a bound at least that of `results/certificate_g3_1.json.gz`. The enclosure radii depend on the
LAPACK run, so a fresh run need not be bit-identical: on a heavily loaded machine (7040 s) it gave
the valid bound 0.636473… against the stored 0.636386…. The scan script writes full spectra to
`results/spectra_g<g>.json` (2.5 MB each, gitignored) for use by `extrapolate.jl` and
`summarize.jl`. `--project=julia` (the core environment) also works for every command.

Measured run times (Julia 1.12.4, one BLAS thread):
- `check.jl` (equivalent to `test/runtests.jl`): 70 tests, 20–35 s including compilation;
- `certify.jl 3//1`: 646 s for the verified spectra of $4$ and $8$ rungs;
- `scripts/replay_audit.jl`: ~450 s; `scripts/check_replay.jl`: ~5 s;
- `scripts/scan.jl`: 90–140 s of dense diagonalization per field.

## What is rigorous

- **Verified enclosures and recurrence:** the exact RP identity; complete dense spectra of all 32 symmetry sectors at $4$ and $8$ rungs (dimension sums $256$ and $65{,}536$), each enclosed by a polar/Weyl whole-spectrum bound that accounts for algebraic matrix-construction error and explicit binary64 BLAS rounding error bounds; outward-rounded MPFR log-partition intervals; exact rational seed defects; and the two-purity doubling recurrence. The enclosures rely on IEEE binary64 round-to-nearest arithmetic with gradual underflow and standard BLAS dot-product error bounds (see [`METHOD.md`](METHOD.md)).
- **Fast vs. full verification:** `check.jl` (fast) uses the stored per-sector eigenvalue radii to recompute $\log Z$, the rational seed defects, and the recurrence bound. `check.jl --full` and `scripts/replay_audit.jl` rebuild every eigenpair and enclosure radius from the Hamiltonian. On the original Linux OpenBLAS machine, `replay_audit.jl` reproduced the stored certificate bit for bit; on an Apple M5, LAPACK tie-breaking produces slightly tighter radii and a valid bound of $0.6364738604216351$ (so `check_replay.jl` flags the bit-identity comparisons while passing all seed, bound, and source-hash checks).
- **Numerical only (not used in the proof):** the seed search ($g \in \{2, 5/2, 3, 4\}$, $n_0 \in \{2, 3, 4\}$, $\beta = 1/8 : 1/8 : 32$), the uncertified $g = 4$ candidate, the negative scan results at $g = 2$ and $5/2$, the extrapolated seed sizes, finite-ring ED reference gaps, and the literature critical field $g_c \approx 1.83214(5)$.

## Layout

- `certify.jl`, `check.jl`: certificate generator and verifier.
- `src/Ladder.jl`: exact RP identity, symmetry sector orbits, block matrices, verified eigenvalue
  enclosures, $\log Z$ intervals, and `seed_defects` (the single implementation of the step from two
  spectra to the rational seeds, shared by `certify.jl`, the tests, and the replay scripts).
- `src/Provenance.jl`: the historical snapshot constants and source-hash verification (the gzip
  certificate is read with the core `Results.read_result` / `read_result_bytes`).
- `scripts/`: `scan.jl` (numerical seed scan), `extrapolate.jl`, `summarize.jl` (updates the
  `## Results` section of `METHOD.md`), `replay_audit.jl`, and `check_replay.jl` (full fixed-seed
  replay and its checker).
- `test/runtests.jl`: fast test suite and certificate checker.
- `results/certificate_g3_1.json.gz`: compressed certificate (gzip of the original 5.3 MB JSON),
  alongside `results/scan_g*.json` (numerical scans), `extrapolation.json`, and `environment.json`
  (Julia, BLAS, and package versions of the certifying run).
- `METHOD.md`: verification method, enclosure theorem, numerical scan tables, and the independent
  audit.

## History

This certificate was generated when the directory was `certificates/two_leg_tfim` and the core
package was named `ThermalMPOGap`. Its `source_sha256` entries reference the audited snapshot
`8a1aad63e1a285392c97bf5cf17e8da03f0b8314`, which the fast checker verifies against the Git object
database when that commit is available (and reports when it is not), along with the SHA-256 hash
(`aba98d98…a74928`) of the uncompressed certificate bytes. The rename changed only the namespace,
and a replay on `CertifiedGaps` 0.2.0 reproduced every number exactly. `METHOD.md` was named
`REPORT.md` until the directory-format standardization. Later the seed step that `certify.jl`,
`test/runtests.jl`, `scripts/replay_audit.jl`, and `scripts/check_replay.jl` each carried moved
into `Ladder.seed_defects` (same operations, same rounding), `certify.jl` gained the optional output
path, and `check.jl --full` stopped writing into `results/`. The stored certificate's provenance
hashes refer to the historical snapshot and are unaffected. Fresh certificates (`certify.jl`,
`replay_audit.jl`) now record only model-local sources under `source_sha256` (compared by the
tests) and the core under `core_sha256` (provenance only). Later, duplicated helpers moved into the shared core `CertifiedGaps`: the gzip
loaders `certificate_bytes`/`read_certificate` became `Results.read_result_bytes`/`read_result`,
and the exact Pauli matrices come from `Spins.exact_pauli`.
