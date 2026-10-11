# Dimerized Heisenberg: verified gaps under general reflected RP

Verified uniform bounds are stored for delta=0.5, 0.3, and 0.25, valid for
**every even L>=4 two-spin cells** (8,12,16,... spins), under the repository's
general internally reflected RP theorem. The original stronger requirement
of real identical factors `C tensor C` is false and is not claimed.

Read [REPORT.md](REPORT.md) for results, scope, dense-sector diagnostics and
scaling comparisons, and [ENCLOSURE_PROOF.md](ENCLOSURE_PROOF.md) for the
residual/polar/Weyl theorem, arithmetic error model and recurrence proof.

The original certificates retain their audited `ThermalMPOGap` provenance at
commit `ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b`. Current main calls the core
`CertifiedGaps`. Use the replay wrapper for historical scripts: it extracts
the pinned core and notes from local git objects into an ignored directory,
without modifying the current checkout. A shallow clone must fetch that
commit first. `replay_pinned.jl --setup` installs its pinned dependencies.

From repository root with Julia 1.12.4 and core dependencies:

```sh
# Historical complete test suite against its original, hashed core.
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/replay_pinned.jl
# Separate reconstruction with the current renamed core.
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/validate_migration.jl
```

The following original commands apply at the historical snapshot. On current
main, insert `certificates/dimerized_heisenberg/replay_pinned.jl` after
`--project=julia` and pass the script path relative to this certificate
directory, e.g. `run_verified.jl 1 2`:

```sh
# Exact local reflection and alternative-convention audit.
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/certify.jl
# Reproduce verified spectra, rational purity seeds and recurrence (about a minute each).
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/run_verified.jl 1 2
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/run_verified.jl 3 10
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/run_verified.jl 1 4
# Fast independent tests plus stored-witness/partition/seed/bound reconstruction.
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/test/runtests.jl
# Optional exploratory grid and heuristic extrapolation.
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/numerical_ed.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/summarize_numerics.jl
```

Local Julia: `certificates/dimerized_heisenberg/julia-1.12.4/bin/julia`.
Set `JULIA_DEPOT_PATH` to the absolute path of
`certificates/dimerized_heisenberg/depot` to reuse installed dependencies.
Installation artifacts are git-ignored.

`verified_delta_*.json` contains exact rational spectral centers, per-block
error witnesses, partition-function intervals and purity seeds. Results use
`ThermalMPOGap.Results.write_result`, source hashes and run metadata. Their
certification scope is explicitly **general reflected RP**, not the original
stronger condition. `rp_audit.json` is a separate algebra audit;
`numerical_*.json` remains explicitly numerical/non-certified.

The fast checker does not repeat dense eigenvector solves. It reconstructs
all scalar enclosure witnesses, complete partition sums, rational seeds and
recurrence bounds from stored data. Repeating `run_verified.jl` reproduces
the eigensolves and residual computations themselves.
