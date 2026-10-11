# Two-leg transverse-field Ising ladder

Certified at **g=3**, with **(n0,β0)=(4,3/2)**: **Δ_L ≥ 0.636386863** for every even L≥4 rungs (N=2L sites, N divisible by four and N≥8). This is the smallest successful field in the tested grid. Exact RP uses the σᶻ basis, C1=Z⊗I and C2=I⊗Z, each with coefficient one. Other field scans and reference gaps are numerical.

See [REPORT.md](REPORT.md) for the mathematics and scope, [AUDIT-2026-10-11.md](AUDIT-2026-10-11.md) for the historical independent replay, and [MIGRATION-2026-10-11.md](MIGRATION-2026-10-11.md) for validation against the current **CertifiedGaps** package. Original certificate and audit JSONs are preserved byte-for-byte.

Run from the repository root with Julia 1.12.4:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
# Exact RP, stored partition enclosures, rational seeds, and current recurrence:
julia --project=julia certificates/two_leg_tfim/test/runtests.jl
# Current-source hashes and independently reconstructed migrated replay comparisons:
julia --project=julia certificates/two_leg_tfim/check_replay.jl results/migration_2026-10-11_replay.json results/my_check.json
# Regenerate every eigenpair and every spectral error bound from the fixed model:
julia --project=julia certificates/two_leg_tfim/replay_audit.jl results/my_replay.json
julia --project=julia certificates/two_leg_tfim/check_replay.jl results/my_replay.json results/my_replay_check.json
```

Choose unused output names: replay and replay-check scripts reject existing outputs.
The fast checker does not regenerate spectral residuals. The full replay does, using
all 32 dense sectors at L=4,8; no Lanczos or ARPACK enters the partition functions.

Historical hashes are checked against Git snapshot
`8a1aad63e1a285392c97bf5cf17e8da03f0b8314`; current replay hashes are checked against
current files. A complete Git checkout is required. If using a shallow clone, fetch
its history (`git fetch --unshallow`) first. A source archive alone cannot verify
historical Git blobs.

The original numerical scans can be reproduced with `scan.jl 2//1`, `scan.jl 5//2`,
`scan.jl 3//1`, and `scan.jl 4//1` using the same project command. Those commands
replace their scan output files; use a disposable checkout to preserve historical
evidence. `certify.jl 3//1` uses the stored winning seed and writes a separate
`certificate_g3_1_certifiedgaps.json`, which the fast checker also discovers.

In the selected cloud environment, Julia is `/tmp/julia-1.12.4/bin/julia` and
`JULIA_DEPOT_PATH=/tmp/tfim-julia-depot`. Prepend that binary directory to PATH and
export the depot path. The home directory is read-only; the approved juliaup installer
returned HTTP 403, while the official versioned Julia binary download succeeded.
