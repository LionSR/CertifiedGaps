# Spin-1 dimerized gap certificates

See [REPORT.md](REPORT.md) for the exact reflection-positive decomposition, the
positive two-site transfer operator, eigenvalue-enclosure proof, scan results,
verified bounds, length ranges, runtimes and limitations.

From the repository root, using Julia 1.12.4:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/spin1_dimerized/scan.jl 7 2,3,4,6
julia --project=julia certificates/spin1_dimerized/certify.jl 7
julia --project=julia certificates/spin1_dimerized/test/runtests.jl
```

The first argument is ten times delta. Scan 9,8,7,6,5 for the requested grid.
The scan caches complete dense spectra and writes `results/scan_d*.json` through
`ThermalMPOGap.Results.write_result`. The verification script reads the best
stored seed and recomputes all spectra needed for that seed with verified errors.
`check.jl` is an alias for the fast checker. `test/verified_blocks.jl` additionally
compares the verified block implementation with direct small-ring spectra.

The certificate checker deliberately treats the stored rational seed enclosures
as inputs: it verifies their conversion from the stored log-partition bounds,
source hashes and the recurrence, but does not repeat the long dense eigensolves.
Full reproducibility of the spectral enclosures requires running `certify.jl`.

`report.jl` regenerates the numerical comparison tables once all parameter scans
are complete. `IDEAS.md` records the initial routes and the positive-transfer
refinement that permits three-cell seeds. All non-certificate spectra and seed
searches are numerical evidence, not rigorous spectral enclosures.

In this selected cloud environment, the installed runtime and depot are temporary:

```sh
export PATH="/tmp/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH=/tmp/julia-depot
```

Set these before the commands above. On another machine use its Julia 1.12.4
installation and instantiate the existing `julia` project normally.
