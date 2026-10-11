# Heisenberg ladder certificate

See [REPORT.md](REPORT.md) for the exact RP identity, enclosure theorem, certified
bound, validity range, numerical threshold search and timings.

From the repository root, using Julia 1.12.4:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/heisenberg_ladder/scan.jl '4//1'
julia --project=julia certificates/heisenberg_ladder/certify.jl
julia --project=julia certificates/heisenberg_ladder/test/runtests.jl
julia --project=julia certificates/heisenberg_ladder/recheck_seeds.jl
```

Run `scan.jl` also for `7//2`, `3//1`, `5//2`, `2//1`, `3//2`.
`smoke.jl` compares symmetry spectra with independent full matrices.
No Lanczos, ARPACK or DMRG is used. `certify.jl` reruns only the best valid even
seed at J_perp=4 through the rigorous path. Scans are numerical.

In the selected environment the installer endpoint returned HTTP 403, so the
official Julia 1.12.4 Linux binary was downloaded to `/tmp/julia-1.12.4`.
The prepared depot is `/tmp/ladder-julia-depot`; prepend its bin directory to
PATH and set `JULIA_DEPOT_PATH` to that depot to reuse this installation.
Dependencies were instantiated and precompiled from the unchanged `julia/` project.

For a lower-memory numerical scan, use `scan_lowmem.jl` with the same coupling
argument. This was used for J_perp=2 and 3/2 and is checked by `test_lowmem.jl`.
After all scans and certification, `summarize.jl` generates the results tables;
`compact_results.jl` retains numerical logZ grids and low levels instead of
millions of unverified individual eigenvalues. The rigorous full eigenvalue
enclosures remain stored. Run `small_summary.jl` before compaction if desired.

This PR targets upstream `22b2759` and imports `CertifiedGaps`, the renamed
core with the same package UUID and unchanged recurrence. Historical numerical
scan metadata was produced before the rename. Fresh certification and checker
results use the current core; no core files are modified by this PR.
