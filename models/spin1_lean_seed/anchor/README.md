# Spin-one comparison anchors

This is an anchor feasibility study for the `(n0,beta0)=(36,49/4)` seed, not a gap
certificate. See [REPORT.md](REPORT.md) for the findings and mathematical qualifications.
The physical Hamiltonian is `sum S_i · S_(i+1)` in the standard Sz basis; the pi twist
negates the two exchange matrix elements on the closing bond. In the seed proof the
reflection-positive representation is the usual sublattice-rotated/Cartesian one;
this study does not reprove reflection positivity or the seed inequalities.

All new files are in this folder. `lean/` and `models/spin1/results/` are unchanged.
Julia 1.12.1 was used, with one BLAS thread per process. No Python is used.

From this directory, using an available Julia 1.12 executable:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
julia --project=. check.jl
```

`Manifest.toml` references the repository core at `../../../julia`. In this cloud session
the executable is `.runtime/julia-1.12.1/bin/julia`, with
`JULIA_DEPOT_PATH="$PWD/.runtime/depot"`. The ignored runtime is not a deliverable.

Reproduce the numerical searches and exact finite checks:

```sh
julia --project=. scan.jl
julia --project=. count_scan.jl
julia --project=. spectra.jl 8 10 12
julia --project=. analyze.jl
julia --project=. tail_bins.jl
julia --project=. certify.jl
julia --project=. check.jl
```

Reproduce the binned-tail continuation for all sectors with nonzero initial counts:

```sh
export ANCHOR_FRESH=1
for n in 8 10 12; do
  for g in I P; do
    julia --project=. homotopy.jl "$n" "$g" binned
  done
done
julia --project=. homotopy_sectors.jl
julia --project=. homotopy_sectors.jl I
julia --project=. homotopy_sectors.jl P
julia --project=. report_tables.jl
```

Without `ANCHOR_FRESH=1`, a completed run is reused and an incomplete run resumes from
its saved stage. The recorded runs include retries and endpoint resumes; a clean run
with the improved eigensolver may use fewer stages and less time.

The continuation scripts are **NON-CERTIFYING Float64 experiments**. Every stage's
retained count, anchor, residual, acceptance margin, and timing is in the corresponding
`results/homotopy_*.json`. Whole runs take seconds to minutes for length 8/10 and minutes
for length 12 on this machine; see the report for actual timings. The unfinished uniform-tail
length-12 pilot is explicitly marked `complete: false`; it was stopped in favor of the
smaller binned-tail calculation.

`certificate.json` contains rational local PSD parameters/pivots, half-chain eigenvalue
enclosures, exact sector counts/anchors, and exact binned-tail arithmetic. `check.jl`
regenerates untrusted floating eigenvector proposals, rounds them to integers, and checks
the resulting finite inequalities exactly. Its checks are not Lean kernel checks. The
eigenvectors are regenerated rather than stored, so the JSON alone is not a standalone
Lean certificate. No certified ring eigenvalue, thermal trace, seed, or gap is claimed here.

The JSON records use `CertifiedGaps.Results`, exact-number strings, source hashes, and
Julia/package metadata. CSV and numerical JSON results are non-certifying unless an exact
scope is explicitly identified. Negative Sz sectors follow by the exact global spin-flip
symmetry and have the same spectra/counts as their positive counterparts.
