# Spin-one large-D certificate

Model: H=sum S_j·S_(j+1)+D sum (Sz_j)^2, one-site cells and periodic rings.
Read REPORT.md for actual results, numerical versus certified status, runtime,
threshold resolution and limitations. METHOD.md contains the exact RP matrices,
symmetry accounting and a self-contained eigenvalue-enclosure proof.
PROMPT.txt preserves the original research prompt byte for byte.

From the repository root, with Julia 1.12.4 and instantiated core dependencies:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/spin1_largeD/test/runtests.jl
julia --project=julia certificates/spin1_largeD/test/small.jl
julia --project=julia certificates/spin1_largeD/test/symmetry.jl
julia --project=julia certificates/spin1_largeD/scan.jl
julia --project=julia certificates/spin1_largeD/scan.jl 25//8
julia --project=julia certificates/spin1_largeD/certify.jl 4 6
julia --project=julia certificates/spin1_largeD/certify.jl 4 12
julia --project=julia certificates/spin1_largeD/seed.jl 4 6 9//4
julia --project=julia certificates/spin1_largeD/certify.jl 1603//500 6
julia --threads=2 --project=julia certificates/spin1_largeD/certify_parallel.jl 1603//500 12
julia --project=julia certificates/spin1_largeD/seed.jl 1603//500 6 381//125
julia --project=julia certificates/spin1_largeD/report.jl
```

The fast test checks exact local RP, provenance, full-spectrum counts, stored
partition enclosures and rational recurrence inputs. Recomputing the dense
residuals requires the separate certify.jl runs. Numerical scan caches are
explicitly noncertified and never used by seed.jl. Each certification invocation
handles one (D,L) pair; large runs should remain below approximately 45 minutes.

Remote setup used the official 1.12.4 tarball under /tmp because the approved
juliaup installer URL returned HTTP 403. Dependencies were installed with
Pkg.instantiate(); Pkg.precompile(), using JULIA_DEPOT_PATH=/tmp/spin1-julia-depot.
The selected runtime is /tmp/julia-1.12.4/bin/julia. No core project or source
files were edited. The certificate Project.toml is an optional path-dependency
manifest of intent; all recorded runs use --project=julia as requested.

Plots use the standalone Vega specification written by report.jl. This run used
Vega 6.4.0 and vega-cli installed under /tmp/spin1-plot. Render with
`vg2svg results/curve.vega.json results/curve.svg` (or vg2png for PNG) from this
directory. The figure uses only the saved numerical and certified result tables.

## Upstream rename and pinned replay

The computation used core revision
`ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b` (`ThermalMPOGap`). During the run,
upstream main advanced to `22b2759e0790197b83e2571b465aa8dc69ebc0be`, renaming
the core to `CertifiedGaps`. Historical code and hashes are deliberately retained
under `provenance/legacy/`; executable imports now use `CertifiedGaps`. The normal
commands above work on the new head, and `results/validation-current.json` records
its separate 134-check replay. To audit the exact original core and original
checker instead, use these entry points (Julia 1.12.4):

```sh
julia --project=julia certificates/spin1_largeD/check_pinned.jl
julia --project=julia certificates/spin1_largeD/check_pinned.jl small.jl
julia --project=julia certificates/spin1_largeD/check_pinned.jl symmetry.jl
```

It uses `git archive` to materialize the exact original core in a temporary
workspace, copies this certificate there, and runs the requested checker. It
requires the original core commit in the clone's Git history and dependencies
installed in the Julia depot; it makes no network request and changes no files
in this checkout. The original source and core hash assertions remain strict.
For a full ED replay, create a separate checkout of the pinned revision, copy
this certificate directory into it, overlay the archived legacy sources, and run
the original commands above. Preserve the recorded results by using a disposable
checkout for dense replay. The current checker additionally tests byte-identical
recurrence source, package-name-only migration of model/driver code, and exact
reproduction of both bounds. This is not a blanket audit of the entire renamed
package or a second full L=12 diagonalization. The explicit user-requested certificate
location is retained despite the upstream move of other models under `models/`.
