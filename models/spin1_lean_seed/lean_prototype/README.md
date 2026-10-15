# Spin-one Lean seed: kernel feasibility prototype

Read [REPORT.md](REPORT.md) for findings and limitations. This is **not** a completed seed
or a proof of the physical low eigenvalue enclosures. It contains kernel-checked arithmetic
witnesses, proved elementary packing lemmas, precise unproved correctness specifications,
and an adversarial audit of the proposed spectral route.

The project pins Lean `v4.35.0-rc3` and Mathlib
`c55e6e786f49471c72fbddbec5415808896aec1e`, matching the original read-only project.
The compiled prototype imports `Std`; Mathlib's source is used for the API audit.
No `sorry`, `native_decide`, or new axioms are used. Open obligations are **definitions of
propositions**, not assumptions imported into any proof.

With Lean/Lake and Julia 1.12 on PATH:

```sh
lake exe cache get
lake build
julia --startup-file=no --project=. scripts/generate.jl
julia --startup-file=no --project=. scripts/experiments.jl
```

`cache get` was attempted here; the network policy denies `cache.mathlib.org` (HTTP 403).
The default build requires no Mathlib object files and succeeds without this cache.
A local tool installation is in ignored `.tools/`. To reproduce these recorded benchmarks:

```sh
export ELAN_HOME="$PWD/.tools/elan"
export PATH="$PWD/.tools/elan/bin:$PATH"
export JULIA_DEPOT_PATH="$PWD/.tools/julia_depot"
lake build
scripts/benchmark.sh
scripts/audit_mathlib.sh
```

`benchmark.sh` compiles each theorem afresh with `decide +kernel`, sequentially. Timings
include imports and elaboration; they are not isolated kernel reduction times. Peak RSS
covers the child process tree. The runner records preexisting child RSS to detect inherited
resource statistics. `measure_direct.jl` provides optional sampled process-only measurements.

`generate.jl` is an untrusted Julia BigInt generator with independent sparse-row checks.
It regenerates small Lean recipe files, not lists of vector entries. Its TSVs are evidence,
not trusted certificates. All Float64 outputs have `noncertifying` in their filename/header.
The original `lean/` and `models/spin1/results/` are untouched.

RP convention: the executable H uses the standard Sz basis; the untwisted even-ring power
filter uses Marshall conjugation. The existing gap proof uses Cartesian spins. Their unitary
bridge, and the twisted-ring extension, remain open here. No new gap bound or seed is claimed.

Additional exact regression and inventory checks:

```sh
julia --startup-file=no --project=. test/runtests.jl
julia --startup-file=no scripts/finalize.jl
julia --startup-file=no check.jl
lake env lean SeedPrototype/AxiomAudit.lean
```

`certificate.json` is an inventory of the limited arithmetic witnesses, with exact rational
strings and source hashes. It explicitly disclaims physical spectral/gap certification.
The report's benchmark TSV is `benchmarks/summary.tsv`. Download the paper excerpts with
`scripts/fetch_references.sh`; downloaded upstream sources are ignored by Git.
