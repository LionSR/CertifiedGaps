# Spin-1 finite seed redesign

Read [REPORT.md](REPORT.md) for the recommendation and the numerical Pareto front.
This project searches and checks **finite scalar seed designs**, not an unconditional Lean proof.
The Hamiltonian is the periodic spin-1 Heisenberg ring with coupling one; the reflection-positive
Cartesian basis and existing all-scale moment inequalities are those of the repository.

Recommended design: nine shifted traces at lengths 6, 8, 10, twists identity, pi, and 2pi/3;
base beta `21/2`; shift `700741/500000`; absolute table error `3/10000`; the existing length-120
integer MPS; lifted beta `105/4`; seed length 60. The exact scalar replay gives gap
`> 0.005966152` for every even length at least 60, **conditional on the physical finite inputs
and analytic sector/filter/lifting lemmas**. Those premises have not been checked in Lean here.

Julia is an untrusted generator. `certificate.json` stores exact rational scalar checks and
`compact_witnesses.json` stores smaller replay recipes. `check.jl` recomputes the scalar
inequalities and the conditions of `ChainCertificate.Valid` with `Rational{BigInt}`.
It audits stored input hashes and rational enclosures; it does not repeat the original thermal
contractions or certify any Lanczos spectrum. The fresh length-144 MPS contraction is exact
Julia/Nemo arithmetic and is likewise not a Lean proof.

From the repository root with Julia 1.12:

```sh
julia --project=models/spin1_lean_seed/seed_design -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/check.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/test/runtests.jl
```

Regenerate the selected exact designs, using the checked-in `trial_scan.json` for length 144:

```sh
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/certify.jl
```

For fresh discovery and a fresh length-144 contraction, set `OPENBLAS_NUM_THREADS=1` and run:

```sh
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/fetch_paper.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/trial_scan.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/spectra.jl 4 5 6 7 8 9 10 11 12
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/search_existing.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/search_new.jl 8
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/search_new.jl 10
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/refine.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/diagnostics.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/shift_search.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/certify.jl
julia --project=models/spin1_lean_seed/seed_design models/spin1_lean_seed/seed_design/summarize.jl
```

All Float64 ED, Lanczos, search results, workload estimates, and `pareto.csv` are explicitly
non-certifying. The spectra cache is regenerable and ignored by Git. Spectra through length 8
are full dense diagonalizations of every nonnegative magnetization block; larger blocks use
Lanczos and have no certified eigenvalue-counting anchor. `summarize.jl` refuses a low-state
count estimate if its requested energy cutoff lies beyond the computed Ritz list.

For this managed environment Julia was installed locally at
`.runtime/julia-1.12.4/bin/julia` with depot `.runtime/depot`, both under this directory.
Use those paths with `JULIA_DEPOT_PATH` if Julia is not on `PATH`. The environment setup and
all generated outputs remain in this task directory. No file under `lean/` or
`models/spin1/results/` was modified.

The project includes the core's `ExactArithmetic.jl` and `Results.jl` directly to reuse their
exact arithmetic and JSON helpers without loading unrelated numerical packages. Original
MPS entry and charge-block routines are included read-only from `models/spin1/src/trials.jl`.
Package versions are pinned in `Manifest.toml`; output metadata records versions and hashes.
