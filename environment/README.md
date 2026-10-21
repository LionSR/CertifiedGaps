# Environment

Run commands from the repository root. Use Julia 1.12 and the Lean version in
[`lean-toolchain`](../lean-toolchain). Set `OPENBLAS_NUM_THREADS=1` for the numerical checks.
The folder migration was checked locally on macOS with Julia 1.12.4; the existing CI runs on
Ubuntu with Julia 1.12.
Use a checkout that preserves relative symlinks; model `results/` paths link to `data/models/`.

## Canonical dependency specifications

| Component | Specification |
|---|---|
| Julia core and most models | [`code/julia/Project.toml`](../code/julia/Project.toml), [`Manifest.toml`](../code/julia/Manifest.toml) |
| Models with pinned environments | `code/models/{spin1,spin2,aklt,lean_seeds}/{Project,Manifest}.toml` |
| Leaderboard validator | [`code/tools/leaderboard/Project.toml`](../code/tools/leaderboard/Project.toml), [`Manifest.toml`](../code/tools/leaderboard/Manifest.toml) |
| Lean proofs and leaderboard | Root [`lakefile.toml`](../lakefile.toml), [`lake-manifest.json`](../lake-manifest.json), [`lean-toolchain`](../lean-toolchain) |
| Website | [`home_page/Gemfile`](../home_page/Gemfile), [`Gemfile.lock`](../home_page/Gemfile.lock) |

These files stay beside their packages or at the root required by their build tools.
Moving the pinned model manifests would change certificate source hashes and their relative
package dependencies. This directory documents the canonical files without duplicating them.

## Julia

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=code/julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=code/julia code/models/heisenberg_ladder/check.jl

julia code/tools/fetch_upstream.jl
julia --project=code/models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=code/models/spin1 code/models/spin1/check.jl

julia --project=code/tools/leaderboard -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=code/tools/leaderboard code/tools/leaderboard/check_records.jl
```

Use `--project=code/julia` for models without their own `Project.toml`. For a model with
its own environment, use `--project=code/models/<slug>`. Check each model's README for
full replay commands and expected runtime.

## Lean and documentation

Install `elan`, then:

```sh
lake exe cache get
lake build
lake exe leaderboard_audit home_page/_data/lean_records.json
julia --project=code/tools/leaderboard code/tools/leaderboard/check_records.jl \
  --lean-records home_page/_data/lean_records.json
```

The root Lake configuration points at `code/lean/` and `code/leaderboard/`.
The standalone configuration in `code/lean/` remains available to the seed-generation scripts.
The [blueprint workflow](../.github/workflows/blueprint.yml) defines the documentation build,
and the [home-page workflow](../.github/workflows/homepage.yml) defines the lightweight site build.

The paper uses LaTeX with REVTeX 4.2 and `latexmk`; the root `.latexmkrc` finds `paper/tnnotes.sty`:

```sh
mkdir -p work/latex-build
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=work/latex-build paper/thermal-mpo-gap-note.tex
```

See the [explainer README](../supplementary/explainer/README.md) for optional browser viewing
and film-rendering dependencies.
