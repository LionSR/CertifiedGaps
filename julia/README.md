# ThermalMPOGap.jl

The shared Julia core of this repository. `spin1/`, `spin2/` and `aklt/` depend on it as a path
dependency; anything two projects need lives here once.

It holds the numerics of the two-purity doubling criterion. A ring with spatial and thermal
purities Pₛ(n, β) = Z₂ₙ(β)/Zₙ(β)² and P_τ(n, β) = Zₙ(2β)/Zₙ(β)² has a uniform gap for every even
L ≥ n₀ when its transfer operator is self-adjoint and the seed defects at (n₀, β₀) are admissible.

## Layout

Each submodule is used on its own, so generic names never collide with a project's:
`using ThermalMPOGap.ExactDiagonalization`. The top level re-exports `Recurrence` only.

| submodule | file | content |
|---|---|---|
| `Recurrence` | `src/Recurrence.jl`, `src/recurrence/` | the exact doubling map in any number type (`map.jl`), the outward-rounded uniform bound (`certified.jl`), the exact-rational chain with the bootstrap tail used by the spin-1 certificates (`rational.jl`) |
| `ExactArithmetic` | `src/ExactArithmetic.jl` | rational helpers with outward rounding: `log_enclosure`, `power_interval`, `integer_root`, `first_true`, `round_up`/`round_down`, exact fraction tokens |
| `Results` | `src/Results.jl` | `write_result` (atomic JSON), `source_hashes`/`file_hashes`, `run_metadata`, `parse_exact_json` (big integers never pass through `Float64`) |
| `Spins` | `src/Spins.jl` | `spin_matrices(s)`, `integer_ladders(j)` and `binomial_weights(j)`, `cartesian_generators(s)` |
| `ReflectionPositivity` | `src/ReflectionPositivity.jl` | `realign`, `rp_test`, `rp_decomposition` |
| `ExactDiagonalization` | `src/ExactDiagonalization.jl` | `Chain` on any bond list (`ring_bonds`, `open_bonds`, `graph_bonds`), charge sectors, `spectrum` (dense or Lanczos, optional spin-flip doubling), `PartialSpectrum`, `logZ_bounds`, `log_partition` |
| `Models` | `src/Models.jl` | `tfim` (with exact `tfim_logZ`), `heisenberg(s)`, `staggered_heisenberg` |
| `Seeds` | `src/Seeds.jl` | `seed_defects` from spectra, `scan_seeds`, `best_seed` |
| `Purification` | `src/Purification.jl`, `ext/ThermalMPOGapITensorsExt.jl` | `purified_log_partition`: thermal purification with ITensors; enabled by `using ITensors, ITensorMPS` |
| `CertifiedMPO` | `src/CertifiedMPO/` | certified MPO arithmetic with Nemo/Arb (norms, squaring, products, operator errors, gauges), moved from spin2 |

## Conventions

* `spin_matrices(s)`: basis m = s, …, −s; `sz`, `sp`, `sm`, `sx` real, `sy` complex.
* `integer_ladders(j)`: basis m = −j, …, j in the binomial gauge, integer entries.
* `ExactDiagonalization`: site 1 is the most significant base-d digit; a two-site term acts on
  `kron(site i, site j)`. Spectra do not depend on this; code that stores basis indices keeps its
  own encoding.

## Use

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=julia julia/scripts/constants.jl    # constants of survey §2
julia --project=julia julia/scripts/tfim.jl         # free-fermion TFIM seeds
julia --project=julia julia/scripts/ed_seeds.jl     # ED seed searches (≈15 min)
```

From another project: `Pkg.develop(path = "../julia")`, then e.g.

```julia
using ThermalMPOGap                                  # doubling map and bounds
using ThermalMPOGap.ExactDiagonalization, ThermalMPOGap.Models
spec = spectrum(heisenberg(1).chain, 8; spin_flip = true)
certified_gap_bound(1 // 100, 1 // 100, 784)         # (bound = 0.0042911…, …)
```

## What is rigorous

`certified_gap_bound` and `rational_gap_chain` are outward rounded: their output is a lower
bound for the stated seed defects. `ExactArithmetic` returns exact rational enclosures. Seed
defects from `spectrum` are floating point; they are exact enclosures when every sector is
diagonalised densely (`PartialSpectrum.exact`), and not certificates when a sector falls back to
Lanczos, which can miss degenerate partners. `purified_log_partition` is a numerical control.

Run Julia with one BLAS thread (`OPENBLAS_NUM_THREADS=1`) when several jobs share the machine.
