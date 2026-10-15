# CertifiedGaps.jl

The shared Julia core of this repository. Every project under `models/` depends on it as a
path dependency (or runs with `--project=julia`); anything two projects need lives here once.

It holds the numerics of the spacetime bootstrap, the two-purity doubling criterion. A ring with spatial and thermal
purities Pₛ(n, β) = Z₂ₙ(β)/Zₙ(β)² and P_τ(n, β) = Zₙ(2β)/Zₙ(β)² has a uniform gap for every even
L ≥ n₀ when its transfer operator is self-adjoint and the seed defects at (n₀, β₀) are admissible.

## Layout

Each submodule is used on its own, so generic names never collide with a project's:
`using CertifiedGaps.ExactDiagonalization`. The top level re-exports `Recurrence` only.

| submodule | file | content |
|---|---|---|
| `Recurrence` | `src/Recurrence.jl`, `src/recurrence/` | the exact doubling map in any number type (`map.jl`), the outward-rounded uniform bound (`certified.jl`), the exact-rational chain with the bootstrap tail used by the spin-1 certificates (`rational.jl`) |
| `ExactArithmetic` | `src/ExactArithmetic.jl` | rational helpers with outward rounding: `log_enclosure`, `power_interval`, `integer_root`, `first_true`, `round_up`/`round_down`, exact fraction tokens |
| `Results` | `src/Results.jl` | `write_result` (atomic JSON), `source_hashes`/`file_hashes`, `run_metadata`, `parse_exact_json` (big integers never pass through `Float64`) |
| `Spins` | `src/Spins.jl` | `spin_matrices(s)`, `integer_ladders(j)` and `binomial_weights(j)`, `cartesian_generators(s)` |
| `ReflectionPositivity` | `src/ReflectionPositivity.jl` | `realign`, `rp_test`, `rp_decomposition` |
| `ExactDiagonalization` | `src/ExactDiagonalization.jl` | `Chain` on any bond list (`ring_bonds`, `open_bonds`, `graph_bonds`), charge sectors, `spectrum` (dense or Lanczos, optional spin-flip doubling), `PartialSpectrum`, `logZ_bounds`, `log_partition` |
| `RingModels` | `src/RingModels.jl` | `RingModel(d, bond, onsite; charge, modulus)` with exact entries, the exact RP certificate `rp_certificate` / `verify_rp` (−h = Σ w C ⊗ C̄, by given factors or exact LDL†), exact JSON encoding |
| `VerifiedSpectra` | `src/VerifiedSpectra.jl` | `verified_ring_spectrum`: complete spectrum of a ring by charge and momentum blocks with rigorous eigenvalue radii (residual, Weyl and Ostrowski bounds with a priori rounding errors), `log_partition_enclosure` |
| `GapCertificates` | `src/GapCertificates.jl` | `certify_gap(model; n0, β0, rp)` and `check_certificate`: RP identity, seed defects from verified spectra, exact `rational_gap_chain` |
| `Models` | `src/Models.jl` | `tfim` (with exact `tfim_logZ`), `heisenberg(s)`, `staggered_heisenberg` |
| `Seeds` | `src/Seeds.jl` | `seed_defects` from spectra, `scan_seeds`, `best_seed` |
| `Purification` | `src/Purification.jl`, `ext/CertifiedGapsITensorsExt.jl` | `purified_log_partition`: thermal purification with ITensors; enabled by `using ITensors, ITensorMPS` |
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
julia --project=julia julia/scripts/constants.jl    # recurrence constants, outward rounded
julia --project=julia julia/scripts/tfim.jl         # free-fermion TFIM seeds (BigFloat)
julia --project=julia julia/scripts/ed_seeds.jl     # float ED seed scans (≈15 min)
```

The scripts reproduce numbers of the local survey (its section numbers appear in their
comments); `ed_seeds.jl` is floating-point screening for the staggered-field and spin-1
Heisenberg rings, and `tfim.jl` pads BigFloat defects heuristically, so neither is a
certificate.

From a project under `models/`: `Pkg.develop(path = "../../julia")`, then e.g.

```julia
using CertifiedGaps                                  # doubling map and bounds
using CertifiedGaps.ExactDiagonalization, CertifiedGaps.Models
spec = spectrum(heisenberg(1).chain, 8; spin_flip = true)
certified_gap_bound(1 // 100, 1 // 100, 784)         # (bound = 0.0042911…, …)
```

## Certify a new model

Any translation-invariant ring whose seed rings n₀ and 2n₀ fit in dense exact diagonalisation
(up to about 2¹⁶ states per ring) is certified end to end:

```julia
using CertifiedGaps.RingModels, CertifiedGaps.GapCertificates, CertifiedGaps.Results
S = (x = [0 1//2; 1//2 0], iy = [0 1//2; -1//2 0], z = [1//2 0; 0 -1//2])
bond = -kron(S.x, S.x) - kron(S.iy, S.iy) - kron(S.z, S.z)     # exact entries only
model = RingModel(2, bond, -S.z; charge = [0, 1], modulus = 2)
rp = rp_certificate(model; factors = [S.x, S.iy, S.z])          # or rp_certificate(model)
cert = certify_gap(model; n0 = 8, β0 = 21 // 4, rp)
write_result("certificate.json", cert)
check_certificate(read_exact_json("certificate.json"))
```

The certificate stores the model, the RP factors, the seed (n₀, β₀, εₛ, εₜ) as exact strings,
the log Z enclosures and the rational chain; these are the hypotheses of the Lean ring theorem.
A unit cell of c sites is one site of dimension dᶜ with the intra-cell terms in `onsite`.
Larger seeds need structure (closed forms, expansions, symmetry counting; see `AGENTS.md`).

## What is rigorous

`certified_gap_bound` and `rational_gap_chain` are outward rounded: their output is a lower
bound for the stated seed defects. `ExactArithmetic` returns exact rational enclosures.
`GapCertificates` is rigorous end to end: the RP identity is exact, every eigenvalue comes
with a proven radius, and log Z is enclosed in interval arithmetic. The radii assume only
IEEE binary64 arithmetic with round-to-nearest in BLAS (the classical a priori error bound for
dot products, valid for any summation order and with FMA). Seed defects from `spectrum`
are floating point (screening only), and when a sector falls back to Lanczos they can miss
degenerate partners. `purified_log_partition` is a numerical control.

Run Julia with one BLAS thread (`OPENBLAS_NUM_THREADS=1`) when several jobs share the machine.
