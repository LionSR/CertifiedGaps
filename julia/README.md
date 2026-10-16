# CertifiedGaps.jl

`CertifiedGaps` is the shared Julia core package of this repository. Every model directory under
`models/` depends on it as a local path dependency (or runs directly with `--project=julia`), so
functionality shared across models lives here once.

The package implements the numerical and algebraic machinery of the **spacetime bootstrap** (the
two-purity doubling criterion): a periodic ring with spatial purity
$P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$ and thermal purity
$P_\tau(n, \beta) = Z_n(2\beta)/Z_n(\beta)^2$ has a unique ground state and a certified uniform lower
bound on its spectral gap $\Delta_L$ for every even ring length $L \ge n_0$ whenever its bond transfer
operator is self-adjoint (reflection positive) and the seed purity defects
$\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ and
$\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$ lie in the admissible basin of the doubling recurrence.

## Layout

Each submodule is imported explicitly (for example, `using CertifiedGaps.ExactDiagonalization`) so that
identifiers do not collide with model-specific code. The top-level module `CertifiedGaps` re-exports
`Recurrence` only.

| submodule | file | content |
|---|---|---|
| `Recurrence` | `src/Recurrence.jl`, `src/recurrence/` | the exact doubling recurrence in any number type (`map.jl`), outward-rounded interval gap bounds (`certified.jl`), and the exact-rational doubling chain with the analytic tail bound used by the certificates (`rational.jl`) |
| `ExactArithmetic` | `src/ExactArithmetic.jl` | rational enclosure helpers with outward rounding: `log_enclosure`, `power_interval`, `integer_root`, `first_true`, `round_up`/`round_down`, and exact fraction parsers |
| `Results` | `src/Results.jl` | `write_result` (atomic JSON output; pretty, or single-line with `compact = true`; gzip-compressed when the path ends in `.gz`), `read_result`/`read_result_bytes` (plain or `.json.gz`), `source_hashes`/`file_hashes`, `model_provenance` (model-local `source_sha256` and core `core_sha256`, keyed by repository path) and `verify_source_hashes` (compares the in-scope paths only), `run_metadata`, and `parse_exact_json` (which parses arbitrary-precision integers and rationals directly without passing through `Float64`) |
| `Spins` | `src/Spins.jl` | `spin_matrices(s)`, `integer_ladders(j)`, `binomial_weights(j)`, and `cartesian_generators(s)`; exact matrices for RP identities: `exact_pauli(T)`, `exact_spin_half(T)` (with the real antisymmetric `iy`) and `integer_cartesian_generators()` (spin 1) |
| `ReflectionPositivity` | `src/ReflectionPositivity.jl` | floating-point screening utilities: `realign`, `rp_test`, and `rp_decomposition` |
| `ExactDiagonalization` | `src/ExactDiagonalization.jl` | `Chain` on arbitrary bond graphs (`ring_bonds`, `open_bonds`, `graph_bonds`), abelian charge sectors, `spectrum` (dense or Lanczos, with optional spin-flip doubling), `PartialSpectrum`, `WeightedSpectrum` (levels with integer multiplicities), `logZ_bounds`, and `log_partition` |
| `RingModels` | `src/RingModels.jl` | `RingModel(d, bond, onsite; charge, modulus)` with exact rational entries, exact reflection-positivity certificates `rp_certificate` and `verify_rp` ($-h = \sum_\alpha w_\alpha C_\alpha \otimes \overline{C_\alpha}$ via explicit factors or exact rational $LDL^\dagger$), and exact JSON serialization |
| `VerifiedSpectra` | `src/VerifiedSpectra.jl` | `verified_ring_spectrum`: complete block-diagonal spectrum of a ring across charge and translation-momentum sectors with rigorous eigenvalue error radii (residual, Weyl, and Ostrowski bounds with a priori IEEE binary64 rounding bounds), plus `log_partition_enclosure` |
| `GapCertificates` | `src/GapCertificates.jl` | end-to-end certification pipeline `certify_gap(model; n0, β0, rp)` and fast/full verifier `check_certificate`: verifies the exact RP identity, encloses the seed defects from verified ring spectra, and checks the exact `rational_gap_chain` |
| `Models` | `src/Models.jl` | standard constructors for screening and testing: `tfim` (including closed-form `tfim_logZ`), `heisenberg(s)`, and `staggered_heisenberg` |
| `Seeds` | `src/Seeds.jl` | floating-point seed screening helpers: `seed_defects(spectra, n, β)`, `scan_seeds`, and `best_seed`; and the rigorous step `seed_defects(logZn, logZ2n, logZ2n_2β)` from three log Z intervals to exact rational upper bounds (εₛ, εₜ) |
| `Purification` | `src/Purification.jl`, `ext/CertifiedGapsITensorsExt.jl` | `purified_log_partition`: thermal MPO purification via ITensors (activated when `using ITensors, ITensorMPS` is loaded) |
| `CertifiedMPO` | `src/CertifiedMPO/` | certified MPO arithmetic over Nemo/Arb ball matrices (norms, squaring, products, operator truncation bounds, and gauge transforms) |

Besides `src/`, the package contains `ext/` (the ITensors extension), `test/` (run by `Pkg.test()`),
and `examples/` (two survey screening drivers that use the API but are not loaded by the package;
see "Usage").

## Conventions

* `spin_matrices(s)`: ordered basis $m = s, s-1, \dots, -s$; `sz`, `sp`, `sm`, and `sx` are real, while `sy` is purely imaginary.
* `integer_ladders(j)`: ordered basis $m = -j, \dots, j$ in the binomial gauge, where the ladder matrices have exact integer entries.
* `ExactDiagonalization`: site $1$ is the most significant base-$d$ digit, and a two-site operator on $(i, j)$ acts as `kron(site i, site j)`. Eigenvalues and partition functions are basis-independent; code that stores explicit basis indices preserves its own digit ordering.

## Usage

From the repository root:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=julia julia/examples/tfim.jl        # free-fermion TFIM seeds (BigFloat, ~1 min)
julia --project=julia julia/examples/ed_seeds.jl    # float ED seed scans (≈15 min)
```

The two examples in `julia/examples/` reproduce reference numbers from the local survey notes (its
section numbers appear in their comments) using the core API; neither produces a rigorous
certificate. `ed_seeds.jl` performs floating-point screening for the staggered-field spin-½ and
spin-1 Heisenberg chains (certified in `models/heisenberg_staggered_field` and reproduced in
`models/spin1`), and `tfim.jl` evaluates free-fermion partition functions in `BigFloat` with
heuristic safety padding for the TFIM open targets of `AGENTS.md`. The recurrence constants that the
former `scripts/constants.jl` printed (the diagonal threshold ε⋆, the note's seeds, odd lengths, and
the Theorem 2.22 closed forms) are now asserted in `test/runtests.jl`.

To use `CertifiedGaps` from a model project under `models/`, link it with
`Pkg.develop(path = "../../julia")`:

```julia
using CertifiedGaps                                  # doubling recurrence and gap bounds
using CertifiedGaps.ExactDiagonalization, CertifiedGaps.Models
spec = spectrum(heisenberg(1).chain, 8; spin_flip = true)
certified_gap_bound(1 // 100, 1 // 100, 784)         # (bound = 0.0042911…, …)
```

## Certifying a new model

Any translation-invariant ring whose seed rings of lengths $n_0$ and $2n_0$ fit within dense exact
diagonalization (up to roughly $2^{16}$ states per ring, reduced by charge and translation symmetry) can
be certified end to end with `GapCertificates`:

```julia
using CertifiedGaps.RingModels, CertifiedGaps.GapCertificates, CertifiedGaps.Results
S = (x = [0 1//2; 1//2 0], iy = [0 1//2; -1//2 0], z = [1//2 0; 0 -1//2])
bond = -kron(S.x, S.x) - kron(S.iy, S.iy) - kron(S.z, S.z)     # exact rational entries
model = RingModel(2, bond, -S.z; charge = [0, 1], modulus = 2)
rp = rp_certificate(model; factors = [S.x, S.iy, S.z])          # or rp_certificate(model)
cert = certify_gap(model; n0 = 8, β0 = 21 // 4, rp)
write_result("certificate.json", cert)
check_certificate(read_exact_json("certificate.json"))
```

The resulting certificate stores the exact model matrices, the reflection-positivity factors, the seed
parameters $(n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ as exact rational strings, the interval
enclosures of $\log Z$, and the exact rational doubling chain (`rational_gap_chain`). These quantities
provide the exact inputs required by the Lean ring theorems in `lean/`. A multi-site unit cell of $c$
spins is represented as a single effective site of local dimension $d^c$ with all intra-cell bonds
collected in `onsite`. Models requiring larger seed lengths $n_0$ need additional analytic or algebraic
structure (such as symmetry-sector trace recursion, polymer expansions, or variational bounds; see
[`../AGENTS.md`](../AGENTS.md)).

## What is rigorous

- **Recurrence and rational arithmetic:** `certified_gap_bound` (interval arithmetic) and
  `rational_gap_chain` (exact rational arithmetic) round outward at every step, guaranteeing that the
  returned gap bound is a strict lower bound for the given seed defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$.
  `ExactArithmetic` produces verified rational enclosures of logarithms and fractional powers.
- **Verified ring spectra (`VerifiedSpectra` and `GapCertificates`):** The reflection-positivity
  decomposition is verified in exact rational arithmetic. In each charge and momentum block, candidate
  eigenvectors are computed in `Float64`, and the residual and orthogonality defects are enclosed using
  explicit a priori IEEE binary64 rounding bounds (assuming round-to-nearest and the standard $\gamma_n$
  BLAS dot-product bound, valid for any summation order and fused multiply-add). Weyl and Ostrowski
  perturbation theorems then enclose every eigenvalue, and $\log Z$ is summed in interval arithmetic.
- **Screening tools (non-rigorous):** Eigenvalues from `ExactDiagonalization.spectrum` are unverified
  floating-point approximations (and when a sector uses Lanczos, degenerate eigenvalues can be missed).
  Similarly, `Purification.purified_log_partition` is an uncertified tensor-network diagnostic.

When running multiple Julia processes concurrently on a shared machine, set `OPENBLAS_NUM_THREADS=1`.
