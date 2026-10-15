# Spin-2 Heisenberg chain

`Spin2Gap` investigates the periodic spin-2 Heisenberg ring with the thermal-MPO and
variational tools of the spin-1 certificate. **No uniform spin-2 gap is certified.** The
reflection-positivity input and the scalar bootstrap carry over unchanged, and exact finite
variational energies are certified up to 1024 sites; the missing ingredient is a rigorous
thermal upper bound (one certified trace or two purities) at a useful seed scale, which no
route here has reached. [`METHOD.md`](METHOD.md) is the dated lab notebook of every
experiment; [`IDEAS.md`](IDEAS.md) is the route memo and its audit.

| parameter | seed (n₀, β₀) | certified gap ≥ | lengths covered |
|---|---|---|---|
| J = 1, S = 2 | none; conditional target n₀ = 256, β₀ = 64 needs a certified \(\hat Z_{256}(64)\le1.08\) | no certificate (0.003019921 if that trace bound were proved) | — |

## Model

\[
H_L=\sum_{j=1}^L \mathbf S_j\cdot\mathbf S_{j+1},\qquad
\mathbf S^2=6I,\quad J=1,\quad \mathbf S_{L+1}=\mathbf S_1,
\]

one spin-2 site per cell (lengths count sites), periodic, no on-site term. The thermal and
variational codes shift energies by a rational \(a\) per bond (e.g. \(a=238061/50000\)); the
shift cancels in both purities. Numerical references: gap ≈ 0.0892 and correlation length
≈ 49.5 (Todo–Kato), so sharp seeds need rings of order 10³ sites.

## Reflection positivity

Real basis of the spin-2 representation on real symmetric traceless 3×3 matrices:
\(S^\alpha=-iG_\alpha\) with real antisymmetric \(G_\alpha\), so
\(-\mathbf S_j\cdot\mathbf S_{j+1}=\sum_\alpha G_\alpha\otimes G_\alpha
=\sum_\alpha G_\alpha\otimes\overline{G_\alpha}\) with \(\|G_\alpha\|=2\). `real_generators()`
(= `CertifiedGaps.Spins.cartesian_generators(2)`) builds the \(G_\alpha\) independently of the
magnetic basis; `test/runtests.jl` checks the identity and the generator relations in
`Float64` (an exact-rational RP certificate is not stored here). The transfer-operator
argument behind the criterion is spelled out in METHOD.md, "Why the exact transfer argument
extends".

## Result

- **No gap certificate.** No admissible seed has been certified or numerically located at
  accessible sizes. DMRG at L = 4, 8, 16 and thermal MPOs at n = 4 (all floating point) show
  that small seeds are hopeless: the necessary β-window screen is empty at n = 4
  (β ≥ 8.42 for thermal purity, β ≤ 0.0136 for spatial purity).
- **Certified finite variational energies** (exact integer or Arb arithmetic, tracked
  certificates): SU(2)-structured periodic MPS trials give trial energy per bond < −4.7595 at
  D = 24 (N = 128, 256, 512; `results/trial-certificate.json`), < −4.76122 at D = 137
  (`results/trial-137-certificate.json`) and < −4.7612463 at D = 295 for N up to 1024
  (`results/trial-295-symmetric-certificate.json`). These are upper bounds on E₀(N)/N only.
- **Conditional one-trace target.** With the D = 295 trial, a certified shifted trace
  \(\hat Z_{256}(64)\le1.08\) would give Δ_L > 0.003019921 for every even L ≥ 512
  (`scripts/one_trace_target.jl`). The analogous n = 128, β = 32 target is ruled out by the
  D = 137 trial. The thermal inequality itself is the open step.
- **Thermal upper bounds (open).** Rigorous building blocks exist and are certified on small
  references: Araki–Lieb–Thirring Trotter upper bounds, local MPO norm and squaring
  certificates, polynomial majorants with a certified starting point, operator-error and
  adjoint trace corrections, and charge-preserving generation (METHOD.md). None yet produces a
  certified bound at n = 256, β = 64: the furthest 256-site trajectory has two independently
  replayed squarings, at β = 2⁻²⁰, with 26 doublings still missing. The October 14 audit
  (IDEAS.md) recommends local gates with final-trace certification next.

## Reproduce

From the repository root, Julia 1.12, `export OPENBLAS_NUM_THREADS=1`; the environment pins
ITensors 0.9.32, ITensorMPS 0.4.2 and Nemo 0.57.0. Times last measured on a heavily
oversubscribed machine (load ≈ 100 on 10 cores), so upper bounds: `test/runtests.jl` 54 min
(dominated by the thermal and DMRG controls), the standalone test files 1.5–18 min each.

```sh
julia --startup-file=no --project=models/spin2 -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/spin2 models/spin2/test/runtests.jl  # core tests
julia --startup-file=no --project=models/spin2 models/spin2/check.jl          # tracked certificates
julia --startup-file=no --project=models/spin2 models/spin2/check.jl --full   # + exact replays
# Standalone test files (each with --project=models/spin2): test/trials.jl, variational.jl,
# certificate.jl, thermal_bounds.jl, polynomial_mpo.jl, l2_bounds.jl, mpo_adjoint_flow.jl,
# one_trace.jl; the CertifiedMPO tests run with the core: julia --project=julia -e 'using Pkg; Pkg.test()'
```

Commands that work on a fresh clone (they read only tracked files or write new ones):

```sh
# DMRG and thermal-MPO pilots (write results/*.json; floating point):
julia --startup-file=no --project=models/spin2 models/spin2/scripts/run.jl \
  models/spin2/scripts/configs/dmrg-pilot.toml models/spin2/results/dmrg-pilot.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/run.jl \
  models/spin2/scripts/configs/thermal-pilot.toml models/spin2/results/thermal-pilot.json
# Exact replay of the D = 24 and D = 137 integer witnesses (as check.jl --full does):
julia --startup-file=no --project=models/spin2 models/spin2/scripts/certify_trial.jl \
  models/spin2/results/trial-replay.json models/spin2/results/trial-certificate.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/trials/certify_reduced_trial.jl \
  models/spin2/results/trial-137-replay.json models/spin2/results/trial-137-certificate.json 137 replay
# Thermal starting-point audit of IDEAS.md (exact/Arb, ≈0.5 s after startup):
julia --startup-file=no --project=models/spin2 models/spin2/scripts/polynomial/audit_thermal_start.jl
# Model-independent scalar targets (the D = 137 witness is tracked):
julia --startup-file=no --project=models/spin2 models/spin2/scripts/one_trace_target.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/polynomial/polynomial_thermal_targets.jl
```

Everything else under `results/` is gitignored and produced by the drivers themselves, so on a
fresh clone each verifier needs its producer to run first. For example,
`scripts/screen.jl` and `scripts/bootstrap_screen.jl` read the DMRG outputs of `run.jl`
(`dmrg-pilot.json`, `dmrg-refinement.json`, then `seed-screen.json`);
`scripts/trials/{search,refine,polish,analytic,large,preconditioned}_trial.jl` form a chain,
each reading the previous output; and the `verify_*.jl` scripts read witnesses written by the
matching probe (`mpo_norm/certify_mpo_norm.jl` → `mpo-local-norm-certificate-n*-witness.json`,
`squaring/mpo_squaring_probe.jl` and `refine_mpo_squaring.jl` → `mpo-squaring-*-witness.json`,
`polynomial/exact_polynomial_seed.jl` → `polynomial-seed-exact-n*-witness.json`,
`polynomial/polynomial_flow_probe.jl` → `polynomial-*-flow-*.json` and its step witnesses,
`adjoint/trace_correction_probe.jl` and `adjoint/adjoint_trace_probe.jl` → their
`*-witness.json`). METHOD.md gives each experiment's exact invocation and parameters; the long
runs take minutes to hours.

## What is rigorous

- **Exact or interval (trusted base: integer/rational arithmetic, Nemo/FLINT Arb balls with
  outward rounding):** the trial-energy certificates (integer contractions, or reduced SU(2)
  blocks in Arb with exact covariance and Schur-lemma checks), the scalar bootstrap and
  one-trace targets (exact rationals), the Araki–Lieb–Thirring Trotter upper bounds on 4- and
  6-site references, and every saved MPO norm, squaring, operator-error and trace-correction
  witness (rational cores, Arb Cholesky or interval residuals; acceptance never uses a
  floating-point eigensolver).
- **Floating point only:** DMRG energies and gaps, thermal-MPO purities, the β-window screens,
  all optimisation and gauge searches, and the RP identity check in `test/runtests.jl`.
- **Not established:** any thermal upper bound at a seed scale, any seed, any gap. Every
  stored result says `uniform_gap_certified = false` or `certified = false`.
- Several results record SHA-256 hashes of the sources that produced them; format only the
  files you change (`.JuliaFormatter.toml`: margin 100, `for … in`).

## Layout

| path | content |
|---|---|
| `src/Spin2Gap.jl`, `model.jl`, `exact.jl`, `dmrg.jl`, `thermal.jl` | module `Spin2Gap`: S = 2 site type, bonds, Hamiltonian MPO, exact sector spectra, DMRG, thermal curves, provenance |
| `src/Trials.jl`, `src/trial_certificate.jl` | `Spin2Gap.Trials`: SU(2) trial templates and exact integer energy certificates |
| `src/Variational.jl` | `Spin2Gap.Variational`: reduced variational energies, analytic gradients, optimisers |
| `src/Certificate.jl` | `Spin2Gap.Certificate`: reduced exact/Arb trial certificates |
| `src/ThermalBounds.jl` | `Spin2Gap.ThermalBounds`: rational Trotter-layer thermal upper bounds |
| `src/PolynomialMPO.jl` | `Spin2Gap.PolynomialMPO`: spin-2 polynomial and quartic SOS seeds |
| `check.jl` | checker of the tracked trial certificates (`--full`: exact replays) |
| `scripts/run.jl`, `scripts/configs/` | DMRG/thermal-MPO runner and its four TOML configurations |
| `scripts/screen.jl`, `bootstrap_screen.jl`, `one_trace_target.jl` | β-window screen, bootstrap targets, one-trace targets (cited by `julia/src/recurrence/rational.jl` and `models/spin1/experiments/optimized-spin1-bounds.tex`, so kept at this path) |
| `scripts/certify_trial.jl` | exact D = 24 trial certificate; kept at this path because `results/trial-certificate.json` records its SHA-256 (`runner_sha256`) and it resolves `../src`, `../results` from here |
| `scripts/trials/` | trial search, refinement, polishing, analytic and large searches, `certify_reduced_trial.jl` |
| `scripts/mpo_norm/` | MPO norm probe, gauge probe, `certify_mpo_norm.jl`, `verify_mpo_norm.jl` |
| `scripts/squaring/` | `mpo_squaring_probe.jl`, `refine_mpo_squaring.jl`, `verify_mpo_squaring.jl` |
| `scripts/compression/` | four-site references: positive, repeated, Hilbert–Schmidt (`repeated_l2_probe.jl`) and Loewner compression, dense diagnostic |
| `scripts/polynomial/` | polynomial seeds, flows and their verifiers, quartic seeds, gauge screens, thermal targets, degree screen, `audit_thermal_start.jl` |
| `scripts/operator_error/` | operator-error probe and refinement, operator seeds and their verifiers, power screen |
| `scripts/adjoint/` | adjoint flow, adjoint trace probe, Schatten backward refinement, trace correction |
| `scripts/benchmarks/` | Arb product and charge-block benchmarks |
| `test/runtests.jl` and standalone `test/*.jl` | tests (see Reproduce) |
| `results/trial-certificate.json`, `trial-137-certificate.json`, `trial-295-symmetric-certificate.json` | the tracked certificates; everything else in `results/` is gitignored |
| `.JuliaFormatter.toml` | formatter settings (editor tooling; kept, not data) |

Exploratory drivers that no document or test invokes by name, kept because METHOD.md reports
their outputs: `scripts/trials/analytic_trial.jl` (`trial-analytic.json`),
`scripts/trials/preconditioned_trial.jl` (`trial-preconditioned.json`),
`scripts/mpo_norm/certify_mpo_norm.jl` (`mpo-local-norm-certificate.json`),
`scripts/squaring/mpo_squaring_probe.jl` (`mpo-squaring-probe.json` and the gauged, capped,
blocked and small-step variants), `scripts/compression/repeated_compression_probe.jl`
(`repeated-positive-compression.json`), `scripts/compression/dense_compression_diagnostic.jl`
(`dense-compression-diagnostic-n4.json`), `scripts/polynomial/exact_polynomial_seed.jl`
(`polynomial-seed-exact-n*.json`), `scripts/polynomial/certify_schmidt_seed_norm.jl`
(`schmidt-seed-norm-n256.json`) and `scripts/operator_error/refine_operator_error.jl`
(refined operator-error records). None duplicates another script.

Hash pinning: `test/trials.jl` compares the `source_sha256` of `results/trial-certificate.json`
with `src/model.jl`, `src/Trials.jl` and `src/trial_certificate.jl`, so those three files must
not change. Hash fields recorded by the probes in gitignored results (`runner_sha256`,
`source_sha256`, module hashes) are provenance and are not compared with the checkout.

## History

- October 7–8, 2026: DMRG and thermal pilots, SU(2) trial search and the first exact trial
  certificates (D = 24, then D = 137 and D = 295).
- October 8–11: thermal upper-bound experiments (positive, Hilbert–Schmidt and spectral
  compression; MPO norms and squaring; polynomial seeds and flows; operator errors; adjoint and
  Schatten trace corrections; charge blocks). The certified MPO toolkit then moved verbatim
  into the core as `CertifiedGaps.CertifiedMPO` (with its tests, `julia/test/certified_mpo/`),
  and the shared spin, ED, purification and recurrence code into the core as well; hash fields
  that scripts record for modules point at the new core files since then.
- October 11 and 14: route memo and the reproducible thermal-starting-point audit
  (`IDEAS.md`, PR #26).
- Layout tidy: the 1400-line lab notebook moved from this README to `METHOD.md` (text
  unchanged apart from script paths); `configs/` moved to `scripts/configs/`; the flat
  `scripts/` was grouped by topic (subdirectories above), with `certify_trial.jl` and the
  externally cited drivers left at the top; `check.jl` and a model `.gitignore` were added.
