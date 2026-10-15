# Spin-1 Heisenberg ring: reproduction of the periodic certificate

`Spin1Certificates` recomputes, in exact integer and rational arithmetic, the finite inputs of
*The periodic spin-one Haldane gap* (OpenAI, September 24, 2026) and checks every integer
against the paper's own evidence files. Using those frozen finite inputs, the scalar and
temperature-lifting analyses in `experiments/` establish significantly stronger certified lower
bounds on the uniform spectral gap without running any new trace or trial-state contractions.
Independent `ITensors` calculations provide numerical physical-MPO, MPS, DMRG, and thermal-MPO
cross-checks (not part of the rigorous certificate).

| Certificate variant | Seed purities / defects | Certified gap $\Delta_L \ge$ | Lengths covered |
| --- | --- | --- | --- |
| Source 1 (published) | $P_{\mathrm{s}}(2304, 784) \ge 1 - 84513/10^7$, $P_\tau(4608, 784) \ge 1 - 27493/(5\cdot 10^6)$; bootstrap $u_0 = 1/100$, $C = 5$ | $\log(20)/784 \approx 0.00382$ | even $L \ge 2304$ |
| Source 2 (published) | $P_{\mathrm{s}}(60, 105/4) > 7/8$, $P_\tau(120, 105/4) > 7/8$; bootstrap $u_0 = 1/8$, $C = 79/10$ | $(4/105)\log(80/79) \approx 0.00048$ | even $L \ge 60$ |
| Scalar optimization (`experiments/`) | Reuses source-1 and source-2 finite inputs with unequal-defect recurrence | $0.016978269$ ($> 1/60$) | even $L \ge 36$ |
| Temperature lifting (`experiments/`) | Reuses source-1 traces ($L \le 12$) and $L=120$ MPS trial; lifted seed at $n_0 = 36, 60$ | $0.072901771$ ($> 7/100$ simple) | even $L \ge 36$ |

## Model

We consider the periodic spin-1 antiferromagnetic Heisenberg ring of length $L$ with unit
coupling:
$$H_L = \sum_{j=1}^{L} \mathbf{S}_j \cdot \mathbf{S}_{j+1}, \qquad \mathbf{S}_j^2 = 2I, \qquad \mathbf{S}_{L+1} \equiv \mathbf{S}_1.$$
Lengths $L$ count sites. In the thermal trace recursion, each bond is shifted by a rational
per-site constant $\kappa = 700741/500000$ (denoted $a$ in the upstream source code), so that
the twisted partition function is evaluated as
$$Z_n(\beta, g) = \operatorname{Tr}\!\left(g\,\exp[-\beta(H_n + n\kappa I)]\right)$$
for a global spin rotation $g \in \mathrm{SO}(3)$ acting on the closing bond. Because the shift
$n\kappa I$ is extensive, it cancels identically in both the spatial purity
$P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$ and the thermal purity
$P_\tau(2n, \beta) = Z_{2n}(2\beta)/Z_{2n}(\beta)^2$ (denoted $S(n, \beta)$ and $T(2n, \beta)$ in
the upstream paper and JSON records).

## Reflection positivity

In the Cartesian spin-1 basis $\{|x\rangle, |y\rangle, |z\rangle\}$, the spin operators
$(S^\alpha)_{\beta\gamma} = -i\epsilon_{\alpha\beta\gamma}$ are purely imaginary and
antisymmetric. Setting $C_\alpha = G_\alpha = i S^\alpha$ (which are real matrices) gives the
exact reflection-positive decomposition on every bond:
$$-h = -\mathbf{S}_i \cdot \mathbf{S}_j = \sum_{\alpha \in \{x,y,z\}} C_\alpha \otimes \overline{C_\alpha}.$$
For this package, the reflection-positivity identity is cited from the upstream paper and is
checked algebraically in the shared core test suite (`julia/test/runtests.jl`) and in the Lean 4
formalization (`lean/CertifiedGaps/Spin1.lean`).

## Seeds and certified bounds

Rather than diagonalizing $H_n$ at $n_0 = 60$ or $n_0 = 2304$ directly, the certificate bounds
the seed purities indirectly by combining two exact finite computations:
1. **Short-ring twisted thermal traces** at lengths $n \le 12$ and inverse temperatures
   $\beta = 49/4$ (source 1) or $\beta = 21/2$ (source 2), evaluated with three boundary twists $g$
   at the closing bond: the identity $I$, a $\pi$-rotation $P = e^{i\pi S^y}$, and a
   $2\pi/3$-rotation $C = e^{i(2\pi/3)S^z}$. Linear combinations of these three twisted traces
   project onto the odd-spin ($O$) and non-singlet ($N$) $\mathrm{SU}(2)$ sectors, bounding excited
   state weights at short lengths.
2. **Periodic MPS trial-state contractions** at lengths $N = 72$ and $N = 120$ using an exact
   bond-dimension-26 integer matrix product state, which provides rigorous variational upper
   bounds on the ground-state energy density ($\langle H_{72}\rangle/72 < -1.401482$ and
   $\langle H_{120}\rangle/120 < -1.401482$).

Polynomial moment filters and spatial doubling inequalities then propagate these short-ring trace
bounds and variational energy bounds to certify the seed defects $\varepsilon_{\mathrm{s}}$ and
$\varepsilon_\tau$ at $n_0 = 60$ and $n_0 = 2304$.

### Component computations

| Component | Fresh computation | Compared with upstream |
| --- | --- | --- |
| `src/thermal.jl`, source 1 | 18 traces: lengths $4\text{--}12$, identity and $\pi$ twists, $\beta = 49/4$ | every full integer total and rational interval endpoint |
| `src/thermal.jl`, source 2 | 9 traces: lengths $6, 8, 10$, all three twists ($I, P, C$), $\beta = 21/2$ | every full integer total, error bound, and interval endpoint |
| `src/trials.jl`, `src/trials_su2.jl` | all 2,028 entries of three $26 \times 26$ integer matrices; periodic MPS contractions at lengths 72 and 120 | full energy numerators, norms, positive variational margins, all length-120 transfer blocks |
| `src/arithmetic.jl` | independent column-action Hamiltonians, all magnetization blocks, twisted symmetries, 27 weighted triad rows | exact equality and integer/rational error bounds |
| `src/scalars.jl` | moment filters, capped masses, both purity initializations, six rounded updates, bootstrap constants | rational inequalities and every published rounded pair |
| `src/itensors_checks.jl` | four-site MPOs for all twists, integer trial imported as an MPS, QN DMRG, thermal MPO | numerical cross-check only |

The integer thermal recurrence uses arithmetic right shifts, which implement exact floor division
for both positive and negative entries. Explicit bounds rule out overflow in `Int64` partial sums
and `Int128` norm chunks; final totals and trial contractions use `BigInt`. Source 2 is an
independent unpacked implementation in the Eisenstein basis $\{1, \omega\}$
($\omega^2 + \omega + 1 = 0$), rather than a translation of the upstream packed-lane arithmetic.
JSON integers are parsed directly from their string tokens so large certificate integers never
pass through `Float64`. Sectors across all lengths and twists run on the available threads in a
largest-first queue, with 32 integer columns per batch.

The default trial verifier works in exact $\mathrm{SU}(2)$ multiplicity blocks of sizes
$26, 56, 48, 24, 7,$ and $1$ instead of the eleven $S^z$ charge blocks (largest dimension 162);
[`METHOD.md`](METHOD.md) derives that reduction, describes the thermal kernel optimizations, and
records the runtime benchmarks.

### Stronger bounds from the same finite inputs (`experiments/`)

Both experiments in `experiments/` reuse the stored, hash-checked finite inputs from `results/`
without running any new thermal or trial contractions. Every inequality is decided with exact
rational logarithm and fractional-power enclosures; floating-point numbers are used only for
human-readable display. Full derivations are given in
[`experiments/optimized-spin1-bounds.tex`](experiments/optimized-spin1-bounds.tex) (with the
compiled PDF alongside). Its later sections also record spin-2 variational certificates and a
conditional one-trace criterion on which [`models/spin2`](../spin2) builds.

| Experiment | Certified gap $\Delta_L >$ | Lengths covered | Record |
| --- | --- | --- | --- |
| Scalar optimization of source 1 (`optimize_bounds.jl`) | $0.016978269$ ($> 1/60$ with six-decimal cutoffs) | even $L \ge 36$ | `bounds.json` |
| Scalar optimization of source 2 (`optimize_bounds.jl`) | $0.011977901$ | even $L \ge 60$ | `bounds.json` |
| Temperature lifting, simple (`temperature_lift.jl`) | $0.070450002 > 7/100$ | even $L \ge 36$ | `temperature-lift.json` |
| Temperature lifting, optimized (`temperature_lift.jl`) | $0.072901771$ | even $L \ge 36$ | `temperature-lift.json` |

The scalar optimization tightens the polynomial filter cutoffs, tracks unequal spatial and
thermal defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ through the coupled recurrence
rather than collapsing to $\max(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$, and covers intermediate
chain lengths via the reflection-positivity interpolation inequality before invoking the uniform
bootstrap tail. Temperature lifting combines the 18 source-1 thermal traces with the length-120
MPS trial, raises the inverse temperature before extracting purity bounds, and retains the
threefold $\mathrm{SU}(2)$ triplet multiplicity when bounding residual moments; a length-72 seed
bridges even chain lengths $36 \le L \le 72$. The script `search_temperature_lift.jl` performs the
exploratory floating-point parameter scan that selected the moment length and temperature
multiplier (`temperature-lift-search.json`, labeled non-certifying), and `temperature_lift.jl`
replays the chosen parameters in exact rational arithmetic.

```sh
julia --startup-file=no --project=models/spin1 models/spin1/experiments/test_bounds.jl            # 549 tests
julia --startup-file=no --project=models/spin1 models/spin1/experiments/test_temperature_lift.jl  # 285 tests
julia --startup-file=no --project=models/spin1 models/spin1/experiments/optimize_bounds.jl        # regenerate
julia --startup-file=no --project=models/spin1 models/spin1/experiments/temperature_lift.jl       # regenerate
```

The generators verify that the package's current source fingerprint matches the hash recorded in
`results/verification.json` and that every component output matches its recorded SHA-256 hash.
The test suites verify the stored script and input hashes and confirm that forged certificates
are rejected.

## Reproduce

Run times of the stored reproduction (4 threads): thermal source 1 takes 2210 s, thermal source 2
takes 35 s, MPS trials take 9.6 s, `ITensors` controls take 278 s, the stored-output audit takes
about 20 s, and each experiment script runs in under 2 s.

From the repository root, using Julia 1.12, first fetch the upstream evidence files (which are
downloaded into the gitignored `references/` directory):

```sh
julia tools/fetch_upstream.jl
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/spin1 models/spin1/test/runtests.jl   # 43 tests, ~1.5 min
julia --startup-file=no --project=models/spin1 models/spin1/check.jl            # stored-output audit
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/check.jl --full  # ~45 min
```

`check.jl` (fast) audits a temporary copy of `results/`, requires the component hashes and the
implementation fingerprint to equal those in the committed `verification.json`, and runs the
corrupted-artifact tests on the copy and the two `experiments/` test files (16 checks; about 13 min
on a heavily loaded machine, a few minutes otherwise). `check.jl --full` runs the complete fresh
reproduction below into a temporary directory and requires every certified integer, interval, and
scalar record to equal the committed one. Neither writes anything into the repository.

The same steps can be run by hand. To verify the stored outputs without repeating the expensive
contractions, run the audit on a **copy** of `results/`. Because the audit rewrites
`verification.json` in the target directory and the scripts in `experiments/` pin the SHA-256 hash
of the committed `results/verification.json`, auditing `results/` in place would invalidate the
experiment hash checks:

```sh
cp -R models/spin1/results /tmp/spin1-audit
julia --startup-file=no --project=models/spin1 models/spin1/scripts/reproduce.jl audit /tmp/spin1-audit
julia --startup-file=no --project=models/spin1 models/spin1/test/audit_failures.jl /tmp/spin1-audit
```

The audit reports `execution_mode: stored-output-audit`. It verifies all full contraction
integers against the upstream evidence, rebuilds their rational intervals, checks sector sums,
and recomputes the arithmetic and scalar assertions without claiming a fresh thermal or trial
contraction. A failed audit overwrites any earlier success record with `status: failed`, and
`audit_failures.jl` (8 tests) verifies that truncated or corrupted artifacts are rejected.

A complete fresh reproduction (about 45 minutes on 4 threads) should likewise be directed to a
separate output directory rather than overwriting the committed `results/`:

```sh
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/scripts/reproduce.jl all /tmp/spin1-fresh
```

This command recomputes every contraction, runs the `ITensors` controls, and writes
`verification.json` as passed only after every component succeeds. Individual components can also
be run separately with `scripts/thermal.jl 1`, `scripts/thermal.jl 2`, `scripts/trials.jl`,
`scripts/scalars.jl`, `scripts/arithmetic.jl`, and `scripts/itensors_crosscheck.jl`; the thermal
scripts take an optional output filename after the source number, while the others take it as
their first argument.

## What is rigorous

- **Rigorous in this package:** every finite inequality reported above is decided in exact integer
  (`Int64`/`Int128`/`BigInt`/FLINT `ZZMatrix`) or exact rational (`Rational{BigInt}`) arithmetic:
  the fixed-point thermal trace recursions with explicit a priori rounding bounds, the integer MPS
  trial contractions, the polynomial moment filters, the purity initializations, and the scalar
  bootstrap updates.
- **Outside this package's scope:** the reflection-positivity identity and the all-scale operator
  induction are formalized separately in [`lean/`](../../lean) and derived in [`notes/`](../../notes).
  The temperature-lifting derivation in `experiments/optimized-spin1-bounds.tex` is proved in
  that note and checked in rational arithmetic by `temperature_lift.jl`, while the `ITensors`
  cross-checks are numerical controls explicitly labeled `certified: false`.

## Layout

The upstream source and evidence files are pinned to commit
`adc7f1241b42e322a6451854ab7e4b4c146bf78a` and read from `references/periodic-gap-source` at the
repository root. `tools/fetch_upstream.jl` downloads the required files and verifies the pinned
manifest hash and every individual file hash; each component re-checks 43 upstream source and
evidence files before running. These are modified Julia implementations, not upstream code; see
[`NOTICE`](NOTICE) and the upstream Apache-2.0 license in [`LICENSE.upstream`](LICENSE.upstream).

- `src/`: the `Spin1Certificates` module (thermal recursions and kernels, MPS trials and their
  $\mathrm{SU}(2)$ reduction, scalar filters and recurrences, independent arithmetic checks,
  stored-output audit, and `ITensors` controls). Uses [`CertifiedGaps`](../../julia) for exact
  rational arithmetic, result I/O, the rational doubling chain with the bootstrap tail, and
  `ITensors` purification.
- `scripts/`: component runners, the reproduction and audit entry point `reproduce.jl`, and three
  benchmark scripts.
- `test/`: fast unit tests (`runtests.jl`) and corrupted-artifact rejection tests
  (`audit_failures.jl`).
- `results/` (frozen): the optimized run's component outputs, benchmarks, and
  `verification.json`, which records component hashes and the implementation fingerprint (SHA-256
  of `Project.toml`, `Manifest.toml`, `.JuliaFormatter.toml`, and every `.jl` file under `src/`,
  `scripts/`, and `test/`). `results/reference/` preserves the initial unoptimized Julia run.
  Because of the fingerprint, modifying any of those files will cause the experiment generators to
  refuse to run until a full reproduction updates `verification.json`.
- `experiments/`: the stronger bounds above (scalar optimization and temperature-lifting scripts,
  JSON records, test suites, and LaTeX derivations).
- `check.jl`: the checker (fast stored-output audit, `--full` fresh reproduction); it lies outside
  the fingerprint, so editing it does not invalidate `verification.json`.
- `IDEAS.md` (route memo, reconstructed) and `METHOD.md` (the $\mathrm{SU}(2)$ reduction, thermal
  kernel optimizations, and benchmarks; formerly `OPTIMIZATION.md`).

Layout exceptions (top-level entries beyond the standard format), each kept on purpose:

- `.JuliaFormatter.toml`: hash-pinned. It is one of the three environment files in the
  implementation fingerprint recorded in `results/verification.json`; moving or editing it would
  invalidate the stored audit and the experiments' generators.
- `experiments/`: the stronger bounds, with records that pin the SHA-256 hashes of their own
  scripts and of `results/verification.json` (`experiments/test_*.jl` check them); moving them
  under `scripts/` would put them inside the fingerprint.
- `NOTICE` and `LICENSE.upstream`: attribution and the Apache-2.0 license of the upstream code
  these Julia files are adapted from, which the license requires to accompany the modified files.

Rationals are stored as `"numerator//denominator"` strings. `Manifest.toml` pins `ITensors`
0.9.32, `ITensorMPS` 0.4.2, `JSON3` 1.14.3, and `Nemo` 0.57.0.

## History

The initial unoptimized Julia port is preserved in `results/reference/`. [`METHOD.md`](METHOD.md)
documents the subsequent $\mathrm{SU}(2)$ multiplicity-block reduction and thermal kernel
optimizations along with their benchmarks. The improved gap bounds in `experiments/` were derived
afterward from the frozen finite inputs.

The format tidy renamed `OPTIMIZATION.md` to `METHOD.md` (no fingerprinted file reads it), added
`IDEAS.md` and the top-level `check.jl`, and documented the layout exceptions; no fingerprinted
file, result, or experiment record changed.
