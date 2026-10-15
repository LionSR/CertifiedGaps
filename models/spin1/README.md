# Spin-1 Heisenberg ring: reproduction of the periodic certificate

`Spin1Certificates` recomputes, in exact integer and rational arithmetic, the finite inputs of
**The periodic spin-one Haldane gap** (OpenAI, September 24, 2026) and checks every integer
against the paper's own evidence files. ITensors supplies independent physical-MPO, MPS, DMRG
and thermal-MPO controls, which are numerical and never certify anything. The experiments in
`experiments/` reuse the same finite inputs to prove stronger constants.

| | |
| --- | --- |
| Source 1 seed | \(S(2304,784)\ge1-84513/10^7\), \(T(4608,784)\ge1-27493/(5\cdot10^6)\); bootstrap \(u_0=1/100\), \(C=5\) |
| Source 2 seed | \(S(60,105/4)>7/8\), \(T(120,105/4)>7/8\); bootstrap \(u_0=1/8\), \(C=79/10\) |
| Published bounds | \(\Delta_L\ge\log(20)/784\approx0.00382\) for even \(L\ge2304\); \((4/105)\log(80/79)\approx0.00048\) for even \(L\ge60\) |
| Improved bounds (`experiments/`) | \(\Delta_L>0.072901771\) (and \(>7/100\) by a simpler choice) for every even \(L\ge36\) |

## Model

\(H_L=\sum_{j=1}^{L}\mathbf S_j\cdot\mathbf S_{j+1}\), spin 1 (\(\mathbf S^2=2I\)), periodic,
coupling 1.

Shift: \(Z_n(\beta,g)=\operatorname{Tr}\exp[-\beta(H_n^g+naI)]\)
with \(a=700741/500000\); it cancels in both purities.

Here \(S(n,\beta)=Z_{2n}(\beta)/Z_n(\beta)^2\) and \(T(2n,\beta)=Z_{2n}(2\beta)/Z_{2n}(\beta)^2\).

## Reflection positivity

Cartesian basis, where every \(S^a\) is purely imaginary, so
\(-\mathbf S_i\cdot\mathbf S_j=\sum_a S^a\otimes\overline{S^a}\).
Taken from the paper; not re-checked by this package.

## Seeds and certified bounds

The seeds and bounds are in the table at the top. The rest of this section says how the finite
inputs are recomputed and how `experiments/` strengthens the bounds.

### What is computed

The twists at the closing bond are the identity, a rotation by \(\pi\), and a rotation by
\(2\pi/3\).

| Component | Fresh computation | Compared with upstream |
| --- | --- | --- |
| `src/thermal.jl`, source 1 | 18 traces: lengths 4–12, identity and \(\pi\) twists, \(\beta=49/4\) | every full integer total and rational interval endpoint |
| `src/thermal.jl`, source 2 | 9 traces: lengths 6, 8, 10, all three twists, \(\beta=21/2\) | every full integer total, error bound and interval endpoint |
| `src/trials.jl`, `src/trials_su2.jl` | all 2,028 entries of three 26×26 integer matrices; periodic MPS contractions at lengths 72 and 120 | full energy numerators, norms, positive variational margins, all length-120 transfer blocks |
| `src/arithmetic.jl` | independent column-action Hamiltonians, all magnetisation blocks, twisted symmetries, 27 weighted triad rows | exact equality and integer/rational error bounds |
| `src/scalars.jl` | moment filters, capped masses, both purity initialisations, six rounded updates, bootstrap constants | rational inequalities and every published rounded pair |
| `src/itensors_checks.jl` | four-site MPOs for all twists, integer trial imported as an MPS, QN DMRG, thermal MPO | numerical control only |

The integer thermal recurrence uses arithmetic right shifts, which implement floor division
also for negative entries. Explicit bounds rule out overflow in `Int64` partial sums and
`Int128` norm chunks; final totals and trial contractions use `BigInt`. Source 2 is an
independent unpacked implementation in the Eisenstein basis \(1,\omega\)
(\(\omega^2+\omega+1=0\)), not a translation of the upstream packed-lane arithmetic. JSON
integers are parsed from their tokens, so large certificate integers never pass through
`Float64`. Sectors across all lengths and twists run on the available threads in a
largest-first queue, with 32 integer columns per batch.

The exact trials establish \(\langle H_{72}\rangle/72<-1.401482\) and
\(\langle H_{120}\rangle/120<-1.401482\). The default trial verifier works in exact SU(2)
multiplicity blocks of sizes 26, 56, 48, 24, 7 and 1 instead of the eleven charge blocks (largest
162); [`METHOD.md`](METHOD.md) derives that reduction, describes the thermal kernel
changes, and records the benchmarks.

### Stronger bounds from the same finite inputs (`experiments/`)

Both experiments reuse the stored, hash-checked finite inputs and run no new thermal or trial
contraction. They decide every inequality with rational logarithm and fractional-power
enclosures; floating-point numbers are display aids only. The derivations are in
[`experiments/optimized-spin1-bounds.tex`](experiments/optimized-spin1-bounds.tex) (compiled PDF
alongside). Its later sections also record spin-2 variational certificates and a conditional
one-trace criterion on which `models/spin2` builds.

| Experiment | Bound | Lengths | Record |
| --- | --- | --- | --- |
| Scalar optimisation of source 1 (`optimize_bounds.jl`) | \(\Delta_L>0.016978269\) (\(>1/60\) with six-decimal cutoffs) | even \(L\ge36\) | `bounds.json` |
| Scalar optimisation of source 2 | \(\Delta_L>0.011977901\) | even \(L\ge60\) | `bounds.json` |
| Temperature lifting, simple (`temperature_lift.jl`) | \(\Delta_L>0.070450002>7/100\) | even \(L\ge36\) | `temperature-lift.json` |
| Temperature lifting, optimised | \(\Delta_L>0.072901771\) | even \(L\ge36\) | `temperature-lift.json` |

The scalar optimisation tightens the polynomial cutoffs, keeps unequal purity defects through
the coupled updates, and covers the early lengths by the interpolation inequality before
invoking the uniform bootstrap theorem. Temperature lifting combines the eighteen source-1
thermal values with the length-120 trial, raises the inverse temperature before extracting
purity bounds, and keeps the threefold sector multiplicity when bounding residual moments; the
length-72 seed bridges lengths 36–72. `search_temperature_lift.jl` is the floating-point
parameter scan that chose the moment length and multiplier (`temperature-lift-search.json`,
labelled noncertifying); the chosen parameters are replayed by the rational verifier.

```sh
julia --startup-file=no --project=models/spin1 models/spin1/experiments/test_bounds.jl            # 549 tests
julia --startup-file=no --project=models/spin1 models/spin1/experiments/test_temperature_lift.jl  # 285 tests
julia --startup-file=no --project=models/spin1 models/spin1/experiments/optimize_bounds.jl        # regenerate
julia --startup-file=no --project=models/spin1 models/spin1/experiments/temperature_lift.jl       # regenerate
```

The generators assert that the package's current source fingerprint equals the one recorded in
`results/verification.json` and that every component output still has its recorded hash. The
tests check the stored script and input hashes and reject forged certificates.

## Reproduce

Run times of the stored run (4 threads): thermal source 1 2210 s, source 2 35 s; trials 9.6 s;
ITensors controls 278 s; stored-output audit about 20 s; each experiment under 2 s.

From the repository root, with Julia 1.12. The upstream evidence is not stored in the
repository; fetch it first (it lands in the gitignored `references/` directory):

```sh
julia tools/fetch_upstream.jl
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/spin1 models/spin1/test/runtests.jl   # 43 tests, ~1.5 min
julia --startup-file=no --project=models/spin1 models/spin1/check.jl            # stored-output audit
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/check.jl --full  # ~45 min
```

`check.jl` (fast) audits a temporary copy of `results/`, requires the component hashes and
the implementation fingerprint to equal those in the committed `verification.json`, runs the
corrupted-artifact tests on the copy and the two `experiments/` test files (16 checks; about
13 min on a heavily loaded machine, a few minutes otherwise). `check.jl --full` runs the
complete fresh reproduction below into a temporary directory and requires every certified
integer, interval and scalar record to equal the committed one. Both write nothing into the
repository. The same steps by hand: to check the stored outputs without repeating the expensive
contractions, audit a **copy** of `results/`. The audit rewrites `verification.json` in the directory it checks, and the
experiments pin the SHA-256 of the committed `results/verification.json`, so an in-place audit
of `results/` would break their tests.

```sh
cp -R models/spin1/results /tmp/spin1-audit
julia --startup-file=no --project=models/spin1 models/spin1/scripts/reproduce.jl audit /tmp/spin1-audit
julia --startup-file=no --project=models/spin1 models/spin1/test/audit_failures.jl /tmp/spin1-audit
```

The audit reports `execution_mode: stored-output-audit`. It checks all full contraction
integers against the upstream evidence, rebuilds their rational intervals, checks sector sums,
and recomputes the inexpensive arithmetic and scalar assertions; it does not claim a fresh
thermal or trial contraction. A failed audit replaces any earlier success record with
`status: failed`, and `audit_failures.jl` (8 tests) checks that truncated or corrupted artifacts
are rejected.

The complete fresh reproduction (about 45 minutes on 4 threads) writes into a directory of
your choice; do not point it at the committed `results/`, which is its default:

```sh
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/scripts/reproduce.jl all /tmp/spin1-fresh
```

It recomputes every contraction, runs the ITensors controls, and writes `verification.json` as
passed only after every component succeeds. The components can also be run separately with
`scripts/thermal.jl 1`, `scripts/thermal.jl 2`, `scripts/trials.jl`, `scripts/scalars.jl`,
`scripts/arithmetic.jl` and `scripts/itensors_crosscheck.jl`; the thermal scripts take an
optional output filename after the source number, the others take it as their first argument.

## What is rigorous

**Rigorous here:** every finite inequality above is decided by integers or rationals: the thermal
trace recursions with explicit rounding bounds, the integer MPS trial contractions, the moment
filters, the purity initialisations and the scalar bootstrap updates.
**Not checked here:** the reflection-positivity input, the paper's analytic transfer
construction and its all-scale induction (cited, not formalised), and the temperature-lifting
derivation in `experiments/`, a new argument written up in a note but not independently
reviewed. The ITensors outputs are labelled `certified: false`.

## Layout

The upstream source and evidence are pinned to commit
`adc7f1241b42e322a6451854ab7e4b4c146bf78a` and read from `references/periodic-gap-source` at the
repository root. `tools/fetch_upstream.jl` downloads exactly the files used here and checks the
pinned manifest hash and every file hash; every component re-checks 43 source and evidence files
before running. These are modified Julia implementations, not upstream code; see
[`NOTICE`](NOTICE) and the upstream Apache-2.0 licence in [`LICENSE.upstream`](LICENSE.upstream).

- `src/`: the package (thermal recursions and kernels, trials and their SU(2) reduction,
  scalars, independent arithmetic checks, audit, ITensors controls). It uses the shared core
  [`CertifiedGaps`](../../julia) for exact arithmetic, result I/O, the rational doubling chain
  with the bootstrap tail, and the ITensors purification.
- `scripts/`: component runners, the reproduction/audit entry point `reproduce.jl`, and three
  benchmarks.
- `test/`: fast tests and the corrupted-artifact audit tests.
- `results/` (frozen): the optimised run's component outputs, the benchmarks, and
  `verification.json` with component hashes and the implementation fingerprint (SHA-256 of
  `Project.toml`, `Manifest.toml`, `.JuliaFormatter.toml` and every `.jl` file under `src/`,
  `scripts/` and `test/`). `results/reference/` keeps the first, unoptimised Julia run. Because of
  the fingerprint, any edit to those files makes the experiments' generators refuse to run until
  the full reproduction is redone.
- `experiments/`: the stronger bounds above, their records, tests and derivation.
- `check.jl`: the checker (fast stored-output audit, `--full` fresh reproduction); outside the
  fingerprint, so editing it does not invalidate `verification.json`.
- `IDEAS.md` (route memo, reconstructed) and `METHOD.md` (the SU(2) reduction, thermal kernel
  changes and benchmarks; formerly `OPTIMIZATION.md`).

Layout exceptions (top-level entries beyond the standard format), each kept on purpose:

- `.JuliaFormatter.toml`: hash-pinned. It is one of the three environment files in the
  implementation fingerprint recorded in `results/verification.json`; moving or editing it
  would invalidate the stored audit and the experiments' generators.
- `experiments/`: the stronger bounds, with records that pin the SHA-256 of their own scripts
  and of `results/verification.json` (`experiments/test_*.jl` check them); moving them under
  `scripts/` would put them inside the fingerprint.
- `NOTICE` and `LICENSE.upstream`: attribution and the Apache-2.0 licence of the upstream
  code these Julia files are adapted from, which the licence requires to accompany the
  modified files.

Rationals are stored as `numerator//denominator` strings. The manifest pins ITensors 0.9.32,
ITensorMPS 0.4.2, JSON3 1.14.3 and Nemo 0.57.0.

## History

The first, unoptimised Julia run is kept in `results/reference/`;
[`METHOD.md`](METHOD.md) records the later kernel and SU(2) optimisations and their
benchmarks. The stronger bounds in `experiments/` were derived afterwards from the frozen
finite inputs.

The format tidy renamed `OPTIMIZATION.md` to `METHOD.md` (no fingerprinted file reads it),
added `IDEAS.md` and the top-level `check.jl`, and documented the layout exceptions; no
fingerprinted file, result or experiment record changed.
